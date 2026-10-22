// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math"
	"math/big"
	"math/bits"
)

// expSplit returns the number s of halvings before the Taylor series of
// exp: Brent's square-and-reduce (MCA §4.3.1) balances about sqrt(w/2)
// squarings against the series length.
func expSplit(w uint) uint {
	return uint(math.Sqrt(float64(w)/2)) + 1
}

// taylorTerms returns the smallest n with 2^(1 + lt·(n+1)) / (n+1)! below
// one ulp of 2^-w, using floor(log2 j) <= log2 j for a rigorous lower bound
// on log2((n+1)!). lt <= -1 is a bound on log2|t|.
func taylorTerms(lt int64, w uint) int64 {
	n := int64(1)
	for f := int64(0); 1+lt*(n+1)-f-int64(bitlenU(n+1)-1) > -int64(w); n++ {
		f += int64(bitlenU(n+1) - 1)
	}
	return n
}

// expCore returns e^r for a fixed-point ball: Taylor series of e^(r/2^s')
// by Horner's rule, then s' squarings, where s' >= s makes |r/2^s'| <= 1/2
// whatever r. The tail Σ_{k>n} |t|^k/k! is at most 2|t|^(n+1)/(n+1)! for
// |t| <= 1/2. The result goes to y (distinct from r; r's scratch is
// written).
func expCore(y, r *fx, s uint) {
	w := r.w
	s += uint(max(0, r.log2bnd()+1))
	ws := getFxs(1)
	defer putFxs(ws)
	t := &(*ws)[0]
	t.shr(r, s)
	n := taylorTerms(min(int64(t.log2bnd()), -1), w)
	expSeries(y, t, n)
	y.e = addUp(y.e, 1) // truncation
	for i := uint(0); i < s; i++ {
		y.mul(y, y)
	}
}

// expSeries returns Σ_{k=0}^{n} t^k/k! for a fixed-point ball t, |t| <= 1/2
// (the tail is the caller's). As oddSeries, the terms go to m accumulators
// by k mod m, acc_i holding them without their factor t^i: the running term
// a = t^(k-i)/k! takes one division per term and one multiplication by t^m
// every m terms (mpmath's exp_basecase runs two such sums, libelefun.py).
// The result goes to dst (distinct from t).
func expSeries(dst, t *fx, n int64) {
	w := t.w
	m := max(2, int(math.Sqrt(float64(n))))
	ws := getFxs(2*m + 3)
	defer putFxs(ws)
	buf := *ws
	pw, acc := buf[:m+1], buf[m+1:2*m+1]
	a, q := &buf[2*m+1], &buf[2*m+2]
	pw[1].set(t)
	for i := 2; i <= m; i++ {
		pw[i].mul(&pw[i-1], t)
	}
	acc[0].setInt(1, w)
	for i := 1; i < m; i++ {
		acc[i].setInt(0, w)
	}
	a.setInt(1, w)
	for k := int64(1); k <= n; k++ {
		i := int(k % int64(m))
		if i == 0 {
			a.mul(a, &pw[m])
		}
		a.divSmall(a, uint64(k))
		acc[i].add(&acc[i], a)
	}
	r := &acc[0]
	for i := 1; i < m; i++ {
		q.mul(&acc[i], &pw[i])
		r.add(r, q)
	}
	dst.set(r)
}

func bitlenU(n int64) int { return bits.Len64(uint64(n)) } // n >= 0

// expReduced encloses e^x = b·2^k, for x within xr of an exact x with
// |x| < 2^34: x = k·ln2 + r with |r| <= 0.36 (plus xr), then expCore on r.
// Absolute error in r becomes relative error in the result, so r is formed
// to wp + 16 bits absolute; the radius xr is carried into the fixed-point
// ball, whose error propagation holds for any point of it.
func expReduced(x *big.Float, xr mag, wp uint) (*ball, int64) {
	xf, _ := x.Float64()
	k := int64(math.Round(xf / math.Ln2))
	s := expSplit(wp)
	w := wp + s + 16
	// r = x - k·ln2 in fixed point at scale wk: the error of k·ln2 is |k|
	// ulps of the cached ln2 (plus one), which the bitlen(k) + 8 extra bits
	// bring below one ulp at scale w.
	ws := getFxs(2)
	defer putFxs(ws)
	r, kl := &(*ws)[0], &(*ws)[1]
	wk := w + uint(bitlenU(abs64(k))) + 8
	r.setFloat(x, wk)
	if k != 0 {
		kl.mulInt(ln2Fx(wk), k)
		r.sub(r, kl)
	}
	r.rescale(r, w)
	if !xr.isZero() {
		// carry the radius of x over, in ulps rounded up (fxFromBall of
		// the ball 0 ± xr: one ulp more than needed)
		var rb ball
		rb.r = xr
		r.e = addUp(r.e, fxFromBall(&rb, w).e)
	}
	b := new(ball)
	expCore(kl, r, s)
	kl.toBall(b)
	return b, k
}

func abs64(k int64) int64 {
	if k < 0 {
		return -k
	}
	return k
}

// Exp sets z to the correctly rounded e^x and returns z. Exp(±0) = 1,
// Exp(+Inf) = +Inf, Exp(-Inf) = +0; results beyond big.Float's exponent
// range overflow to +Inf or underflow to +0 as big.Float's arithmetic does.
func Exp(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch {
	case x.IsInf() && x.Signbit():
		return z.SetInt64(0)
	case x.IsInf():
		return z.SetInf(false)
	case x.Sign() == 0:
		return z.SetInt64(1)
	}
	p := z.Prec()
	ex := int64(x.MantExp(nil))
	one := big.NewFloat(1)
	switch {
	case ex+1 <= nearBound(one, p):
		// |e^x - 1| <= 2|x| < 2^(ex+1) for |x| <= ln 2.
		return finish(z, proxyNear(one, x.Sign(), p), 0)
	case ex > 34 && x.Sign() > 0:
		return overflow(z, false)
	case ex > 34:
		return underflow(z, false)
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		b, k := expReduced(x, mag{}, wp)
		return b, k, true
	})
}

// expm1Core returns e^x - 1 for a fixed-point ball with |x| < 1/2:
// Taylor series of expm1(x/2^s), then s doublings E <- E(E+2).
func expm1Core(x *fx, s uint) *fx {
	w := x.w
	var t fx
	t.shr(x, s)
	n := taylorTerms(min(int64(t.log2bnd()), -1), w)
	var y, one, two fx
	one.setInt(1, w)
	two.setInt(2, w)
	y.setInt(1, w)
	for k := n; k >= 2; k-- {
		y.mul(&y, &t)
		y.divSmall(&y, uint64(k))
		y.add(&y, &one)
	}
	y.mul(&y, &t)
	y.e = addUp(y.e, 1)
	var u fx
	for i := uint(0); i < s; i++ {
		u.add(&y, &two)
		y.mul(&y, &u)
	}
	return &y
}

// expm1Split is the number of doublings for expm1 of an argument of
// exponent ex < 0: fewer for smaller arguments, none for tiny ones.
func expm1Split(wp uint, ex int64) uint {
	return uint(max(0, int64(expSplit(wp))+ex+1))
}

// Expm1 sets z to the correctly rounded e^x - 1 and returns z. Expm1(±0) =
// ±0, Expm1(+Inf) = +Inf, Expm1(-Inf) = -1.
func Expm1(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch {
	case x.IsInf() && x.Signbit():
		return z.SetInt64(-1)
	case x.IsInf():
		return z.SetInf(false)
	case x.Sign() == 0:
		return z.Set(x)
	}
	p := z.Prec()
	ex := int64(x.MantExp(nil))
	m1 := big.NewFloat(-1)
	xf, _ := x.Float64() // ±Inf beyond float64's range: only used when ex <= 34
	switch {
	case 2*ex <= nearBound(x, p):
		// e^x - 1 - x is in (0, x²) for 0 < |x| <= 1/2.
		return finish(z, proxyNear(x, 1, p), 0)
	case ex > 34 && x.Sign() > 0:
		return overflow(z, false)
	case x.Sign() < 0 && (ex > 34 || xf*math.Log2E < float64(nearBound(m1, p))-1):
		// e^x < 2^(nearBound(-1)): -1 proxies, from above.
		return finish(z, proxyNear(m1, 1, p), 0)
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		if ex >= 0 {
			// |x| >= 1/2: e^x - 1 = 2^k(b - 2^-k), |e^x - 1| > 0.39.
			b, k := expReduced(x, mag{}, wp+4)
			if k > int64(wp)+64 {
				// 2^-k is far below the rounding of b (b >= 1/2): a radius.
				// (int(-k) would also wrap on 32-bit targets for k >= 2^31.)
				b.addRad(magPow2(-k))
				return b, k, true
			}
			var c ball
			c.m.SetMantExp(big.NewFloat(1), int(-k))
			b.sub(b, &c, wp+4)
			return b, k, true
		}
		s := expm1Split(wp, ex)
		w := wp + uint(-ex) + 2*s + 16
		var b ball
		expm1Core(fxFromFloat(x, w), s).toBall(&b)
		return &b, 0, true
	})
}
