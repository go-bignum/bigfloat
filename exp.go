// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math"
	"math/big"
)

// expSplit returns the number s of halvings before the Taylor series of
// exp: Brent's square-and-reduce (MCA §4.3.1) balances about sqrt(w/2)
// squarings against the series length.
func expSplit(w uint) uint {
	return uint(math.Sqrt(float64(w)/2)) + 1
}

// expCore returns e^r for a fixed-point ball with |r| < 1: Taylor series of
// e^(r/2^s) by Horner's rule, then s squarings.
func expCore(r *fx, s uint) *fx {
	w := r.w
	var t fx
	t.shr(r, s)
	// |t| <= 2^lt, lt <= -s. Tail Σ_{k>n} |t|^k/k! <= 2|t|^(n+1)/(n+1)!
	// <= 2^(1 + lt(n+1) - log2fact(n+1)): stop when below one ulp.
	lt := int64(t.log2bnd())
	if lt > -1 {
		lt = -1
	}
	n := int64(1)
	for f := int64(0); 1+lt*(n+1)-f-int64(bitlenU(n+1)-1) > -int64(w); n++ {
		f += int64(bitlenU(n+1) - 1)
	}
	var y fx
	y.setInt(1, w)
	var one fx
	one.setInt(1, w)
	for k := n; k >= 1; k-- {
		y.mul(&y, &t)
		y.divSmall(&y, uint64(k))
		y.add(&y, &one)
	}
	y.e = addUp(y.e, 1) // truncation
	for i := uint(0); i < s; i++ {
		y.mul(&y, &y)
	}
	return &y
}

func bitlenU(n int64) int { return big.NewInt(n).BitLen() }

// expReduced encloses e^x = b·2^k for an exact x with |x| < 2^34: x = k·ln2
// + r with |r| <= 0.36, then expCore on r. Absolute error in r becomes
// relative error in the result, so r is formed to wp + 10 bits absolute.
func expReduced(x *big.Float, xr mag, wp uint) (*ball, int64, bool) {
	xf, _ := x.Float64()
	k := int64(math.Round(xf / math.Ln2))
	s := expSplit(wp)
	w := wp + s + 16
	var r ball
	if k != 0 {
		l := ln2Ball(w + uint(bitlenU(abs64(k))) + 8)
		var kb, t ball
		kb.setInt64(k)
		t.mul(l, &kb, w+uint(bitlenU(abs64(k)))+8)
		var xb ball
		xb.setFloat(x, w+48)
		r.sub(&xb, &t, w+8)
	} else {
		r.setFloat(x, w+8)
	}
	r.addRad(xr)
	if u := r.upper(); u.e > 0 {
		return nil, 0, false // |r| may reach 1: cannot happen with a sane k
	}
	y := expCore(fxFromBall(&r, w), s)
	var b ball
	y.toBall(&b)
	return &b, k, true
}

func abs64(k int64) int64 {
	if k < 0 {
		return -k
	}
	return k
}

// expBall encloses e^t = b·2^K for a ball t with |t| < 2^33 and radius at
// most 1/2: e^T is within e^t·(e^r - 1) <= 2r·e^t of e^t... which the
// reduced evaluation already carries, since the radius is folded into r.
func expBall(t *ball, wp uint) (*ball, int64, bool) {
	if t.r.e > 0 || t.r.e == 0 && t.r.m > 0.5 {
		return nil, 0, false
	}
	return expReduced(&t.m, t.r, wp)
}

// Exp sets z to the correctly rounded e^x and returns z. Exp(±0) = 1,
// Exp(+Inf) = +Inf, Exp(-Inf) = +0; results beyond big.Float's exponent
// range overflow to +Inf or underflow to +0 as big.Float's arithmetic does.
func Exp(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch {
	case x.IsInf():
		if x.Signbit() {
			return z.SetInt64(0)
		}
		return z.SetInf(false)
	case x.Sign() == 0:
		return z.SetInt64(1)
	}
	p := z.Prec()
	ex := int64(x.MantExp(nil))
	one := big.NewFloat(1)
	// |e^x - 1| <= 2|x| < 2^(ex+1) for |x| <= ln 2.
	if ex+1 <= nearBound(one, p) {
		return finish(z, proxyNear(one, x.Sign(), p), 0)
	}
	if ex > 34 {
		if x.Sign() > 0 {
			return overflow(z, false)
		}
		return underflow(z, false)
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		return expReduced(x, mag{}, wp)
	})
}

// expm1Core returns e^x - 1 for a fixed-point ball with |x| < 1/2:
// Taylor series of expm1(x/2^s), then s doublings E <- E(E+2).
func expm1Core(x *fx, s uint) *fx {
	w := x.w
	var t fx
	t.shr(x, s)
	lt := int64(t.log2bnd())
	if lt > -1 {
		lt = -1
	}
	n := int64(1)
	for f := int64(0); 1+lt*(n+1)-f-int64(bitlenU(n+1)-1) > -int64(w); n++ {
		f += int64(bitlenU(n+1) - 1)
	}
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

// Expm1 sets z to the correctly rounded e^x - 1 and returns z. Expm1(±0) =
// ±0, Expm1(+Inf) = +Inf, Expm1(-Inf) = -1.
func Expm1(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch {
	case x.IsInf():
		if x.Signbit() {
			return z.SetInt64(-1)
		}
		return z.SetInf(false)
	case x.Sign() == 0:
		return z.Set(x)
	}
	p := z.Prec()
	ex := int64(x.MantExp(nil))
	// e^x - 1 - x is in (0, x²) for 0 < |x| <= 1/2.
	if 2*ex <= nearBound(x, p) {
		return finish(z, proxyNear(x, 1, p), 0)
	}
	if ex > 34 && x.Sign() > 0 {
		return overflow(z, false)
	}
	m1 := big.NewFloat(-1)
	// x < 0 with e^x < 2^(nearBound(-1)): -1 proxies, from above.
	if x.Sign() < 0 {
		xf, _ := x.Float64() // finite enough: |x| below 2^34 or -Inf-like
		if ex > 34 || xf*math.Log2E < float64(nearBound(m1, p))-1 {
			return finish(z, proxyNear(m1, 1, p), 0)
		}
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		if ex >= 0 {
			// |x| >= 1/2: e^x - 1 = 2^k(b - 2^-k), |e^x - 1| > 0.39.
			b, k, ok := expReduced(x, mag{}, wp+4)
			if !ok {
				return nil, 0, false
			}
			var c ball
			c.m.SetMantExp(big.NewFloat(1), int(-k))
			b.sub(b, &c, wp+4)
			return b, k, true
		}
		s := uint(max(0, int64(expSplit(wp))+ex+1))
		w := wp + uint(-ex) + 2*s + 16
		y := expm1Core(fxFromFloat(x, w), s)
		var b ball
		y.toBall(&b)
		return &b, 0, true
	})
}
