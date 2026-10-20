// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"fmt"
	"math"
	"math/big"
)

// trigSplit returns the number of argument halvings before the sin/cos and
// atan series.
func trigSplit(w uint) uint { return uint(math.Sqrt(float64(w))/3) + 1 }

// sincosCore returns sin r and cos r for a fixed-point ball with |r| < 1:
// both Taylor series at t = r/2^j (alternating, decreasing terms: the tail
// is below the first omitted term), then j doublings sin 2a = 2 sin a cos a,
// cos 2a = 1 - 2 sin² a.
func sincosCore(r *fx, j uint) (*fx, *fx) {
	w := r.w
	var t fx // log2bnd writes its scratch: a local, not the caller's r
	t.shr(r, j)
	lt := min(int64(t.log2bnd()), -1)
	// cos tail <= |t|^(2n+2)/(2n+2)!, sin tail smaller.
	n := int64(1)
	for lt*(2*n+2)-log2FactLower(2*n+2) > -int64(w)-1 {
		n++
	}
	ws := getFxs(2)
	defer putFxs(ws)
	one, s2 := &(*ws)[0], &(*ws)[1]
	s, c := sinCosSeries(&t, n)
	s.e = addUp(s.e, 1)
	c.e = addUp(c.e, 1)
	one.setInt(1, w)
	for i := uint(0); i < j; i++ {
		s2.mul(s, s)
		s.mul(s, c)
		s.shl(s, 1)
		c.shl(s2, 1)
		c.sub(one, c)
	}
	return s, c
}

// sinCosSeries returns the Taylor polynomials of sin t and cos t to the
// terms k = n (t^(2n+1)/(2n+1)! and t^(2n)/(2n)!), tails left to the
// caller. Both run on m accumulators of u = t² as oddSeries does, sharing
// the powers of u: the running cos term a = u^(k-i)/(2k)! takes one
// division, the sin term is a/(2k+1), and a advances by u^m every m terms.
func sinCosSeries(t *fx, n int64) (*fx, *fx) {
	w := t.w
	m := max(2, int(math.Sqrt(float64(n)))&^1)
	ws := getFxs(3*m + 4)
	defer putFxs(ws)
	buf := *ws
	pw, ac, as := buf[:m+1], buf[m+1:2*m+1], buf[2*m+1:3*m+1]
	a, b, q := &buf[3*m+1], &buf[3*m+2], &buf[3*m+3]
	pw[1].mul(t, t)
	for i := 2; i <= m; i++ {
		pw[i].mul(&pw[i-1], &pw[1])
	}
	ac[0].setInt(1, w)
	as[0].setInt(1, w)
	for i := 1; i < m; i++ {
		ac[i].setInt(0, w)
		as[i].setInt(0, w)
	}
	a.setInt(1, w)
	for k := int64(1); k <= n; k++ {
		i := int(k % int64(m))
		if i == 0 {
			a.mul(a, &pw[m])
		}
		a.divSmall(a, uint64((2*k-1)*(2*k)))
		ac[i].add(&ac[i], a)
		b.divSmall(a, uint64(2*k+1))
		as[i].add(&as[i], b)
	}
	c, s := &ac[0], &as[0]
	for i := 1; i < m; i++ {
		// m even: the sign of term jm+i is (-1)^i
		q.mul(&ac[i], &pw[i])
		if i%2 == 1 {
			c.sub(c, q)
		} else {
			c.add(c, q)
		}
		q.mul(&as[i], &pw[i])
		if i%2 == 1 {
			s.sub(s, q)
		} else {
			s.add(s, q)
		}
	}
	S := new(fx).mul(s, t)
	return S, new(fx).set(c)
}

var quarterPi = big.NewFloat(0.78)

// MaxReductionExponent bounds the binary exponent of the arguments of Sin,
// Cos and Tan: an argument near 2^e is reduced modulo π/2 with about e
// bits of π, which costs time and memory growing with e (the first such
// call at 2^22 computes four million bits of π, in seconds, later calls
// reuse the cache). Beyond the bound these functions panic with ErrLimit
// instead of consuming unbounded resources on an untrusted input. Set it
// before use, not concurrently with calls.
var MaxReductionExponent = 1 << 22

// ErrLimit is the panic value of a function asked for work beyond one of
// this package's documented resource bounds (MaxReductionExponent).
type ErrLimit struct {
	Msg string
}

func (e ErrLimit) Error() string { return e.Msg }

// reduceHalfPi writes x = q·π/2 + r with |r| <= π/4 + ε and returns r as
// a ball of absolute accuracy about 2^-w, and q mod 4. π is taken to as
// many bits as x's exponent requires, so that the reduction stays exact
// for any argument (MPFR's mpfr_sin does the same); the first call at a
// huge exponent pays for computing those bits of π, later ones reuse the
// cache.
//
// r is a fixed-point ball at scale w, within two ulps: the reduction runs
// at scale W = w + bitlen(k) + 9, where the |k| ulps of error of k·π/2 are
// below one ulp of 2^-w.
func reduceHalfPi(x *big.Float, w uint) (*fx, int) {
	r := new(fx)
	var ax big.Float
	if ax.Abs(x).Cmp(quarterPi) <= 0 {
		return r.setFloat(x, w), 0
	}
	if e := x.MantExp(nil); e > MaxReductionExponent {
		panic(ErrLimit{Msg: fmt.Sprintf("bigfloat: argument 2^%d beyond MaxReductionExponent (2^%d): reducing it needs that many bits of π", e, MaxReductionExponent)})
	}
	ex := uint(max(x.MantExp(nil), 0))
	var kint *big.Int
	if ex <= 40 {
		// float64 quotient: absolute error below 2^-11, far below 1/2
		xf, _ := x.Float64()
		kint = big.NewInt(int64(math.Round(xf / (math.Pi / 2))))
	} else {
		// k = round(x / (π/2)) from a rough π: error far below 1/2.
		pr := piBall(ex + 64)
		var hp, kf big.Float
		hp.SetPrec(ex+64).SetMantExp(&pr.m, -1)
		kf.SetPrec(ex+64).Quo(x, &hp)
		kint = roundInt(&kf)
	}
	var q4 big.Int
	q := int(q4.And(kint, big.NewInt(3)).Int64())
	W := w + uint(kint.BitLen()) + 8
	// the integer of π at scale W is π/2 at scale W+1
	var kl fx
	kl.mulInt(piFx(W), kint)
	kl.w = W + 1
	r.setFloat(x, W+1)
	r.sub(r, &kl)
	return r.rescale(r, w), q
}

// roundInt rounds x to an integer near it (ties away; any rounding will
// do for a reduction quotient).
func roundInt(x *big.Float) *big.Int {
	var h big.Float
	h.SetPrec(x.Prec() + 2)
	if x.Sign() >= 0 {
		h.Add(x, big.NewFloat(0.5))
	} else {
		h.Sub(x, big.NewFloat(0.5))
	}
	n, _ := h.Int(nil)
	return n
}

// sincosBall returns balls for sin x and cos x (exact x), each with about
// wp bits of relative accuracy. The reduced argument r needs as many more
// bits of absolute accuracy as it has leading zeros; ok = false when the
// reduction at this precision cannot even separate r from zero (x within
// 2^-wp of a multiple of π/2), so that the Ziv loop raises the precision.
func sincosBall(x *big.Float, wp uint) (*ball, *ball, bool) {
	j := trigSplit(wp)
	base := wp + 2*j + 20
	w := base
	var r *fx
	var q int
	for {
		r, q = reduceHalfPi(x, w+8)
		er, ok := r.expLower()
		if !ok {
			return nil, nil, false
		}
		need := base + uint(max(0, -er))
		if need <= w {
			break
		}
		w = need
	}
	s, c := sincosCore(r.rescale(r, w), j)
	var sb, cb ball
	s.toBall(&sb)
	c.toBall(&cb)
	switch q {
	case 1:
		sb, cb = cb, sb
		cb.neg(&cb)
	case 2:
		sb.neg(&sb)
		cb.neg(&cb)
	case 3:
		sb, cb = cb, sb
		sb.neg(&sb)
	}
	return &sb, &cb, true
}

// Sin sets z to the correctly rounded sin x and returns z. Sin(±0) = ±0;
// Sin(±Inf) panics with ErrNaN.
func Sin(z, x *big.Float) *big.Float { return trig(z, x, 0) }

// Cos sets z to the correctly rounded cos x and returns z. Cos(±0) = 1;
// Cos(±Inf) panics with ErrNaN.
func Cos(z, x *big.Float) *big.Float { return trig(z, x, 1) }

// Tan sets z to the correctly rounded tan x and returns z. Tan(±0) = ±0;
// Tan(±Inf) panics with ErrNaN. (No binary number is an odd multiple of
// π/2, so the result is always finite.)
func Tan(z, x *big.Float) *big.Float { return trig(z, x, 2) }

func trig(z, x *big.Float, which int) *big.Float {
	precOf(z, x.Prec())
	if x.IsInf() {
		nan([]string{"Sin", "Cos", "Tan"}[which] + " of an infinity")
	}
	if x.Sign() == 0 {
		if which == 1 {
			return z.SetInt64(1)
		}
		return z.Set(x)
	}
	p := z.Prec()
	ex := int64(x.MantExp(nil))
	one := big.NewFloat(1)
	switch {
	case which == 0 && 3*ex <= nearBound(x, p):
		// x - sin x in (0, x³/6) for x > 0
		return finish(z, proxyNear(x, -x.Sign(), p), 0)
	case which == 1 && 2*ex-1 <= nearBound(one, p):
		// 1 - cos x in (0, x²/2)
		return finish(z, proxyNear(one, -1, p), 0)
	case which == 2 && 3*ex <= nearBound(x, p):
		// tan x - x in (0, |x|³) for 0 < x < 1/2
		return finish(z, proxyNear(x, x.Sign(), p), 0)
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		s, c, ok := sincosBall(x, wp)
		switch {
		case !ok:
			return nil, 0, false
		case which == 0:
			return s, 0, true
		case which == 1:
			return c, 0, true
		}
		// c has relative accuracy: it excludes zero
		var t ball
		must(t.quo(s, c, wp))
		return &t, 0, true
	})
}

// atanCore returns atan t for a fixed-point ball with |t| <= 1: j
// reductions t <- t/(1+sqrt(1+t²)) (atan t = 2 atan of the result), then
// the alternating series Σ (-1)^k t^(2k+1)/(2k+1), whose tail is below
// |t|^(2n+3).
func atanCore(t *fx, j uint) *fx {
	w := t.w
	var y, one, u fx
	one.setInt(1, w)
	y.set(t)
	for i := uint(0); i < j; i++ {
		u.mul(&y, &y)
		u.add(&u, &one)
		must(u.sqrt(&u))
		u.add(&u, &one)
		must(y.quo(&y, &u))
	}
	// lt >= 1-w: log2bnd counts the error bound
	lt := min(int64(y.log2bnd()), -1)
	n := (int64(w)+1+(-lt)-1)/(-lt)/2 + 1
	s := oddSeries(&y, n, true)
	s.e = addUp(s.e, 1) // tail
	return s.shl(s, j)
}

// atanShift sets the grid of the atan table: atan(n/2^atanShift), as
// mpmath's ATAN_TAYLOR_SHIFT.
const atanShift = 7

// atanTable caches atan(n/2^7), n = 1 … 2^7, computed by atanCore.
var atanTable = &fxTable{compute: func(n int64, W uint) *fx {
	var a fx
	a.setInt(n, W)
	a.v.Rsh(&a.v, atanShift) // exact: n/2^7 at scale W > 7
	return atanCore(&a, trigSplit(W))
}}

// atanFx returns atan t for a fixed-point ball with |t| <= 1 (up to its
// error bound). Up to tableMaxPrec it uses mpmath's table reduction
// (atan_taylor, libelefun.py): with a = n/2^7 the grid point below |t|,
// atan t = atan a + atan d, d = (t-a)/(1+a·t) in [0, 2^-7), so that the
// series gains 14 bits per term with no square root; atan a comes from the
// table. Above, or for |t| < 2^-7, atanCore with j halvings.
func atanFx(t *fx, j uint) *fx {
	w := t.w
	if w > tableMaxPrec {
		return atanCore(t, j)
	}
	neg := t.v.Sign() < 0
	var at fx
	at.set(t)
	if neg {
		at.v.Neg(&at.v)
	}
	var nn fx
	nn.v.Rsh(&at.v, w-atanShift)
	n := min(nn.v.Int64(), 1<<atanShift)
	if n == 0 {
		return atanCore(t, 0)
	}
	var a, one, num, den fx
	a.setInt(n, w)
	a.v.Rsh(&a.v, atanShift)
	one.setInt(1, w)
	num.sub(&at, &a)
	den.mul(&a, &at)
	den.add(&den, &one)
	must(num.quo(&num, &den)) // den >= 1 - ulps
	s := atanCore(&num, 0)
	s.add(s, atanTable.get(n, w))
	if neg {
		s.v.Neg(&s.v)
	}
	return s
}

// atanBall encloses atan t for a ball t that excludes zero, with about wp
// bits of relative accuracy: |atan T - atan t| <= |T - t|.
func atanBall(t *ball, wp uint) *ball {
	el, ok := t.expLower()
	must(ok)
	z := new(ball)
	eu := t.expUpper()
	j := trigSplit(wp)
	w := wp + 2*j + 16
	switch {
	case eu < -int64(wp/2)-4:
		// atan t = t(1 - t²/3 + ...): within |t|³ of t.
		z.set(t)
		z.addRad(magMul(magMul(t.upper(), t.upper()), t.upper()))
	case eu <= 0 || cmpAbs(&t.m, big.NewFloat(1)) <= 0:
		w += uint(max(0, -el))
		atanFx(fxFromBall(t, w), j).toBall(z)
	default:
		// |t| > 1: sign(t)·π/2 - atan(1/t), absolute accuracy suffices.
		var one, u, a ball
		one.setInt64(1)
		must(u.quo(&one, t, w))
		atanFx(fxFromBall(&u, w), j).toBall(&a)
		h := halfPiBall(w, t.m.Sign() < 0)
		z.sub(h, &a, wp+8)
	}
	return z
}

// halfPiBall encloses ±π/2.
func halfPiBall(w uint, neg bool) *ball {
	h := piBall(w)
	h.mul2exp(h, -1)
	if neg {
		h.neg(h)
	}
	return h
}

// atan2Ball encloses atan2(y, x) for exact non-zero finite y, x.
func atan2Ball(y, x *big.Float, wp uint) *ball {
	var yb, xb, q ball
	yb.setExact(y)
	xb.setExact(x)
	w := wp + 16
	if cmpAbs(y, x) <= 0 {
		must(q.quo(&yb, &xb, w))
		a := atanBall(&q, w)
		if x.Sign() < 0 {
			pb := piBall(w)
			if y.Sign() < 0 {
				pb.neg(pb)
			}
			a.add(a, pb, w)
		}
		return a
	}
	must(q.quo(&xb, &yb, w))
	a := atanBall(&q, w)
	return a.sub(halfPiBall(w, y.Sign() < 0), a, w)
}

// Atan sets z to the correctly rounded arctangent of x and returns z.
// Atan(±0) = ±0, Atan(±Inf) = ±π/2.
func Atan(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch {
	case x.IsInf():
		return piMultiple(z, x.Signbit(), 1, 1)
	case x.Sign() == 0:
		return z.Set(x)
	}
	p := z.Prec()
	// x - atan x in (0, x³/3) for x > 0
	if ex := int64(x.MantExp(nil)); 3*ex <= nearBound(x, p) {
		return finish(z, proxyNear(x, -x.Sign(), p), 0)
	}
	return atan2Exact(z, x, big.NewFloat(1))
}

// piMultiple sets z to ±k·π/2^d (correctly rounded): the special values
// of Atan, Atan2, Asin and Acos.
func piMultiple(z *big.Float, neg bool, k int64, d int) *big.Float {
	precOf(z, 64)
	return ziv(z, 0, func(wp uint) (*ball, int64, bool) {
		b := piBall(wp + 8)
		var kb ball
		kb.setInt64(k)
		b.mul(b, &kb, wp+8)
		b.mul2exp(b, -d)
		if neg {
			b.neg(b)
		}
		return b, 0, true
	})
}

func atan2Exact(z, y, x *big.Float) *big.Float {
	y = new(big.Float).Set(y)
	x = new(big.Float).Set(x)
	return ziv(z, y.Prec()+x.Prec(), func(wp uint) (*ball, int64, bool) {
		return atan2Ball(y, x, wp), 0, true
	})
}

// Atan2 sets z to the correctly rounded angle of the point (x, y), in
// [-π, π], and returns z. If z's precision is 0 it is set to the larger of
// x's and y's. The special values follow C99 Annex F (atan2(±0, +0) = ±0,
// atan2(±0, -0) = ±π, atan2(±Inf, ±Inf) = ±π/4 or ±3π/4, …).
func Atan2(z, y, x *big.Float) *big.Float {
	precOf(z, max(x.Prec(), y.Prec()))
	ny := y.Signbit()
	switch {
	case y.IsInf() && x.IsInf() && x.Signbit():
		return piMultiple(z, ny, 3, 2)
	case y.IsInf() && x.IsInf():
		return piMultiple(z, ny, 1, 2)
	case y.IsInf() || x.Sign() == 0 && y.Sign() != 0:
		return piMultiple(z, ny, 1, 1)
	case (x.IsInf() || y.Sign() == 0) && x.Signbit():
		return piMultiple(z, ny, 1, 0)
	case x.IsInf() || y.Sign() == 0:
		return zeroSign(z.SetInt64(0), ny)
	}
	p := z.Prec()
	if x.Sign() > 0 && y.MantExp(nil) < x.MantExp(nil) {
		// atan(y/x) with y/x dyadic and tiny: proxy as Atan does.
		q := new(big.Float).SetPrec(y.Prec() + x.Prec() + p + 8)
		q.Quo(y, x) // exact when y/x is dyadic: it then has at most prec(y) bits
		exact := q.Acc() == big.Exact
		q.SetPrec(q.MinPrec()) // the shortest base allows the widest proxy
		eq := int64(q.MantExp(nil))
		if exact && 3*eq <= nearBound(q, p) {
			return finish(z, proxyNear(q, -q.Sign(), p), 0)
		}
	}
	return atan2Exact(z, y, x)
}

// zeroSign sets z (zero) to -0 when neg.
func zeroSign(z *big.Float, neg bool) *big.Float {
	if neg {
		z.Neg(z)
	}
	return z
}

// sqrt1mx2 encloses sqrt(1 - x²) for exact |x| < 1, from the exact
// product (1-x)(1+x) > 0.
func sqrt1mx2(x *big.Float, wp uint) *ball {
	pr := x.Prec() + uint(max(-x.MantExp(nil), 0)) + 4
	a := new(big.Float).SetPrec(pr).Sub(big.NewFloat(1), x)
	b := new(big.Float).SetPrec(pr).Add(big.NewFloat(1), x)
	c := new(big.Float).SetPrec(2*pr).Mul(a, b)
	var cb, s ball
	cb.setExact(c)
	must(s.sqrt(&cb, wp))
	return &s
}

// Asin sets z to the correctly rounded arcsine of x and returns z.
// Asin(±0) = ±0, Asin(±1) = ±π/2; |x| > 1 panics with ErrNaN.
func Asin(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch c := cmpAbs(x, big.NewFloat(1)); {
	case c > 0:
		nan("Asin of a number beyond [-1, 1]")
	case c == 0:
		return piMultiple(z, x.Signbit(), 1, 1)
	case x.Sign() == 0:
		return z.Set(x)
	}
	p := z.Prec()
	// asin x - x in (0, |x|³/4) for 0 < x <= 1/2
	if ex := int64(x.MantExp(nil)); 3*ex <= nearBound(x, p) {
		return finish(z, proxyNear(x, x.Sign(), p), 0)
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		s := sqrt1mx2(x, wp+16)
		var xb, q ball
		xb.setExact(x)
		must(q.quo(&xb, s, wp+16))
		return atanBall(&q, wp), 0, true
	})
}

// Acos sets z to the correctly rounded arccosine of x, in [0, π], and
// returns z. Acos(1) = +0, Acos(-1) = π; |x| > 1 panics with ErrNaN.
func Acos(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch c := cmpAbs(x, big.NewFloat(1)); {
	case c > 0:
		nan("Acos of a number beyond [-1, 1]")
	case c == 0 && x.Sign() > 0:
		return z.SetInt64(0)
	case c == 0:
		return piMultiple(z, false, 1, 0)
	case x.Sign() == 0:
		return piMultiple(z, false, 1, 1)
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		s := sqrt1mx2(x, wp+16)
		var xb, q ball
		xb.setExact(x)
		must(q.quo(s, &xb, wp+16))
		b := atanBall(&q, wp+8)
		if x.Sign() < 0 {
			b.add(b, piBall(wp+8), wp+8)
		}
		return b, 0, true
	})
}
