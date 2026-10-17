// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
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
	var t fx
	t.shr(r, j)
	lt := int64(t.log2bnd())
	if lt > -1 {
		lt = -1
	}
	// cos tail <= |t|^(2n+2)/(2n+2)!, sin tail smaller.
	n := int64(1)
	for lt*(2*n+2)-log2FactLower(2*n+2) > -int64(w)-1 {
		n++
	}
	var u, one, s, c fx
	u.mul(&t, &t)
	one.setInt(1, w)
	s.setInt(1, w)
	c.setInt(1, w)
	for k := n; k >= 1; k-- {
		s.mul(&s, &u)
		s.divSmall(&s, uint64(2*k*(2*k+1)))
		s.sub(&one, &s)
		c.mul(&c, &u)
		c.divSmall(&c, uint64((2*k-1)*(2*k)))
		c.sub(&one, &c)
	}
	s.mul(&s, &t)
	s.e = addUp(s.e, 1)
	c.e = addUp(c.e, 1)
	var s2 fx
	for i := uint(0); i < j; i++ {
		s2.mul(&s, &s)
		s.mul(&s, &c)
		s.shl(&s, 1)
		c.shl(&s2, 1)
		c.sub(&one, &c)
	}
	return &s, &c
}

var quarterPi = big.NewFloat(0.78)

// reduceHalfPi writes x = q·π/2 + r with |r| <= π/4 + ε and returns r as
// a ball of absolute accuracy about 2^-w, and q mod 4. π is taken to as
// many bits as x's exponent requires, so that the reduction stays exact
// for any argument (MPFR's mpfr_sin does the same); the first call at a
// huge exponent pays for computing those bits of π, later ones reuse the
// cache.
func reduceHalfPi(x *big.Float, w uint) (*ball, int) {
	r := new(ball)
	var ax big.Float
	if ax.Abs(x).Cmp(quarterPi) <= 0 {
		r.setExact(x)
		return r, 0
	}
	ex := uint(max(x.MantExp(nil), 0))
	// k = round(x / (π/2)) from a rough π: error far below 1/2.
	pr := piBall(ex + 64)
	var hp, kf big.Float
	hp.SetPrec(ex+64).SetMantExp(&pr.m, -1)
	kf.SetPrec(ex+64).Quo(x, &hp)
	kint, _ := roundInt(&kf)
	q := int(new(big.Int).And(kint, big.NewInt(3)).Int64())
	wp := ex + w + 16
	pb := piBall(wp)
	var kb, t, xb ball
	kb.m.SetPrec(0).SetInt(kint)
	pb.mul2exp(pb, -1)
	t.mul(pb, &kb, wp)
	xb.setExact(x)
	r.sub(&xb, &t, w+16)
	return r, q
}

// roundInt rounds x to the nearest integer (ties away; any rounding will do
// for a reduction quotient).
func roundInt(x *big.Float) (*big.Int, big.Accuracy) {
	var h big.Float
	h.SetPrec(x.Prec() + 2)
	if x.Sign() >= 0 {
		h.Add(x, big.NewFloat(0.5))
	} else {
		h.Sub(x, big.NewFloat(0.5))
	}
	return h.Int(nil)
}

// sincosBall returns balls for sin x and cos x (exact x), each with about
// wp bits of relative accuracy.
func sincosBall(x *big.Float, wp uint) (*ball, *ball, bool) {
	r, q := reduceHalfPi(x, wp+32)
	er, ok := r.expLower()
	if !ok {
		return nil, nil, false
	}
	j := trigSplit(wp)
	w := wp + 2*j + 20
	if er < 0 {
		w += uint(-er)
		r, q = reduceHalfPi(x, w+8)
		if _, ok := r.expLower(); !ok {
			return nil, nil, false
		}
	}
	s, c := sincosCore(fxFromBall(r, w), j)
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
	switch which {
	case 0: // x - sin x in (0, x³/6) for x > 0
		if 3*ex <= nearBound(x, p) {
			return finish(z, proxyNear(x, -x.Sign(), p), 0)
		}
	case 1: // 1 - cos x in (0, x²/2)
		if one := big.NewFloat(1); 2*ex-1 <= nearBound(one, p) {
			return finish(z, proxyNear(one, -1, p), 0)
		}
	case 2: // tan x - x in (0, |x|³) for 0 < x < 1/2
		if 3*ex <= nearBound(x, p) {
			return finish(z, proxyNear(x, x.Sign(), p), 0)
		}
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		s, c, ok := sincosBall(x, wp)
		if !ok {
			return nil, 0, false
		}
		switch which {
		case 0:
			return s, 0, true
		case 1:
			return c, 0, true
		}
		var t ball
		if !t.quo(s, c, wp) {
			return nil, 0, false
		}
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
		u.sqrt(&u)
		u.add(&u, &one)
		y.quo(&y, &u)
	}
	lt := int64(y.log2bnd())
	if lt > -1 {
		lt = -1
	}
	var s fx
	if lt >= -int64(w) {
		n := (int64(w)+1+(-lt)-1)/(-lt)/2 + 1
		u.mul(&y, &y)
		var c fx
		c.w = w
		coef := func(k int64) *fx {
			c.v.Lsh(big.NewInt(1), w)
			c.v.Quo(&c.v, big.NewInt(2*k+1))
			c.e = 1
			return &c
		}
		s.set(coef(n))
		for k := n - 1; k >= 0; k-- {
			s.mul(&s, &u)
			s.sub(coef(k), &s)
		}
		s.mul(&s, &y)
	} else {
		s.set(&y)
	}
	s.e = addUp(s.e, 1)
	return s.shl(&s, j)
}

// atanBall encloses atan t for a ball t, with about wp bits of relative
// accuracy: |atan T - atan t| <= |T - t|.
func atanBall(t *ball, wp uint) (*ball, bool) {
	el, ok := t.expLower()
	if !ok {
		return nil, false
	}
	z := new(ball)
	eu := t.expUpper()
	switch {
	case eu < -int64(wp/2)-4:
		// atan t = t(1 - t²/3 + ...): within |t|³ of t.
		z.set(t)
		z.addRad(magMul(magMul(t.upper(), t.upper()), t.upper()))
	case eu <= 0 || cmpAbs(&t.m, big.NewFloat(1)) <= 0:
		j := trigSplit(wp)
		w := wp + 2*j + 16
		if el < 0 {
			w += uint(-el)
		}
		atanCore(fxFromBall(t, w), j).toBall(z)
	default:
		// |t| > 1: sign(t)·π/2 - atan(1/t), absolute accuracy suffices.
		j := trigSplit(wp)
		w := wp + 2*j + 16
		var one, u ball
		one.setInt64(1)
		if !u.quo(&one, t, w) {
			return nil, false
		}
		var a ball
		atanCore(fxFromBall(&u, w), j).toBall(&a)
		h := piBall(w)
		h.mul2exp(h, -1)
		if t.m.Sign() < 0 {
			h.neg(h)
		}
		z.sub(h, &a, wp+8)
	}
	return z, true
}

// atan2Ball encloses atan2(y, x) for exact non-zero finite y, x.
func atan2Ball(y, x *big.Float, wp uint) (*ball, bool) {
	var yb, xb, q ball
	yb.setExact(y)
	xb.setExact(x)
	w := wp + 16
	if cmpAbs(y, x) <= 0 {
		q.quo(&yb, &xb, w)
		a, ok := atanBall(&q, w)
		if !ok {
			return nil, false
		}
		if x.Sign() < 0 {
			pb := piBall(w)
			if y.Sign() < 0 {
				pb.neg(pb)
			}
			a.add(a, pb, w)
		}
		return a, true
	}
	q.quo(&xb, &yb, w)
	a, ok := atanBall(&q, w)
	if !ok {
		return nil, false
	}
	h := piBall(w)
	h.mul2exp(h, -1)
	if y.Sign() < 0 {
		h.neg(h)
	}
	a.sub(h, a, w)
	return a, true
}

// Atan sets z to the correctly rounded arctangent of x and returns z.
// Atan(±0) = ±0, Atan(±Inf) = ±π/2.
func Atan(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	if x.IsInf() {
		return halfPiSigned(z, x.Signbit(), 1)
	}
	if x.Sign() == 0 {
		return z.Set(x)
	}
	p := z.Prec()
	// x - atan x in (0, x³/3) for x > 0
	if ex := int64(x.MantExp(nil)); 3*ex <= nearBound(x, p) {
		return finish(z, proxyNear(x, -x.Sign(), p), 0)
	}
	return atan2Exact(z, x, big.NewFloat(1))
}

// halfPiSigned sets z to ±k·π/2 (k = 1: ±π/2, k = 2: ±π); used for the
// special values of Atan, Atan2, Asin and Acos.
func halfPiSigned(z *big.Float, neg bool, k int64) *big.Float {
	precOf(z, 64)
	return ziv(z, 0, func(wp uint) (*ball, int64, bool) {
		b := piBall(wp + 8)
		var kb ball
		kb.setInt64(k)
		b.mul(b, &kb, wp+8)
		b.mul2exp(b, -1)
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
		b, ok := atan2Ball(y, x, wp)
		return b, 0, ok
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
	case y.IsInf() && x.IsInf():
		if x.Signbit() {
			return quarterPiTimes(z, ny, 3)
		}
		return quarterPiTimes(z, ny, 1)
	case y.IsInf():
		return halfPiSigned(z, ny, 1)
	case x.IsInf():
		if x.Signbit() {
			return halfPiSigned(z, ny, 2)
		}
		return zeroSign(z.SetInt64(0), ny)
	case y.Sign() == 0:
		if x.Signbit() {
			return halfPiSigned(z, ny, 2)
		}
		return zeroSign(z.SetInt64(0), ny)
	case x.Sign() == 0:
		return halfPiSigned(z, ny, 1)
	}
	p := z.Prec()
	if x.Sign() > 0 {
		// atan(y/x) with y/x dyadic and tiny: proxy as Atan does.
		ey, ex := int64(y.MantExp(nil)), int64(x.MantExp(nil))
		if d := ey - ex; d < 0 {
			q := new(big.Float).SetPrec(y.Prec() + x.Prec() + p + 8)
			q.Quo(y, x)
			eq := int64(q.MantExp(nil))
			if q.Acc() == big.Exact && 3*eq <= nearBound(q, p) {
				return finish(z, proxyNear(q, -q.Sign(), p), 0)
			}
		}
	}
	return atan2Exact(z, y, x)
}

// zeroSign sets z (zero) to -0 when neg.
func zeroSign(z *big.Float, neg bool) *big.Float {
	if neg != z.Signbit() {
		z.Neg(z)
	}
	return z
}

func quarterPiTimes(z *big.Float, neg bool, k int64) *big.Float {
	precOf(z, 64)
	return ziv(z, 0, func(wp uint) (*ball, int64, bool) {
		b := piBall(wp + 8)
		var kb ball
		kb.setInt64(k)
		b.mul(b, &kb, wp+8)
		b.mul2exp(b, -2)
		if neg {
			b.neg(b)
		}
		return b, 0, true
	})
}

// sqrt1mx2 encloses sqrt(1 - x²) for exact |x| < 1, from the exact
// product (1-x)(1+x).
func sqrt1mx2(x *big.Float, wp uint) *ball {
	pr := x.Prec() + uint(max(-x.MantExp(nil), 0)) + 4
	a := new(big.Float).SetPrec(pr).Sub(big.NewFloat(1), x)
	b := new(big.Float).SetPrec(pr).Add(big.NewFloat(1), x)
	c := new(big.Float).SetPrec(2*pr).Mul(a, b)
	var cb, s ball
	cb.setExact(c)
	s.sqrt(&cb, wp)
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
		return halfPiSigned(z, x.Signbit(), 1)
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
		if !q.quo(&xb, s, wp+16) {
			return nil, 0, false
		}
		b, ok := atanBall(&q, wp)
		return b, 0, ok
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
		return halfPiSigned(z, false, 2)
	case x.Sign() == 0:
		return halfPiSigned(z, false, 1)
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		s := sqrt1mx2(x, wp+16)
		var xb, q ball
		xb.setExact(x)
		if !q.quo(s, &xb, wp+16) {
			return nil, 0, false
		}
		b, ok := atanBall(&q, wp+8)
		if !ok {
			return nil, 0, false
		}
		if x.Sign() < 0 {
			b.add(b, piBall(wp+8), wp+8)
		}
		return b, 0, true
	})
}
