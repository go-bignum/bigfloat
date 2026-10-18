// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math/big"
)

// expPair returns e^|x| = b·2^k as a ball and c = 2^(-2k)/b, so that
// e^-|x| = c·2^k: sinh and cosh are then (b ∓ c)·2^(k-1), free of
// overflow before the final scaling.
func expPair(x *big.Float, wp uint) (b, c *ball, k int64) {
	ax := new(big.Float).Abs(x)
	b, k = expReduced(ax, mag{}, wp)
	c = new(ball)
	if -2*k < big.MinExp+64 {
		// 2^(-2k)/b would underflow: below 2^(-2k+1), it becomes a radius.
		c.setInt64(0)
		c.r = magPow2(-2*k + 1)
		return b, c, k
	}
	var t ball
	t.m.SetMantExp(big.NewFloat(1), int(-2*k))
	must(c.quo(&t, b, wp))
	return b, c, k
}

// Sinh sets z to the correctly rounded hyperbolic sine of x and returns z.
// Sinh(±0) = ±0, Sinh(±Inf) = ±Inf.
func Sinh(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	if x.IsInf() || x.Sign() == 0 {
		return z.Set(x)
	}
	p := z.Prec()
	ex := int64(x.MantExp(nil))
	switch {
	case 3*ex <= nearBound(x, p):
		// sinh x - x in (0, |x|³/5) for 0 < x <= 1/2
		return finish(z, proxyNear(x, x.Sign(), p), 0)
	case ex > 34:
		return overflow(z, x.Signbit())
	}
	x = new(big.Float).Set(x)
	neg := x.Sign() < 0
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		var r ball
		var K int64
		if ex < 0 {
			// |x| < 1/2: E = expm1|x|, sinh = E(E+2)/(2(E+1)), E+1 >= 1.
			w := wp + uint(-ex) + 24
			ax := new(big.Float).Abs(x)
			E := expm1Core(fxFromFloat(ax, w), expm1Split(wp, ex))
			var e1, e2, one, two fx
			one.setInt(1, w)
			two.setInt(2, w)
			e1.add(E, &one)
			e2.add(E, &two)
			e2.mul(E, &e2)
			must(e2.quo(&e2, &e1))
			e2.toBall(&r)
			K = -1
		} else {
			b, c, k := expPair(x, wp+8)
			r.sub(b, c, wp+8)
			K = k - 1
		}
		if neg {
			r.neg(&r)
		}
		return &r, K, true
	})
}

// Cosh sets z to the correctly rounded hyperbolic cosine of x and returns
// z. Cosh(±0) = 1, Cosh(±Inf) = +Inf.
func Cosh(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch {
	case x.IsInf():
		return z.SetInf(false)
	case x.Sign() == 0:
		return z.SetInt64(1)
	}
	p := z.Prec()
	ex := int64(x.MantExp(nil))
	one := big.NewFloat(1)
	switch {
	case 2*ex <= nearBound(one, p):
		// cosh x - 1 in (0, x²) for |x| <= 1
		return finish(z, proxyNear(one, 1, p), 0)
	case ex > 34:
		return overflow(z, false)
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		b, c, k := expPair(x, wp+8)
		var r ball
		r.add(b, c, wp+8)
		return &r, k - 1, true
	})
}

// Tanh sets z to the correctly rounded hyperbolic tangent of x and returns
// z. Tanh(±0) = ±0, Tanh(±Inf) = ±1.
func Tanh(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch {
	case x.IsInf() && x.Signbit():
		return z.SetInt64(-1)
	case x.IsInf():
		return z.SetInt64(1)
	case x.Sign() == 0:
		return z.Set(x)
	}
	p := z.Prec()
	ex := int64(x.MantExp(nil))
	one := big.NewFloat(float64(x.Sign()))
	// 1 - tanh|x| = 2/(e^(2|x|)+1) < 2^(1 - 2|x|·log2 e) <= 2^nearBound(±1)
	// when |x| >= (1 - nearBound)·0.35.
	lim := big.NewFloat(float64(1-nearBound(one, p)) * 0.35)
	switch {
	case 3*ex <= nearBound(x, p):
		// x - tanh x in (0, x³/3) for x > 0
		return finish(z, proxyNear(x, -x.Sign(), p), 0)
	case cmpAbs(x, lim) >= 0:
		return finish(z, proxyNear(one, -x.Sign(), p), 0)
	}
	x = new(big.Float).Set(x)
	neg := x.Sign() < 0
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		var r ball
		ax := new(big.Float).Abs(x)
		ax.SetMantExp(ax, 1) // 2|x|, exact
		if ex < -1 {
			// 2|x| < 1/2: E = expm1(2|x|), tanh = E/(E+2), E+2 > 2.
			w := wp + uint(-ex) + 24
			E := expm1Core(fxFromFloat(ax, w), expm1Split(wp, ex+1))
			var d, two fx
			two.setInt(2, w)
			d.add(E, &two)
			must(E.quo(E, &d))
			E.toBall(&r)
		} else {
			// e^(2|x|) = b·2^k, tanh = (b - 2^-k)/(b + 2^-k)
			b, k := expReduced(ax, mag{}, wp+8)
			var c, n, d ball
			c.m.SetMantExp(big.NewFloat(1), int(-k))
			n.sub(b, &c, wp+8)
			d.add(b, &c, wp+8)
			must(r.quo(&n, &d, wp+8))
		}
		if neg {
			r.neg(&r)
		}
		return &r, 0, true
	})
}

// log1pBall encloses log(1+a) for a ball a > 0 (every point), with
// relative accuracy near a = 0: 2·atanh(a/(2+a)) when |a| < 1/2.
func log1pBall(a *ball, wp uint) *ball {
	z := new(ball)
	if a.expUpper() < 0 {
		el, ok := a.expLower()
		must(ok)
		var two, d, t ball
		two.setInt64(2)
		d.add(&two, a, wp+16)
		must(t.quo(a, &d, wp+16))
		w := wp + 16 + uint(max(0, -el))
		y := atanhSeries(fxFromBall(&t, w))
		y.shl(y, 1).toBall(z)
		return z
	}
	var one, u ball
	one.setInt64(1)
	u.add(a, &one, wp+16)
	must(logBallB(z, &u, wp))
	return z
}

// Asinh sets z to the correctly rounded inverse hyperbolic sine of x and
// returns z. Asinh(±0) = ±0, Asinh(±Inf) = ±Inf.
func Asinh(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	if x.IsInf() || x.Sign() == 0 {
		return z.Set(x)
	}
	p := z.Prec()
	ex := int64(x.MantExp(nil))
	// x - asinh x in (0, x³/6) for x > 0
	if 3*ex <= nearBound(x, p) {
		return finish(z, proxyNear(x, -x.Sign(), p), 0)
	}
	x = new(big.Float).Set(x)
	neg := x.Sign() < 0
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		ax := new(big.Float).Abs(x)
		r := new(ball)
		switch {
		case ex > int64(wp/2)+8:
			// asinh x - log 2x in (0, 1/(4x²))
			t := new(big.Float).SetMantExp(ax, 1)
			logBall(r, t, wp+8)
			r.addRad(magPow2(-2*ex + 2))
		case ex < 0:
			// log1p(|x| + x²/(1 + sqrt(1+x²)))
			var xb, x2, s, one, a ball
			xb.setExact(ax)
			one.setInt64(1)
			x2.mul(&xb, &xb, wp+24)
			s.add(&x2, &one, wp+24)
			must(s.sqrt(&s, wp+24))
			s.add(&s, &one, wp+24)
			must(x2.quo(&x2, &s, wp+24))
			a.add(&xb, &x2, wp+24)
			r = log1pBall(&a, wp+8)
		default:
			// log(|x| + sqrt(x²+1))
			var xb, s, one ball
			xb.setExact(ax)
			one.setInt64(1)
			s.mul(&xb, &xb, wp+16)
			s.add(&s, &one, wp+16)
			must(s.sqrt(&s, wp+16))
			s.add(&s, &xb, wp+16)
			must(logBallB(r, &s, wp+8))
		}
		if neg {
			r.neg(r)
		}
		return r, 0, true
	})
}

// Acosh sets z to the correctly rounded inverse hyperbolic cosine of x and
// returns z. Acosh(1) = +0, Acosh(+Inf) = +Inf; x < 1 panics with ErrNaN.
func Acosh(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	one := big.NewFloat(1)
	switch c := x.Cmp(one); {
	case c < 0:
		nan("Acosh of a number below 1")
	case c == 0:
		return z.SetInt64(0)
	case x.IsInf():
		return z.SetInf(false)
	}
	ex := int64(x.MantExp(nil))
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		if ex > int64(wp/2)+8 {
			// acosh x - log 2x in (-1/(2x²), 0)
			t := new(big.Float).SetMantExp(x, 1)
			r := new(ball)
			logBall(r, t, wp+8)
			r.addRad(magPow2(-2*ex + 4))
			return r, 0, true
		}
		// log1p(d + sqrt(d(2+d))), d = x-1 exactly
		pr := x.Prec() + uint(max(ex, 0)) + 4
		d := new(big.Float).SetPrec(pr).Sub(x, one)
		e := new(big.Float).SetPrec(pr).Add(d, big.NewFloat(2))
		e.SetPrec(2*pr).Mul(e, d)
		var db, s, a ball
		db.setExact(d)
		s.setExact(e)
		must(s.sqrt(&s, wp+24))
		a.add(&db, &s, wp+24)
		return log1pBall(&a, wp+8), 0, true
	})
}

// Atanh sets z to the correctly rounded inverse hyperbolic tangent of x and
// returns z. Atanh(±0) = ±0, Atanh(±1) = ±Inf; |x| > 1 panics with ErrNaN.
func Atanh(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch c := cmpAbs(x, big.NewFloat(1)); {
	case c > 0:
		nan("Atanh of a number beyond [-1, 1]")
	case c == 0:
		return z.SetInf(x.Signbit())
	case x.Sign() == 0:
		return z.Set(x)
	}
	p := z.Prec()
	ex := int64(x.MantExp(nil))
	// atanh x - x in (0, |x|³/2) for 0 < x <= 1/2
	if 3*ex <= nearBound(x, p) {
		return finish(z, proxyNear(x, x.Sign(), p), 0)
	}
	x = new(big.Float).Set(x)
	neg := x.Sign() < 0
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		ax := new(big.Float).Abs(x)
		r := new(ball)
		if ex < 0 {
			w := wp + uint(-ex) + 16
			atanhSeries(fxFromFloat(ax, w)).toBall(r)
		} else {
			// ½·log1p(2|x|/(1-|x|)), 1-|x| exact
			d := new(big.Float).SetPrec(x.Prec()+2).Sub(big.NewFloat(1), ax)
			var n, db, a ball
			n.setExact(ax)
			n.mul2exp(&n, 1)
			db.setExact(d)
			must(a.quo(&n, &db, wp+16))
			r = log1pBall(&a, wp+8)
			r.mul2exp(r, -1)
		}
		if neg {
			r.neg(r)
		}
		return r, 0, true
	})
}
