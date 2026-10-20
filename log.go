// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math/big"
)

// atanhSeries returns atanh(t) = Σ t^(2k+1)/(2k+1) for a fixed-point ball
// with |t| <= 1/2. The tail after the k = n term is at most
// |t|^(2n+3)/(1-t²) <= 2|t|^(2n+3).
func atanhSeries(t *fx) *fx {
	w := t.w
	// every caller has |t| <= 1/3; log2bnd counts the error bound, so
	// lt >= 1 - w and n stays positive
	lt := min(int64(t.log2bnd()), -1)
	// 1 + lt(2n+3) <= -w
	n := (int64(w)+1+(-lt)-1)/(-lt)/2 + 1
	y := oddSeries(t, n, false)
	y.e = addUp(y.e, 1) // tail
	return y
}

// logReduced encloses log m for an exact m in [1/sqrt2, sqrt2) as
// 2·atanh((m-1)/(m+1)), |(m-1)/(m+1)| <= 0.1716, at scale w: the quotient
// is formed from the exact integer mantissa, truncated once.
//
// Up to tableMaxPrec it uses mpmath's table reduction instead
// (log_taylor_cached, libelefun.py): with a = n/2^9 the grid point below m,
// log m = log a + 2·atanh((m-a)/(m+a)), |(m-a)/(m+a)| < 2^-9, so that the
// series gains 18 bits per term; log a comes from the table. Within 2^-8
// of 1 the direct quotient is already that small (and keeps the relative
// accuracy of log m near 0, which a sum log a + … would cancel), so it is
// used there.
func logReduced(m *big.Float, w uint) *fx {
	M, E := mantInt(m) // m = M·2^E, E <= 0
	if w <= tableMaxPrec {
		c := min(E, -logShift)
		Mi := new(big.Int).Lsh(M, uint(E-c)) // m = Mi·2^c
		n := new(big.Int).Rsh(Mi, uint(-logShift-c)).Int64()
		if n < 1<<logShift-2 || n > 1<<logShift+1 {
			// a = n·2^-9 = Ai·2^c
			Ai := new(big.Int).Lsh(big.NewInt(n), uint(-logShift-c))
			num := new(big.Int).Sub(Mi, Ai)
			den := Mi.Add(Mi, Ai)
			num.Lsh(num, w)
			t := &fx{w: w, e: 1}
			t.v.Quo(num, den)
			y := atanhSeries(t)
			y.shl(y, 1)
			return y.add(y, logTable.get(n, w))
		}
	}
	one := new(big.Int).Lsh(big.NewInt(1), uint(-E))
	num := new(big.Int).Sub(M, one)
	den := new(big.Int).Add(M, one)
	num.Lsh(num, w)
	t := &fx{w: w, e: 1}
	t.v.Quo(num, den)
	y := atanhSeries(t)
	return y.shl(y, 1)
}

// logShift sets the grid of the log table: log(n/2^logShift), as mpmath's
// LOG_TAYLOR_SHIFT.
const logShift = 9

// logTable caches log(n/2^9) = 2·atanh((n-2^9)/(n+2^9)) for the grid points
// n of [1/sqrt2, sqrt2), |(n-2^9)/(n+2^9)| <= 0.172.
var logTable = &fxTable{compute: func(n int64, W uint) *fx {
	num := big.NewInt(n - 1<<logShift)
	num.Lsh(num, W)
	t := &fx{w: W, e: 1}
	t.v.Quo(num, big.NewInt(n+1<<logShift))
	y := atanhSeries(t)
	return y.shl(y, 1)
}}

var sqrtHalf = big.NewFloat(0.70710678118654752440)

// logBall sets z to an enclosure of log x, x > 0 exact finite, with about
// wp bits of relative accuracy (the reduction keeps it near x = 1). x = 1
// gives the exact 0 (t = 0 in the series).
func logBall(z *ball, x *big.Float, wp uint) {
	m := new(big.Float)
	e := x.MantExp(m)
	if m.Cmp(sqrtHalf) < 0 {
		m.SetMantExp(m, 1)
		e--
	}
	w := wp + 12
	if e == 0 {
		// near 1: log m ~ m-1 needs as many more bits as m-1 has leading
		// zeros (m-1 is exact at prec(m)+2 bits)
		var d big.Float
		d.SetPrec(m.Prec()+2).Sub(m, big.NewFloat(1))
		w += uint(max(0, -d.MantExp(nil)))
	}
	var rad mag
	// a long m is rounded first, to w+64 bits; since m >= 1/sqrt2,
	// |log m' - log m| <= 2|m'-m|. (m - 1 keeps its leading bits.)
	if m.Prec() > w+96 {
		var mb ball
		mb.setFloat(m, w+64)
		rad = magMul(mb.r, magPow2(1))
		m = &mb.m
	}
	lm := logReduced(m, w)
	var lb ball
	lm.toBall(&lb)
	if e != 0 {
		ee := int64(e)
		bl := uint(bitlenU(abs64(ee)))
		l := ln2Ball(w + bl + 4)
		var eb, t ball
		eb.setInt64(ee)
		t.mul(l, &eb, w+bl+4)
		z.add(&t, &lb, w)
	} else {
		z.set(&lb)
	}
	z.addRad(rad)
}

// logBallB encloses log of a ball: |log X - log x| <= r/(x - r).
func logBallB(z *ball, b *ball, wp uint) bool {
	lb := b.lower()
	if b.m.Sign() <= 0 || lb.isZero() {
		return false
	}
	logBall(z, &b.m, wp)
	z.addRad(magDiv(b.r, lb))
	return true
}

// Log sets z to the correctly rounded natural logarithm of x and returns z.
// Log(+Inf) = +Inf, Log(±0) = -Inf, Log(1) = +0; x < 0 panics with
// big.ErrNaN.
func Log(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch {
	case x.Sign() < 0:
		nan("Log of a negative number")
	case x.Sign() == 0:
		return z.SetInf(true)
	case x.IsInf():
		return z.SetInf(false)
	case x.Cmp(big.NewFloat(1)) == 0:
		return z.SetInt64(0)
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		var b ball
		logBall(&b, x, wp)
		return &b, 0, true
	})
}

// Log2 sets z to the correctly rounded base-2 logarithm of x and returns z.
// It is exact for powers of two.
func Log2(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch {
	case x.Sign() < 0:
		nan("Log2 of a negative number")
	case x.Sign() == 0:
		return z.SetInf(true)
	case x.IsInf():
		return z.SetInf(false)
	}
	m := new(big.Float)
	e := x.MantExp(m)
	if m.Cmp(big.NewFloat(0.5)) == 0 {
		return z.SetInt64(int64(e) - 1)
	}
	return logBase(z, x, ln2Ball)
}

// Log10 sets z to the correctly rounded base-10 logarithm of x and returns
// z. It is exact for x = 10^k, k >= 0 (the only powers of ten a binary
// number can hold exactly).
func Log10(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch {
	case x.Sign() < 0:
		nan("Log10 of a negative number")
	case x.Sign() == 0:
		return z.SetInf(true)
	case x.IsInf():
		return z.SetInf(false)
	}
	if k, ok := powerOfTen(x); ok {
		return z.SetInt64(k)
	}
	return logBase(z, x, ln10Ball)
}

// powerOfTen reports whether x = 10^k for an integer k >= 0. x = 10^k is
// an integer whose odd part 5^k has about 2.32k bits, so only k up to
// prec(x)/2.32 need be tried.
func powerOfTen(x *big.Float) (int64, bool) {
	if !x.IsInt() {
		return 0, false
	}
	M, E := mantInt(x)
	tz := M.TrailingZeroBits()
	M.Rsh(M, tz)
	k := E + int64(tz) // x = odd·2^k, need odd = 5^k
	if k < 0 || k > int64(M.BitLen()) {
		return 0, false
	}
	f := new(big.Int).Exp(big.NewInt(5), big.NewInt(k), nil)
	return k, f.Cmp(M) == 0
}

func logBase(z, x *big.Float, base func(uint) *ball) *big.Float {
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		var b ball
		logBall(&b, x, wp+4)
		c := base(wp + 8)
		b.quo(&b, c, wp+4)
		return &b, 0, true
	})
}

// Log1p sets z to the correctly rounded log(1+x) and returns z.
// Log1p(±0) = ±0, Log1p(-1) = -Inf, Log1p(+Inf) = +Inf; x < -1 panics
// with big.ErrNaN.
func Log1p(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	m1 := big.NewFloat(-1)
	switch c := x.Cmp(m1); {
	case c < 0:
		nan("Log1p of a number below -1")
	case c == 0:
		return z.SetInf(true)
	case x.Sign() == 0:
		return z.Set(x)
	case x.IsInf():
		return z.SetInf(false)
	}
	p := z.Prec()
	ex := int64(x.MantExp(nil))
	// x - log1p(x) is in (0, x²) for 0 < |x| <= 1/2.
	if 2*ex <= nearBound(x, p) {
		return finish(z, proxyNear(x, -1, p), 0)
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		var b ball
		switch {
		case ex < 0:
			// |x| < 1/2: 2·atanh(x/(2+x)), formed exactly in integers.
			w := wp + uint(-ex) + 12
			M, E := mantInt(x) // |x| = M·2^E
			if x.Sign() < 0 {
				M.Neg(M)
			}
			two := new(big.Int).Lsh(big.NewInt(1), uint(1-E))
			den := new(big.Int).Add(two, M)
			M.Lsh(M, w)
			t := &fx{w: w, e: 1}
			t.v.Quo(M, den)
			y := atanhSeries(t)
			y.shl(y, 1).toBall(&b)
		case x.Sign() < 0 || ex <= int64(wp):
			// 1+x exactly (x >= -1/2 has at most a few more bits).
			u := new(big.Float).SetPrec(x.Prec() + uint(max(ex, 0)) + 2)
			u.Add(x, big.NewFloat(1))
			logBall(&b, u, wp)
		default:
			var u, xb, one ball
			xb.setFloat(x, wp+16)
			one.setInt64(1)
			u.add(&xb, &one, wp+16)
			logBallB(&b, &u, wp)
		}
		return &b, 0, true
	})
}
