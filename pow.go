// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math/big"
)

// powExactLimit bounds the size, in bits, of an exact integer power that
// Pow computes outright to settle an exact case.
const powExactLimit = 1 << 20

// isOddInt reports whether y is an odd integer.
func isOddInt(y *big.Float) bool {
	if !y.IsInt() || y.Sign() == 0 {
		return false
	}
	// y = M·2^E with M of prec bits; y is even once E > 0.
	M, E := mantInt(y)
	tz := int64(M.TrailingZeroBits())
	return E+tz == 0
}

// Pow sets z to the correctly rounded x**y and returns z. If z's precision
// is 0 it is set to the larger of x's and y's. Special values follow C99
// Annex F (F.9.4.4) and IEEE 754 pow: Pow(x, ±0) = 1 for any x,
// Pow(1, y) = 1 for any y, Pow(-1, ±Inf) = 1, Pow(±0, y) and Pow(±Inf, y)
// by the sign and parity of y, Pow(x, ±Inf) by |x| against 1. x < 0 with a
// finite non-integer y panics with ErrNaN. Exact results (integer powers
// that fit, perfect powers under dyadic exponents) are recognised and
// rounded exactly.
func Pow(z, x, y *big.Float) *big.Float {
	precOf(z, max(x.Prec(), y.Prec()))
	one := big.NewFloat(1)
	switch {
	case y.Sign() == 0, x.Cmp(one) == 0:
		return z.SetInt64(1)
	case y.IsInf():
		c := cmpAbs(x, one)
		switch {
		case c == 0: // x = -1
			return z.SetInt64(1)
		case (c < 0) == y.Signbit():
			return z.SetInf(false)
		default:
			return z.SetInt64(0)
		}
	case x.Sign() == 0 || x.IsInf():
		// |x|^y is 0 or Inf; the sign is x's when y is an odd integer.
		neg := x.Signbit() && isOddInt(y)
		big0 := (x.Sign() == 0) == (y.Sign() > 0)
		if big0 {
			return zeroSign(z.SetInt64(0), neg)
		}
		return z.SetInf(neg)
	case x.Sign() < 0 && !y.IsInt():
		nan("Pow of a negative number to a non-integer power")
	}
	neg := x.Sign() < 0 && isOddInt(y)
	ax := new(big.Float).Abs(x)
	y = new(big.Float).Set(y)
	p := z.Prec()
	if v, K, ok := powExact(ax, y, p); ok {
		if neg {
			v.Neg(v)
		}
		return finish(z, v, K)
	}
	// t = y·log|x|: its magnitude decides overflow, underflow and the
	// near-1 results first, from a 64-bit enclosure.
	t := powLog(ax, y, 64)
	if t.lower().e > 34 {
		if t.m.Sign() > 0 {
			return overflow(z, neg)
		}
		return underflow(z, neg)
	}
	// |x^y - 1| <= 2|t| for |t| <= ln 2.
	if t.expUpper()+1 <= nearBound(one, p) {
		v := proxyNear(one, t.sign(), p)
		if neg {
			v.Neg(v)
		}
		return finish(z, v, 0)
	}
	return ziv(z, x.Prec()+y.Prec(), func(wp uint) (*ball, int64, bool) {
		eu := max(t.expUpper(), 0)
		t := powLog(ax, y, wp+uint(eu)+16)
		b, K := expReduced(&t.m, t.r, wp)
		if neg {
			b.neg(b)
		}
		return b, K, true
	})
}

// powLog encloses y·log x (x > 0, x != 1) to about wp bits of relative
// accuracy: logBall keeps its accuracy relative near x = 1, so the
// enclosure of log x never contains zero and t has the sign of
// y·(x - 1).
func powLog(x, y *big.Float, wp uint) *ball {
	var l, yb ball
	logBall(&l, x, wp)
	yb.setExact(y)
	t := new(ball)
	return t.mul(&l, &yb, wp)
}

// powExact settles the exact cases of |x|^y (x > 0 finite, y finite
// non-zero, x != 1): it returns v and K with |x|^y = v·2^K exactly when
// |x|^y is a dyadic number of at most powExactLimit bits, a power of two
// of any size, or else reports false — the result is then either
// irrational or a dyadic number too long to be a precision-p breakpoint,
// and the Ziv loop terminates on it.
func powExact(x, y *big.Float, p uint) (*big.Float, int64, bool) {
	// y = k/2^j, k odd: x^y is dyadic iff x is a perfect 2^j-th power
	// (if x^(k/2^j) is rational, so is x^(1/2^j): a·k + b·2^j = 1).
	b := new(big.Float).Set(x)
	k := new(big.Float).Set(y)
	for !k.IsInt() {
		M, E := mantInt(b)
		tz := M.TrailingZeroBits()
		M.Rsh(M, tz)
		E += int64(tz)
		if E%2 != 0 {
			return nil, 0, false
		}
		r, exact := iroot(M, 2)
		if !exact {
			return nil, 0, false
		}
		b.SetPrec(0).SetInt(r)
		b.SetMantExp(b, int(E/2))
		k.SetMantExp(k, 1)
	}
	// b^k for an integer k: b = m·2^e, m odd.
	M, E := mantInt(b)
	tz := M.TrailingZeroBits()
	M.Rsh(M, tz)
	E += int64(tz)
	if k.MantExp(nil) > 62 {
		// |k| >= 2^62: only m = 1 can be exact, and then 2^(E·k)
		// overflows or underflows unless E = 0 (b = 1, excluded).
		if M.BitLen() != 1 || E == 0 {
			return nil, 0, false
		}
		v := big.NewFloat(1)
		if (E > 0) == (k.Sign() > 0) {
			return v, big.MaxExp + 2, true
		}
		return v, big.MinExp - 2, true
	}
	kk, _ := k.Int64()
	exp := mulClamp(E, kk)
	if M.BitLen() == 1 {
		return big.NewFloat(1), exp, true
	}
	if kk < 0 {
		return nil, 0, false // 1/m^|k| with m > 1 odd is not dyadic
	}
	if uint64(kk) > powExactLimit/uint64(M.BitLen()-1) {
		return nil, 0, false // more than powExactLimit bits: not a breakpoint
	}
	// Exact and small enough to compute outright: round it directly.
	r := new(big.Int).Exp(M, big.NewInt(kk), nil)
	v := new(big.Float).SetPrec(0).SetInt(r)
	return v, exp, true
}

// mulClamp returns a·b clamped to ±2^40 (beyond which any power-of-two
// scaling overflows or underflows anyway).
func mulClamp(a, b int64) int64 {
	const lim = 1 << 40
	if a == 0 || b == 0 {
		return 0
	}
	if abs64(a) > lim/abs64(b) {
		if (a > 0) == (b > 0) {
			return lim
		}
		return -lim
	}
	return a * b
}
