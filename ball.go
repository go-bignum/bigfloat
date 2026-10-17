// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math/big"
)

// ball is a midpoint-radius enclosure [m-r, m+r] of a real number: the
// rigorous error model of every computation in this package (Johansson,
// "Arb", IEEE TC 66(8), 2017). Every operation returns a ball guaranteed to
// contain the exact result of the operation applied to any points of its
// operand balls; rounding of the midpoint is added to the radius. The
// midpoint is always finite.
type ball struct {
	m big.Float
	r mag
}

// roundErr bounds the rounding error of a midpoint just computed at
// precision wp to nearest: half an ulp, bounded by 2^(exp(m)-wp).
func roundErr(m *big.Float, wp uint) mag {
	if m.Acc() == big.Exact || m.Sign() == 0 {
		return mag{}
	}
	return magPow2(int64(m.MantExp(nil)) - int64(wp))
}

// prep returns a fresh midpoint at precision wp: operations compute into it
// and then move it into z, so that z may alias an operand (setting the
// precision of z.m first would round that operand).
func prep(wp uint) *big.Float {
	return new(big.Float).SetPrec(wp).SetMode(big.ToNearestEven)
}

// setMid moves m (a fresh value from prep, not used afterwards) into z.
func (z *ball) setMid(m *big.Float, r mag, wp uint) *ball {
	z.m = *m
	z.r = magAdd(r, roundErr(&z.m, wp))
	return z
}

// setFloat sets z to x rounded to wp bits (exact when x fits).
func (z *ball) setFloat(x *big.Float, wp uint) *ball {
	return z.setMid(prep(wp).Set(x), mag{}, wp)
}

// setExact sets z to x exactly, whatever its precision.
func (z *ball) setExact(x *big.Float) *ball {
	z.m.SetPrec(x.Prec()).SetMode(big.ToNearestEven).Set(x)
	z.r = mag{}
	return z
}

// setInt64 sets z to the exact integer v.
func (z *ball) setInt64(v int64) *ball {
	z.m.SetPrec(64).SetMode(big.ToNearestEven).SetInt64(v)
	z.r = mag{}
	return z
}

func (z *ball) set(a *ball) *ball {
	z.m.SetPrec(a.m.Prec()).SetMode(big.ToNearestEven).Set(&a.m)
	z.r = a.r
	return z
}

func (z *ball) add(a, b *ball, wp uint) *ball {
	return z.setMid(prep(wp).Add(&a.m, &b.m), magAdd(a.r, b.r), wp)
}

func (z *ball) sub(a, b *ball, wp uint) *ball {
	return z.setMid(prep(wp).Sub(&a.m, &b.m), magAdd(a.r, b.r), wp)
}

func (z *ball) neg(a *ball) *ball {
	z.set(a)
	z.m.Neg(&z.m)
	return z
}

func (z *ball) abs(a *ball) *ball {
	z.set(a)
	z.m.Abs(&z.m)
	return z
}

// mul: |AB - ab| <= |a|·rb + |b|·ra + ra·rb.
func (z *ball) mul(a, b *ball, wp uint) *ball {
	r := magAdd(magAdd(magMul(magFloat(&a.m), b.r), magMul(magFloat(&b.m), a.r)), magMul(a.r, b.r))
	return z.setMid(prep(wp).Mul(&a.m, &b.m), r, wp)
}

// mul2exp multiplies by 2^k exactly (k small enough not to leave
// big.Float's exponent range).
func (z *ball) mul2exp(a *ball, k int) *ball {
	z.set(a)
	if z.m.Sign() != 0 {
		z.m.SetMantExp(&z.m, k)
	}
	if !z.r.isZero() {
		z.r.e += int64(k)
	}
	return z
}

// lower bounds |x| from below for every x in the ball (0 if it may vanish).
func (a *ball) lower() mag { return magSubLower(magFloatLower(&a.m), a.r) }

// upper bounds |x| from above for every x in the ball.
func (a *ball) upper() mag { return magAdd(magFloat(&a.m), a.r) }

// quo: |A/B - a/b| <= (ra + |a/b|·rb) / (|b| - rb). It reports false when
// the divisor ball contains zero.
func (z *ball) quo(a, b *ball, wp uint) bool {
	lb := b.lower()
	if lb.isZero() {
		return false
	}
	q := magDiv(magFloat(&a.m), magFloatLower(&b.m))
	r := magDiv(magAdd(a.r, magMul(q, b.r)), lb)
	z.setMid(prep(wp).Quo(&a.m, &b.m), r, wp)
	return true
}

// sqrt: |sqrt(X) - sqrt(x)| <= r/(sqrt(x-r)+sqrt(x)) <= r/sqrt(x-r). It
// reports false unless the ball lies in x > 0 (or is the exact point 0).
func (z *ball) sqrt(a *ball, wp uint) bool {
	if a.m.Sign() < 0 {
		return false
	}
	var r mag
	if !a.r.isZero() {
		lb := a.lower()
		if lb.isZero() || a.m.Sign() == 0 {
			return false
		}
		r = magDiv(a.r, magSqrtLower(lb))
	}
	z.setMid(sqrtRound(prep(wp), &a.m), r, wp)
	return true
}

// sign reports the sign of every point of the ball, or 0 when the ball
// contains zero.
func (a *ball) sign() int {
	if a.lower().isZero() {
		return 0
	}
	return a.m.Sign()
}

// exponent returns e with |x| < 2^e for every point (the ball must not be
// zero-free for the bound to be tight; it is always an upper bound).
func (a *ball) expUpper() int64 {
	u := a.upper()
	if u.isZero() {
		return -1 << 40
	}
	return u.e
}

// expLower returns e with |x| >= 2^(e-1) for every point, or ok=false when
// the ball contains zero.
func (a *ball) expLower() (int64, bool) {
	l := a.lower()
	if l.isZero() {
		return 0, false
	}
	return l.e, true
}

// addRad widens the ball by r.
func (z *ball) addRad(r mag) *ball {
	z.r = magAdd(z.r, r)
	return z
}
