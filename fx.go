// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math"
	"math/big"
)

// fx is a fixed-point ball: the real interval [(v-e)·2^-w, (v+e)·2^-w].
// The series cores (exp, log, sin/cos, atan) run in fixed point, as Brent &
// Zimmermann ("Modern Computer Arithmetic", §4.4) and Arb do, because a
// big.Int multiply and shift costs far less than a normalising big.Float
// operation. The error count e is in units of 2^-w, a float64 kept as an
// upper bound by rounding every operation on it upwards. Values handled in
// fixed point are small (|x| < 2^60), which keeps the magnitude bounds
// representable.
//
// The operations do not allocate once the operands have reached their size:
// big.Int reallocates a destination that aliases an operand of Mul, and
// allocates a remainder for Quo, so each fx carries scratch integers that
// take those roles (t, the full product or the remainder; k, a small
// divisor) and keep their capacity from one series term to the next.
type fx struct {
	v    big.Int
	e    float64
	w    uint
	t, k big.Int
	u    big.Int   // isqrt scratch
	f    big.Float // conversion scratch
}

func addUp(a, b float64) float64 { return up(a + b) }
func mulUp(a, b float64) float64 {
	if a == 0 || b == 0 {
		return 0
	}
	return up(a * b)
}

// bnd returns an upper bound on |x| (real units) for every point x of a.
func (a *fx) bnd() float64 {
	n := a.v.BitLen()
	// Both terms are kept at or above the smallest subnormal, where Ldexp
	// would otherwise round them to zero or downwards at scales w > 1074.
	b := math.Ldexp(1, max(n-int(a.w), -1074))
	if n == 0 {
		b = 0
	}
	return addUp(b, up(math.Ldexp(a.e, -int(a.w))))
}

// log2bnd returns an integer b with |x| <= 2^b for every point of a.
func (a *fx) log2bnd() int {
	a.t.Abs(&a.v)
	a.k.SetUint64(uint64(min(math.Ceil(a.e), 1<<62))) // e is finite in practice
	a.t.Add(&a.t, &a.k)
	return a.t.BitLen() - int(a.w)
}

// expLower returns e with |x| >= 2^(e-1) for every point of a, or false
// when a may contain zero (as ball.expLower; exact integer arithmetic, so
// that no scale is too large for it). It writes a's scratch.
func (a *fx) expLower() (int64, bool) {
	a.t.Abs(&a.v)
	a.k.SetUint64(uint64(min(math.Ceil(a.e), 1<<62)))
	a.t.Sub(&a.t, &a.k)
	if a.t.Sign() <= 0 {
		return 0, false
	}
	return int64(a.t.BitLen()) - int64(a.w), true
}

// lowbnd returns a lower bound on |x| for every point of a (0 if it may
// vanish).
func (a *fx) lowbnd() float64 {
	n := a.v.BitLen()
	if n == 0 {
		return 0
	}
	l := math.Ldexp(1, n-1-int(a.w)) - math.Ldexp(a.e, -int(a.w))
	if l <= 0 {
		return 0
	}
	return down(l)
}

func (z *fx) setInt(c int64, w uint) *fx {
	z.w = w
	z.v.SetInt64(c)
	z.v.Lsh(&z.v, w)
	z.e = 0
	return z
}

func (z *fx) set(a *fx) *fx {
	z.w = a.w
	z.v.Set(&a.v)
	z.e = a.e
	return z
}

// mul: |ab - AB| <= |a|·eb + (|b|+eb)·ea, plus one ulp for the floor.
func (z *fx) mul(a, b *fx) *fx {
	ma, mb := a.bnd(), b.bnd()
	e := addUp(addUp(mulUp(ma, b.e), mulUp(mb, a.e)), 1)
	z.w = a.w
	z.t.Mul(&a.v, &b.v) // z.t is never an operand: no aliasing copy
	z.v.Rsh(&z.t, a.w)
	z.e = e
	return z
}

func (z *fx) add(a, b *fx) *fx {
	e := addUp(a.e, b.e)
	z.w = a.w
	z.v.Add(&a.v, &b.v)
	z.e = e
	return z
}

func (z *fx) sub(a, b *fx) *fx {
	e := addUp(a.e, b.e)
	z.w = a.w
	z.v.Sub(&a.v, &b.v)
	z.e = e
	return z
}

// divSmall divides by the positive integer k (truncation: one ulp).
func (z *fx) divSmall(a *fx, k uint64) *fx {
	e := addUp(up(a.e/float64(k)), 1)
	z.w = a.w
	z.k.SetUint64(k)
	z.v.QuoRem(&a.v, &z.k, &z.t) // QuoRem reuses z.t; Quo would allocate a remainder
	z.e = e
	return z
}

// shr divides by 2^s (floor: one ulp).
func (z *fx) shr(a *fx, s uint) *fx {
	e := addUp(up(math.Ldexp(a.e, -int(s))), 1)
	z.w = a.w
	z.v.Rsh(&a.v, s)
	z.e = e
	return z
}

// shl multiplies by 2^s exactly.
func (z *fx) shl(a *fx, s uint) *fx {
	e := math.Ldexp(a.e, int(s))
	z.w = a.w
	z.v.Lsh(&a.v, s)
	z.e = e
	return z
}

// quo divides a by b: |A/B - a/b| <= (ea + |a/b|·eb)/|B|, with |B| bounded
// below; one more ulp for the truncation. It reports false when b may
// vanish.
func (z *fx) quo(a, b *fx) bool {
	lb := b.lowbnd()
	if lb == 0 {
		return false
	}
	q := up(a.bnd() / lb)
	e := addUp(up(addUp(a.e, mulUp(q, b.e))/lb), 1)
	z.t.Lsh(&a.v, a.w)
	z.w = a.w
	z.v.QuoRem(&z.t, &b.v, &z.k)
	z.e = e
	return true
}

// sqrt of a ball whose points are all >= 1/16 (lowbnd, a power of two
// below the value less its error, must show it): |sqrt(A)-sqrt(a)| <=
// |A-a|/(2·sqrt(min)), at most twice the error, plus one ulp for the
// integer square root.
func (z *fx) sqrt(a *fx) bool {
	lb := a.lowbnd()
	if lb < 1.0/16 || a.v.Sign() < 0 {
		return false
	}
	e := addUp(up(a.e/(2*down(math.Sqrt(lb)))), 1)
	z.t.Lsh(&a.v, a.w)
	z.w = a.w
	isqrt(&z.v, &z.t, &z.k, &z.u)
	z.e = e
	return true
}

// isqrt sets z to floor(sqrt(n)) for n > 0, with s and r as scratch (all
// four distinct). big.Int.Sqrt allocates on every call; this is the same
// Newton iteration from above (z <- (z + n/z)/2, stopping when it no
// longer decreases: Brent & Zimmermann, MCA §1.5.1, SqrtInt), started from
// a float64 estimate of the top 62 bits raised past the root, so that the
// iteration is quadratic from the first step and reuses its integers.
func isqrt(z, n, s, r *big.Int) {
	L := n.BitLen()
	shift := uint(max(L-62, 0)) &^ 1 // even, top keeps 61 or 62 bits
	s.Rsh(n, shift)
	// sqrt(n) < sqrt(top+1)·2^(shift/2) <= (floor(sqrt(top)) + 2)·2^(shift/2):
	// float64(top) is within 2^-53 relative, sqrt(top+1) - sqrt(top) < 1.
	z.SetUint64(uint64(math.Sqrt(float64(s.Uint64()))) + 2)
	z.Lsh(z, shift/2)
	for {
		s.QuoRem(n, z, r)
		s.Add(s, z)
		s.Rsh(s, 1)
		if s.Cmp(z) >= 0 {
			return
		}
		z.Set(s)
	}
}

// fxFromBall converts a ball to fixed point at scale w (truncation adds one
// ulp, the radius is carried over rounded up).
func fxFromBall(b *ball, w uint) *fx { return new(fx).setBall(b, w) }

// setBall sets z to the ball b at scale w, as fxFromBall, reusing z's
// storage.
func (z *fx) setBall(b *ball, w uint) *fx {
	z.w = w
	z.f.SetPrec(0) // Set below takes m's precision instead of rounding
	z.f.SetMantExp(&b.m, int(w))
	z.f.Int(&z.v)
	e := 1.0
	if !b.r.isZero() {
		x := b.r.e + int64(w)
		if x > 1000 {
			e = math.Inf(1)
		} else if x > -1000 {
			e = addUp(e, up(math.Ldexp(b.r.m, int(x))))
		} else {
			e = addUp(e, math.Ldexp(1, -1000))
		}
	}
	z.e = e
	return z
}

// fxFromFloat converts an exact big.Float.
func fxFromFloat(x *big.Float, w uint) *fx { return new(fx).setFloat(x, w) }

// setFloat sets z to the exact finite x at scale w: truncated towards
// zero, one ulp of error unless x fits.
func (z *fx) setFloat(x *big.Float, w uint) *fx {
	z.w = w
	z.f.SetPrec(0) // Set below takes x's precision instead of rounding
	z.f.SetMantExp(x, int(w))
	z.e = 0
	if _, acc := z.f.Int(&z.v); acc != big.Exact {
		z.e = 1
	}
	return z
}

// mulInt sets z to a·k for an integer k, exactly: the error scales by |k|
// (float64(k) is within half an ulp, up covers it).
func (z *fx) mulInt(a *fx, k int64) *fx {
	e := mulUp(a.e, up(math.Abs(float64(k))))
	z.w = a.w
	z.k.SetInt64(k)
	z.v.Mul(&a.v, &z.k)
	z.e = e
	return z
}

// rescale sets z to a at the smaller scale w <= a.w (floor: one ulp).
func (z *fx) rescale(a *fx, w uint) *fx {
	d := a.w - w
	e := addUp(up(math.Ldexp(a.e, -int(d))), 1)
	z.v.Rsh(&a.v, d)
	z.w = w
	z.e = e
	return z
}

// toBall converts back to a ball (exactly: the midpoint is the integer v
// scaled by 2^-w).
func (a *fx) toBall(z *ball) *ball {
	z.m.SetPrec(0)
	z.m.SetMode(big.ToNearestEven)
	z.m.SetInt(&a.v)
	if z.m.Sign() != 0 {
		z.m.SetMantExp(&z.m, -int(a.w))
	}
	z.r = magMul(magNorm(up(a.e), 0), magPow2(-int64(a.w)))
	if math.IsInf(a.e, 1) {
		z.r = magPow2(1 << 40)
	}
	return z
}
