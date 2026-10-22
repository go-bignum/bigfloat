// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math"
	"math/big"
)

// mag is an upper (or, where a function says so, lower) bound on a
// non-negative real number, held as m·2^e with m = 0 or m in [0.5, 1).
// It is the radius type of [ball]: cheap float64 arithmetic with every
// operation rounded in the safe direction, over an exponent range wide
// enough for any big.Float (Arb's mag_t plays the same role).
type mag struct {
	m float64
	e int64
}

// up and down move a non-negative float64 at least one ulp towards +Inf or
// towards 0 (at most two). Applied to a result rounded to nearest, they
// give a bound in that direction. For a normal x, x·2^-52 is at least one
// ulp of x, so x(1 ± 2^-52) rounds to nearest at least one ulp away (the
// rounding cannot cross back past a representable neighbour); for zero
// and subnormals, where that product rounds back to x, the smallest
// subnormal 2^-1074 added or removed is the ulp, exactly. Two flops in
// place of math.Nextafter, which the fixed-point bounds call on every
// operation.
func up(x float64) float64 { return x*(1+0x1p-52) + 0x1p-1074 }

func down(x float64) float64 { return max(x*(1-0x1p-52)-0x1p-1074, 0) }

// magNorm builds a normalised mag from m·2^e (m >= 0, finite) without
// changing its value.
func magNorm(m float64, e int64) mag {
	if m == 0 {
		return mag{}
	}
	f, x := math.Frexp(m)
	return mag{f, e + int64(x)}
}

// magPow2 is exactly 2^k.
func magPow2(k int64) mag { return mag{0.5, k + 1} }

func (a mag) isZero() bool { return a.m == 0 }

// magMul bounds a·b from above.
func magMul(a, b mag) mag {
	if a.m == 0 || b.m == 0 {
		return mag{}
	}
	return magNorm(up(a.m*b.m), a.e+b.e)
}

// magAdd bounds a+b from above.
func magAdd(a, b mag) mag {
	if a.m == 0 {
		return b
	}
	if b.m == 0 {
		return a
	}
	if a.e < b.e {
		a, b = b, a
	}
	d := a.e - b.e
	if d > 900 {
		// b < 2^(a.e-900): one ulp of a.m (2^(a.e-53) at least) covers it.
		return magNorm(up(a.m), a.e)
	}
	return magNorm(up(a.m+math.Ldexp(b.m, -int(d))), a.e)
}

// magDiv bounds a/b from above; b must be non-zero.
func magDiv(a, b mag) mag {
	if a.m == 0 {
		return mag{}
	}
	return magNorm(up(a.m/b.m), a.e-b.e)
}

// magSubLower bounds a-b from below, for a lower bound a and an upper bound
// b. The result is zero when the difference cannot be shown positive.
func magSubLower(a, b mag) mag {
	if b.m == 0 {
		return a
	}
	if a.m == 0 || a.e < b.e {
		return mag{}
	}
	d := a.e - b.e
	var bb float64
	if d > 900 {
		bb = math.Ldexp(1, -900) // b < 2^(a.e-d) <= this, rounded up
	} else {
		bb = math.Ldexp(b.m, -int(d))
	}
	s := a.m - bb
	if s <= 0 {
		return mag{}
	}
	return magNorm(down(s), a.e)
}

// magSqrtLower bounds sqrt(a) from below (a non-zero).
func magSqrtLower(a mag) mag {
	m, e := a.m, a.e
	if e%2 != 0 {
		m *= 2
		e--
	}
	return magNorm(down(math.Sqrt(m)), e/2)
}

// magFloat bounds |x| from above (finite x).
func magFloat(x *big.Float) mag {
	if x.Sign() == 0 {
		return mag{}
	}
	e := x.MantExp(nil)
	if e > -1000 && e < 1000 {
		f, _ := x.Float64()
		return magNorm(up(math.Abs(f)), 0)
	}
	// |x| < 2^e.
	return magPow2(int64(e))
}

// magFloatLower bounds |x| from below (finite x).
func magFloatLower(x *big.Float) mag {
	if x.Sign() == 0 {
		return mag{}
	}
	e := x.MantExp(nil)
	if e > -1000 && e < 1000 {
		f, acc := x.Float64()
		f = math.Abs(f)
		if acc != big.Exact {
			f = down(f)
		}
		return magNorm(f, 0)
	}
	// |x| >= 2^(e-1).
	return magPow2(int64(e) - 1)
}

// magUint bounds the integer u from above.
func magUint(u uint64) mag { return magNorm(up(float64(u)), 0) }

// float sets z to the exact value of a (z's precision is set to 53).
// Exponents beyond big.Float's range give ±Inf or 0, which only ever makes a
// radius useless (the Ziv loop then retries), never wrong in the unsafe
// direction for the radii this package builds.
func (a mag) float(z *big.Float) *big.Float {
	z.SetPrec(53)
	z.SetFloat64(a.m)
	if a.m == 0 {
		return z
	}
	e := a.e
	if e > big.MaxExp {
		return z.SetInf(false)
	}
	if e < big.MinExp {
		e = big.MinExp
	}
	return z.SetMantExp(z, int(e))
}
