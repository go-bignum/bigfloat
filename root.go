// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math"
	"math/big"
)

// rootDirectLimit bounds the size, in bits, of the integer whose n-th root
// Root extracts directly; beyond it Root goes through exp(log(x)/n).
const rootDirectLimit = 1 << 20

// mantInt returns the integer M and exponent E with x = M·2^E (x finite,
// non-zero, sign ignored).
func mantInt(x *big.Float) (*big.Int, int64) {
	e := x.MantExp(nil)
	q := x.Prec()
	var t big.Float
	t.SetPrec(0)
	t.SetMantExp(x, int(q)-e)
	t.Abs(&t)
	M, _ := t.Int(nil)
	return M, int64(e) - int64(q)
}

// iroot returns floor(N^(1/n)) for N > 0, n >= 2, and whether it is exact.
func iroot(N *big.Int, n uint) (*big.Int, bool) {
	if n == 2 {
		r := new(big.Int).Sqrt(N)
		var s big.Int
		return r, s.Mul(r, r).Cmp(N) == 0
	}
	// Newton from above: x <- ((n-1)x + N/x^(n-1)) / n decreases
	// monotonically to floor(N^(1/n)) from any starting x >= the root, and
	// never goes below it. From a power of two the first steps only shrink
	// x by (n-1)/n each, hundreds of steps for large n; rootStart starts
	// within half the root's bits instead (precision doubling, as mpmath's
	// nthroot_fixed, libelefun.py), so that one step reaches the root. Once
	// a step has decreased x, x^n <= N proves x = floor(root): one power
	// instead of the further Newton step that would confirm it.
	x := rootStart(N, n)
	nn := big.NewInt(int64(n))
	n1 := big.NewInt(int64(n - 1))
	var t, y, r big.Int
	for {
		t.Exp(x, n1, nil)
		y.QuoRem(N, &t, &r)
		t.Mul(x, n1)
		y.Add(&y, &t)
		y.QuoRem(&y, nn, &r)
		if y.Cmp(x) >= 0 {
			break
		}
		x.Set(&y)
		if t.Exp(x, nn, nil).Cmp(N) <= 0 {
			return x, t.Cmp(N) == 0
		}
	}
	t.Exp(x, nn, nil)
	return x, t.Cmp(N) == 0
}

// rootStart returns an integer at or above N^(1/n) (n >= 3); a variable
// for tests. (Set in init: rootAbove calls iroot, which reads it.)
var rootStart func(N *big.Int, n uint) *big.Int

func init() { rootStart = rootAbove }

// rootAbove returns an integer above N^(1/n) within a relative 2^-(b/2)
// for a root of b bits: rootGuess for roots of up to 64 bits, else
// (floor((N >> ns)^(1/n)) + 1)·2^s with s half the root's bits, the exact
// root of the top part from iroot (recursively: the precision doubles from
// one level to the next). N < ((N >> ns) + 1)·2^(ns) gives
// root(N) < root((N >> ns) + 1)·2^s <= (floor(root(N >> ns)) + 1)·2^s.
func rootAbove(N *big.Int, n uint) *big.Int {
	b := (N.BitLen() + int(n) - 1) / int(n)
	if b <= 64 {
		return rootGuess(N, n)
	}
	s := uint(b / 2)
	top := new(big.Int).Rsh(N, n*s)
	r, _ := iroot(top, n)
	r.Add(r, big.NewInt(1))
	return r.Lsh(r, s)
}

// rootGuess returns an integer above N^(1/n) (N > 0, n >= 3), within a
// relative 2^-19 of it: log2 N from the top 64 bits, divided by n, read
// back as 53 bits and a shift, then raised by 2^-20 and 1. The float64
// error of log2 N / n is below 2^-40 relative to the root even for N of
// 2^40 bits, well inside the margin (Newton from below would stop at once).
func rootGuess(N *big.Int, n uint) *big.Int {
	L := N.BitLen()
	s := max(L-64, 0)
	top := new(big.Int).Rsh(N, uint(s))
	m, _ := new(big.Float).SetInt(top).Float64()
	l2 := (math.Log2(m) + float64(s)) / float64(n) // log2 of the root
	e := max(int(math.Floor(l2))-52, 0)
	g, _ := big.NewFloat(math.Exp2(l2 - float64(e))).Int(nil)
	g.Add(g, new(big.Int).Rsh(g, 20))
	g.Add(g, big.NewInt(2))
	return g.Lsh(g, uint(e))
}

// rootProxy returns a proxy (see proxyNear) for |x|^(1/n) at precision p,
// computed exactly with integers: the root r of a scaled mantissa with at
// least p+2 bits, with a sticky bit appended when it is not exact. The true
// root lies in (r, r+1) and no precision-p breakpoint lies strictly inside
// an interval of integers that wide, so r + 1/2 rounds like the root.
func rootProxy(x *big.Float, n uint, p uint) *big.Float {
	M, E := mantInt(x)
	t := int64(n)*int64(p+2) - int64(M.BitLen())
	if t < 0 {
		t = 0
	}
	// make E - t divisible by n
	if r := (E - t) % int64(n); r != 0 {
		if r < 0 {
			r += int64(n)
		}
		t += r
	}
	N := new(big.Int).Lsh(M, uint(t))
	r, exact := iroot(N, n)
	exp := (E - t) / int64(n)
	if !exact {
		r.Lsh(r, 1)
		r.SetBit(r, 0, 1)
		exp--
	}
	v := new(big.Float).SetPrec(0).SetInt(r)
	return v.SetMantExp(v, int(exp))
}

// sqrtRound sets z to the correctly rounded square root of x > 0, at z's
// precision and mode.
func sqrtRound(z, x *big.Float) *big.Float {
	return z.Set(rootProxy(x, 2, z.Prec()))
}

// Sqrt sets z to the correctly rounded square root of x and returns z.
// math/big's own Float.Sqrt is not correctly rounded in every mode (it
// computes an inverse square root by Newton's method and rounds a product);
// this one extracts an integer square root and is exact in every mode.
// Sqrt(±0) = ±0, Sqrt(+Inf) = +Inf; x < 0 panics with big.ErrNaN.
func Sqrt(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch {
	case x.Sign() < 0:
		nan("Sqrt of a negative number")
	case x.Sign() == 0 || x.IsInf():
		return z.Set(x)
	}
	return finish(z, rootProxy(x, 2, z.Prec()), 0)
}

// Cbrt sets z to the correctly rounded cube root of x and returns z.
// Cbrt(±0) = ±0, Cbrt(±Inf) = ±Inf.
func Cbrt(z, x *big.Float) *big.Float { return Root(z, x, 3) }

// Root sets z to the correctly rounded n-th root of x and returns z.
// Following C23 rootn: Root(±0, n) is ±0 for odd n and +0 for even n,
// Root(+Inf, n) = +Inf, Root(-Inf, n) = -Inf for odd n. n = 0, and x < 0
// with n even, panic with big.ErrNaN.
func Root(z, x *big.Float, n uint) *big.Float {
	precOf(z, x.Prec())
	switch {
	case n == 0:
		nan("Root of order 0")
	case x.Sign() < 0 && n%2 == 0:
		nan("even Root of a negative number")
	case n == 1:
		return z.Set(x)
	case x.Sign() == 0:
		if n%2 == 0 {
			return z.SetInt64(0)
		}
		return z.Set(x)
	case x.IsInf():
		return z.Set(x)
	}
	p := z.Prec()
	neg := x.Sign() < 0
	M, _ := mantInt(x)
	if uint64(n)*uint64(p+2)+uint64(M.BitLen()) <= rootDirectLimit {
		v := rootProxy(x, n, p)
		if neg {
			v.Neg(v)
		}
		return finish(z, v, 0)
	}
	// Large n: exact only when x is a perfect n-th power of a dyadic number.
	if v := exactRoot(x, n); v != nil {
		if neg {
			v.Neg(v)
		}
		return finish(z, v, 0)
	}
	ax := new(big.Float).Abs(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		var l, nb, t ball
		logBall(&l, ax, wp+64)
		nb.m.SetUint64(uint64(n))
		must(t.quo(&l, &nb, wp+64))
		b, K := expReduced(&t.m, t.r, wp)
		if neg {
			b.neg(b) // odd n: the root of -|x| is -(|x|^(1/n))
		}
		return b, K, true
	})
}

// exactRoot returns x^(1/n) when it is a dyadic number, else nil
// (x != 0 finite; sign ignored). x = m·2^e with m odd: the root is dyadic
// exactly when m is a perfect n-th power and n divides e.
func exactRoot(x *big.Float, n uint) *big.Float {
	M, E := mantInt(x)
	tz := M.TrailingZeroBits()
	M.Rsh(M, tz)
	E += int64(tz)
	if E%int64(n) != 0 {
		return nil
	}
	// m = 1: the root is 2^(e/n) (no integer root to extract); m > 1 an
	// n-th power needs m >= 3^n, more than 1.5n bits.
	if M.BitLen() == 1 {
		return new(big.Float).SetMantExp(big.NewFloat(1), int(E/int64(n)))
	}
	if uint(M.BitLen()) < n+n/2 {
		return nil
	}
	r, exact := iroot(M, n)
	if !exact {
		return nil
	}
	v := new(big.Float).SetPrec(0).SetInt(r)
	return v.SetMantExp(v, int(E/int64(n)))
}

// Hypot sets z to the correctly rounded sqrt(x²+y²) and returns z. If z's
// precision is 0 it is set to the larger of x's and y's. Hypot(±Inf, y) =
// Hypot(x, ±Inf) = +Inf.
func Hypot(z, x, y *big.Float) *big.Float {
	precOf(z, max(x.Prec(), y.Prec()))
	if x.IsInf() || y.IsInf() {
		return z.SetInf(false)
	}
	a := new(big.Float).Abs(x)
	b := new(big.Float).Abs(y)
	if a.Cmp(b) < 0 {
		a, b = b, a
	}
	if a.Sign() == 0 {
		return z.SetInt64(0)
	}
	if b.Sign() == 0 {
		return z.Set(a)
	}
	p := z.Prec()
	ea, eb := a.MantExp(nil), b.MantExp(nil)
	// hypot = a·sqrt(1+(b/a)²) in (a, a + a(b/a)²/2) and a(b/a)²/2 <
	// 2^(2eb-ea+1): when that is below nearBound(a), a proxies the result.
	if int64(2*eb-ea+1) <= nearBound(a, p) {
		return finish(z, proxyNear(a, 1, p), 0)
	}
	// Scale by 2^-ea (exact), then x²+y² exactly: its terms differ in
	// exponent by at most about p+prec(a), so the exact sum fits.
	a.SetMantExp(a, -ea)
	b.SetMantExp(b, -ea)
	s := new(big.Float).SetPrec(2*a.Prec() + 2*b.Prec() + 2*uint(ea-eb) + 8)
	t := new(big.Float).SetPrec(2 * b.Prec())
	s.Mul(a, a)
	t.Mul(b, b)
	s.Add(s, t)
	return finish(z, rootProxy(s, 2, p), int64(ea))
}
