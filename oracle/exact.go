// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package oracle

import (
	"math/big"
)

// Exact reports the exact value of fn(args) when it is a dyadic number
// (representable by a big.Float of some precision), the only case where an
// interval from an approximate oracle can never decide a rounding. It
// covers the exact cases of every function in [Functions] for finite
// arguments: zeros and ones (exp(0), log(1), sin(0), …), integer
// logarithms (log2 of powers of two, log10 of 10^k), perfect powers under
// roots and hypot, and dyadic powers x^y. Everything else these functions
// return at a dyadic argument is irrational (Lindemann–Weierstrass and
// Gelfond–Schneider), which is why the list is complete.
//
// The arithmetic here is deliberately independent of the bigfloat package
// it is used to test: plain big.Int and big.Rat, binary search for roots.
func Exact(fn string, args ...*big.Float) (*big.Float, bool) {
	one := big.NewFloat(1)
	x := args0(args)
	switch fn {
	case "exp", "cos", "cosh":
		if x.Sign() == 0 {
			return one, true
		}
	case "expm1", "log1p", "sin", "tan", "asin", "atan", "sinh", "tanh", "asinh", "atanh":
		if x.Sign() == 0 {
			return x, true
		}
	case "log", "acos", "acosh":
		if x.Cmp(one) == 0 {
			return new(big.Float), true
		}
	case "log2":
		if m, e := oddPart(x); x.Sign() > 0 && m.Cmp(big.NewInt(1)) == 0 {
			return new(big.Float).SetInt64(e), true
		}
	case "log10":
		if k, ok := log10Exact(x); ok {
			return new(big.Float).SetInt64(k), true
		}
	case "sqrt":
		return rootExact(x, 2)
	case "cbrt":
		return rootExact(x, 3)
	case "root":
		n, _ := args[1].Uint64()
		return rootExact(x, n)
	case "hypot":
		y := args[1]
		s := new(big.Float).SetPrec(4*(x.Prec()+y.Prec()) + 2*expGap(x, y) + 16)
		t := new(big.Float).SetPrec(s.Prec())
		s.Mul(x, x)
		t.Mul(y, y)
		s.Add(s, t)
		return rootExact(s, 2)
	case "pow":
		return powExact(x, args[1])
	}
	return nil, false
}

func args0(args []*big.Float) *big.Float {
	if len(args) == 0 {
		return new(big.Float)
	}
	return args[0]
}

// expGap is the difference of the binary exponents of two non-zero numbers
// (0 if either is zero).
func expGap(x, y *big.Float) uint {
	if x.Sign() == 0 || y.Sign() == 0 {
		return 0
	}
	d := x.MantExp(nil) - y.MantExp(nil)
	if d < 0 {
		d = -d
	}
	return uint(d)
}

// oddPart returns m odd and e with |x| = m·2^e (x finite non-zero), or
// (0, 0) for x = 0.
func oddPart(x *big.Float) (*big.Int, int64) {
	if x.Sign() == 0 {
		return new(big.Int), 0
	}
	var r big.Rat
	x.Rat(&r)
	num := new(big.Int).Abs(r.Num())
	den := r.Denom()
	e := int64(num.TrailingZeroBits())
	num.Rsh(num, uint(e))
	e -= int64(den.BitLen() - 1) // den is a power of two
	return num, e
}

// log10Exact: x = 10^k exactly, k >= 0 (10^-k is not dyadic).
func log10Exact(x *big.Float) (int64, bool) {
	if x.Sign() <= 0 || !x.IsInt() {
		return 0, false
	}
	n, _ := x.Int(nil)
	ten := big.NewInt(10)
	k := int64(0)
	var q, r big.Int
	for n.Cmp(big.NewInt(1)) != 0 {
		q.QuoRem(n, ten, &r)
		if r.Sign() != 0 {
			return 0, false
		}
		n.Set(&q)
		k++
	}
	return k, true
}

// nthRoot returns r with r^n = m exactly, by binary search.
func nthRoot(m *big.Int, n uint64) (*big.Int, bool) {
	lo, hi := big.NewInt(0), new(big.Int).Lsh(big.NewInt(1), uint(m.BitLen())/uint(n)+1)
	N := new(big.Int).SetUint64(n)
	var mid, p big.Int
	one := big.NewInt(1)
	for lo.Cmp(hi) <= 0 {
		mid.Add(lo, hi)
		mid.Rsh(&mid, 1)
		p.Exp(&mid, N, nil)
		switch p.Cmp(m) {
		case 0:
			return new(big.Int).Set(&mid), true
		case -1:
			lo.Add(&mid, one)
		default:
			hi.Sub(&mid, one)
		}
	}
	return nil, false
}

// rootExact: |x|^(1/n) is dyadic iff the odd part is a perfect n-th power
// and n divides the exponent.
func rootExact(x *big.Float, n uint64) (*big.Float, bool) {
	if x.Sign() == 0 {
		return new(big.Float), true
	}
	m, e := oddPart(x)
	if e%int64(n) != 0 {
		return nil, false
	}
	r, ok := nthRoot(m, n)
	if !ok {
		return nil, false
	}
	v := new(big.Float).SetInt(r)
	v.SetMantExp(v, int(e/int64(n)))
	if x.Sign() < 0 {
		v.Neg(v)
	}
	return v, true
}

// powExact: x^y for x > 0 dyadic, y = k/2^j: dyadic iff x is a perfect
// 2^j-th power b and b^k is dyadic (k < 0 needs b a power of two). Only
// |k| <= 1<<16 is computed; larger exact powers are not used in fixtures.
func powExact(x, y *big.Float) (*big.Float, bool) {
	if y.Sign() == 0 || x.Cmp(big.NewFloat(1)) == 0 {
		return big.NewFloat(1), true
	}
	if x.Sign() <= 0 {
		if x.Sign() < 0 && y.IsInt() {
			v, ok := powExact(new(big.Float).Neg(x), y)
			if ok && isOdd(y) {
				v.Neg(v)
			}
			return v, ok
		}
		return nil, false
	}
	var yr big.Rat
	y.Rat(&yr)
	j := yr.Denom().BitLen() - 1
	k := yr.Num()
	if j > 64 || k.BitLen() > 17 {
		return nil, false
	}
	b := new(big.Float).Set(x)
	for i := 0; i < j; i++ {
		s, ok := rootExact(b, 2)
		if !ok {
			return nil, false
		}
		b = s
	}
	var br, p big.Rat
	b.Rat(&br)
	kk := k.Int64()
	neg := kk < 0
	if neg {
		kk = -kk
	}
	num := new(big.Int).Exp(br.Num(), big.NewInt(kk), nil)
	den := new(big.Int).Exp(br.Denom(), big.NewInt(kk), nil)
	if neg {
		num, den = den, num
	}
	p.SetFrac(num, den)
	d := p.Denom()
	if d.BitLen()-1 != int(d.TrailingZeroBits()) {
		return nil, false // denominator not a power of two
	}
	v := new(big.Float).SetPrec(uint(num.BitLen()) + 64)
	v.SetRat(&p)
	return v, true
}

func isOdd(y *big.Float) bool {
	n, _ := y.Int(nil)
	return n.Bit(0) == 1
}
