// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat_test

import (
	"math/big"
	"strings"
	"testing"

	"github.com/go-bignum/bigfloat"
	"github.com/go-bignum/bigfloat/oracle"
)

// f parses a decimal or hexadecimal ("0x.8p+1") number exactly, at 200
// bits, or ±Inf / -0.
func f(s string) *big.Float {
	switch s {
	case "+Inf", "-Inf":
		return new(big.Float).SetInf(s[0] == '-')
	case "-0":
		return new(big.Float).Neg(new(big.Float))
	}
	x, _, err := big.ParseFloat(s, 0, 200, big.ToNearestEven)
	if err != nil {
		panic(err)
	}
	return x
}

func z(p uint, m big.RoundingMode) *big.Float { return new(big.Float).SetPrec(p).SetMode(m) }

// piCase expects k·π/2^d, rounded to nearest at 53 bits: the rounding of a
// 400-bit π times k, which is far from any 53-bit rounding boundary.
func piCase(name string, got func() *big.Float, k int64, d int) special {
	v := bigfloat.Pi(z(400, big.ToNearestEven))
	v.Mul(v, big.NewFloat(float64(k)))
	v.SetMantExp(v, -d)
	r := z(53, big.ToNearestEven).Set(v)
	return special{name, got, oracle.Format(r), r.Acc()}
}

func intS(n int64) string { return oracle.Format(new(big.Float).SetInt64(n)) }

type special struct {
	name string
	got  func() *big.Float
	want string
	acc  big.Accuracy
}

func checkSpecials(t *testing.T, cases []special) {
	t.Helper()
	for _, c := range cases {
		g := c.got()
		if s := oracle.Format(g); s != c.want || g.Acc() != c.acc {
			t.Errorf("%s = %s (%v), want %s (%v)", c.name, s, g.Acc(), c.want, c.acc)
		}
	}
}

const (
	ne = big.ToNearestEven
	up = big.ToPositiveInf
	dn = big.ToNegativeInf
	ex = big.Exact
	ab = big.Above
	be = big.Below
)

func TestSpecialValues(t *testing.T) {
	one := func(fn func(z, x *big.Float) *big.Float, x string, m big.RoundingMode) func() *big.Float {
		return func() *big.Float { return fn(z(53, m), f(x)) }
	}
	two := func(fn func(z, x, y *big.Float) *big.Float, x, y string, m big.RoundingMode) func() *big.Float {
		return func() *big.Float { return fn(z(53, m), f(x), f(y)) }
	}
	minHalf := "0x.8p-2147483648" // 2^(MinExp-1), the smallest positive big.Float
	checkSpecials(t, []special{
		{"Exp(+Inf)", one(bigfloat.Exp, "+Inf", ne), "+Inf", ex},
		{"Exp(-Inf)", one(bigfloat.Exp, "-Inf", ne), "0", ex},
		{"Exp(-0)", one(bigfloat.Exp, "-0", ne), "0x.8p+1", ex},
		{"Exp(1e11)", one(bigfloat.Exp, "1e11", ne), "+Inf", ab},
		{"Exp(-1e11)", one(bigfloat.Exp, "-1e11", ne), "0", be},
		{"Exp(1.6e9)", one(bigfloat.Exp, "1.6e9", ne), "+Inf", ab},
		{"Exp(-1.6e9)", one(bigfloat.Exp, "-1.6e9", ne), "0", be},
		{"Exp(near MinExp-1)", func() *big.Float {
			// (MinExp-1)·ln2 - 2^-80: e^x is 2^(MinExp-1)(1-2^-80·…), which
			// rounds up to the smallest positive number.
			l := bigfloat.Ln2(z(400, ne))
			l.Mul(l, new(big.Float).SetInt64(big.MinExp-1))
			l.Sub(l, f("0x1p-80"))
			return bigfloat.Exp(z(24, ne), l)
		}, minHalf, ex},
		{"Expm1(+Inf)", one(bigfloat.Expm1, "+Inf", ne), "+Inf", ex},
		{"Expm1(-Inf)", one(bigfloat.Expm1, "-Inf", ne), "-0x.8p+1", ex},
		{"Expm1(-0)", one(bigfloat.Expm1, "-0", ne), "-0", ex},
		{"Expm1(1e11)", one(bigfloat.Expm1, "1e11", ne), "+Inf", ab},
		{"Expm1(-1e11)", one(bigfloat.Expm1, "-1e11", ne), "-0x.8p+1", be},
		{"Expm1(-1e11) up", one(bigfloat.Expm1, "-1e11", up), "-0x.fffffffffffff8p+0", ab},
		{"Log(0)", one(bigfloat.Log, "0", ne), "-Inf", ex},
		{"Log(+Inf)", one(bigfloat.Log, "+Inf", ne), "+Inf", ex},
		{"Log(1)", one(bigfloat.Log, "1", ne), "0", ex},
		{"Log2(-0)", one(bigfloat.Log2, "-0", ne), "-Inf", ex},
		{"Log2(+Inf)", one(bigfloat.Log2, "+Inf", ne), "+Inf", ex},
		{"Log2(2^-1074)", one(bigfloat.Log2, "0x1p-1074", ne), intS(-1074), ex},
		{"Log10(0)", one(bigfloat.Log10, "0", ne), "-Inf", ex},
		{"Log10(+Inf)", one(bigfloat.Log10, "+Inf", ne), "+Inf", ex},
		{"Log10(1000)", one(bigfloat.Log10, "1000", ne), "0x.cp+2", ex},
		{"Log1p(-1)", one(bigfloat.Log1p, "-1", ne), "-Inf", ex},
		{"Log1p(+Inf)", one(bigfloat.Log1p, "+Inf", ne), "+Inf", ex},
		{"Log1p(-0)", one(bigfloat.Log1p, "-0", ne), "-0", ex},
		{"Sqrt(-0)", one(bigfloat.Sqrt, "-0", ne), "-0", ex},
		{"Sqrt(+Inf)", one(bigfloat.Sqrt, "+Inf", ne), "+Inf", ex},
		{"Sqrt(4)", one(bigfloat.Sqrt, "4", dn), "0x.8p+2", ex},
		{"Cbrt(-Inf)", one(bigfloat.Cbrt, "-Inf", ne), "-Inf", ex},
		{"Cbrt(-0)", one(bigfloat.Cbrt, "-0", ne), "-0", ex},
		{"Root(-0, 2)", func() *big.Float { return bigfloat.Root(z(53, ne), f("-0"), 2) }, "0", ex},
		{"Root(x, 1)", func() *big.Float { return bigfloat.Root(z(53, ne), f("-3"), 1) }, "-0x.cp+2", ex},
		{"Root(2^(3·2^20), 2^20)", func() *big.Float {
			x := new(big.Float).SetMantExp(big.NewFloat(1), 3<<20)
			return bigfloat.Root(z(53, ne), x, 1<<20)
		}, "0x.8p+4", ex},
		{"Root(-3^1001, 1001) at 2000 bits", func() *big.Float {
			n := new(big.Int).Exp(big.NewInt(3), big.NewInt(1001), nil)
			x := new(big.Float).SetInt(n)
			return bigfloat.Root(z(2000, ne), x.Neg(x), 1001)
		}, "-0x.cp+2", ex},
		{"Hypot(-Inf, 1)", two(bigfloat.Hypot, "-Inf", "1", ne), "+Inf", ex},
		{"Hypot(1, +Inf)", two(bigfloat.Hypot, "1", "+Inf", ne), "+Inf", ex},
		{"Hypot(0, -0)", two(bigfloat.Hypot, "0", "-0", ne), "0", ex},
		{"Hypot(-0, -3)", two(bigfloat.Hypot, "-0", "-3", ne), "0x.cp+2", ex},
		{"Hypot(-3, 4)", two(bigfloat.Hypot, "-3", "4", ne), "0x.ap+3", ex},
		{"Pow(x, -0)", two(bigfloat.Pow, "-Inf", "-0", ne), "0x.8p+1", ex},
		{"Pow(1, NaN-ish Inf)", two(bigfloat.Pow, "1", "+Inf", ne), "0x.8p+1", ex},
		{"Pow(-1, +Inf)", two(bigfloat.Pow, "-1", "+Inf", ne), "0x.8p+1", ex},
		{"Pow(0.5, +Inf)", two(bigfloat.Pow, "0.5", "+Inf", ne), "0", ex},
		{"Pow(0.5, -Inf)", two(bigfloat.Pow, "0.5", "-Inf", ne), "+Inf", ex},
		{"Pow(-2, +Inf)", two(bigfloat.Pow, "-2", "+Inf", ne), "+Inf", ex},
		{"Pow(2, -Inf)", two(bigfloat.Pow, "2", "-Inf", ne), "0", ex},
		{"Pow(-0, -3)", two(bigfloat.Pow, "-0", "-3", ne), "-Inf", ex},
		{"Pow(-0, -2)", two(bigfloat.Pow, "-0", "-2", ne), "+Inf", ex},
		{"Pow(-0, 3)", two(bigfloat.Pow, "-0", "3", ne), "-0", ex},
		{"Pow(-0, 2)", two(bigfloat.Pow, "-0", "2", ne), "0", ex},
		{"Pow(-0, 0.5)", two(bigfloat.Pow, "-0", "0.5", ne), "0", ex},
		{"Pow(-Inf, 3)", two(bigfloat.Pow, "-Inf", "3", ne), "-Inf", ex},
		{"Pow(-Inf, -3)", two(bigfloat.Pow, "-Inf", "-3", ne), "-0", ex},
		{"Pow(-Inf, 2)", two(bigfloat.Pow, "-Inf", "2", ne), "+Inf", ex},
		{"Pow(+Inf, -0.5)", two(bigfloat.Pow, "+Inf", "-0.5", ne), "0", ex},
		{"Pow(3, 1e11)", two(bigfloat.Pow, "3", "1e11", ne), "+Inf", ab},
		{"Pow(3, -1e11)", two(bigfloat.Pow, "3", "-1e11", ne), "0", be},
		{"Pow(-3, 2^40+1)", two(bigfloat.Pow, "-3", "1099511627777", ne), "-Inf", be},
		{"Pow(2, 2^70)", two(bigfloat.Pow, "2", "0x1p70", ne), "+Inf", ab},
		{"Pow(2, -2^70)", two(bigfloat.Pow, "2", "-0x1p70", ne), "0", be},
		{"Pow(0.5, 2^70)", two(bigfloat.Pow, "0.5", "0x1p70", ne), "0", be},
		{"Pow(3, 2^70)", two(bigfloat.Pow, "3", "0x1p70", ne), "+Inf", ab},
		{"Pow(1+2^-50, 2^-60)", two(bigfloat.Pow, "0x1.0000000000004p0", "0x1p-60", ne), "0x.8p+1", be},
		{"Pow(1+2^-50, 2^-60) up", two(bigfloat.Pow, "0x1.0000000000004p0", "0x1p-60", up), "0x.80000000000008p+1", ab},
		{"Pow(-(1+2^-50), 1)", two(bigfloat.Pow, "-0x1.0000000000004p0", "1", ne), "-0x.8000000000002p+1", ex},
		{"Pow(-1.5, 4)", two(bigfloat.Pow, "-1.5", "4", ne), "0x.a2p+3", ex},
		{"Pow(-3, -(2^40+1))", two(bigfloat.Pow, "-3", "-1099511627777", ne), "-0", ab},
		// -(1+2^-150)^(2^63+1) = -(1 + 2^-87·…): just below -1
		{"Pow(-(1+2^-150), 2^63+1)", two(bigfloat.Pow, "-0x1.000000000000000000000000000000000000004p0", "9223372036854775809", ne), "-0x.8p+1", ab}, {"Sin(-0)", one(bigfloat.Sin, "-0", ne), "-0", ex},
		{"Cos(-0)", one(bigfloat.Cos, "-0", ne), "0x.8p+1", ex},
		{"Tan(-0)", one(bigfloat.Tan, "-0", ne), "-0", ex},
		piCase("Atan(+Inf)", one(bigfloat.Atan, "+Inf", ne), 1, 1),
		piCase("Atan(-Inf)", one(bigfloat.Atan, "-Inf", ne), -1, 1),
		{"Atan(-0)", one(bigfloat.Atan, "-0", ne), "-0", ex},
		{"Atan2(+0, +0)", two(bigfloat.Atan2, "0", "0", ne), "0", ex},
		{"Atan2(-0, +0)", two(bigfloat.Atan2, "-0", "0", ne), "-0", ex},
		piCase("Atan2(+0, -0)", two(bigfloat.Atan2, "0", "-0", ne), 1, 0),
		piCase("Atan2(-0, -0)", two(bigfloat.Atan2, "-0", "-0", ne), -1, 0),
		piCase("Atan2(-0, -1)", two(bigfloat.Atan2, "-0", "-1", ne), -1, 0),
		{"Atan2(0, 1)", two(bigfloat.Atan2, "0", "1", ne), "0", ex},
		piCase("Atan2(1, -0)", two(bigfloat.Atan2, "1", "-0", ne), 1, 1),
		piCase("Atan2(-1, 0)", two(bigfloat.Atan2, "-1", "0", ne), -1, 1),
		piCase("Atan2(+Inf, +Inf)", two(bigfloat.Atan2, "+Inf", "+Inf", ne), 1, 2),
		piCase("Atan2(-Inf, -Inf)", two(bigfloat.Atan2, "-Inf", "-Inf", ne), -3, 2),
		piCase("Atan2(+Inf, 1)", two(bigfloat.Atan2, "+Inf", "1", ne), 1, 1),
		{"Atan2(1, +Inf)", two(bigfloat.Atan2, "1", "+Inf", ne), "0", ex},
		{"Atan2(-1, +Inf)", two(bigfloat.Atan2, "-1", "+Inf", ne), "-0", ex},
		piCase("Atan2(1, -Inf)", two(bigfloat.Atan2, "1", "-Inf", ne), 1, 0),
		piCase("Atan2(-1, -Inf)", two(bigfloat.Atan2, "-1", "-Inf", ne), -1, 0),
		{"Atan2(2^-100, 1)", two(bigfloat.Atan2, "0x1p-100", "1", ne), "0x.8p-99", ab},
		{"Atan2(2^-100, 1) down", two(bigfloat.Atan2, "0x1p-100", "1", dn), "0x.fffffffffffff8p-100", be},
		piCase("Asin(1)", one(bigfloat.Asin, "1", ne), 1, 1),
		piCase("Asin(-1)", one(bigfloat.Asin, "-1", ne), -1, 1),
		{"Asin(-0)", one(bigfloat.Asin, "-0", ne), "-0", ex},
		{"Acos(1)", one(bigfloat.Acos, "1", ne), "0", ex},
		piCase("Acos(-1)", one(bigfloat.Acos, "-1", ne), 1, 0),
		piCase("Acos(-0)", one(bigfloat.Acos, "-0", ne), 1, 1),
		{"Sinh(-Inf)", one(bigfloat.Sinh, "-Inf", ne), "-Inf", ex},
		{"Sinh(-0)", one(bigfloat.Sinh, "-0", ne), "-0", ex},
		{"Sinh(-1e11)", one(bigfloat.Sinh, "-1e11", ne), "-Inf", be},
		{"Cosh(-Inf)", one(bigfloat.Cosh, "-Inf", ne), "+Inf", ex},
		{"Cosh(-0)", one(bigfloat.Cosh, "-0", ne), "0x.8p+1", ex},
		{"Cosh(1e11)", one(bigfloat.Cosh, "1e11", ne), "+Inf", ab},
		{"Tanh(+Inf)", one(bigfloat.Tanh, "+Inf", ne), "0x.8p+1", ex},
		{"Tanh(-Inf)", one(bigfloat.Tanh, "-Inf", ne), "-0x.8p+1", ex},
		{"Tanh(-0)", one(bigfloat.Tanh, "-0", ne), "-0", ex},
		{"Asinh(-Inf)", one(bigfloat.Asinh, "-Inf", ne), "-Inf", ex},
		{"Asinh(-0)", one(bigfloat.Asinh, "-0", ne), "-0", ex},
		{"Acosh(1)", one(bigfloat.Acosh, "1", ne), "0", ex},
		{"Acosh(+Inf)", one(bigfloat.Acosh, "+Inf", ne), "+Inf", ex},
		{"Atanh(1)", one(bigfloat.Atanh, "1", ne), "+Inf", ex},
		{"Atanh(-1)", one(bigfloat.Atanh, "-1", ne), "-Inf", ex},
		{"Atanh(-0)", one(bigfloat.Atanh, "-0", ne), "-0", ex},
	})
}

// TestNaN checks that every domain error panics with bigfloat.ErrNaN.
func TestNaN(t *testing.T) {
	cases := map[string]func(){
		"Sqrt(-1)":     func() { bigfloat.Sqrt(z(53, ne), f("-1")) },
		"Root(1, 0)":   func() { bigfloat.Root(z(53, ne), f("1"), 0) },
		"Root(-8, 2)":  func() { bigfloat.Root(z(53, ne), f("-8"), 2) },
		"Log(-1)":      func() { bigfloat.Log(z(53, ne), f("-1")) },
		"Log2(-1)":     func() { bigfloat.Log2(z(53, ne), f("-1")) },
		"Log10(-1)":    func() { bigfloat.Log10(z(53, ne), f("-1")) },
		"Log1p(-2)":    func() { bigfloat.Log1p(z(53, ne), f("-2")) },
		"Pow(-2, 0.5)": func() { bigfloat.Pow(z(53, ne), f("-2"), f("0.5")) },
		"Sin(+Inf)":    func() { bigfloat.Sin(z(53, ne), f("+Inf")) },
		"Cos(-Inf)":    func() { bigfloat.Cos(z(53, ne), f("-Inf")) },
		"Tan(+Inf)":    func() { bigfloat.Tan(z(53, ne), f("+Inf")) },
		"Asin(2)":      func() { bigfloat.Asin(z(53, ne), f("2")) },
		"Acos(-2)":     func() { bigfloat.Acos(z(53, ne), f("-2")) },
		"Acosh(0.5)":   func() { bigfloat.Acosh(z(53, ne), f("0.5")) },
		"Atanh(-1.5)":  func() { bigfloat.Atanh(z(53, ne), f("-1.5")) },
	}
	for name, fn := range cases {
		func() {
			defer func() {
				e, ok := recover().(bigfloat.ErrNaN)
				if !ok || !strings.HasPrefix(e.Error(), "bigfloat: ") {
					t.Errorf("%s: recovered %v, want a bigfloat.ErrNaN", name, e)
				}
			}()
			fn()
		}()
	}
}

// TestLimit: a trigonometric argument beyond MaxReductionExponent panics
// with ErrLimit instead of computing millions of bits of π.
func TestLimit(t *testing.T) {
	x := new(big.Float).SetMantExp(big.NewFloat(1), bigfloat.MaxReductionExponent+1)
	defer func() {
		e, ok := recover().(bigfloat.ErrLimit)
		if !ok || !strings.Contains(e.Error(), "MaxReductionExponent") {
			t.Errorf("recovered %v, want an ErrLimit", e)
		}
	}()
	bigfloat.Cos(z(53, ne), x)
}

// TestDefaultPrecision checks the big.Float conventions: a receiver of
// precision 0 takes the operands' (the larger of two), constants 64 bits;
// the receiver may alias an operand.
func TestDefaultPrecision(t *testing.T) {
	x := new(big.Float).SetPrec(100).SetInt64(2)
	y := new(big.Float).SetPrec(30).SetInt64(3)
	if p := bigfloat.Exp(new(big.Float), x).Prec(); p != 100 {
		t.Errorf("Exp: prec %d, want 100", p)
	}
	if p := bigfloat.Pow(new(big.Float), y, x).Prec(); p != 100 {
		t.Errorf("Pow: prec %d, want 100", p)
	}
	if p := bigfloat.Pi(new(big.Float)).Prec(); p != 64 {
		t.Errorf("Pi: prec %d, want 64", p)
	}
	want := bigfloat.Log(new(big.Float).SetPrec(100), x)
	if bigfloat.Log(x, x).Cmp(want) != 0 {
		t.Errorf("Log(x, x) = %v, want %v", x, want)
	}
}

// TestReductionNearMultiple checks a sine whose argument is within 2^-2900
// of a multiple of π/2: the reduction must raise its precision (the first
// enclosures of the reduced argument contain zero).
func TestReductionNearMultiple(t *testing.T) {
	x := bigfloat.Pi(z(3000, ne))
	x.Mul(x, big.NewFloat(1e6))
	// sin(x) for x ≈ 1e6·π is ±(x - 1e6·π) (to first order)
	pi := bigfloat.Pi(z(6000, ne))
	pi.Mul(pi, big.NewFloat(1e6))
	d := new(big.Float).SetPrec(6000).Sub(x, pi)
	want := z(24, ne).Set(d) // sin(1e6π + d) = sin(d) ≈ d
	got := bigfloat.Sin(z(24, ne), x)
	if got.Cmp(want) != 0 {
		t.Errorf("Sin(1e6·π rounded to 3000 bits) = %v, want %v", got, want)
	}
}
