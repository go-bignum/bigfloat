// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat_test

import (
	"errors"
	"math/big"
	"os/exec"
	"slices"
	"testing"

	"github.com/go-bignum/bigfloat/oracle"
)

// hostileCase is an argument tuple chosen at a boundary of the
// implementation: a grid point of a table reduction (exactly on it, or
// within an ulp), the edge between a table and the direct series, the
// fixed-point asin bound, huge and tiny exponents.
type hostileCase struct {
	fn   string
	args []string
}

func hf(s string) *big.Float {
	x, _, err := big.ParseFloat(s, 0, 4096, big.ToNearestEven)
	if err != nil {
		panic(err)
	}
	return x.SetPrec(max(x.MinPrec(), 1))
}

// hostileCases lists them; an argument "a+u" or "a-u" is a perturbed by
// 2^-u, at full precision.
var hostileCases = func() []hostileCase {
	var cs []hostileCase
	add := func(fn string, args ...string) { cs = append(cs, hostileCase{fn, args}) }
	perturb := []string{"", "+0x1p-60", "-0x1p-60", "+0x1p-130", "-0x1p-1100"}
	for _, n := range []string{"1", "64", "127", "128"} { // atan grid n/2^7
		for _, d := range perturb {
			add("atan", n+"/128"+d)
			add("atan", "-"+n+"/128"+d)
		}
	}
	add("atan", "0x1p-7")
	add("atan", "1e300")
	add("atan", "-0x1p-1000")
	add("atan", "3")
	for _, n := range []string{"362", "509", "510", "511", "513", "514", "724"} { // log grid n/2^9
		for _, d := range perturb {
			for _, fn := range []string{"log", "log2", "log10"} {
				add(fn, n+"/512"+d)
			}
			add("log1p", n+"/512-1"+d)
		}
	}
	for _, x := range []string{"1+0x1p-8", "1-0x1p-8", "1+0x1p-9", "1-0x1p-9", "1+0x1p-200", "0x1.3p100000", "0x1.3p-100000"} {
		add("log", x)
		add("log2", x)
	}
	for _, n := range []string{"1", "100", "201", "202"} { // sin/cos grid n/2^8
		for _, d := range perturb {
			for _, fn := range []string{"sin", "cos", "tan"} {
				add(fn, n+"/256"+d)
			}
		}
	}
	for _, x := range []string{"0x1p-8", "0x1p100", "0x1.8p1000", "1e22", "-1e22", "0x1p-500"} {
		for _, fn := range []string{"sin", "cos", "tan"} {
			add(fn, x)
		}
	}
	for _, x := range []string{"0.86", "0.86+0x1p-60", "0.86-0x1p-60", "-0.86", "0.5", "0x1p-500", "1-0x1p-100", "-1+0x1p-100"} {
		add("asin", x)
		add("acos", x)
	}
	for _, x := range []string{"0.7", "-0.7", "1000.5", "17179869183", "-17179869183", "0x1p-100", "709.78"} {
		for _, fn := range []string{"exp", "expm1", "sinh", "cosh", "tanh"} {
			add(fn, x)
		}
	}
	for _, x := range []string{"0.7", "1e100", "0x1p-300", "-3.5"} {
		add("asinh", x)
	}
	for _, x := range []string{"1.7", "1+0x1p-90", "1e100"} {
		add("acosh", x)
	}
	for _, x := range []string{"0.7", "-0.5", "1-0x1p-90", "0x1p-300"} {
		add("atanh", x)
	}
	add("pow", "1.7", "2.3")
	add("pow", "2", "0.5")
	add("pow", "1+0x1p-60", "0x1p40")
	add("pow", "0.5", "0x1p33")
	add("pow", "10", "-400")
	add("pow", "1.000001", "1e-30")
	add("root", "2", "7")
	add("root", "3", "5")
	add("cbrt", "2")
	add("cbrt", "-0x1.5p-3000")
	add("atan2", "0.7", "1.3")
	add("atan2", "-1", "-0x1p-200")
	return cs
}()

// hostilePrecs straddle the table limits (sin/cos tables up to 640 bits).
var hostilePrecs = []uint{2, 24, 53, 113, 200, 640, 641, 1000, 2100}

// parseArg evaluates a sum of terms "a" or "a/b" (b a power of two), such
// as "362/512-1+0x1p-60", exactly.
func parseArg(s string) *big.Float {
	sum := new(big.Float).SetPrec(8192)
	start := 0
	for i := 1; i <= len(s); i++ {
		if i < len(s) && !((s[i] == '+' || s[i] == '-') && s[i-1] != 'p' && s[i-1] != 'e') {
			continue
		}
		term := s[start:i]
		neg := term[0] == '-'
		if term[0] == '+' || term[0] == '-' {
			term = term[1:]
		}
		var x *big.Float
		if j := slices.Index([]byte(term), '/'); j > 0 {
			x = new(big.Float).SetPrec(8192).Quo(hf(term[:j]), hf(term[j+1:]))
		} else {
			x = hf(term)
		}
		if neg {
			sum.Sub(sum, x)
		} else {
			sum.Add(sum, x)
		}
		start = i
	}
	return sum.SetPrec(max(sum.MinPrec(), 1))
}

func hostileRun(t *testing.T, o oracle.Oracle) {
	t.Helper()
	rd := &oracle.Rounder{Oracle: o}
	bad := 0
	for _, c := range hostileCases {
		args := make([]*big.Float, len(c.args))
		for i, a := range c.args {
			args[i] = parseArg(a)
		}
		for _, p := range hostilePrecs {
			want, err := rd.Round(c.fn, p, args...)
			if errors.Is(err, oracle.ErrUndecided) {
				continue
			}
			if err != nil {
				t.Fatalf("%s%v at %d bits: %v", c.fn, c.args, p, err)
			}
			for i, m := range oracle.Modes {
				z := new(big.Float).SetPrec(p).SetMode(m)
				oracle.Call(c.fn, z, args...)
				if z.Cmp(want[i].Value) != 0 || z.Acc() != want[i].Acc {
					if bad++; bad <= 10 {
						t.Errorf("%s%v at %d bits, %v: got %v (%v), want %v (%v) [%s]", c.fn, c.args, p, m, z, z.Acc(), want[i].Value, want[i].Acc, o.Name())
					}
				}
			}
		}
	}
	if bad > 10 {
		t.Errorf("... %d mismatches in all", bad)
	}
}

// TestHostileSelf checks the hostile points against the package itself at
// more than twice the precision (through the Rounder, as the fuzz test).
func TestHostileSelf(t *testing.T) {
	if testing.Short() {
		t.Skip("short")
	}
	hostileRun(t, oracle.Self{})
}

// TestHostileMpmath checks them against mpmath when uv is installed (CI
// need not have it; run locally before changing an evaluation path).
func TestHostileMpmath(t *testing.T) {
	if _, err := exec.LookPath("uv"); err != nil || testing.Short() {
		t.Skip("uv not installed")
	}
	o, err := oracle.StartMpmath(nil)
	if err != nil {
		t.Skip(err)
	}
	defer o.Close()
	hostileRun(t, o)
}
