// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math"
	"math/big"
	"strings"
	"testing"
)

func panics(t *testing.T, name, want string, fn func()) {
	t.Helper()
	defer func() {
		e := recover()
		if s, _ := e.(string); !strings.Contains(s, want) {
			t.Errorf("%s: recovered %v, want a panic containing %q", name, e, want)
		}
	}()
	fn()
}

func TestMust(t *testing.T) {
	must(true)
	panics(t, "must(false)", "invariant", func() { must(false) })
}

// TestZivLimit: an enclosure of an exact value never decides; the loop
// gives up with a panic past its limit instead of running forever.
func TestZivLimit(t *testing.T) {
	saved := zivLimit
	defer func() { zivLimit = saved }()
	zivLimit = func(p, in uint) uint { return 200 }
	panics(t, "ziv", "did not converge", func() {
		ziv(new(big.Float).SetPrec(24), 0, func(wp uint) (*ball, int64, bool) {
			var b ball
			b.setInt64(1)
			b.r = magPow2(-int64(wp))
			return &b, 0, true
		})
	})
}

// TestBallFailures exercises the refusals of ball operations whose
// operands may vanish.
func TestBallFailures(t *testing.T) {
	var zero, one, z ball
	zero.setInt64(0)
	zero.r = magPow2(-10)
	one.setInt64(1)
	if z.quo(&one, &zero, 64) {
		t.Error("quo by a ball containing 0 succeeded")
	}
	if z.sqrt(&zero, 64) {
		t.Error("sqrt of a ball containing 0 succeeded")
	}
	var neg ball
	neg.setInt64(-4)
	if z.sqrt(&neg, 64) {
		t.Error("sqrt of a negative ball succeeded")
	}
	if zero.sign() != 0 || one.sign() != 1 {
		t.Error("sign of a ball containing 0 is not 0, or of 1 not 1")
	}
	if _, ok := zero.expLower(); ok {
		t.Error("expLower of a ball containing 0 succeeded")
	}
	if logBallB(&z, &zero, 64) {
		t.Error("log of a ball containing 0 succeeded")
	}
}

// TestFixedFailures exercises the fixed-point refusals and bounds.
func TestFixedFailures(t *testing.T) {
	var zero, tiny, one, z fx
	zero.setInt(0, 64)
	one.setInt(1, 64)
	if zero.lowbnd() != 0 {
		t.Error("lowbnd of 0 is not 0")
	}
	tiny.setInt(0, 64)
	tiny.v.SetInt64(3)
	tiny.e = 100 // [3-100, 3+100] ulps contains 0
	if tiny.lowbnd() != 0 {
		t.Error("lowbnd of a ball containing 0 is not 0")
	}
	if z.quo(&one, &tiny) {
		t.Error("fx quo by a ball containing 0 succeeded")
	}
	if z.sqrt(&tiny) {
		t.Error("fx sqrt of a small ball succeeded")
	}
	// radii far beyond a fixed-point ulp, and far below it
	var b ball
	b.setInt64(1)
	b.r = magPow2(2000)
	if f := fxFromBall(&b, 64); !math.IsInf(f.e, 1) {
		t.Errorf("huge radius: e = %v, want +Inf", f.e)
	}
	b.r = magPow2(-3000)
	if f := fxFromBall(&b, 64); f.e <= 1 || f.e > 2 {
		t.Errorf("tiny radius: e = %v, want just above 1", f.e)
	}
	f := fxFromBall(one.toBallOf(), 64)
	f.e = math.Inf(1)
	var back ball
	f.toBall(&back)
	if back.r.e < 1<<39 {
		t.Errorf("infinite fixed-point error gave radius %v", back.r)
	}
}

func (a *fx) toBallOf() *ball {
	var b ball
	a.toBall(&b)
	return &b
}

// TestMagEdges exercises the zero and out-of-range paths of mag.
func TestMagEdges(t *testing.T) {
	if !magNorm(0, 5).isZero() {
		t.Error("magNorm(0) is not zero")
	}
	if !magSubLower(magPow2(0), magPow2(1)).isZero() {
		t.Error("1 - 2 is not bounded below by 0")
	}
	if !magSubLower(magNorm(0.75, 0), magNorm(0.75, 0)).isZero() {
		t.Error("a - a is not bounded below by 0")
	}
	var z big.Float
	if (mag{}).float(&z).Sign() != 0 {
		t.Error("zero mag is not 0")
	}
	if !magPow2(1 << 40).float(&z).IsInf() {
		t.Error("huge mag is not +Inf")
	}
	if magPow2(-1<<40).float(&z).Sign() <= 0 {
		t.Error("tiny mag is not clamped to the smallest positive value")
	}
	huge := new(big.Float).SetMantExp(big.NewFloat(1), 5000)
	if m := magFloat(huge); m.e != 5002 { // 2^5001 = 0.5·2^5002 bounds 2^5000
		t.Errorf("magFloat(2^5000) = %v", m)
	}
	if m := magFloatLower(huge); m.e > 5001 {
		t.Errorf("magFloatLower(2^5000) = %v", m)
	}
	if !magFloatLower(new(big.Float)).isZero() || !magFloat(new(big.Float)).isZero() {
		t.Error("magFloat(0) is not zero")
	}
}

// TestRootStart: iroot recovers from a starting point below the root.
func TestRootStart(t *testing.T) {
	saved := rootStart
	defer func() { rootStart = saved }()
	rootStart = func(*big.Int, uint) *big.Int { return big.NewInt(1) }
	N := new(big.Int).Exp(big.NewInt(12345), big.NewInt(7), nil)
	if r, exact := iroot(N, 7); !exact || r.Int64() != 12345 {
		t.Errorf("iroot(12345^7, 7) = %v, %v", r, exact)
	}
	rootStart = saved
	N.Add(N, big.NewInt(1))
	if r, exact := iroot(N, 7); exact || r.Int64() != 12345 {
		t.Errorf("iroot(12345^7+1, 7) = %v, %v", r, exact)
	}
}

// TestExactRoot covers the exactness test of Root's large-n path.
func TestExactRoot(t *testing.T) {
	n := new(big.Int).Exp(big.NewInt(3), big.NewInt(1000), nil)
	x := new(big.Float).SetInt(n)
	if v := exactRoot(x, 1000); v == nil || v.Cmp(big.NewFloat(3)) != 0 {
		t.Errorf("exactRoot(3^1000, 1000) = %v, want 3", v)
	}
	// 3^1000 + 2: odd, 1585 bits (enough to be a 1000th power), and not one
	m := new(big.Int).Exp(big.NewInt(3), big.NewInt(1000), nil)
	m.Add(m, big.NewInt(2))
	x = new(big.Float).SetInt(m)
	if v := exactRoot(x, 1000); v != nil {
		t.Errorf("exactRoot(3^1000+2, 1000) = %v, want nil", v)
	}
	if v := exactRoot(big.NewFloat(3), 1000); v != nil {
		t.Errorf("exactRoot(3, 1000) = %v, want nil", v)
	}
	if v := exactRoot(big.NewFloat(2), 1000); v != nil {
		t.Errorf("exactRoot(2, 1000) = %v, want nil", v)
	}
}
