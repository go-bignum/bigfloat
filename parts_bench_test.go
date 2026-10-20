// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math/big"
	"testing"
)

// BenchmarkParts times the stages of a 53-bit Sin call separately (the
// macOS CPU profiler attributes most samples to runtime.kevent).
func BenchmarkParts(b *testing.B) {
	x := big.NewFloat(0.7)
	const wp = 77
	b.Run("reduce", func(b *testing.B) {
		for b.Loop() {
			reduceHalfPi(x, 111)
		}
	})
	r, _ := reduceHalfPi(x, 103)
	s, c := new(fx), new(fx)
	b.Run("sincosFx", func(b *testing.B) {
		for b.Loop() {
			sincosFx(s, c, r, 3)
		}
	})
	b.Run("sincosCore", func(b *testing.B) {
		for b.Loop() {
			sincosCore(s, c, r, 3)
		}
	})
	sincosFx(s, c, r, 3)
	var sb ball
	b.Run("toBall", func(b *testing.B) {
		for b.Loop() {
			s.toBall(&sb)
		}
	})
	b.Run("sincosBall", func(b *testing.B) {
		for b.Loop() {
			sincosBall(x, wp)
		}
	})
	ws := new(zivWS)
	b.Run("decide", func(b *testing.B) {
		for b.Loop() {
			decide(&sb, 53, big.ToNearestEven, ws)
		}
	})
	z := new(big.Float).SetPrec(53)
	b.Run("Sin", func(b *testing.B) {
		for b.Loop() {
			Sin(z, x)
		}
	})
	var p, q fx
	q.set(r)
	b.Run("fxmul", func(b *testing.B) {
		for b.Loop() {
			p.mul(&q, &q)
		}
	})
	b.Run("fxdiv", func(b *testing.B) {
		for b.Loop() {
			p.divSmall(&q, 12345)
		}
	})
}
