// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat_test

import (
	"fmt"
	"math/big"
	"testing"

	"github.com/go-bignum/bigfloat/oracle"
)

// BenchmarkFunctions times the calls of cmd/bench (the docs/perf.md table)
// without mpmath, for profiling: go test -bench Functions/atan -cpuprofile.
func BenchmarkFunctions(b *testing.B) {
	cases := []struct {
		fn   string
		args []float64
	}{
		{"sqrt", []float64{2}}, {"cbrt", []float64{2}}, {"root", []float64{2, 7}}, {"hypot", []float64{0.75, 1.25}},
		{"exp", []float64{0.7}}, {"expm1", []float64{0.7}}, {"log", []float64{1.7}}, {"log1p", []float64{0.7}},
		{"log2", []float64{1.7}}, {"log10", []float64{1.7}}, {"pow", []float64{1.7, 2.3}},
		{"sin", []float64{0.7}}, {"cos", []float64{0.7}}, {"tan", []float64{0.7}},
		{"asin", []float64{0.7}}, {"acos", []float64{0.7}}, {"atan", []float64{0.7}}, {"atan2", []float64{0.7, 1.3}},
		{"sinh", []float64{0.7}}, {"cosh", []float64{0.7}}, {"tanh", []float64{0.7}},
		{"asinh", []float64{0.7}}, {"acosh", []float64{1.7}}, {"atanh", []float64{0.7}},
		{"sin", []float64{1e22}}, {"exp", []float64{1000.5}},
	}
	for _, c := range cases {
		args := make([]*big.Float, len(c.args))
		for i, a := range c.args {
			args[i] = new(big.Float).SetFloat64(a)
		}
		for _, p := range []uint{53, 113, 200, 1000} {
			b.Run(fmt.Sprintf("%s/%g/%d", c.fn, c.args[0], p), func(b *testing.B) {
				z := new(big.Float)
				for b.Loop() {
					oracle.Call(c.fn, z.SetPrec(p), args...)
				}
			})
		}
	}
}
