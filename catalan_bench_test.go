// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"fmt"
	"math/big"
	"testing"
)

// BenchmarkCatalan times G from scratch, as MPFR's mpfr_const_catalan is
// timed after mpfr_free_cache: the caches of G and of π are emptied on every
// iteration.
func BenchmarkCatalan(b *testing.B) {
	for _, p := range []uint{53, 200, 1000, 10000, 100000} {
		b.Run(fmt.Sprint(p), func(b *testing.B) {
			z := new(big.Float)
			for b.Loop() {
				piCache = &constCache{compute: piBallCompute}
				catalanCache = &constCache{compute: catalanBallCompute}
				Catalan(z.SetPrec(p))
			}
		})
	}
}
