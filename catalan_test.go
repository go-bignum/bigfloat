// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat_test

import (
	"math/big"
	"testing"

	"github.com/go-bignum/bigfloat"
)

// catalan60 is G to 60 decimals (OEIS A006752), an independent check of the
// fixtures' bits.
const catalan60 = "0.915965594177219015054603514932384110774149374281672134266498"

func TestCatalanDecimals(t *testing.T) {
	// 60 decimals are worth 199 bits: at 150 bits the nearest rounding of
	// the truncated decimal is G's.
	want, _, err := big.ParseFloat(catalan60, 10, 400, big.ToNearestEven)
	if err != nil {
		t.Fatal(err)
	}
	z := bigfloat.Catalan(new(big.Float).SetPrec(150))
	if want.SetPrec(150); z.Cmp(want) != 0 {
		t.Errorf("Catalan at 150 bits = %v, want %v", z, want)
	}
	if z.Acc() == big.Exact {
		t.Error("G is not representable: accuracy must be Below or Above")
	}
}

func TestCatalanDefaultPrecision(t *testing.T) {
	z := bigfloat.Catalan(new(big.Float))
	if z.Prec() != 64 || z.Text('p', 0) != "0x.ea7cb89f409ae845p+0" {
		t.Errorf("default: prec %d, %s", z.Prec(), z.Text('p', 0))
	}
	y := bigfloat.Catalan(new(big.Float).SetPrec(10).SetMode(big.ToZero))
	if y.Prec() != 10 || y.Text('p', 0) != "0x.ea4p+0" || y.Acc() != big.Below {
		t.Errorf("10 bits toward zero: prec %d, %s, %v", y.Prec(), y.Text('p', 0), y.Acc())
	}
}
