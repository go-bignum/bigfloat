// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat_test

import (
	"fmt"
	"math/big"

	"github.com/go-bignum/bigfloat"
)

// sin(10^22) at 53 bits, rounded down and up: the two adjacent doubles
// around sin(10^22) = -0.852200849767188801772705893753…, whatever the
// size of the argument (the reduction modulo π/2 is exact).
func ExampleSin() {
	x := new(big.Float).SetFloat64(1e22)
	lo := bigfloat.Sin(new(big.Float).SetPrec(53).SetMode(big.ToNegativeInf), x)
	hi := bigfloat.Sin(new(big.Float).SetPrec(53).SetMode(big.ToPositiveInf), x)
	fmt.Println(lo.Text('g', 20), lo.Acc())
	fmt.Println(hi.Text('g', 20), hi.Acc())
	// Output:
	// -0.85220084976718890601 Below
	// -0.85220084976718879499 Above
}

// Pi to 200 bits.
func ExamplePi() {
	fmt.Println(bigfloat.Pi(new(big.Float).SetPrec(200)).Text('g', 60))
	// Output: 3.14159265358979323846264338327950288419716939937510582097494
}

// An exact case is recognised and rounded exactly: 4^1.5 = 8.
func ExamplePow() {
	z := bigfloat.Pow(new(big.Float).SetPrec(24), big.NewFloat(4), big.NewFloat(1.5))
	fmt.Println(z, z.Acc())
	// Output: 8 Exact
}
