// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat_test

import (
	"errors"
	"maps"
	"math/big"
	"slices"
	"testing"

	"github.com/go-bignum/bigfloat"
	"github.com/go-bignum/bigfloat/oracle"
)

var fuzzNames = slices.Sorted(maps.Keys(oracle.Functions))

// fuzzFloat builds a number from fuzz bytes: up to 32 bytes of mantissa, a
// bounded exponent, a sign.
func fuzzFloat(mant []byte, exp int16, neg bool) *big.Float {
	if len(mant) > 32 {
		mant = mant[:32]
	}
	m := new(big.Int).SetBytes(mant)
	x := new(big.Float).SetPrec(uint(max(m.BitLen(), 1))).SetInt(m)
	x.SetMantExp(x, int(exp)%4096)
	if neg {
		x.Neg(x)
	}
	return x
}

// FuzzFunctions checks, for arbitrary arguments, precisions and modes,
// that every function returns what the oracle Rounder derives from the
// same package evaluated at more than twice the precision (the rounding
// is then settled by the Rounder's own interval test), that the accuracy
// it reports agrees, and that the only panics are ErrNaN (domain errors)
// and ErrLimit (resource bounds).
func FuzzFunctions(f *testing.F) {
	f.Add(uint8(0), []byte{3}, int16(0), false, []byte{1}, int16(0), false, uint16(53), uint8(0))
	f.Add(uint8(7), []byte{0xff, 0x13}, int16(-40), true, []byte{5}, int16(3), false, uint16(24), uint8(4))
	f.Add(uint8(20), []byte{0x80, 0, 0, 1}, int16(1000), false, []byte{9}, int16(-2), true, uint16(113), uint8(5))
	f.Fuzz(func(t *testing.T, fi uint8, xm []byte, xe int16, xn bool, ym []byte, ye int16, yn bool, prec uint16, mode uint8) {
		fn := fuzzNames[int(fi)%len(fuzzNames)]
		p := uint(prec)%300 + 2
		m := oracle.Modes[int(mode)%len(oracle.Modes)]
		args := []*big.Float{fuzzFloat(xm, xe, xn), fuzzFloat(ym, ye, yn)}[:oracle.Functions[fn]]
		if fn == "root" {
			k, _ := args[1].Abs(args[1]).Uint64()
			args[1].SetUint64(k%50 + 1)
		}
		z := new(big.Float).SetPrec(p).SetMode(m)
		if !evalOK(t, func() { oracle.Call(fn, z, args...) }) {
			return
		}
		if z.IsInf() || z.Sign() == 0 {
			return // special values and over/underflow: no finite rounding to compare
		}
		res, err := (&oracle.Rounder{Oracle: oracle.Self{}, MaxPrec: 32 * (p + 300)}).Round(fn, p, args...)
		if errors.Is(err, oracle.ErrUndecided) {
			// within 2^-MaxPrec of a breakpoint (expm1 of a huge negative
			// number is -1 + e^x): this property cannot judge it
			t.Skip("undecidable by the Rounder")
		}
		if err != nil {
			t.Fatalf("%s%v at %d bits: oracle: %v", fn, args, p, err)
		}
		i := slices.Index(oracle.Modes, m)
		if z.Cmp(res[i].Value) != 0 || z.Acc() != res[i].Acc {
			t.Fatalf("%s%v at %d bits, %v: got %v (%v), want %v (%v)", fn, args, p, m, z, z.Acc(), res[i].Value, res[i].Acc)
		}
	})
}

// evalOK runs fn and reports false when it panicked with a documented
// panic value; any other panic fails the test.
func evalOK(t *testing.T, fn func()) (ok bool) {
	defer func() {
		switch e := recover().(type) {
		case nil:
		case bigfloat.ErrNaN, bigfloat.ErrLimit:
			ok = false
		default:
			t.Fatalf("unexpected panic: %v", e)
		}
	}()
	fn()
	return true
}
