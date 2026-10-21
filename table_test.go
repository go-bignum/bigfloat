// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math/big"
	"testing"
)

// TestTableBounds: the memoised truncations of a table are bounded (the
// map starts afresh past maxScaled) and a recomputed entry is the same
// ball; an entry encloses the function at its grid point.
func TestTableBounds(t *testing.T) {
	saved := maxScaled
	defer func() { maxScaled = saved }()
	maxScaled = 2
	first := atanTable.get(5, 100)
	for _, w := range []uint{101, 102, 103} {
		atanTable.get(5, w)
	}
	if n := len(atanTable.scaled); n > 2 {
		t.Errorf("scaled copies: %d, want at most 2", n)
	}
	again := atanTable.get(5, 100)
	if first == again || first.v.Cmp(&again.v) != 0 || first.e != again.e {
		t.Errorf("recomputed entry differs or was not recomputed: %v±%v vs %v±%v", &first.v, first.e, &again.v, again.e)
	}
	// atan(5/128) against the package's own Atan at the same scale.
	var b ball
	again.toBall(&b)
	want := Atan(new(big.Float).SetPrec(90), big.NewFloat(5.0/128))
	var d, rad big.Float
	d.SetPrec(200).Sub(&b.m, want)
	b.r.float(&rad)
	rad.Add(&rad, new(big.Float).SetMantExp(big.NewFloat(1), -90))
	if d.Abs(&d).Cmp(&rad) > 0 {
		t.Errorf("atan table entry %v ± %v does not enclose atan(5/128) = %v", &b.m, &rad, want)
	}
}

// TestFxPools: workspaces come back with the requested length, from the
// pools when they fit and freshly allocated beyond the largest pool.
func TestFxPools(t *testing.T) {
	for _, n := range []int{1, 2, 3, 1 << 11, 1<<11 + 1, 5000} {
		p := getFxs(n)
		if len(*p) != n {
			t.Errorf("getFxs(%d) has length %d", n, len(*p))
		}
		putFxs(p)
	}
}

// TestConstFxBound: the fixed-point copies of a constant are bounded and
// the uncached huge scale still works.
func TestConstFxBound(t *testing.T) {
	for w := uint(64); w < 64+maxConstFx+4; w++ {
		ln2Fx(w)
	}
	if n := len(ln2Cache.fxs); n > maxConstFx {
		t.Errorf("constant copies: %d, want at most %d", n, maxConstFx)
	}
	if f := piFx(tableMaxPrec + 1); f.v.BitLen() != tableMaxPrec+3 { // π·2^w: w+2 bits
		t.Errorf("π at scale %d has %d bits", tableMaxPrec+1, f.v.BitLen())
	}
}
