// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math"
	"math/bits"
	"sync"
)

// tableMaxPrec is the largest fixed-point scale served from a value table;
// above it the functions fall back to their table-free reductions. It
// bounds the memory of every table: one entry per grid point and per
// power-of-two precision class up to this scale (mpmath stops its log and
// atan tables at 2500 and 3000 bits, libelefun.py).
const tableMaxPrec = 1 << 13

// tableClass returns the scale at which table entries serving scale w are
// computed: the next power of two at or above w, plus guard bits, as
// mpmath's cache_prec_steps does, so that a table is computed at most once
// per doubling of the precision.
func tableClass(w uint) uint { return 1<<bits.Len(w-1) + 32 }

// fxTable caches f(n) for the grid points n of a function's table reduction
// (atan(n/2^7), log(n/2^9)), as fixed-point balls at the class scale,
// computed on first use. Entries are never modified once stored.
type fxTable struct {
	mu      sync.Mutex
	entries map[[2]uint64]*fx // by (class scale, n)
	scaled  map[[2]uint64]*fx // by (scale, n), at most maxScaled
	compute func(n int64, W uint) *fx
}

// maxScaled bounds the truncated copies an fxTable keeps; a new one beyond
// it starts that map afresh (the class entries stay). A variable for tests.
var maxScaled = 1 << 12

// get returns f(n) at scale w <= tableMaxPrec: the class entry truncated to
// w bits (one more ulp of error). The value is shared: an operand only,
// never a destination (nor the receiver of log2bnd or expLower, which
// write scratch).
func (c *fxTable) get(n int64, w uint) *fx {
	key := [2]uint64{uint64(w), uint64(n)}
	c.mu.Lock()
	defer c.mu.Unlock()
	if z := c.scaled[key]; z != nil {
		return z
	}
	W := tableClass(w)
	ckey := [2]uint64{uint64(W), uint64(n)}
	e := c.entries[ckey]
	if e == nil {
		if c.entries == nil {
			c.entries = map[[2]uint64]*fx{}
		}
		e = c.compute(n, W)
		c.entries[ckey] = e
	}
	z := &fx{w: w}
	z.v.Rsh(&e.v, W-w)
	z.e = addUp(up(math.Ldexp(e.e, -int(W-w))), 1) // up: Ldexp may round in the subnormal range
	if c.scaled == nil || len(c.scaled) >= maxScaled {
		c.scaled = map[[2]uint64]*fx{}
	}
	c.scaled[key] = z
	return z
}
