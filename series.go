// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math"
	"sync"
)

// oddSeries returns Σ_{k=0}^{n} s_k·y^(2k+1)/(2k+1) for a fixed-point ball
// y, with s_k = (-1)^k when alt (the atan series) and 1 otherwise (atanh).
// The truncation tail is the caller's to add.
//
// The terms are spread over m accumulators (m even), acc_i collecting the
// terms k ≡ i (mod m) without their factor u^i = y^(2i): a term costs one
// division by a small integer, and the power advances by u^m once every m
// terms; the accumulators are then recombined with the powers u^i. That is
// mpmath's two-sum loop of log_taylor and atan_taylor (libelefun.py)
// generalised to m sums, the baby-step giant-step of Paterson & Stockmeyer
// (SIAM J. Comput. 2(1), 1973) applied to the multiplications: about
// n/m + 2m multiplications instead of n, with m near sqrt(n). As m is even,
// the sign of term k = jm + i is (-1)^i, which the recombination applies.
// The fixed-point operations carry every rounding error into the result's
// error bound.
func oddSeries(y *fx, n int64, alt bool) *fx {
	w := y.w
	m := max(2, int(math.Sqrt(float64(n+1)))&^1)
	// pw[1..m] = u^1..u^m, acc[0..m-1], then P and q, from a pooled
	// workspace whose integers keep their capacity between calls.
	ws := getFxs(2*m + 3)
	defer putFxs(ws)
	buf := *ws
	pw, acc := buf[:m+1], buf[m+1:2*m+1]
	P, q := &buf[2*m+1], &buf[2*m+2]
	pw[1].mul(y, y)
	for i := 2; i <= m; i++ {
		pw[i].mul(&pw[i-1], &pw[1])
	}
	for i := range acc {
		acc[i].setInt(0, w)
	}
	P.set(y)
	for k := int64(0); k <= n; {
		for i := 0; i < m && k <= n; i, k = i+1, k+1 {
			q.divSmall(P, uint64(2*k+1))
			acc[i].add(&acc[i], q)
		}
		if k <= n {
			P.mul(P, &pw[m])
		}
	}
	r := &acc[0]
	for i := 1; i < m; i++ {
		q.mul(&acc[i], &pw[i])
		if alt && i%2 == 1 {
			r.sub(r, q)
		} else {
			r.add(r, q)
		}
	}
	return new(fx).set(r)
}

// fxPool recycles fixed-point workspaces (slices of fx whose integers have
// grown to the size of a series' operands), so that a series costs no
// allocation once the pool is warm.
var fxPool sync.Pool

// getFxs returns a workspace of at least n fixed-point values.
func getFxs(n int) *[]fx {
	if p, _ := fxPool.Get().(*[]fx); p != nil && cap(*p) >= n {
		*p = (*p)[:n]
		return p
	}
	s := make([]fx, n)
	return &s
}

func putFxs(p *[]fx) { fxPool.Put(p) }
