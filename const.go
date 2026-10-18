// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"math/big"
	"sync"
)

// bsResult holds the products of a binary-splitting range [n1, n2) of a
// series S = Σ a(k)/b(k) · Π_{j<=k} p(j)/q(j)  (Haible & Papanikolaou,
// "Fast multiprecision evaluation of series of rational numbers", ANTS-III,
// LNCS 1423, 1998; Brent & Zimmermann, MCA §4.9):
// P = Π p, Q = Π q, B = Π b, T = B·Q·S over the range.
type bsResult struct{ P, Q, B, T *big.Int }

type bsSeries interface {
	term(k int64) (a, b, p, q *big.Int)
}

func binarySplit(s bsSeries, n1, n2 int64) bsResult {
	if n2-n1 == 1 {
		a, b, p, q := s.term(n1)
		return bsResult{p, q, b, new(big.Int).Mul(a, p)}
	}
	m := (n1 + n2) / 2
	l := binarySplit(s, n1, m)
	r := binarySplit(s, m, n2)
	var x, y big.Int
	// T = Br·Qr·Tl + Bl·Pl·Tr
	x.Mul(r.B, r.Q)
	x.Mul(&x, l.T)
	y.Mul(l.B, l.P)
	y.Mul(&y, r.T)
	x.Add(&x, &y)
	return bsResult{
		P: l.P.Mul(l.P, r.P),
		Q: l.Q.Mul(l.Q, r.Q),
		B: l.B.Mul(l.B, r.B),
		T: &x,
	}
}

// bsBall returns a ball enclosing S = T/(B·Q) of a split series widened by
// the relative truncation bound 2^-relBits.
func bsBall(r bsResult, wp uint, relBits int64) *ball {
	var num, den ball
	num.m.SetPrec(0).SetInt(r.T)
	d := new(big.Int).Mul(r.B, r.Q)
	den.m.SetPrec(0).SetInt(d)
	z := new(ball)
	z.quo(&num, &den, wp)
	z.addRad(magMul(z.upper(), magPow2(-relBits)))
	return z
}

// log2FactLower returns Σ_{j=2}^{k} floor(log2 j) <= log2(k!).
func log2FactLower(k int64) int64 {
	var s int64
	for j := int64(2); j <= k; j++ {
		s += int64(big.NewInt(j).BitLen() - 1)
	}
	return s
}

// --- e: Σ 1/k! -----------------------------------------------------------

type eSeries struct{}

func (eSeries) term(k int64) (a, b, p, q *big.Int) {
	return big.NewInt(1), big.NewInt(1), big.NewInt(1), big.NewInt(max(k, 1))
}

func eBallCompute(wp uint) *ball {
	// Σ_{k>=N} 1/k! <= 2/N!; relative to e > 2 that is below 1/N!.
	N := int64(2)
	for log2FactLower(N) < int64(wp)+4 {
		N += 8
	}
	return bsBall(binarySplit(eSeries{}, 0, N), wp, log2FactLower(N))
}

// --- atanh(1/x) = Σ 1/((2k+1) x^(2k+1)) ----------------------------------

type atanhInvSeries struct{ x int64 }

func (s atanhInvSeries) term(k int64) (a, b, p, q *big.Int) {
	if k == 0 {
		return big.NewInt(1), big.NewInt(1), big.NewInt(1), big.NewInt(s.x)
	}
	return big.NewInt(1), big.NewInt(2*k + 1), big.NewInt(1), big.NewInt(s.x * s.x)
}

// atanhInvBall encloses atanh(1/x), x >= 2. The tail after N terms is at
// most x^-(2N+1)/(1-x^-2) <= 2·x^-(2N+1), below 2·x^-2N relative to the sum
// (which exceeds 1/x).
func atanhInvBall(x int64, wp uint) *ball {
	lx := int64(big.NewInt(x).BitLen() - 1) // 2^lx <= x
	N := (int64(wp)+8)/(2*lx) + 1
	return bsBall(binarySplit(atanhInvSeries{x}, 0, N), wp, 2*N*lx-1)
}

// --- π: Chudnovsky -------------------------------------------------------

type chudnovsky struct{}

var chudC3over24 = big.NewInt(10939058860032000) // 640320³/24

func (chudnovsky) term(k int64) (a, b, p, q *big.Int) {
	a = big.NewInt(545140134)
	a.Mul(a, big.NewInt(k))
	a.Add(a, big.NewInt(13591409))
	if k == 0 {
		return a, big.NewInt(1), big.NewInt(1), big.NewInt(1)
	}
	p = big.NewInt(6*k - 5)
	p.Mul(p, big.NewInt(2*k-1))
	p.Mul(p, big.NewInt(6*k-1))
	p.Neg(p)
	q = big.NewInt(k)
	q.Mul(q, q)
	q.Mul(q, big.NewInt(k))
	q.Mul(q, chudC3over24)
	return a, big.NewInt(1), p, q
}

// piBallCompute: π = 426880·sqrt(10005)/S with S the Chudnovsky sum
// (Chudnovsky & Chudnovsky 1988; binary splitting as in MCA §4.9.1). Every
// term ratio is below 2^-45 in absolute value (|p(k)/q(k)| <= 1728/640320³
// and the linear factor grows by at most 41.1 at k = 1, 1.98 after), so
// the tail after N terms is below 2^(2-45N) relative to S.
func piBallCompute(wp uint) *ball {
	N := (int64(wp)+8)/45 + 1
	s := bsBall(binarySplit(chudnovsky{}, 0, N), wp+8, 45*N-2)
	var c, r, z ball
	c.setInt64(10005)
	r.sqrt(&c, wp+8)
	c.setInt64(426880)
	r.mul(&r, &c, wp+8)
	z.quo(&r, s, wp)
	return &z
}

// ln2 = 18 atanh(1/26) - 2 atanh(1/4801) + 8 atanh(1/8749).
func ln2BallCompute(wp uint) *ball {
	w := wp + 10
	a := atanhInvBall(26, w)
	b := atanhInvBall(4801, w)
	c := atanhInvBall(8749, w)
	var k, t, z ball
	k.setInt64(18)
	z.mul(a, &k, w)
	k.setInt64(2)
	t.mul(b, &k, w)
	z.sub(&z, &t, w)
	k.setInt64(8)
	t.mul(c, &k, w)
	z.add(&z, &t, wp)
	return &z
}

// ln10 = 3 ln2 + 2 atanh(1/9)  (log(5/4) = 2 atanh(1/9)).
func ln10BallCompute(wp uint) *ball {
	w := wp + 10
	l := ln2Ball(w)
	a := atanhInvBall(9, w)
	var k, z, t ball
	k.setInt64(3)
	z.mul(l, &k, w)
	t.mul2exp(a, 1)
	z.add(&z, &t, wp)
	return &z
}

// eulerBallCompute uses Brent & McMillan's algorithm B1 (Math. Comp. 34,
// 1980; MCA §4.12): with B_0 = 1, A_0 = -log n, B_k = B_{k-1}n²/k²,
// A_k = (A_{k-1}n²/k + B_k)/k, U = ΣA_k, V = ΣB_k, one has
// U/V - γ = K0(2n)/I0(2n), between 0 and π·e^(-4n). The sums are truncated
// once K >= 2n (successive B then shrink by 4 at least) and B_K is
// negligible: the tail of V is below B_K and that of U below
// 2(K+1)·B_K (|A_j| <= B_j(H_j + log n) <= 2j·B_j), both bounded from the
// computed balls.
func eulerBallCompute(wp uint) *ball {
	w := wp + 32
	// 2π·e^(-4n) <= 2^-(wp+6)  <=  4n·log2(e) >= wp+10
	n := int64(wp+10)*1000/5770 + 1
	var A, B, U, V, t, nn, kb ball
	logInt(&A, n, w)
	A.neg(&A)
	B.setInt64(1)
	U.set(&A)
	V.set(&B)
	nn.setInt64(n * n)
	var k int64
	for k = 1; ; k++ {
		kb.setInt64(k)
		t.setInt64(k * k)
		B.mul(&B, &nn, w)
		B.quo(&B, &t, w)
		A.mul(&A, &nn, w)
		A.quo(&A, &kb, w)
		A.add(&A, &B, w)
		A.quo(&A, &kb, w)
		U.add(&U, &A, w)
		V.add(&V, &B, w)
		if k >= 2*n && B.expUpper() < -int64(w)-int64(big.NewInt(k).BitLen())-8 {
			break
		}
	}
	tail := magMul(B.upper(), magUint(uint64(4*(k+1))))
	U.addRad(tail)
	V.addRad(tail)
	var g ball
	g.quo(&U, &V, wp+8)
	// γ is within π·e^(-4n) of U/V (on one side; taken symmetric here).
	return g.addRad(magPow2(-int64(wp) - 6))
}

// logInt encloses log n for an integer n >= 1.
func logInt(z *ball, n int64, wp uint) {
	logBall(z, new(big.Float).SetInt64(n), wp)
}

// constCache memoises a constant at the highest precision asked so far;
// lower precisions are rounded from it (with the rounding added to the
// radius).
type constCache struct {
	mu      sync.Mutex
	prec    uint
	b       *ball
	compute func(uint) *ball
}

func (c *constCache) get(wp uint) *ball {
	c.mu.Lock()
	if c.b == nil || c.prec < wp {
		p := max(wp, c.prec+c.prec/2)
		c.b = c.compute(p)
		c.prec = p
	}
	src := c.b
	c.mu.Unlock()
	z := new(ball)
	z.setFloat(&src.m, wp)
	z.addRad(src.r)
	return z
}

var (
	piCache    = &constCache{compute: piBallCompute}
	eCache     = &constCache{compute: eBallCompute}
	ln2Cache   = &constCache{compute: ln2BallCompute}
	ln10Cache  = &constCache{compute: ln10BallCompute}
	eulerCache = &constCache{compute: eulerBallCompute}
)

func piBall(wp uint) *ball    { return piCache.get(wp) }
func ln2Ball(wp uint) *ball   { return ln2Cache.get(wp) }
func ln10Ball(wp uint) *ball  { return ln10Cache.get(wp) }
func eBall(wp uint) *ball     { return eCache.get(wp) }
func eulerBall(wp uint) *ball { return eulerCache.get(wp) }

func constant(z *big.Float, f func(uint) *ball) *big.Float {
	precOf(z, 64)
	return ziv(z, 0, func(wp uint) (*ball, int64, bool) { return f(wp), 0, true })
}

// Pi sets z to π correctly rounded at z's precision and mode (64 bits when
// z has none) and returns z.
func Pi(z *big.Float) *big.Float { return constant(z, piBall) }

// E sets z to e = exp(1), correctly rounded.
func E(z *big.Float) *big.Float { return constant(z, eBall) }

// Ln2 sets z to log 2, correctly rounded.
func Ln2(z *big.Float) *big.Float { return constant(z, ln2Ball) }

// Ln10 sets z to log 10, correctly rounded.
func Ln10(z *big.Float) *big.Float { return constant(z, ln10Ball) }

// EulerGamma sets z to the Euler–Mascheroni constant γ, correctly rounded.
func EulerGamma(z *big.Float) *big.Float { return constant(z, eulerBall) }
