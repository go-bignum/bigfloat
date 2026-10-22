// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"fmt"
	"math"
	"math/big"
)

// trigSplit returns the number of argument halvings before the sin/cos and
// atan series.
func trigSplit(w uint) uint { return uint(math.Sqrt(float64(w))/3) + 1 }

// sincosCore returns sin r and cos r for a fixed-point ball with |r| < 1:
// both Taylor series at t = r/2^j (alternating, decreasing terms: the tail
// is below the first omitted term), then j doublings sin 2a = 2 sin a cos a,
// cos 2a = 1 - 2 sin² a. The results go to s and c (distinct from r and
// from each other).
func sincosCore(s, c, r *fx, j uint) {
	w := r.w
	ws := getFxs(3)
	defer putFxs(ws)
	t, one, s2 := &(*ws)[0], &(*ws)[1], &(*ws)[2]
	t.shr(r, j)
	lt := min(int64(t.log2bnd()), -1)
	// cos tail <= |t|^(2n+2)/(2n+2)!, sin tail smaller.
	n := int64(1)
	for lt*(2*n+2)-log2FactLower(2*n+2) > -int64(w)-1 {
		n++
	}
	sinCosSeries(s, c, t, n)
	s.e = addUp(s.e, 1)
	c.e = addUp(c.e, 1)
	one.setInt(1, w)
	for i := uint(0); i < j; i++ {
		s2.mul(s, s)
		s.mul(s, c)
		s.shl(s, 1)
		c.shl(s2, 1)
		c.sub(one, c)
	}
}

// sinCosShift sets the grid of the sin/cos tables, n/2^8 (mpmath's
// COS_SIN_CACHE_STEP), and sinCosTableMax the largest scale they serve
// (mpmath's COS_SIN_CACHE_PREC is 400 bits): above it the halvings of
// sincosCore cost about as much as the table saves.
const (
	sinCosShift    = 8
	sinCosTableMax = 640
)

// sinTable and cosTable cache sin(n/2^8) and cos(n/2^8) for 1 <= n <= 202
// (|r| <= π/4 + ε), computed by sincosCore.
var (
	sinTable = &fxTable{compute: func(n int64, W uint) *fx { s, _ := sinCosAt(n, W); return s }}
	cosTable = &fxTable{compute: func(n int64, W uint) *fx { _, c := sinCosAt(n, W); return c }}
)

func sinCosAt(n int64, W uint) (*fx, *fx) {
	var a fx
	a.setInt(n, W)
	a.v.Rsh(&a.v, sinCosShift) // exact: n/2^8 at scale W > 8
	s, c := new(fx), new(fx)
	sincosCore(s, c, &a, trigSplit(W))
	return s, c
}

// sincosFx returns sin r and cos r for a fixed-point ball with |r| < 1. Up
// to sinCosTableMax it follows mpmath's cos_sin_basecase (libelefun.py):
// with a = n/2^8 the grid point below |r|, d = |r| - a in [0, 2^-8), the
// series of d need no halving, and sin(a+d) = sin a cos d + cos a sin d,
// cos(a+d) = cos a cos d - sin a sin d with sin a, cos a from the tables.
// Otherwise, or for |r| < 2^-8, sincosCore with j halvings. The results go
// to s and c, as for sincosCore.
func sincosFx(s, c, r *fx, j uint) {
	w := r.w
	if w > sinCosTableMax {
		sincosCore(s, c, r, j)
		return
	}
	neg := r.v.Sign() < 0
	ws := getFxs(5)
	defer putFxs(ws)
	ar, sd, cd, p, q := &(*ws)[0], &(*ws)[1], &(*ws)[2], &(*ws)[3], &(*ws)[4]
	ar.set(r)
	if neg {
		ar.v.Neg(&ar.v)
	}
	p.v.Rsh(&ar.v, w-sinCosShift)
	n := p.v.Int64()
	if n == 0 {
		sincosCore(s, c, r, 0)
		return
	}
	p.setInt(n, w)
	p.v.Rsh(&p.v, sinCosShift) // a = n/2^8, exact
	ar.sub(ar, p)
	sincosCore(sd, cd, ar, 0)
	sa, ca := sinTable.get(n, w), cosTable.get(n, w)
	p.mul(sa, cd)
	q.mul(ca, sd)
	s.add(p, q)
	p.mul(ca, cd)
	q.mul(sa, sd)
	c.sub(p, q)
	if neg {
		s.v.Neg(&s.v)
	}
}

// sinCosSeries returns the Taylor polynomials of sin t and cos t to the
// terms k = n (t^(2n+1)/(2n+1)! and t^(2n)/(2n)!), tails left to the
// caller. Both run on m accumulators of u = t² as oddSeries does, sharing
// the powers of u: the running cos term a = u^(k-i)/(2k)! takes one
// division, the sin term is a/(2k+1), and a advances by u^m every m terms.
// The results go to s and c (distinct from t).
func sinCosSeries(S, C, t *fx, n int64) {
	w := t.w
	m := max(2, int(math.Sqrt(float64(n)))&^1)
	ws := getFxs(3*m + 4)
	defer putFxs(ws)
	buf := *ws
	pw, ac, as := buf[:m+1], buf[m+1:2*m+1], buf[2*m+1:3*m+1]
	a, b, q := &buf[3*m+1], &buf[3*m+2], &buf[3*m+3]
	pw[1].mul(t, t)
	for i := 2; i <= m; i++ {
		pw[i].mul(&pw[i-1], &pw[1])
	}
	ac[0].setInt(1, w)
	as[0].setInt(1, w)
	for i := 1; i < m; i++ {
		ac[i].setInt(0, w)
		as[i].setInt(0, w)
	}
	a.setInt(1, w)
	for k := int64(1); k <= n; k++ {
		i := int(k % int64(m))
		if i == 0 {
			a.mul(a, &pw[m])
		}
		a.divSmall(a, uint64((2*k-1)*(2*k)))
		ac[i].add(&ac[i], a)
		b.divSmall(a, uint64(2*k+1))
		as[i].add(&as[i], b)
	}
	c, s := &ac[0], &as[0]
	for i := 1; i < m; i++ {
		// m even: the sign of term jm+i is (-1)^i
		q.mul(&ac[i], &pw[i])
		if i%2 == 1 {
			c.sub(c, q)
		} else {
			c.add(c, q)
		}
		q.mul(&as[i], &pw[i])
		if i%2 == 1 {
			s.sub(s, q)
		} else {
			s.add(s, q)
		}
	}
	S.mul(s, t)
	C.set(c)
}

var quarterPi = big.NewFloat(0.78)

// MaxReductionExponent bounds the binary exponent of the arguments of Sin,
// Cos and Tan: an argument near 2^e is reduced modulo π/2 with about e
// bits of π, which costs time and memory growing with e (the first such
// call at 2^22 computes four million bits of π, in seconds, later calls
// reuse the cache). Beyond the bound these functions panic with ErrLimit
// instead of consuming unbounded resources on an untrusted input. Set it
// before use, not concurrently with calls.
var MaxReductionExponent = 1 << 22

// ErrLimit is the panic value of a function asked for work beyond one of
// this package's documented resource bounds (MaxReductionExponent).
type ErrLimit struct {
	Msg string
}

func (e ErrLimit) Error() string { return e.Msg }

// reduceHalfPi writes x = q·π/2 + r with |r| <= π/4 + ε and returns r as
// a ball of absolute accuracy about 2^-w, and q mod 4. π is taken to as
// many bits as x's exponent requires, so that the reduction stays exact
// for any argument (MPFR's mpfr_sin does the same); the first call at a
// huge exponent pays for computing those bits of π, later ones reuse the
// cache.
//
// r is a fixed-point ball at scale w, within two ulps: the reduction runs
// at scale W = w + bitlen(k) + 9, where the |k| ulps of error of k·π/2 are
// below one ulp of 2^-w. r and kl (scratch) are the caller's.
func reduceHalfPi(r, kl *fx, x *big.Float, w uint) int {
	e := x.MantExp(nil)
	xf, _ := x.Float64() // ±Inf beyond float64's range: only used when e <= 40
	if e <= 0 && math.Abs(xf) <= 0.78 {
		// |x| <= 0.78·(1 + 2^-53) < π/4
		r.setFloat(x, w)
		return 0
	}
	if e > MaxReductionExponent {
		panic(ErrLimit{Msg: fmt.Sprintf("bigfloat: argument 2^%d beyond MaxReductionExponent (2^%d): reducing it needs that many bits of π", e, MaxReductionExponent)})
	}
	ex := uint(max(e, 0))
	kint := &kl.k
	if ex <= 40 {
		// float64 quotient: absolute error below 2^-11, far below 1/2
		kint.SetInt64(int64(math.Round(xf / (math.Pi / 2))))
	} else {
		// k = round(x/(π/2)) = floor((2X + P)/(2P)) with P = π·2^V ± e
		// (V = ex + 64; e of a few units) and X = x·2^(V+1) (truncated when
		// x has bits below 2^-(V+1), possible only for a long x): the
		// quotient is off by at most about |k|·2^-(V-2) + 2^-V < 1/4 of x/(π/2),
		// far below the 1/2 that would move |r| past π/4 + ε.
		V := ex + 64
		pf := piFx(V)
		kl.setFloat(x, V+1) // X in kl.v
		kl.v.Lsh(&kl.v, 1)
		kl.v.Add(&kl.v, &pf.v)
		kl.t.Lsh(&pf.v, 1)
		kint.Div(&kl.v, &kl.t) // Euclidean: floor for the positive divisor
	}
	// k mod 4 from the low bits of |k| (big.Int keeps sign and magnitude)
	q := 0
	if kint.Sign() != 0 { // k = 0 for 0.78 < |x| < π/4
		q = int(kint.Bits()[0] & 3)
	}
	if kint.Sign() < 0 {
		q = (4 - q) & 3
	}
	W := w + uint(kint.BitLen()) + 8
	// The integer of π at scale W is π/2 at scale W+1. k·π/2 is formed
	// exactly from it; its error, |k|·e_π < 2^bitlen(k)·e_π ulps at scale
	// W+1, is added after the rescaling to w, where it is below e_π/2^9
	// (|k| may exceed float64's range: |k|·e_π itself is not formed).
	pf := piFx(W)
	kl.v.Mul(&pf.v, kint)
	kl.w, kl.e = W+1, 0
	r.setFloat(x, W+1)
	r.sub(r, kl)
	r.rescale(r, w)
	r.e = addUp(r.e, up(pf.e/512))
	return q
}

// sincosBall returns balls for sin x and cos x (exact x), each with about
// wp bits of relative accuracy. The reduced argument r needs as many more
// bits of absolute accuracy as it has leading zeros; ok = false when the
// reduction at this precision cannot even separate r from zero (x within
// 2^-wp of a multiple of π/2), so that the Ziv loop raises the precision.
// want selects the balls built: 0 sin only, 1 cos only, 2 both (the other
// is nil).
func sincosBall(x *big.Float, wp uint, want int) (*ball, *ball, bool) {
	j := trigSplit(wp)
	// The table path of sincosFx (scales up to sinCosTableMax) needs no
	// guard bits for doublings; sincosCore's j doublings lose about 2j bits.
	base := wp + 20
	w := base
	ws := getFxs(4)
	defer putFxs(ws)
	r, s, c := &(*ws)[0], &(*ws)[1], &(*ws)[2]
	var q int
	for {
		q = reduceHalfPi(r, &(*ws)[3], x, w+8)
		er, ok := r.expLower()
		if !ok {
			return nil, nil, false
		}
		need := base + uint(max(0, -er))
		if need > sinCosTableMax && base == wp+20 {
			base += 2 * j
			need += 2 * j
		}
		if need <= w {
			break
		}
		w = need
	}
	sincosFx(s, c, r.rescale(r, w), j)
	// x = qπ/2 + r: sin x, cos x = (sin r, cos r), (cos r, -sin r),
	// (-sin r, -cos r), (-cos r, sin r) for q = 0 … 3
	sf, cf := s, c
	sneg, cneg := q == 2 || q == 3, q == 1 || q == 2
	if q%2 == 1 {
		sf, cf = c, s
	}
	var sb, cb *ball
	if want != 1 {
		sb = sf.toBall(new(ball))
		if sneg {
			sb.neg(sb)
		}
	}
	if want != 0 {
		cb = cf.toBall(new(ball))
		if cneg {
			cb.neg(cb)
		}
	}
	return sb, cb, true
}

// Sin sets z to the correctly rounded sin x and returns z. Sin(±0) = ±0;
// Sin(±Inf) panics with ErrNaN.
func Sin(z, x *big.Float) *big.Float { return trig(z, x, 0) }

// Cos sets z to the correctly rounded cos x and returns z. Cos(±0) = 1;
// Cos(±Inf) panics with ErrNaN.
func Cos(z, x *big.Float) *big.Float { return trig(z, x, 1) }

// Tan sets z to the correctly rounded tan x and returns z. Tan(±0) = ±0;
// Tan(±Inf) panics with ErrNaN. (No binary number is an odd multiple of
// π/2, so the result is always finite.)
func Tan(z, x *big.Float) *big.Float { return trig(z, x, 2) }

func trig(z, x *big.Float, which int) *big.Float {
	precOf(z, x.Prec())
	if x.IsInf() {
		nan([]string{"Sin", "Cos", "Tan"}[which] + " of an infinity")
	}
	if x.Sign() == 0 {
		if which == 1 {
			return z.SetInt64(1)
		}
		return z.Set(x)
	}
	p := z.Prec()
	ex := int64(x.MantExp(nil))
	one := big.NewFloat(1)
	switch {
	case which == 0 && 3*ex <= nearBound(x, p):
		// x - sin x in (0, x³/6) for x > 0
		return finish(z, proxyNear(x, -x.Sign(), p), 0)
	case which == 1 && 2*ex-1 <= nearBound(one, p):
		// 1 - cos x in (0, x²/2)
		return finish(z, proxyNear(one, -1, p), 0)
	case which == 2 && 3*ex <= nearBound(x, p):
		// tan x - x in (0, |x|³) for 0 < x < 1/2
		return finish(z, proxyNear(x, x.Sign(), p), 0)
	}
	x = new(big.Float).Set(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		s, c, ok := sincosBall(x, wp, which)
		switch {
		case !ok:
			return nil, 0, false
		case which == 0:
			return s, 0, true
		case which == 1:
			return c, 0, true
		}
		// c has relative accuracy: it excludes zero
		var t ball
		must(t.quo(s, c, wp))
		return &t, 0, true
	})
}

// atanCore returns atan t for a fixed-point ball with |t| <= 1: j
// reductions t <- t/(1+sqrt(1+t²)) (atan t = 2 atan of the result), then
// the alternating series Σ (-1)^k t^(2k+1)/(2k+1), whose tail is below
// |t|^(2n+3). The result goes to dst (distinct from t).
func atanCore(dst, t *fx, j uint) {
	w := t.w
	ws := getFxs(3)
	defer putFxs(ws)
	y, one, u := &(*ws)[0], &(*ws)[1], &(*ws)[2]
	one.setInt(1, w)
	y.set(t)
	for i := uint(0); i < j; i++ {
		u.mul(y, y)
		u.add(u, one)
		must(u.sqrt(u))
		u.add(u, one)
		must(y.quo(y, u))
	}
	// lt >= 1-w: log2bnd counts the error bound
	lt := min(int64(y.log2bnd()), -1)
	n := (int64(w)+1+(-lt)-1)/(-lt)/2 + 1
	oddSeries(dst, y, n, true)
	dst.e = addUp(dst.e, 1) // tail
	dst.shl(dst, j)
}

// atanShift sets the grid of the atan table: atan(n/2^atanShift), as
// mpmath's ATAN_TAYLOR_SHIFT.
const atanShift = 7

// atanTable caches atan(n/2^7), n = 1 … 2^7, computed by atanCore.
var atanTable = &fxTable{compute: func(n int64, W uint) *fx {
	var a fx
	a.setInt(n, W)
	a.v.Rsh(&a.v, atanShift) // exact: n/2^7 at scale W > 7
	z := new(fx)
	atanCore(z, &a, trigSplit(W))
	return z
}}

// atanFx returns atan t for a fixed-point ball with |t| <= 1 (up to its
// error bound). Up to tableMaxPrec it uses mpmath's table reduction
// (atan_taylor, libelefun.py): with a = n/2^7 the grid point below |t|,
// atan t = atan a + atan d, d = (t-a)/(1+a·t) in [0, 2^-7), so that the
// series gains 14 bits per term with no square root; atan a comes from the
// table. Above, or for |t| < 2^-7, atanCore with j halvings. The result
// goes to dst (distinct from t).
func atanFx(dst, t *fx, j uint) {
	w := t.w
	if w > tableMaxPrec {
		atanCore(dst, t, j)
		return
	}
	neg := t.v.Sign() < 0
	ws := getFxs(4)
	defer putFxs(ws)
	at, a, one, den := &(*ws)[0], &(*ws)[1], &(*ws)[2], &(*ws)[3]
	at.set(t)
	if neg {
		at.v.Neg(&at.v)
	}
	a.v.Rsh(&at.v, w-atanShift)
	n := min(a.v.Int64(), 1<<atanShift)
	if n == 0 {
		atanCore(dst, t, 0)
		return
	}
	a.setInt(n, w)
	a.v.Rsh(&a.v, atanShift)
	one.setInt(1, w)
	den.mul(a, at)
	den.add(den, one)
	at.sub(at, a)
	must(at.quo(at, den)) // den >= 1 - ulps
	atanCore(dst, at, 0)
	dst.add(dst, atanTable.get(n, w))
	if neg {
		dst.v.Neg(&dst.v)
	}
}

// atanBall encloses atan t for a ball t that excludes zero, with about wp
// bits of relative accuracy: |atan T - atan t| <= |T - t|.
func atanBall(t *ball, wp uint) *ball {
	el, ok := t.expLower()
	must(ok)
	z := new(ball)
	eu := t.expUpper()
	j := trigSplit(wp)
	w := wp + 2*j + 16
	switch {
	case eu < -int64(wp/2)-4:
		// atan t = t(1 - t²/3 + ...): within |t|³ of t.
		z.set(t)
		z.addRad(magMul(magMul(t.upper(), t.upper()), t.upper()))
	case eu <= 0 || cmpAbs(&t.m, big.NewFloat(1)) <= 0:
		w += uint(max(0, -el))
		ws := getFxs(2)
		defer putFxs(ws)
		atanFx(&(*ws)[1], (*ws)[0].setBall(t, w), j)
		(*ws)[1].toBall(z)
	default:
		// |t| > 1: sign(t)·π/2 - atan(1/t), absolute accuracy suffices.
		var one, u, a ball
		one.setInt64(1)
		must(u.quo(&one, t, w))
		ws := getFxs(2)
		defer putFxs(ws)
		atanFx(&(*ws)[1], (*ws)[0].setBall(&u, w), j)
		(*ws)[1].toBall(&a)
		h := halfPiBall(w, t.m.Sign() < 0)
		z.sub(h, &a, wp+8)
	}
	return z
}

// halfPiBall encloses ±π/2.
func halfPiBall(w uint, neg bool) *ball {
	h := piBall(w)
	h.mul2exp(h, -1)
	if neg {
		h.neg(h)
	}
	return h
}

// atan2Ball encloses atan2(y, x) for exact non-zero finite y, x.
func atan2Ball(y, x *big.Float, wp uint) *ball {
	var yb, xb, q ball
	yb.setExact(y)
	xb.setExact(x)
	w := wp + 16
	if cmpAbs(y, x) <= 0 {
		must(q.quo(&yb, &xb, w))
		a := atanBall(&q, w)
		if x.Sign() < 0 {
			pb := piBall(w)
			if y.Sign() < 0 {
				pb.neg(pb)
			}
			a.add(a, pb, w)
		}
		return a
	}
	must(q.quo(&xb, &yb, w))
	a := atanBall(&q, w)
	return a.sub(halfPiBall(w, y.Sign() < 0), a, w)
}

// Atan sets z to the correctly rounded arctangent of x and returns z.
// Atan(±0) = ±0, Atan(±Inf) = ±π/2.
func Atan(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch {
	case x.IsInf():
		return piMultiple(z, x.Signbit(), 1, 1)
	case x.Sign() == 0:
		return z.Set(x)
	}
	p := z.Prec()
	// x - atan x in (0, x³/3) for x > 0
	if ex := int64(x.MantExp(nil)); 3*ex <= nearBound(x, p) {
		return finish(z, proxyNear(x, -x.Sign(), p), 0)
	}
	// atan of the exact x itself: atan2(x, 1) would divide by one first
	var xb ball
	xb.setExact(x)
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		return atanBall(&xb, wp), 0, true
	})
}

// piMultiple sets z to ±k·π/2^d (correctly rounded): the special values
// of Atan, Atan2, Asin and Acos.
func piMultiple(z *big.Float, neg bool, k int64, d int) *big.Float {
	precOf(z, 64)
	return ziv(z, 0, func(wp uint) (*ball, int64, bool) {
		b := piBall(wp + 8)
		var kb ball
		kb.setInt64(k)
		b.mul(b, &kb, wp+8)
		b.mul2exp(b, -d)
		if neg {
			b.neg(b)
		}
		return b, 0, true
	})
}

func atan2Exact(z, y, x *big.Float) *big.Float {
	y = new(big.Float).Set(y)
	x = new(big.Float).Set(x)
	return ziv(z, y.Prec()+x.Prec(), func(wp uint) (*ball, int64, bool) {
		return atan2Ball(y, x, wp), 0, true
	})
}

// Atan2 sets z to the correctly rounded angle of the point (x, y), in
// [-π, π], and returns z. If z's precision is 0 it is set to the larger of
// x's and y's. The special values follow C99 Annex F (atan2(±0, +0) = ±0,
// atan2(±0, -0) = ±π, atan2(±Inf, ±Inf) = ±π/4 or ±3π/4, …).
func Atan2(z, y, x *big.Float) *big.Float {
	precOf(z, max(x.Prec(), y.Prec()))
	ny := y.Signbit()
	switch {
	case y.IsInf() && x.IsInf() && x.Signbit():
		return piMultiple(z, ny, 3, 2)
	case y.IsInf() && x.IsInf():
		return piMultiple(z, ny, 1, 2)
	case y.IsInf() || x.Sign() == 0 && y.Sign() != 0:
		return piMultiple(z, ny, 1, 1)
	case (x.IsInf() || y.Sign() == 0) && x.Signbit():
		return piMultiple(z, ny, 1, 0)
	case x.IsInf() || y.Sign() == 0:
		return zeroSign(z.SetInt64(0), ny)
	}
	p := z.Prec()
	if x.Sign() > 0 && y.MantExp(nil) < x.MantExp(nil) {
		// atan(y/x) with y/x dyadic and tiny: proxy as Atan does.
		q := new(big.Float).SetPrec(y.Prec() + x.Prec() + p + 8)
		q.Quo(y, x) // exact when y/x is dyadic: it then has at most prec(y) bits
		exact := q.Acc() == big.Exact
		q.SetPrec(q.MinPrec()) // the shortest base allows the widest proxy
		eq := int64(q.MantExp(nil))
		if exact && 3*eq <= nearBound(q, p) {
			return finish(z, proxyNear(q, -q.Sign(), p), 0)
		}
	}
	return atan2Exact(z, y, x)
}

// zeroSign sets z (zero) to -0 when neg.
func zeroSign(z *big.Float, neg bool) *big.Float {
	if neg {
		z.Neg(z)
	}
	return z
}

// sqrt1mx2 encloses sqrt(1 - x²) for exact |x| < 1, from the exact
// product (1-x)(1+x) > 0.
func sqrt1mx2(x *big.Float, wp uint) *ball {
	pr := x.Prec() + uint(max(-x.MantExp(nil), 0)) + 4
	a := new(big.Float).SetPrec(pr).Sub(big.NewFloat(1), x)
	b := new(big.Float).SetPrec(pr).Add(big.NewFloat(1), x)
	c := new(big.Float).SetPrec(2*pr).Mul(a, b)
	var cb, s ball
	cb.setExact(c)
	must(s.sqrt(&cb, wp))
	return &s
}

// Asin sets z to the correctly rounded arcsine of x and returns z.
// Asin(±0) = ±0, Asin(±1) = ±π/2; |x| > 1 panics with ErrNaN.
func Asin(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch c := cmpAbs(x, big.NewFloat(1)); {
	case c > 0:
		nan("Asin of a number beyond [-1, 1]")
	case c == 0:
		return piMultiple(z, x.Signbit(), 1, 1)
	case x.Sign() == 0:
		return z.Set(x)
	}
	p := z.Prec()
	// asin x - x in (0, |x|³/4) for 0 < x <= 1/2
	if ex := int64(x.MantExp(nil)); 3*ex <= nearBound(x, p) {
		return finish(z, proxyNear(x, x.Sign(), p), 0)
	}
	x = new(big.Float).Set(x)
	small := cmpAbs(x, asinFxLimit) <= 0
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		if small {
			// relative accuracy near 0: as many more bits as x has
			// leading zeros
			w := wp + 16 + uint(max(0, -x.MantExp(nil)))
			b := new(ball)
			asinFx(x, w, trigSplit(wp), func(a *fx) { a.toBall(b) })
			return b, 0, true
		}
		s := sqrt1mx2(x, wp+16)
		var xb, q ball
		xb.setExact(x)
		must(q.quo(&xb, s, wp+16))
		return atanBall(&q, wp), 0, true
	})
}

// asinFxLimit bounds the arguments of asinFx: 1 - x² >= 1/4 for the
// fixed-point square root.
var asinFxLimit = big.NewFloat(0.86)

// asinFx evaluates asin x for an exact |x| <= 0.86 at scale w by mpmath's
// formula (mpf_asin, libelefun.py) asin x = 2·atan(x/(1+sqrt(1-x²))),
// whose atan argument stays below 0.6, all in fixed point, and hands the
// result to use before the pooled workspace is released.
func asinFx(x *big.Float, w, j uint, use func(*fx)) {
	ws := getFxs(5)
	defer putFxs(ws)
	X, C, one, T, A := &(*ws)[0], &(*ws)[1], &(*ws)[2], &(*ws)[3], &(*ws)[4]
	X.setFloat(x, w)
	one.setInt(1, w)
	C.mul(X, X)
	C.sub(one, C)
	must(C.sqrt(C)) // 1 - x² >= 1 - 0.86² > 1/4
	C.add(C, one)
	must(T.quo(X, C))
	atanFx(A, T, j)
	A.shl(A, 1)
	use(A)
}

// Acos sets z to the correctly rounded arccosine of x, in [0, π], and
// returns z. Acos(1) = +0, Acos(-1) = π; |x| > 1 panics with ErrNaN.
func Acos(z, x *big.Float) *big.Float {
	precOf(z, x.Prec())
	switch c := cmpAbs(x, big.NewFloat(1)); {
	case c > 0:
		nan("Acos of a number beyond [-1, 1]")
	case c == 0 && x.Sign() > 0:
		return z.SetInt64(0)
	case c == 0:
		return piMultiple(z, false, 1, 0)
	case x.Sign() == 0:
		return piMultiple(z, false, 1, 1)
	}
	x = new(big.Float).Set(x)
	small := cmpAbs(x, asinFxLimit) <= 0
	return ziv(z, x.Prec(), func(wp uint) (*ball, int64, bool) {
		if small {
			// acos x = π/2 - asin x, in [0.52, 2.62]: absolute accuracy
			// suffices. π at scale w-1 is π/2 at scale w.
			w := wp + 16
			b := new(ball)
			asinFx(x, w, trigSplit(wp), func(a *fx) {
				var h fx
				h.set(piFx(w - 1))
				h.w = w
				h.sub(&h, a)
				h.toBall(b)
			})
			return b, 0, true
		}
		s := sqrt1mx2(x, wp+16)
		var xb, q ball
		xb.setExact(x)
		must(q.quo(s, &xb, wp+16))
		b := atanBall(&q, wp+8)
		if x.Sign() < 0 {
			b.add(b, piBall(wp+8), wp+8)
		}
		return b, 0, true
	})
}
