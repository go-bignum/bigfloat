// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat

import (
	"fmt"
	"math/big"
	"sync"
)

// A proxy is an exact big.Float v standing in for a real result y: rounding
// v to the target precision, in the target mode, gives the correctly rounded
// y, and v lies on the same side of that rounded value as y. Setting z from
// a proxy therefore yields both the right value and the right z.Acc().
//
// proxyNear returns base + dir·2^(exp(base)-P-1), P = max(p+2, prec(base))+2.
// It is a proxy for any y with 0 < (y-base)·dir < 2^(exp(base)-P+1): base
// and every rounding breakpoint of precision p near it (representable
// numbers and midpoints) lie on a common grid of spacing 2^(exp(base)-P+2)
// or half of it below a power of two, so neither y nor v crosses one.
// base is finite and non-zero, dir is ±1.
func proxyNear(base *big.Float, dir int, p uint) *big.Float {
	var d big.Float
	return proxyNearTo(new(big.Float), &d, base, dir, p)
}

// proxyNearTo is proxyNear computed into z, with d as scratch (both
// distinct from base), so that a caller holding them across calls does not
// allocate.
func proxyNearTo(z, d, base *big.Float, dir int, p uint) *big.Float {
	P := max(p+2, base.Prec()) + 2
	e := base.MantExp(nil)
	d.SetInt64(int64(dir))
	d.SetMantExp(d, e-int(P)-1)
	return z.SetPrec(P+2).SetMode(big.ToNearestEven).Add(base, d) // exact
}

// must panics on an internal invariant that does not hold: a ball operation
// that cannot fail for the arguments its caller builds (a divisor or a
// logarithm's argument bounded away from zero by construction).
func must(ok bool) {
	if !ok {
		panic("bigfloat: internal error: enclosure invariant violated")
	}
}

// nearBound returns the exponent b such that proxyNear(base, ±1, p) may be
// used when |y - base| < 2^b.
func nearBound(base *big.Float, p uint) int64 {
	P := max(p+2, base.Prec()) + 2
	return int64(base.MantExp(nil)) - int64(P) + 1
}

// finish sets z to the rounding of v·2^K at z's precision and mode, with the
// accuracy and the overflow/underflow semantics of big.Float itself: the
// result is ±Inf when the rounded exponent exceeds MaxExp, ±0 when it falls
// below MinExp.
//
// v is scaled exactly first (at its own precision) and z.Set rounds it once:
// that keeps z's precision and mode and sets z.Acc(). (SetMantExp would
// copy the operand's precision, and resets the accuracy to Exact even when
// it only moves an exponent.) A value whose rounded exponent is MaxExp+1
// overflows inside z.Set, as in big.Float's own arithmetic. v is finite and
// non-zero (zero and infinite results are special values, set directly).
// v is the caller's to give up: it is scaled in place.
func finish(z, v *big.Float, K int64) *big.Float {
	e := int64(v.MantExp(nil)) + K
	switch {
	case e > big.MaxExp:
		return overflow(z, v.Signbit())
	case e < big.MinExp-1:
		return underflow(z, v.Signbit())
	case e == big.MinExp-1:
		// Just below the smallest normal binade: rounding may carry it
		// up into range. Round one binade higher, then move down; the
		// value and the underflow are big.Float's, the accuracy of such a
		// result is reported as Exact.
		z.Set(scale2(v, K+1))
		return z.SetMantExp(z, -1)
	}
	return z.Set(scale2(v, K))
}

// scale2 sets v to v·2^K exactly, at v's precision, and returns it. |K|
// may reach 2^32: it is applied in two steps so that each fits an int on
// 32-bit targets (the intermediate exponent lies between the first and the
// last, both in range).
func scale2(v *big.Float, K int64) *big.Float {
	k1 := K / 2
	v.SetMantExp(v, int(k1))
	return v.SetMantExp(v, int(K-k1))
}

// overflow sets z to ±Inf the way big.Float's own arithmetic does on
// exponent overflow (with z.Acc() reporting the direction).
func overflow(z *big.Float, neg bool) *big.Float {
	z.SetInt64(1)
	if neg {
		z.Neg(z)
	}
	return z.SetMantExp(z, big.MaxExp)
}

// underflow sets z to ±0 the way big.Float does on exponent underflow.
func underflow(z *big.Float, neg bool) *big.Float {
	z.SetFloat64(0.25)
	if neg {
		z.Neg(z)
	}
	return z.SetMantExp(z, big.MinExp)
}

// evalFn evaluates a function at working precision wp and returns a ball
// containing f(x)·2^-K. ok = false asks for more precision (for instance
// when an argument reduction could not separate the result from zero).
type evalFn func(wp uint) (b *ball, K int64, ok bool)

// zivLimit is the working precision beyond which the Ziv loop gives up with
// a panic: no correctly rounded result of this package needs it (worst
// cases of the table maker's dilemma need about 2p bits plus the
// argument's own precision), so reaching it reveals a missed exact case.
var zivLimit = func(p, inPrec uint) uint { return 64*(p+inPrec) + 1<<16 }

// ziv implements Ziv's strategy (ACM TOMS 17(3), 1991): evaluate an
// enclosure at a working precision, and when both of its ends round to the
// same value, that value is the correctly rounded result; otherwise raise
// the precision and try again. Exact results never reach it: each function
// detects them first, since an enclosure of a representable value or of a
// rounding midpoint never rounds to a single value.
func ziv(z *big.Float, inPrec uint, eval evalFn) *big.Float {
	p := z.Prec()
	mode := z.Mode()
	limit := zivLimit(p, inPrec)
	ws := zivPool.Get().(*zivWS)
	defer zivPool.Put(ws)
	for wp := p + 24; ; wp += wp/2 + 16 {
		if wp > limit {
			panic(fmt.Sprintf("bigfloat: internal error: Ziv loop did not converge (precision %d, working precision %d)", p, wp))
		}
		b, K, ok := eval(wp)
		if !ok {
			continue
		}
		if v := decide(b, p, mode, ws); v != nil {
			return finish(z, v, K)
		}
	}
}

// zivWS holds the scratch values of decide, recycled through zivPool so
// that their mantissas keep their capacity from one call to the next.
type zivWS struct{ rad, lo, hi, d, pd, pv big.Float }

var zivPool = sync.Pool{New: func() any { return new(zivWS) }}

// decide returns a proxy of the correctly rounded value of every point of b,
// or nil when the points do not all round alike or when the rounded value
// lies inside the ball (its side is then unknown). The proxy is ws.pv.
func decide(b *ball, p uint, mode big.RoundingMode, ws *zivWS) *big.Float {
	rad, lo, hi, d := &ws.rad, &ws.lo, &ws.hi, &ws.d
	b.r.float(rad)
	lo.SetPrec(p).SetMode(mode).Sub(&b.m, rad)
	hi.SetPrec(p).SetMode(mode).Add(&b.m, rad)
	if lo.Cmp(hi) != 0 || lo.Signbit() != hi.Signbit() || lo.Sign() == 0 || lo.IsInf() {
		return nil
	}
	// d = lo - m, exactly: the bits of lo and m together span at most
	// max(p, prec(m)) + |exp(lo) - exp(m)| positions, plus a carry.
	em, el := b.m.MantExp(nil), lo.MantExp(nil)
	ed := uint(max(em-el, el-em))
	d.SetPrec(max(p, b.m.Prec()) + ed + 4).SetMode(big.ToNearestEven)
	d.Sub(lo, &b.m)
	sd := d.Sign()
	if d.Abs(d).Cmp(rad) <= 0 {
		return nil
	}
	return proxyNearTo(&ws.pv, &ws.pd, lo, -sd, p)
}

// cmpAbs compares |a| and |b|.
func cmpAbs(a, b *big.Float) int {
	var x, y big.Float
	return x.Abs(a).Cmp(y.Abs(b))
}

// precOf returns the precision of z, set from def when z has none (the
// big.Float convention for operations: a zero-precision receiver takes its
// operands' precision).
func precOf(z *big.Float, def uint) uint {
	if z.Prec() == 0 {
		z.SetPrec(def)
	}
	return z.Prec()
}

// ErrNaN is the panic value of a function whose result would be NaN, which
// big.Float cannot represent. It mirrors big.ErrNaN, the panic value of
// big.Float's own operations in that case (Sqrt of a negative number,
// ∞-∞, 0·∞, …); big.ErrNaN's message field is unexported, hence a type of
// the same shape here. Recover it as big.ErrNaN is recovered:
//
//	defer func() {
//		if e, ok := recover().(bigfloat.ErrNaN); ok { ... }
//	}()
type ErrNaN struct {
	Msg string
}

func (e ErrNaN) Error() string { return e.Msg }

func nan(msg string) {
	panic(ErrNaN{Msg: "bigfloat: " + msg})
}
