// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

// Package oracle produces reference values for validating numerical code:
// correctly rounded results of mathematical functions, at any precision
// and in all six rounding modes of math/big, obtained from an [Oracle]
// that can evaluate a function to a requested number of bits.
//
// Two oracles are provided:
//
//   - [Mpmath] drives mpmath (Python, BSD-3) as a subprocess over a
//     JSON-lines pipe. Inputs and outputs travel exactly, as an integer
//     mantissa and a binary exponent, never as decimal strings.
//   - [Self] evaluates with this module's own correctly rounded functions,
//     so that once validated against mpmath the package can stand in for it
//     (the plan for go-numerics/special and go-statistics/stats).
//
// A [Rounder] turns approximations into certified roundings: it asks the
// oracle for far more bits than the target precision, widens the answer by
// the oracle's declared slack, and raises the precision (Ziv's strategy)
// until both ends of that interval round to the same value in every mode.
// A decision is then confirmed by a second evaluation at twice the
// precision, which must decide the same way: an oracle that is less
// accurate than it claims at one precision (a cancellation it did not
// compensate) is caught rather than trusted. Exact results — whose interval always straddles a representable value —
// are recognised beforehand by [Exact], with integer arithmetic written
// independently of the package under test.
package oracle

import (
	"errors"
	"fmt"
	"math/big"
)

// An Oracle evaluates fn at args with about prec bits of relative accuracy:
// the returned value must lie within 2^(Slack()-prec) of the true value,
// relatively. Function names are the lower-case names of the bigfloat
// functions ("exp", "atan2", "pi", …; see [Functions]).
type Oracle interface {
	Eval(fn string, prec uint, args ...*big.Float) (*big.Float, error)
	// Slack is the number of low-order bits of an Eval result that are not
	// trusted.
	Slack() uint
	// Name identifies the oracle in reports ("mpmath 1.4.1", "self").
	Name() string
}

// Functions lists the functions an Oracle must know, with their arity.
var Functions = map[string]int{
	"sqrt": 1, "cbrt": 1, "root": 2, "hypot": 2,
	"exp": 1, "expm1": 1, "log": 1, "log1p": 1, "log2": 1, "log10": 1, "pow": 2,
	"sin": 1, "cos": 1, "tan": 1, "asin": 1, "acos": 1, "atan": 1, "atan2": 2,
	"sinh": 1, "cosh": 1, "tanh": 1, "asinh": 1, "acosh": 1, "atanh": 1,
	"pi": 0, "e": 0, "ln2": 0, "ln10": 0, "euler": 0,
}

// Modes lists big.Float's six rounding modes, in the order fixtures store
// them.
var Modes = []big.RoundingMode{
	big.ToNearestEven, big.ToNearestAway, big.ToZero,
	big.AwayFromZero, big.ToNegativeInf, big.ToPositiveInf,
}

// ErrUndecided is returned by Round when the oracle's precision limit is
// reached before the rounding is settled.
var ErrUndecided = errors.New("oracle: rounding not decided within the precision limit")

// Result is a correctly rounded value together with the side of the true
// value it lies on, as big.Float.Acc reports it (Below: the result is below
// the true value; Exact: the true value is representable).
type Result struct {
	Value *big.Float
	Acc   big.Accuracy
}

// Rounder computes correctly rounded results from an Oracle.
type Rounder struct {
	Oracle Oracle
	// MaxPrec bounds the working precision asked of the oracle; 0 means
	// 64 times the target precision plus 65536.
	MaxPrec uint
}

// Round returns fn(args) correctly rounded to prec bits in every mode of
// Modes (one result per mode, in that order). The arguments must lie in
// the function's domain (no NaN results).
func (r *Rounder) Round(fn string, prec uint, args ...*big.Float) ([]Result, error) {
	n, ok := Functions[fn]
	if !ok {
		return nil, fmt.Errorf("oracle: unknown function %q", fn)
	}
	if len(args) != n {
		return nil, fmt.Errorf("oracle: %s takes %d arguments, got %d", fn, n, len(args))
	}
	if v, ok := Exact(fn, args...); ok {
		return roundAll(v, prec), nil
	}
	limit := r.MaxPrec
	if limit == 0 {
		limit = 64*prec + 1<<16
	}
	slack := r.Oracle.Slack()
	var prev []Result
	for q := 2*prec + 64; q <= limit; q *= 2 {
		v, err := r.Oracle.Eval(fn, q, args...)
		if err != nil {
			return nil, err
		}
		res, ok := decideAll(v, q-slack, prec)
		switch {
		case ok && prev != nil && sameResults(res, prev):
			return res, nil
		case ok:
			prev = res
		default:
			prev = nil
		}
	}
	return nil, ErrUndecided
}

func sameResults(a, b []Result) bool {
	for i := range a {
		if a[i].Value.Cmp(b[i].Value) != 0 || a[i].Acc != b[i].Acc {
			return false
		}
	}
	return true
}

// roundAll rounds an exact value in every mode.
func roundAll(v *big.Float, prec uint) []Result {
	out := make([]Result, len(Modes))
	for i, m := range Modes {
		z := new(big.Float).SetPrec(prec).SetMode(m).Set(v)
		out[i] = Result{z, z.Acc()}
	}
	return out
}

// decideAll rounds the interval v·(1 ± 2^-rel) in every mode; ok when both
// ends round alike in all of them and the rounded value lies outside the
// interval (so that its side is known). An infinite or zero v is never
// decided: the oracle returns those only for exact results, which Exact
// handles, or on overflow, which the caller must not ask about.
func decideAll(v *big.Float, rel, prec uint) ([]Result, bool) {
	if v.IsInf() || v.Sign() == 0 {
		return nil, false
	}
	e := v.MantExp(nil)
	rad := new(big.Float).SetMantExp(big.NewFloat(1), e-int(rel))
	lo := new(big.Float).SetPrec(v.Prec()+rel+8).Sub(v, rad)
	hi := new(big.Float).SetPrec(v.Prec()+rel+8).Add(v, rad)
	out := make([]Result, len(Modes))
	for i, m := range Modes {
		a := new(big.Float).SetPrec(prec).SetMode(m).Set(lo)
		b := new(big.Float).SetPrec(prec).SetMode(m).Set(hi)
		if a.Cmp(b) != 0 {
			return nil, false
		}
		var acc big.Accuracy
		switch {
		case a.Cmp(lo) < 0:
			acc = big.Below
		case a.Cmp(hi) > 0:
			acc = big.Above
		default:
			return nil, false
		}
		out[i] = Result{a, acc}
	}
	return out, true
}
