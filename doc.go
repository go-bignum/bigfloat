// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

// Package bigfloat provides correctly rounded elementary functions and
// constants on math/big's Float, at any precision and in all six rounding
// modes.
//
// Every function sets its receiver z to the exact result rounded to z's
// precision in z's rounding mode — the value an infinitely precise
// computation would round to — and z.Acc() reports whether that value is
// below, equal to or above the exact result:
//
//	z := new(big.Float).SetPrec(53).SetMode(big.ToNegativeInf)
//	bigfloat.Sin(z, x) // sin x rounded down to 53 bits, for any x
//
// The functions follow math/big's conventions: a receiver of precision 0
// takes the operand's precision (the larger of two; 64 bits for the
// constants), the receiver may alias an operand, and results beyond
// big.Float's exponent range overflow to ±Inf or underflow to ±0 exactly as
// big.Float's own arithmetic does. Special values follow C99 Annex F. An
// operation whose result would be NaN panics with ErrNaN, as big.Float.Sqrt
// panics with big.ErrNaN; a trigonometric argument beyond
// MaxReductionExponent panics with ErrLimit.
//
// Functions: Sqrt, Cbrt, Root, Hypot; Exp, Expm1, Log, Log1p, Log2, Log10,
// Pow; Sin, Cos, Tan, Asin, Acos, Atan, Atan2; Sinh, Cosh, Tanh, Asinh,
// Acosh, Atanh. Constants: Pi, E, Ln2, Ln10, EulerGamma.
//
// Correct rounding comes from rigorous enclosures (midpoint-radius balls,
// fixed-point series with bounded errors) evaluated by Ziv's strategy, exact
// cases being detected beforehand; docs/design.md describes the algorithms
// and their references. Package oracle validates the results against
// mpmath, bit for bit, and can stand in for mpmath elsewhere.
package bigfloat
