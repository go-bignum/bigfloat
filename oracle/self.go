// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package oracle

import (
	"fmt"
	"math/big"

	"github.com/go-bignum/bigfloat"
)

// Self is an Oracle backed by this module's own correctly rounded
// functions: an Eval at prec bits is the round-to-nearest value, within
// half an ulp, so a single bit of slack covers it. It is what downstream
// packages use once bigfloat has been validated against mpmath (by this
// package's fixtures), without a Python toolchain.
type Self struct{}

// Name implements Oracle.
func (Self) Name() string { return "self" }

// Slack implements Oracle.
func (Self) Slack() uint { return 1 }

// Unary lists the bigfloat functions of one argument by oracle name.
var Unary = map[string]func(z, x *big.Float) *big.Float{
	"sqrt": bigfloat.Sqrt, "cbrt": bigfloat.Cbrt,
	"exp": bigfloat.Exp, "expm1": bigfloat.Expm1,
	"log": bigfloat.Log, "log1p": bigfloat.Log1p, "log2": bigfloat.Log2, "log10": bigfloat.Log10,
	"sin": bigfloat.Sin, "cos": bigfloat.Cos, "tan": bigfloat.Tan,
	"asin": bigfloat.Asin, "acos": bigfloat.Acos, "atan": bigfloat.Atan,
	"sinh": bigfloat.Sinh, "cosh": bigfloat.Cosh, "tanh": bigfloat.Tanh,
	"asinh": bigfloat.Asinh, "acosh": bigfloat.Acosh, "atanh": bigfloat.Atanh,
}

// Binary lists the bigfloat functions of two arguments by oracle name.
var Binary = map[string]func(z, x, y *big.Float) *big.Float{
	"hypot": bigfloat.Hypot, "pow": bigfloat.Pow, "atan2": bigfloat.Atan2,
	"root": func(z, x, n *big.Float) *big.Float {
		k, _ := n.Uint64()
		return bigfloat.Root(z, x, uint(k))
	},
}

// Constants lists the bigfloat constants by oracle name.
var Constants = map[string]func(z *big.Float) *big.Float{
	"pi": bigfloat.Pi, "e": bigfloat.E, "ln2": bigfloat.Ln2, "ln10": bigfloat.Ln10, "euler": bigfloat.EulerGamma,
}

// Call evaluates the named bigfloat function into z (at z's precision and
// mode), panicking on an unknown name or a wrong number of arguments.
func Call(fn string, z *big.Float, args ...*big.Float) *big.Float {
	switch n := len(args); {
	case n == 0 && Constants[fn] != nil:
		return Constants[fn](z)
	case n == 1 && Unary[fn] != nil:
		return Unary[fn](z, args[0])
	case n == 2 && Binary[fn] != nil:
		return Binary[fn](z, args[0], args[1])
	}
	panic(fmt.Sprintf("oracle: no function %s of %d arguments", fn, len(args)))
}

// Eval implements Oracle.
func (Self) Eval(fn string, prec uint, args ...*big.Float) (v *big.Float, err error) {
	defer func() {
		if e := recover(); e != nil {
			v, err = nil, fmt.Errorf("oracle: self %s: %v", fn, e)
		}
	}()
	return Call(fn, new(big.Float).SetPrec(prec), args...), nil
}
