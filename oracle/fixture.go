// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package oracle

import (
	"compress/gzip"
	"encoding/json"
	"fmt"
	"io"
	"math/big"
	"os"
	"strings"
)

// A Fixture is the golden data of one function: correctly rounded results
// in every mode for a set of arguments, as produced by a Rounder.
type Fixture struct {
	Fn     string `json:"fn"`
	Oracle string `json:"oracle"`
	Seed   uint64 `json:"seed"`
	Cases  []Case `json:"cases"`
}

// A Case is one evaluation: arguments, target precision, and the result in
// each mode of Modes. Numbers are written by [Format] (exact hexadecimal),
// accuracies as one character per mode: '<' Below, '=' Exact, '>' Above.
type Case struct {
	Region string   `json:"r"`
	Prec   uint     `json:"p"`
	Args   []string `json:"a"`
	Want   []string `json:"w"`
	Acc    string   `json:"c"`
}

// Format writes x exactly: hexadecimal mantissa and binary exponent
// ("0x.8p+1"), "+Inf"/"-Inf", "0"/"-0".
func Format(x *big.Float) string {
	if x.Sign() == 0 {
		if x.Signbit() {
			return "-0"
		}
		return "0"
	}
	return x.Text('p', 0)
}

// Parse reads a number written by Format, at the smallest precision that
// holds it exactly but at least minPrec.
func Parse(s string, minPrec uint) (*big.Float, error) {
	switch s {
	case "+Inf", "-Inf":
		return new(big.Float).SetPrec(minPrec).SetInf(s[0] == '-'), nil
	case "-0":
		return new(big.Float).SetPrec(minPrec).Neg(new(big.Float)), nil
	}
	// a hex mantissa of n digits needs at most 4n bits
	x, _, err := big.ParseFloat(s, 0, uint(4*len(s)+64), big.ToNearestEven)
	if err != nil {
		return nil, err
	}
	return x.SetPrec(max(x.MinPrec(), minPrec, 1)), nil
}

// AccChar maps an accuracy to its fixture character.
func AccChar(a big.Accuracy) byte {
	switch a {
	case big.Below:
		return '<'
	case big.Above:
		return '>'
	}
	return '='
}

// NewCase builds a Case from a Rounder's results.
func NewCase(region string, prec uint, args []*big.Float, res []Result) Case {
	c := Case{Region: region, Prec: prec}
	var acc strings.Builder
	for _, a := range args {
		c.Args = append(c.Args, Format(a))
	}
	for _, r := range res {
		c.Want = append(c.Want, Format(r.Value))
		acc.WriteByte(AccChar(r.Acc))
	}
	c.Acc = acc.String()
	return c
}

// ParseArgs returns the case's arguments; each argument keeps at least the
// case's precision, as when the fixture was generated.
func (c *Case) ParseArgs() ([]*big.Float, error) {
	out := make([]*big.Float, len(c.Args))
	for i, s := range c.Args {
		x, err := Parse(s, c.Prec)
		if err != nil {
			return nil, err
		}
		out[i] = x
	}
	return out, nil
}

// WriteFixture writes f gzip-compressed JSON to w.
func WriteFixture(w io.Writer, f *Fixture) error {
	zw := gzip.NewWriter(w)
	enc := json.NewEncoder(zw)
	if err := enc.Encode(f); err != nil {
		return err
	}
	return zw.Close()
}

// ReadFixture reads a fixture written by WriteFixture.
func ReadFixture(r io.Reader) (*Fixture, error) {
	zr, err := gzip.NewReader(r)
	if err != nil {
		return nil, err
	}
	var f Fixture
	if err := json.NewDecoder(zr).Decode(&f); err != nil {
		return nil, fmt.Errorf("oracle: decoding fixture: %w", err)
	}
	return &f, nil
}

// LoadFixture reads a fixture file.
func LoadFixture(path string) (*Fixture, error) {
	fh, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer fh.Close()
	return ReadFixture(fh)
}

// Check evaluates fn with eval (into a z of the case's precision and each
// mode in turn) and returns a description of every disagreement with the
// fixture — value, sign of zero or accuracy — or nil.
func (c *Case) Check(fn string, eval func(z *big.Float, args ...*big.Float) *big.Float) []string {
	args, err := c.ParseArgs()
	if err != nil {
		return []string{err.Error()}
	}
	var bad []string
	for i, m := range Modes {
		z := new(big.Float).SetPrec(c.Prec).SetMode(m)
		eval(z, args...)
		got := Format(z)
		if got != c.Want[i] || AccChar(z.Acc()) != c.Acc[i] {
			bad = append(bad, fmt.Sprintf("%s(%s) prec %d %v: got %s (%c), want %s (%c) [%s]",
				fn, strings.Join(c.Args, ", "), c.Prec, m, got, AccChar(z.Acc()), c.Want[i], c.Acc[i], c.Region))
		}
	}
	return bad
}
