// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

// Command oracle generates the golden fixtures of package bigfloat: for
// every function, random arguments per region (tiny, near 1, huge, near
// singularities, exact cases …) at precisions 24, 53, 64, 113, 200 and
// 1000 bits, rounded correctly in all six modes by an oracle (mpmath by
// default) at a precision well above the target.
//
//	go run ./cmd/oracle -out testdata            # all functions
//	go run ./cmd/oracle -fn exp,log -scale 0.1   # a few, fewer points
//	go run ./cmd/oracle -check                   # compare bigfloat, write nothing
//
// The fixtures are deterministic for a given seed.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"maps"
	"math/big"
	"math/rand/v2"
	"os"
	"path/filepath"
	"slices"
	"strings"

	"github.com/go-bignum/bigfloat/oracle"
)

// Precs are the fixture precisions, with the number of cases per function
// at each.
var Precs = []struct {
	Prec  uint
	Cases int
}{{24, 240}, {53, 240}, {64, 240}, {113, 120}, {200, 60}, {1000, 16}}

// constPrecs are the precisions at which constants are checked: every
// precision up to 300 bits, then a few larger ones.
func constPrecs() []uint {
	var ps []uint
	for p := uint(2); p <= 300; p++ {
		ps = append(ps, p)
	}
	return append(ps, 500, 1000, 2000)
}

func main() { os.Exit(run(os.Args[1:], os.Stdout, os.Stderr)) }

func run(args []string, stdout, stderr io.Writer) int {
	fs := flag.NewFlagSet("oracle", flag.ContinueOnError)
	fs.SetOutput(stderr)
	out := fs.String("out", "testdata", "output directory")
	fns := fs.String("fn", "", "comma-separated functions (default: all)")
	seed := fs.Uint64("seed", 20261017, "random seed")
	scale := fs.Float64("scale", 1, "multiply the number of cases")
	check := fs.Bool("check", false, "compare bigfloat with the oracle; write nothing")
	self := fs.Bool("self", false, "use bigfloat itself as the oracle (no Python)")
	cmdJSON := fs.String("cmd", "", "oracle command as a JSON array (default: uv + mpmath)")
	if err := fs.Parse(args); err != nil {
		return 2
	}
	var o oracle.Oracle
	switch {
	case *self:
		o = oracle.Self{}
	default:
		var command []string
		if *cmdJSON != "" {
			if err := json.Unmarshal([]byte(*cmdJSON), &command); err != nil || len(command) == 0 {
				fmt.Fprintf(stderr, "oracle: -cmd: want a JSON array of strings\n")
				return 2
			}
		}
		m, err := oracle.StartMpmath(command)
		if err != nil {
			fmt.Fprintln(stderr, err)
			return 1
		}
		defer m.Close()
		o = m
	}
	names := slices.Sorted(maps.Keys(oracle.Functions))
	if *fns != "" {
		names = strings.Split(*fns, ",")
	}
	rd := &oracle.Rounder{Oracle: o}
	status := 0
	for _, fn := range names {
		if _, ok := oracle.Functions[fn]; !ok {
			fmt.Fprintf(stderr, "oracle: unknown function %q\n", fn)
			return 2
		}
		f, err := generate(rd, fn, *seed, *scale)
		if err != nil {
			fmt.Fprintf(stderr, "%s: %v\n", fn, err)
			return 1
		}
		if *check {
			bad := 0
			for i := range f.Cases {
				for _, msg := range f.Cases[i].Check(fn, func(z *big.Float, a ...*big.Float) *big.Float { return oracle.Call(fn, z, a...) }) {
					bad++
					if bad <= 20 {
						fmt.Fprintln(stdout, msg)
					}
				}
			}
			fmt.Fprintf(stdout, "%-6s %5d cases × %d modes: %d mismatches\n", fn, len(f.Cases), len(oracle.Modes), bad)
			if bad > 0 {
				status = 1
			}
			continue
		}
		if err := writeFile(filepath.Join(*out, fn+".json.gz"), f); err != nil {
			fmt.Fprintln(stderr, err)
			return 1
		}
		fmt.Fprintf(stdout, "%-6s %5d cases\n", fn, len(f.Cases))
	}
	return status
}

// generate draws and rounds the cases of one function.
func generate(rd *oracle.Rounder, fn string, seed uint64, scale float64) (*oracle.Fixture, error) {
	f := &oracle.Fixture{Fn: fn, Oracle: rd.Oracle.Name(), Seed: seed}
	h := uint64(0)
	for _, c := range fn {
		h = h*131 + uint64(c)
	}
	r := rand.New(rand.NewPCG(seed, h))
	add := func(name string, p uint, args []*big.Float) error {
		res, err := rd.Round(fn, p, args...)
		if err != nil {
			return fmt.Errorf("%s at %d bits, %v: %w", name, p, args, err)
		}
		f.Cases = append(f.Cases, oracle.NewCase(name, p, args, res))
		return nil
	}
	if oracle.Functions[fn] == 0 {
		for _, p := range constPrecs() {
			if err := add("const", p, nil); err != nil {
				return nil, err
			}
		}
		return f, nil
	}
	regions := regionsOf(fn)
	for _, pc := range Precs {
		per := max(int(float64(pc.Cases)*scale)/len(regions), 2)
		for _, rg := range regions {
			for i := 0; i < per; i++ {
				if err := add(rg.name, pc.Prec, rg.gen(r, pc.Prec)); err != nil {
					return nil, err
				}
			}
		}
	}
	return f, nil
}

func writeFile(path string, f *oracle.Fixture) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	fh, err := os.Create(path)
	if err != nil {
		return err
	}
	if err := oracle.WriteFixture(fh, f); err != nil {
		fh.Close()
		return err
	}
	return fh.Close()
}
