// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

// Command mpfrcatalan writes the golden fixture of bigfloat.Catalan from
// MPFR itself: julia's Base.MathConstants.catalan calls mpfr_const_catalan,
// and catalan.jl prints its result at each precision in each of MPFR's five
// rounding modes (RNDN, RNDZ, RNDA, RNDD, RNDU).
//
//	go run ./cmd/mpfrcatalan -out testdata
//
// big.ToNearestAway has no MPFR counterpart (RNDNA is not supported by
// mpfr_const_catalan). It differs from ToNearestEven only on a tie, which
// needs G to be a number of p+1 bits; the command refuses to write a case
// unless the p+1-bit roundings toward and away from zero differ, i.e. G is
// not such a number, so ToNearestAway is the ToNearestEven result.
package main

import (
	"bufio"
	"bytes"
	_ "embed"
	"errors"
	"flag"
	"fmt"
	"io"
	"math/big"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/go-bignum/bigfloat/oracle"
)

//go:embed catalan.jl
var script string

// pkgx is the launcher of julia, replaced by tests.
var pkgx = "pkgx"

// julia runs the script at the given precisions and returns its output. It
// is replaced by tests.
var julia = func(precs []uint) ([]byte, error) {
	args := []string{"julia@1.12.6", "--startup-file=no", "-e", script}
	for _, p := range precs {
		args = append(args, strconv.FormatUint(uint64(p), 10))
	}
	return exec.Command(pkgx, args...).Output()
}

// exit is os.Exit, replaceable by tests.
var exit = os.Exit

func main() { exit(run(os.Args[1:], os.Stdout, os.Stderr)) }

// precisions returns the precisions of the fixture: 1 to 1000, then 4096
// and 10000.
func precisions() []uint {
	var ps []uint
	for p := uint(1); p <= 1000; p++ {
		ps = append(ps, p)
	}
	return append(ps, 4096, 10000)
}

func run(args []string, stdout, stderr io.Writer) int {
	fs := flag.NewFlagSet("mpfrcatalan", flag.ContinueOnError)
	fs.SetOutput(stderr)
	out := fs.String("out", "testdata", "output directory")
	if err := fs.Parse(args); err != nil {
		return 2
	}
	f, err := generate(precisions())
	if err != nil {
		fmt.Fprintln(stderr, "mpfrcatalan:", err)
		return 1
	}
	path := filepath.Join(*out, "catalan.json.gz")
	if err := os.MkdirAll(*out, 0o755); err != nil {
		fmt.Fprintln(stderr, "mpfrcatalan:", err)
		return 1
	}
	var buf bytes.Buffer
	_ = oracle.WriteFixture(&buf, f) // a bytes.Buffer never refuses a write
	if err := os.WriteFile(path, buf.Bytes(), 0o644); err != nil {
		fmt.Fprintln(stderr, "mpfrcatalan:", err)
		return 1
	}
	fmt.Fprintf(stdout, "catalan %5d cases (%s)\n", len(f.Cases), f.Oracle)
	return 0
}

// generate runs julia at the precisions and their successors (the p+1-bit
// results exclude ties) and builds the fixture.
func generate(precs []uint) (*oracle.Fixture, error) {
	ask := map[uint]bool{}
	for _, p := range precs {
		ask[p], ask[p+1] = true, true
	}
	var list []uint
	for p := range ask {
		list = append(list, p)
	}
	raw, err := julia(list)
	if err != nil {
		return nil, fmt.Errorf("julia: %w", err)
	}
	rows, version, err := parse(bytes.NewReader(raw))
	if err != nil {
		return nil, err
	}
	f := &oracle.Fixture{Fn: "catalan", Oracle: "MPFR " + version + " mpfr_const_catalan via julia 1.12.6"}
	for _, p := range precs {
		c, err := buildCase(p, rows)
		if err != nil {
			return nil, err
		}
		f.Cases = append(f.Cases, c)
	}
	return f, nil
}

// row is MPFR's result at one precision, by mode letter (N, Z, A, D, U).
type row map[string]*big.Float

// parse reads the script's output: lines "prec mode mantissa exponent" and
// a final "# version".
func parse(r io.Reader) (map[uint]row, string, error) {
	rows := map[uint]row{}
	version := ""
	sc := bufio.NewScanner(r)
	for sc.Scan() {
		line := sc.Text()
		if v, ok := strings.CutPrefix(line, "# "); ok {
			version = v
			continue
		}
		fld := strings.Fields(line)
		if len(fld) != 4 {
			return nil, "", fmt.Errorf("bad line %q", line)
		}
		p, err := strconv.ParseUint(fld[0], 10, 32)
		mant, ok := new(big.Int).SetString(fld[2], 16)
		e, err2 := strconv.Atoi(fld[3])
		if err != nil || !ok || err2 != nil || p == 0 {
			return nil, "", fmt.Errorf("bad line %q", line)
		}
		x := new(big.Float).SetPrec(uint(p)).SetInt(mant)
		if x.Acc() != big.Exact || mant.BitLen() != int(p) {
			return nil, "", fmt.Errorf("mantissa of %q is not %d bits", line, p)
		}
		x.SetMantExp(x, e-int(p))
		if rows[uint(p)] == nil {
			rows[uint(p)] = row{}
		}
		rows[uint(p)][fld[1]] = x
	}
	if version == "" {
		return nil, "", errors.New("no MPFR version line: output truncated")
	}
	return rows, version, nil
}

// buildCase assembles the case at precision p, checking what the fixture
// takes for granted: G is inexact (toward zero differs from away), nearest
// is one of the two, down is toward zero and up is away (G > 0), and G is
// no tie (the p+1-bit roundings differ).
func buildCase(p uint, rows map[uint]row) (oracle.Case, error) {
	r, next := rows[p], rows[p+1]
	for _, m := range []string{"N", "Z", "A", "D", "U"} {
		if r[m] == nil || next[m] == nil {
			return oracle.Case{}, fmt.Errorf("precision %d: no %s result", p, m)
		}
	}
	z, a := r["Z"], r["A"]
	switch {
	case z.Cmp(a) >= 0:
		return oracle.Case{}, fmt.Errorf("precision %d: toward zero %v >= away %v", p, z, a)
	case r["D"].Cmp(z) != 0 || r["U"].Cmp(a) != 0:
		return oracle.Case{}, fmt.Errorf("precision %d: down/up differ from toward/away for G > 0", p)
	case r["N"].Cmp(z) != 0 && r["N"].Cmp(a) != 0:
		return oracle.Case{}, fmt.Errorf("precision %d: nearest %v is neither neighbour", p, r["N"])
	case next["Z"].Cmp(next["A"]) == 0:
		return oracle.Case{}, fmt.Errorf("precision %d: G has %d bits, a tie is possible", p, p+1)
	}
	res := make([]oracle.Result, len(oracle.Modes))
	for i, m := range oracle.Modes {
		var v *big.Float
		switch m {
		case big.ToNearestEven, big.ToNearestAway:
			v = r["N"]
		case big.ToZero:
			v = z
		case big.AwayFromZero:
			v = a
		case big.ToNegativeInf:
			v = r["D"]
		default:
			v = r["U"]
		}
		acc := big.Below
		if v.Cmp(a) == 0 {
			acc = big.Above
		}
		res[i] = oracle.Result{Value: v, Acc: acc}
	}
	return oracle.NewCase("const", p, nil, res), nil
}
