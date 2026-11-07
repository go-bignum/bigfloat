// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package main

import (
	"bytes"
	"errors"
	"fmt"
	"math/big"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/go-bignum/bigfloat"
	"github.com/go-bignum/bigfloat/oracle"
)

// With FAKE_JULIA set, the test binary stands for pkgx: it prints a version
// line ("ok") or fails ("fail").
func TestMain(m *testing.M) {
	switch os.Getenv("FAKE_JULIA") {
	case "ok":
		fmt.Println("# 4.2.2")
		os.Exit(0)
	case "fail":
		os.Exit(1)
	}
	os.Exit(m.Run())
}

// fakeJulia prints what catalan.jl prints, computed with bigfloat.Catalan.
func fakeJulia(precs []uint) ([]byte, error) {
	var b strings.Builder
	for _, p := range precs {
		for _, m := range []struct {
			name string
			mode big.RoundingMode
		}{{"N", big.ToNearestEven}, {"Z", big.ToZero}, {"A", big.AwayFromZero}, {"D", big.ToNegativeInf}, {"U", big.ToPositiveInf}} {
			x := bigfloat.Catalan(new(big.Float).SetPrec(p).SetMode(m.mode))
			frac := new(big.Float)
			e := x.MantExp(frac)
			mant, _ := new(big.Float).SetMantExp(frac, int(p)).Int(nil)
			fmt.Fprintf(&b, "%d %s %s %d\n", p, m.name, mant.Text(16), e)
		}
	}
	b.WriteString("# 4.2.2\n")
	return []byte(b.String()), nil
}

func withJulia(t *testing.T, f func([]uint) ([]byte, error)) {
	old := julia
	julia = f
	t.Cleanup(func() { julia = old })
}

func TestRun(t *testing.T) {
	withJulia(t, fakeJulia)
	dir := t.TempDir()
	var out, errb bytes.Buffer
	if st := run([]string{"-out", dir}, &out, &errb); st != 0 {
		t.Fatalf("status %d: %s", st, errb.String())
	}
	if !strings.Contains(out.String(), "1002 cases (MPFR 4.2.2") {
		t.Errorf("output %q", out.String())
	}
	f, err := oracle.LoadFixture(filepath.Join(dir, "catalan.json.gz"))
	if err != nil {
		t.Fatal(err)
	}
	if f.Fn != "catalan" || len(f.Cases) != 1002 || f.Cases[0].Prec != 1 || f.Cases[1001].Prec != 10000 {
		t.Errorf("fixture %s: %d cases", f.Fn, len(f.Cases))
	}
	for i := range f.Cases {
		if bad := f.Cases[i].Check("catalan", func(z *big.Float, a ...*big.Float) *big.Float { return bigfloat.Catalan(z) }); bad != nil {
			t.Fatal(bad)
		}
	}
}

func TestRunFailures(t *testing.T) {
	var out, errb bytes.Buffer
	if st := run([]string{"-nope"}, &out, &errb); st != 2 {
		t.Errorf("bad flag: status %d", st)
	}
	withJulia(t, func([]uint) ([]byte, error) { return nil, errors.New("no julia") })
	if st := run([]string{"-out", t.TempDir()}, &out, &errb); st != 1 || !strings.Contains(errb.String(), "no julia") {
		t.Errorf("julia failure: status %d: %s", st, errb.String())
	}
	withJulia(t, fakeJulia)
	file := filepath.Join(t.TempDir(), "file")
	if err := os.WriteFile(file, nil, 0o644); err != nil {
		t.Fatal(err)
	}
	if st := run([]string{"-out", filepath.Join(file, "sub")}, &out, &errb); st != 1 {
		t.Errorf("output below a file: status %d", st)
	}
	dir := t.TempDir()
	if err := os.Mkdir(filepath.Join(dir, "catalan.json.gz"), 0o755); err != nil {
		t.Fatal(err)
	}
	if st := run([]string{"-out", dir}, &out, &errb); st != 1 {
		t.Errorf("output file is a directory: status %d", st)
	}
}

func TestMain_(t *testing.T) {
	var got int
	old := exit
	exit = func(c int) { got = c }
	defer func() { exit = old }()
	main() // the test binary's own flags are not ours: status 2
	if got != 2 {
		t.Errorf("status %d", got)
	}
}

func TestJulia(t *testing.T) {
	old := pkgx
	pkgx = os.Args[0]
	defer func() { pkgx = old }()
	t.Setenv("FAKE_JULIA", "ok")
	out, err := julia([]uint{1, 2})
	if err != nil || string(out) != "# 4.2.2\n" {
		t.Errorf("got %q, %v", out, err)
	}
	t.Setenv("FAKE_JULIA", "fail")
	if _, err := julia([]uint{1}); err == nil {
		t.Error("a failing julia must be an error")
	}
}

func TestParseErrors(t *testing.T) {
	for _, in := range []string{
		"1 N 1 1\n",            // no version line
		"1 N 1\n# 4.2.2\n",     // too few fields
		"x N 1 1\n# 4.2.2\n",   // precision
		"0 N 1 1\n# 4.2.2\n",   // zero precision
		"1 N zz 1\n# 4.2.2\n",  // mantissa
		"1 N 1 x\n# 4.2.2\n",   // exponent
		"1 N 3 1\n# 4.2.2\n",   // mantissa longer than the precision
		"2 N 1 1\n# 4.2.2\n",   // mantissa shorter than the precision
		"# \n1 N 1 1\n",        // empty version
		"1 N 1 1 extra\n# v\n", // too many fields
	} {
		if _, _, err := parse(strings.NewReader(in)); err == nil {
			t.Errorf("%q: no error", in)
		}
	}
	rows, v, err := parse(strings.NewReader("2 N 3 0\n# 4.2.2\n"))
	if err != nil || v != "4.2.2" || rows[2]["N"].Cmp(big.NewFloat(0.75)) != 0 {
		t.Errorf("valid input: %v %q %v", rows, v, err)
	}
}

func TestBuildCaseErrors(t *testing.T) {
	good := func() map[uint]row {
		out, _ := fakeJulia([]uint{5, 6})
		rows, _, err := parse(bytes.NewReader(out))
		if err != nil {
			t.Fatal(err)
		}
		return rows
	}
	if _, err := buildCase(5, good()); err != nil {
		t.Fatal(err)
	}
	for name, mutate := range map[string]func(map[uint]row){
		"missing mode":      func(r map[uint]row) { delete(r[5], "U") },
		"missing next":      func(r map[uint]row) { delete(r, 6) },
		"zero >= away":      func(r map[uint]row) { r[5]["A"] = r[5]["Z"] },
		"down != zero":      func(r map[uint]row) { r[5]["D"] = r[5]["A"] },
		"up != away":        func(r map[uint]row) { r[5]["U"] = r[5]["Z"] },
		"nearest elsewhere": func(r map[uint]row) { r[5]["N"] = big.NewFloat(0.5) },
		"tie possible":      func(r map[uint]row) { r[6]["A"] = r[6]["Z"] },
	} {
		rows := good()
		mutate(rows)
		if _, err := buildCase(5, rows); err == nil {
			t.Errorf("%s: no error", name)
		}
	}
}

// Rows from julia that lack a precision the fixture asks for.
func TestGenerateMissing(t *testing.T) {
	withJulia(t, func([]uint) ([]byte, error) { return []byte("# 4.2.2\n"), nil })
	if _, err := generate([]uint{3}); err == nil {
		t.Error("no rows: no error")
	}
	withJulia(t, func([]uint) ([]byte, error) { return []byte("garbage\n"), nil })
	if _, err := generate([]uint{3}); err == nil {
		t.Error("unparsable output: no error")
	}
}
