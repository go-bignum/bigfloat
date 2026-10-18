// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"fmt"
	"math/big"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/go-bignum/bigfloat/oracle"
)

// With FAKE_ORACLE set, the test binary is an oracle process that answers
// every request with 1.5 + 2^-2999: wrong for every function, but
// decidable at every fixture precision, so that -check finds mismatches;
// with FAKE_ORACLE=error it reports an error for every request.
func TestMain(m *testing.M) {
	if mode := os.Getenv("FAKE_ORACLE"); mode != "" {
		in := bufio.NewScanner(os.Stdin)
		in.Buffer(make([]byte, 1<<20), 1<<26)
		man := new(big.Int).Lsh(big.NewInt(3), 2998)
		man.Add(man, big.NewInt(1))
		for in.Scan() {
			if mode == "error" {
				fmt.Println(`{"x":"no"}`)
				continue
			}
			fmt.Printf("{\"s\":0,\"m\":%q,\"e\":-2999}\n", man.String())
		}
		os.Exit(0)
	}
	os.Exit(m.Run())
}

func fakeCmd(t *testing.T, mode string) string {
	t.Setenv("FAKE_ORACLE", mode)
	b, _ := json.Marshal([]string{os.Args[0], "-test.run=^$"})
	return string(b)
}

func TestGenerateSelf(t *testing.T) {
	dir := t.TempDir()
	var out, errb bytes.Buffer
	if st := run([]string{"-self", "-out", dir, "-fn", "exp,pi,pow,root", "-scale", "0.02"}, &out, &errb); st != 0 {
		t.Fatalf("status %d: %s", st, errb.String())
	}
	for _, fn := range []string{"exp", "pi", "pow", "root"} {
		f, err := oracle.LoadFixture(filepath.Join(dir, fn+".json.gz"))
		if err != nil || f.Fn != fn || f.Oracle != "self" || len(f.Cases) == 0 {
			t.Errorf("%s: %v %v", fn, f, err)
		}
	}
	out.Reset()
	if st := run([]string{"-self", "-check", "-fn", "sin", "-scale", "0.02"}, &out, &errb); st != 0 || !strings.Contains(out.String(), "0 mismatches") {
		t.Errorf("check: %d %q", st, out.String())
	}
}

// TestEveryRegion draws one case from every region of every function, so
// that each generator is exercised (Self as the oracle).
func TestEveryRegion(t *testing.T) {
	rd := &oracle.Rounder{Oracle: oracle.Self{}}
	for fn := range oracle.Functions {
		if testing.Short() && fn != "exp" {
			continue
		}
		if _, err := generate(rd, fn, 1, 0.001); err != nil {
			t.Errorf("%s: %v", fn, err)
		}
	}
}

func TestCheckMismatches(t *testing.T) {
	var out, errb bytes.Buffer
	st := run([]string{"-cmd", fakeCmd(t, "lie"), "-check", "-fn", "atan", "-scale", "0.05"}, &out, &errb)
	if st != 1 || !strings.Contains(out.String(), "mismatches") || strings.Contains(out.String(), " 0 mismatches") {
		t.Errorf("status %d, output %q, stderr %q", st, out.String(), errb.String())
	}
}

func TestRunErrors(t *testing.T) {
	file := filepath.Join(t.TempDir(), "file")
	os.WriteFile(file, nil, 0o644)
	dirAsFile := t.TempDir()
	os.Mkdir(filepath.Join(dirAsFile, "pi.json.gz"), 0o755)
	cases := []struct {
		args []string
		want int
		msg  string
	}{
		{[]string{"-nope"}, 2, "flag provided but not defined"},
		{[]string{"-cmd", "{"}, 2, "JSON array"},
		{[]string{"-cmd", "[]"}, 2, "JSON array"},
		{[]string{"-cmd", `["` + filepath.Join(t.TempDir(), "missing") + `"]`}, 1, "starting mpmath"},
		{[]string{"-self", "-fn", "nope"}, 2, "unknown function"},
		{[]string{"-cmd", fakeCmd(t, "error"), "-fn", "exp", "-scale", "0.01"}, 1, "no"},
		{[]string{"-self", "-fn", "pi", "-out", filepath.Join(file, "sub")}, 1, ""},
		{[]string{"-self", "-fn", "pi", "-out", dirAsFile}, 1, ""},
	}
	for _, c := range cases {
		var out, errb bytes.Buffer
		if st := run(c.args, &out, &errb); st != c.want || !strings.Contains(errb.String(), c.msg) {
			t.Errorf("%v: status %d, stderr %q; want %d, %q", c.args, st, errb.String(), c.want, c.msg)
		}
	}
}

func TestMain_(t *testing.T) {
	saved, savedArgs := exit, os.Args
	defer func() { exit, os.Args = saved, savedArgs }()
	got := -1
	exit = func(c int) { got = c }
	os.Args = []string{"oracle", "-nope"}
	main()
	if got != 2 {
		t.Errorf("main exit status %d, want 2", got)
	}
}
