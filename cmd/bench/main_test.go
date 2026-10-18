// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// With FAKE_TIMER set the test binary answers timing requests with 2 µs
// (or an error, with FAKE_TIMER=error), until its standard input ends.
func TestMain(m *testing.M) {
	if mode := os.Getenv("FAKE_TIMER"); mode != "" {
		in := bufio.NewScanner(os.Stdin)
		for in.Scan() {
			if mode == "error" {
				fmt.Println(`{"x":"no"}`)
				continue
			}
			fmt.Println(`{"t":2e-06}`)
		}
		os.Exit(0)
	}
	os.Exit(m.Run())
}

func fake(t *testing.T, mode string) string {
	t.Setenv("FAKE_TIMER", mode)
	return jsonArray(os.Args[0], "-test.run=^$")
}

// jsonArray encodes a command line for -cmd (paths with backslashes on
// Windows must be escaped).
func jsonArray(s ...string) string {
	b, _ := json.Marshal(s)
	return string(b)
}

func TestTable(t *testing.T) {
	path := filepath.Join(t.TempDir(), "perf.md")
	var out, errb bytes.Buffer
	args := []string{"-cmd", fake(t, "ok"), "-prec", "53", "-rounds", "1", "-target", "1ms"}
	if st := run(append(args, "-out", path), &out, &errb); st != 0 {
		t.Fatalf("status %d: %s", st, errb.String())
	}
	data, _ := os.ReadFile(path)
	if !strings.Contains(string(data), "| exp | 0.7 | 53 |") || strings.Count(string(data), "\n") != len(cases)+2 {
		t.Errorf("table:\n%s", data)
	}
	if st := run(args, &out, &errb); st != 0 || !strings.Contains(out.String(), "| sqrt | 2 | 53 |") {
		t.Errorf("stdout: %d %q", st, out.String())
	}
}

func TestErrors(t *testing.T) {
	for _, c := range []struct {
		args []string
		want int
	}{
		{[]string{"-nope"}, 2},
		{[]string{"-prec", "1"}, 2},
		{[]string{"-prec", "x"}, 2},
		{[]string{"-cmd", "{"}, 2},
		{[]string{"-cmd", jsonArray(filepath.Join(t.TempDir(), "missing"))}, 1},
		{[]string{"-cmd", fake(t, "ok"), "-prec", "24", "-rounds", "1", "-target", "1ms", "-out", filepath.Join(t.TempDir(), "no", "dir")}, 1},
	} {
		var out, errb bytes.Buffer
		if st := run(c.args, &out, &errb); st != c.want {
			t.Errorf("%v: status %d, want %d (%s)", c.args, st, c.want, errb.String())
		}
	}
	// the fake's mode is read from the environment when it starts
	var out, errb bytes.Buffer
	if st := run([]string{"-cmd", fake(t, "error"), "-prec", "24", "-rounds", "1", "-target", "1ms"}, &out, &errb); st != 1 {
		t.Errorf("oracle error: status %d", st)
	}
}

func TestDur(t *testing.T) {
	for s, want := range map[float64]string{5e-8: "50 ns", 2.5e-5: "25.0 µs", 0.0125: "12.50 ms"} {
		if got := dur(s); got != want {
			t.Errorf("dur(%g) = %q, want %q", s, got, want)
		}
	}
}

func TestMainExit(t *testing.T) {
	saved, savedArgs := exit, os.Args
	defer func() { exit, os.Args = saved, savedArgs }()
	got := -1
	exit = func(c int) { got = c }
	os.Args = []string{"bench", "-nope"}
	main()
	if got != 2 {
		t.Errorf("exit %d", got)
	}
}
