// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

// Command bench times every bigfloat function against mpmath at the same
// precision and writes the table of docs/perf.md.
//
//	go run ./cmd/bench -out docs/perf.md
//
// Go runs on one thread (GOMAXPROCS=1); mpmath is timed inside its own
// process (the pipe is not counted), the best of five loops. Rounds
// alternate Go and mpmath, and each side keeps its best round, so that a
// burst of load on a shared machine hits both and is filtered out. Run it
// on a quiet machine: the numbers are only as good as the host's idleness.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"math/big"
	"os"
	"runtime"
	"strconv"
	"strings"
	"time"

	"github.com/go-bignum/bigfloat/oracle"
)

// cases are the timed calls: representative arguments, exact in binary.
var cases = []struct {
	fn   string
	args []float64
}{
	{"sqrt", []float64{2}}, {"cbrt", []float64{2}}, {"root", []float64{2, 7}}, {"hypot", []float64{0.75, 1.25}},
	{"exp", []float64{0.7}}, {"expm1", []float64{0.7}}, {"log", []float64{1.7}}, {"log1p", []float64{0.7}},
	{"log2", []float64{1.7}}, {"log10", []float64{1.7}}, {"pow", []float64{1.7, 2.3}},
	{"sin", []float64{0.7}}, {"cos", []float64{0.7}}, {"tan", []float64{0.7}},
	{"asin", []float64{0.7}}, {"acos", []float64{0.7}}, {"atan", []float64{0.7}}, {"atan2", []float64{0.7, 1.3}},
	{"sinh", []float64{0.7}}, {"cosh", []float64{0.7}}, {"tanh", []float64{0.7}},
	{"asinh", []float64{0.7}}, {"acosh", []float64{1.7}}, {"atanh", []float64{0.7}},
	{"sin", []float64{1e22}}, {"exp", []float64{1000.5}},
}

// exit is os.Exit, replaceable by tests.
var exit = os.Exit

func main() { exit(run(os.Args[1:], os.Stdout, os.Stderr)) }

// timer is what bench needs of an mpmath oracle.
type timer interface {
	Time(fn string, prec uint, n int, args ...*big.Float) (float64, error)
}

func run(args []string, stdout, stderr io.Writer) int {
	fs := flag.NewFlagSet("bench", flag.ContinueOnError)
	fs.SetOutput(stderr)
	out := fs.String("out", "-", "output file (- for standard output)")
	precList := fs.String("prec", "53,113,200,1000", "comma-separated precisions")
	rounds := fs.Int("rounds", 3, "interleaved rounds")
	target := fs.Duration("target", 50*time.Millisecond, "duration of one timed loop")
	cmdJSON := fs.String("cmd", "", "oracle command as a JSON array (default: uv + mpmath)")
	if err := fs.Parse(args); err != nil {
		return 2
	}
	var precs []uint
	for _, s := range strings.Split(*precList, ",") {
		p, err := strconv.ParseUint(s, 10, 32)
		if err != nil || p < 2 {
			fmt.Fprintf(stderr, "bench: bad precision %q\n", s)
			return 2
		}
		precs = append(precs, uint(p))
	}
	var command []string
	if *cmdJSON != "" {
		if err := json.Unmarshal([]byte(*cmdJSON), &command); err != nil || len(command) == 0 {
			fmt.Fprintln(stderr, "bench: -cmd: want a JSON array of strings")
			return 2
		}
	}
	m, err := oracle.StartMpmath(command)
	if err != nil {
		fmt.Fprintln(stderr, err)
		return 1
	}
	defer m.Close()
	runtime.GOMAXPROCS(1)
	var b strings.Builder
	if err := table(&b, m, precs, *rounds, *target); err != nil {
		fmt.Fprintln(stderr, err)
		return 1
	}
	if *out == "-" {
		io.WriteString(stdout, b.String())
		return 0
	}
	if err := os.WriteFile(*out, []byte(b.String()), 0o644); err != nil {
		fmt.Fprintln(stderr, err)
		return 1
	}
	return 0
}

// table writes the Markdown comparison.
func table(w io.Writer, m timer, precs []uint, rounds int, target time.Duration) error {
	fmt.Fprintf(w, "| function | argument | bits | bigfloat | mpmath %s | ratio |\n|---|---|---:|---:|---:|---:|\n", oracle.MpmathVersion)
	for _, c := range cases {
		for _, p := range precs {
			args := make([]*big.Float, len(c.args))
			text := make([]string, len(c.args))
			for i, a := range c.args {
				args[i] = new(big.Float).SetFloat64(a)
				text[i] = strconv.FormatFloat(a, 'g', -1, 64)
			}
			g, py, err := measure(m, c.fn, p, args, rounds, target)
			if err != nil {
				return err
			}
			fmt.Fprintf(w, "| %s | %s | %d | %s | %s | %.2f |\n", c.fn, strings.Join(text, ", "), p, dur(g), dur(py), g/py)
		}
	}
	return nil
}

// measure returns the best Go and mpmath times per call over interleaved
// rounds.
func measure(m timer, fn string, p uint, args []*big.Float, rounds int, target time.Duration) (float64, float64, error) {
	call := func() { oracle.Call(fn, new(big.Float).SetPrec(p), args...) }
	n := 1
	for {
		t0 := time.Now()
		for i := 0; i < n; i++ {
			call()
		}
		if time.Since(t0) >= target/4 || n >= 1<<24 {
			n = max(1, int(float64(n)*float64(target)/float64(time.Since(t0)+1)))
			break
		}
		n *= 4
	}
	bestGo, bestPy := 1e9, 1e9
	for r := 0; r < rounds; r++ {
		t0 := time.Now()
		for i := 0; i < n; i++ {
			call()
		}
		bestGo = min(bestGo, time.Since(t0).Seconds()/float64(n))
		py, err := m.Time(fn, p, max(1, n/4), args...)
		if err != nil {
			return 0, 0, err
		}
		bestPy = min(bestPy, py)
	}
	return bestGo, bestPy, nil
}

// dur formats seconds per call.
func dur(s float64) string {
	switch {
	case s < 1e-6:
		return fmt.Sprintf("%.0f ns", s*1e9)
	case s < 1e-3:
		return fmt.Sprintf("%.1f µs", s*1e6)
	}
	return fmt.Sprintf("%.2f ms", s*1e3)
}
