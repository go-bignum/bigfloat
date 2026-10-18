// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package oracle

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/json"
	"errors"
	"fmt"
	"math/big"
	"math/bits"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// The test binary doubles as a fake mpmath: with FAKE_MPMATH set it speaks
// MpmathScript's protocol, answering from Self (or misbehaving, as the
// variable's value asks), so that the pipe client is tested without Python.
func TestMain(m *testing.M) {
	if mode := os.Getenv("FAKE_MPMATH"); mode != "" {
		fakeMpmath(mode)
		os.Exit(0)
	}
	os.Exit(m.Run())
}

func fakeMpmath(mode string) {
	in := bufio.NewScanner(os.Stdin)
	in.Buffer(make([]byte, 1<<20), 1<<26)
	for in.Scan() {
		var q struct {
			F string
			P uint
			A [][2]json.Number
		}
		if err := json.Unmarshal(in.Bytes(), &q); err != nil {
			fmt.Println(`{"x":"bad request"}`)
			continue
		}
		args := make([]*big.Float, len(q.A))
		for i, a := range q.A {
			m, _ := new(big.Int).SetString(a[0].String(), 10)
			e, _ := a[1].Int64()
			args[i] = FromMantExp(m, e)
		}
		switch mode {
		case "die":
			return
		case "garbage":
			fmt.Println("not json")
			continue
		case "badmant":
			fmt.Println(`{"s":0,"m":"12x","e":0}`)
			continue
		case "error":
			fmt.Println(`{"x":"ValueError('math domain error')"}`)
			continue
		case "inf":
			fmt.Println(`{"i":-1}`)
			continue
		}
		v, err := Self{}.Eval(q.F, q.P, args...)
		if err != nil {
			fmt.Printf("{\"x\":%q}\n", err.Error())
			continue
		}
		m, e := ToMantExp(v)
		s := 0
		if m.Sign() < 0 {
			s = 1
			m.Neg(m)
		}
		fmt.Printf("{\"s\":%d,\"m\":%q,\"e\":%d}\n", s, m.String(), e)
	}
}

func startFake(t *testing.T, mode string) *Mpmath {
	t.Helper()
	t.Setenv("FAKE_MPMATH", mode)
	o, err := StartMpmath([]string{os.Args[0], "-test.run=^$"})
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { o.Close() })
	return o
}

func fl(s string) *big.Float {
	x, _, err := big.ParseFloat(s, 0, 200, big.ToNearestEven)
	if err != nil {
		panic(err)
	}
	return x
}

func TestMpmathPipe(t *testing.T) {
	o := startFake(t, "self")
	if o.Name() != "mpmath "+MpmathVersion || o.Slack() != 8 {
		t.Errorf("name %q slack %d", o.Name(), o.Slack())
	}
	v, err := o.Eval("exp", 100, fl("-1.5"))
	if err != nil {
		t.Fatal(err)
	}
	want, _ := Self{}.Eval("exp", 100, fl("-1.5"))
	if v.Cmp(want) != 0 {
		t.Errorf("exp(-1.5) over the pipe = %v, want %v", v, want)
	}
	if v, err := o.Eval("sin", 100, fl("-2")); err != nil || v.Sign() >= 0 {
		t.Errorf("sin(-2) = %v, %v", v, err)
	}
	if _, err := o.Eval("exp", 100, new(big.Float).SetInf(false)); err == nil {
		t.Error("infinite argument accepted")
	}
	// Self reports a NaN panic as an error, relayed as {"x": ...}
	if _, err := o.Eval("log", 100, fl("-1")); err == nil || !strings.Contains(err.Error(), "Log of a negative") {
		t.Errorf("log(-1): %v", err)
	}
	r := &Rounder{Oracle: o}
	res, err := r.Round("atan", 53, fl("0.75"))
	if err != nil || len(res) != len(Modes) {
		t.Fatalf("Round over the pipe: %v %v", res, err)
	}
}

func TestMpmathFailures(t *testing.T) {
	for mode, want := range map[string]string{
		"garbage": "bad reply",
		"badmant": "bad mantissa",
		"error":   "math domain error",
		"die":     "reading from mpmath",
	} {
		o := startFake(t, mode)
		_, err := o.Eval("exp", 64, fl("1"))
		if err == nil || !strings.Contains(err.Error(), want) {
			t.Errorf("%s: err %v, want %q", mode, err, want)
		}
	}
	o := startFake(t, "inf")
	if v, err := o.Eval("exp", 64, fl("1")); err != nil || !v.IsInf() || !v.Signbit() {
		t.Errorf("inf reply: %v %v", v, err)
	}
	// writing after the process has gone
	o = startFake(t, "die")
	o.Eval("exp", 64, fl("1"))
	o.in.Close()
	if _, err := o.Eval("exp", 64, fl("1")); err == nil || !strings.Contains(err.Error(), "writing to mpmath") {
		t.Errorf("write to closed pipe: %v", err)
	}
	if _, err := StartMpmath([]string{filepath.Join(t.TempDir(), "no-such-binary")}); err == nil {
		t.Error("missing binary started")
	}
	c := exec.Command("true")
	c.Stdin = strings.NewReader("")
	if _, err := start(c); err == nil {
		t.Error("start with Stdin already set succeeded")
	}
	c = exec.Command("true")
	c.Stdout = new(bytes.Buffer)
	if _, err := start(c); err == nil {
		t.Error("start with Stdout already set succeeded")
	}
	if c := DefaultMpmathCommand(); c[0] != "uv" || !strings.Contains(strings.Join(c, " "), "mpmath=="+MpmathVersion) {
		t.Errorf("default command %v", c)
	}
}

// TestMpmathReal runs the real oracle when uv is installed (CI need not).
func TestMpmathReal(t *testing.T) {
	if _, err := exec.LookPath("uv"); err != nil || testing.Short() {
		t.Skip("uv not installed")
	}
	o, err := StartMpmath(nil)
	if err != nil {
		t.Skip(err)
	}
	defer o.Close()
	r := &Rounder{Oracle: o}
	res, err := r.Round("sin", 53, fl("1e300"))
	if err != nil {
		t.Fatal(err)
	}
	self, _ := (&Rounder{Oracle: Self{}}).Round("sin", 53, fl("1e300"))
	for i := range res {
		if res[i].Value.Cmp(self[i].Value) != 0 || res[i].Acc != self[i].Acc {
			t.Errorf("sin(1e300) mode %v: mpmath %v, self %v", Modes[i], res[i], self[i])
		}
	}
}

// scripted answers a or b by turns (by the parity of the precision's bit
// length), whatever it is asked.
type scripted struct{ a, b float64 }

func (scripted) Name() string { return "scripted" }
func (scripted) Slack() uint  { return 1 }
func (s scripted) Eval(fn string, prec uint, args ...*big.Float) (*big.Float, error) {
	if bits.Len(prec)%2 == 0 {
		return big.NewFloat(s.a), nil
	}
	return big.NewFloat(s.b), nil
}

// liar's answers each decide a rounding, but no two successive ones agree.
var liar = scripted{1.25 + 0x1p-40, 1.5 + 0x1p-40}

// midpoint always answers a 24-bit rounding midpoint, which never decides.
var midpoint = scripted{1 + 0x1p-24, 1 + 0x1p-24}

type failing struct{ scripted }

func (failing) Eval(string, uint, ...*big.Float) (*big.Float, error) {
	return nil, errors.New("boom")
}

func TestRounder(t *testing.T) {
	r := &Rounder{Oracle: Self{}}
	if Self.Name(Self{}) != "self" || (Self{}).Slack() != 1 {
		t.Error("Self identity")
	}
	res, err := r.Round("pi", 24)
	if err != nil || Format(res[0].Value) != "0x.c90fdbp+2" || res[0].Acc != big.Above {
		t.Errorf("pi at 24 bits: %v %v", res, err)
	}
	if res, err := r.Round("sqrt", 24, fl("4")); err != nil || res[0].Acc != big.Exact {
		t.Errorf("exact sqrt: %v %v", res, err)
	}
	if _, err := r.Round("nope", 24); err == nil {
		t.Error("unknown function accepted")
	}
	if _, err := r.Round("exp", 24); err == nil {
		t.Error("wrong arity accepted")
	}
	if _, err := (&Rounder{Oracle: failing{}}).Round("exp", 24, fl("1")); err == nil {
		t.Error("oracle error not returned")
	}
	// a liar decides differently at every precision: never confirmed
	if _, err := (&Rounder{Oracle: liar, MaxPrec: 4000}).Round("exp", 24, fl("1")); !errors.Is(err, ErrUndecided) {
		t.Errorf("liar: %v, want ErrUndecided", err)
	}
	if _, err := (&Rounder{Oracle: midpoint, MaxPrec: 4000}).Round("exp", 24, fl("1")); !errors.Is(err, ErrUndecided) {
		t.Errorf("midpoint: %v, want ErrUndecided", err)
	}
	// an infinite answer (overflow) never decides
	if _, err := (&Rounder{Oracle: Self{}, MaxPrec: 500}).Round("exp", 24, fl("1e10")); !errors.Is(err, ErrUndecided) {
		t.Errorf("overflow: %v, want ErrUndecided", err)
	}
	// a value straddling a breakpoint at the oracle's precision: exp of a
	// tiny argument is 1 + x and cannot be decided below about 2·|log2 x|
	if _, err := (&Rounder{Oracle: Self{}, MaxPrec: 200}).Round("exp", 24, fl("0x1p-150")); !errors.Is(err, ErrUndecided) {
		t.Errorf("tiny exp at a low limit: %v, want ErrUndecided", err)
	}
}

func TestSelfEvalAndCall(t *testing.T) {
	if _, err := (Self{}).Eval("sqrt", 53, fl("-1")); err == nil {
		t.Error("Self.Eval did not report the NaN panic")
	}
	for _, fn := range []string{"root", "hypot", "pow", "atan2"} {
		v, err := (Self{}).Eval(fn, 53, fl("2"), fl("3"))
		if err != nil || v.Sign() <= 0 {
			t.Errorf("%s(2, 3) = %v %v", fn, v, err)
		}
	}
	defer func() {
		if recover() == nil {
			t.Error("Call of an unknown function did not panic")
		}
	}()
	Call("nope", new(big.Float))
}

func TestExact(t *testing.T) {
	cases := []struct {
		fn   string
		args []string
		want string // "" for not exact
	}{
		{"exp", []string{"0"}, "0x.8p+1"}, {"cosh", []string{"1"}, ""},
		{"sin", []string{"0"}, "0"}, {"atanh", []string{"0.5"}, ""},
		{"log", []string{"1"}, "0"}, {"acos", []string{"0.5"}, ""},
		{"log2", []string{"0x1p-30"}, "-0x.fp+5"}, {"log2", []string{"3"}, ""}, {"log2", []string{"-4"}, ""},
		{"log10", []string{"1000"}, "0x.cp+2"}, {"log10", []string{"1001"}, ""}, {"log10", []string{"0.5"}, ""}, {"log10", []string{"20"}, ""},
		{"sqrt", []string{"0x1.9p+10"}, "0x.ap+6"}, {"sqrt", []string{"2"}, ""}, {"sqrt", []string{"3"}, ""}, {"sqrt", []string{"0x1p-3"}, ""}, {"sqrt", []string{"0"}, "0"},
		{"cbrt", []string{"-27"}, "-0x.cp+2"},
		{"root", []string{"0x1p100", "4"}, "0x.8p+26"},
		{"hypot", []string{"3", "-4"}, "0x.ap+3"}, {"hypot", []string{"0", "0"}, "0"}, {"hypot", []string{"1", "1"}, ""},
		{"pow", []string{"5", "0"}, "0x.8p+1"}, {"pow", []string{"1", "7.5"}, "0x.8p+1"},
		{"pow", []string{"-2", "3"}, "-0x.8p+4"}, {"pow", []string{"-2", "2"}, "0x.8p+3"}, {"pow", []string{"-2", "0.5"}, ""},
		{"pow", []string{"0", "2"}, ""}, {"pow", []string{"9", "1.5"}, "0x.d8p+5"}, {"pow", []string{"8", "0.5"}, ""},
		{"pow", []string{"4", "-1.5"}, "0x.8p-2"}, {"pow", []string{"3", "-2"}, ""},
		{"pow", []string{"2", "0x1p-70"}, ""}, {"pow", []string{"2", "0x1p20"}, ""},
		{"pi", nil, ""},
	}
	for _, c := range cases {
		args := make([]*big.Float, len(c.args))
		for i, a := range c.args {
			args[i] = fl(a)
		}
		v, ok := Exact(c.fn, args...)
		got := ""
		if ok {
			got = Format(v)
		}
		if got != c.want {
			t.Errorf("Exact(%s, %v) = %q, want %q", c.fn, c.args, got, c.want)
		}
	}
	if m, e := oddPart(new(big.Float)); m.Sign() != 0 || e != 0 {
		t.Error("oddPart(0)")
	}
	if expGap(new(big.Float), fl("8")) != 0 || expGap(fl("1"), fl("8")) != 3 {
		t.Error("expGap")
	}
}

func TestFixtureFormat(t *testing.T) {
	for _, s := range []string{"0", "-0", "+Inf", "-Inf", "0x.8p+1", "-0x.c90fdbp+2"} {
		x, err := Parse(s, 24)
		if err != nil || Format(x) != s {
			t.Errorf("round trip %q: %v %v", s, x, err)
		}
	}
	if _, err := Parse("0xg", 24); err == nil {
		t.Error("bad number parsed")
	}
	for a, c := range map[big.Accuracy]byte{big.Below: '<', big.Exact: '=', big.Above: '>'} {
		if AccChar(a) != c {
			t.Errorf("AccChar(%v)", a)
		}
	}
	res, _ := (&Rounder{Oracle: Self{}}).Round("exp", 24, fl("1"))
	c := NewCase("r", 24, []*big.Float{fl("1")}, res)
	f := &Fixture{Fn: "exp", Oracle: "self", Seed: 1, Cases: []Case{c}}
	var buf bytes.Buffer
	if err := WriteFixture(&buf, f); err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(t.TempDir(), "exp.json.gz")
	if err := os.WriteFile(path, buf.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
	g, err := LoadFixture(path)
	if err != nil || len(g.Cases) != 1 || g.Cases[0].Acc != c.Acc {
		t.Fatalf("load: %v %v", g, err)
	}
	if bad := g.Cases[0].Check("exp", func(z *big.Float, a ...*big.Float) *big.Float { return Call("exp", z, a...) }); bad != nil {
		t.Errorf("check of a correct case: %v", bad)
	}
	if bad := g.Cases[0].Check("exp", func(z *big.Float, a ...*big.Float) *big.Float { return z.SetInt64(3) }); len(bad) != len(Modes) {
		t.Errorf("check of a wrong evaluation: %d complaints", len(bad))
	}
	g.Cases[0].Args[0] = "0xg"
	if bad := g.Cases[0].Check("exp", nil); len(bad) != 1 {
		t.Errorf("check of a bad argument: %v", bad)
	}
	if _, err := LoadFixture(filepath.Join(t.TempDir(), "missing")); err == nil {
		t.Error("missing fixture loaded")
	}
	if _, err := ReadFixture(strings.NewReader("not gzip")); err == nil {
		t.Error("non-gzip fixture read")
	}
	buf.Reset()
	zw := gzip.NewWriter(&buf)
	zw.Write([]byte("{not json"))
	zw.Close()
	if _, err := ReadFixture(&buf); err == nil {
		t.Error("bad JSON fixture read")
	}
	if err := WriteFixture(failWriter{}, f); err == nil {
		t.Error("write to a failing writer succeeded")
	}
}

type failWriter struct{}

func (failWriter) Write([]byte) (int, error) { return 0, errors.New("full") }

func TestMantExp(t *testing.T) {
	for _, s := range []string{"0", "1", "-0x.c90fdbp+2", "0x1p-3000"} {
		x := fl(s)
		m, e := ToMantExp(x)
		if y := FromMantExp(m, e); y.Cmp(x) != 0 {
			t.Errorf("%s: round trip gave %v", s, y)
		}
	}
	if !FromMantExp(big.NewInt(-1), 1<<40).IsInf() {
		t.Error("overflow is not Inf")
	}
	if FromMantExp(big.NewInt(1), -1<<40).Sign() != 0 {
		t.Error("underflow is not 0")
	}
}
