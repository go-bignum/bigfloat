// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package oracle

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"math/big"
	"os/exec"
	"sync"
)

// MpmathVersion is the mpmath release the default command pins.
const MpmathVersion = "1.4.1"

// MpmathScript is the whole of the Python side: it reads one JSON request
// per line — {"f": name, "p": bits, "a": [[mantissa, exponent], ...]} —
// sets mpmath's working precision to the larger of "p" and the arguments'
// own lengths (an mpf is rounded to the working precision when it is
// built, so a shorter one would round the arguments), evaluates, and
// writes {"s": sign, "m": mantissa, "e": exponent} (or {"i": ±1} for an
// infinity, {"x": message} for an error). Mantissas are decimal integer
// strings, so values cross the pipe exactly. Odd roots of negative numbers
// are taken as real (mpmath returns the principal complex root).
const MpmathScript = `import sys,json,mpmath as M
from mpmath import mp,mpf
sys.set_int_max_str_digits(0)
r1=lambda f:(lambda x,*a:-f(-x,*a) if x<0 else f(x,*a))
F=dict(sqrt=M.sqrt,cbrt=r1(M.cbrt),root=r1(lambda x,n:M.root(x,int(n))),hypot=M.hypot,exp=M.exp,expm1=M.expm1,log=M.log,log1p=M.log1p,log2=lambda x:M.log(x,2),log10=M.log10,pow=M.power,sin=M.sin,cos=M.cos,tan=M.tan,asin=M.asin,acos=M.acos,atan=M.atan,atan2=M.atan2,sinh=M.sinh,cosh=M.cosh,tanh=M.tanh,asinh=M.asinh,acosh=M.acosh,atanh=M.atanh,pi=lambda:+mp.pi,e=lambda:+mp.e,ln2=lambda:+mp.ln2,ln10=lambda:+mp.ln10,euler=lambda:+mp.euler)
for l in sys.stdin:
  q=json.loads(l);a=[(int(m),e) for m,e in q["a"]]
  mp.prec=max([q["p"]]+[abs(m).bit_length()+1 for m,e in a])
  try:
    v=F[q["f"]](*[mpf(t) for t in a])
    if M.isinf(v):o={"i":1 if v>0 else -1}
    else:s,m,e,b=v._mpf_;o={"s":s,"m":str(m),"e":e}
  except Exception as x:o={"x":repr(x)}
  print(json.dumps(o),flush=True)
`

// DefaultMpmathCommand runs MpmathScript under uv with mpmath pinned.
func DefaultMpmathCommand() []string {
	return []string{"uv", "run", "--quiet", "--with", "mpmath==" + MpmathVersion, "python", "-c", MpmathScript}
}

// Mpmath is an Oracle backed by an mpmath subprocess. It is safe for
// concurrent use (requests are serialised on the pipe).
type Mpmath struct {
	mu   sync.Mutex
	cmd  *exec.Cmd
	in   io.WriteCloser
	out  *bufio.Reader
	name string
}

// StartMpmath starts the oracle process with the given command (nil means
// DefaultMpmathCommand). The command must speak MpmathScript's protocol.
func StartMpmath(command []string) (*Mpmath, error) {
	if command == nil {
		command = DefaultMpmathCommand()
	}
	return start(exec.Command(command[0], command[1:]...))
}

func start(cmd *exec.Cmd) (*Mpmath, error) {
	in, err := cmd.StdinPipe()
	if err != nil {
		return nil, err
	}
	out, err := cmd.StdoutPipe()
	if err != nil {
		return nil, err
	}
	if err := cmd.Start(); err != nil {
		return nil, fmt.Errorf("oracle: starting mpmath: %w", err)
	}
	return &Mpmath{cmd: cmd, in: in, out: bufio.NewReaderSize(out, 1<<20), name: "mpmath " + MpmathVersion}, nil
}

// Name implements Oracle.
func (o *Mpmath) Name() string { return o.name }

// Slack implements Oracle: mpmath does not promise correct rounding, but
// its elementary functions carry guard bits internally; 8 untrusted bits
// leave a wide margin, which only costs an occasional retry at a higher
// precision.
func (o *Mpmath) Slack() uint { return 8 }

type request struct {
	F string   `json:"f"`
	P uint     `json:"p"`
	A [][2]any `json:"a"`
}

type response struct {
	S int    `json:"s"`
	M string `json:"m"`
	E int64  `json:"e"`
	I int    `json:"i"`
	X string `json:"x"`
}

// Eval implements Oracle.
func (o *Mpmath) Eval(fn string, prec uint, args ...*big.Float) (*big.Float, error) {
	req := request{F: fn, P: prec, A: make([][2]any, len(args))}
	for i, a := range args {
		if a.IsInf() {
			return nil, fmt.Errorf("oracle: infinite argument to %s", fn)
		}
		m, e := ToMantExp(a)
		req.A[i] = [2]any{json.Number(m.String()), e}
	}
	// strings, an integer and json.Numbers built from big.Int.String:
	// encoding cannot fail
	line, _ := json.Marshal(req)
	o.mu.Lock()
	defer o.mu.Unlock()
	if _, err := o.in.Write(append(line, '\n')); err != nil {
		return nil, fmt.Errorf("oracle: writing to mpmath: %w", err)
	}
	reply, err := o.out.ReadBytes('\n')
	if err != nil {
		return nil, fmt.Errorf("oracle: reading from mpmath: %w", err)
	}
	var r response
	if err := json.Unmarshal(reply, &r); err != nil {
		return nil, fmt.Errorf("oracle: bad reply %q: %w", reply, err)
	}
	return r.value(fn)
}

func (r *response) value(fn string) (*big.Float, error) {
	switch {
	case r.X != "":
		return nil, fmt.Errorf("oracle: mpmath %s: %s", fn, r.X)
	case r.I != 0:
		return new(big.Float).SetInf(r.I < 0), nil
	}
	m, ok := new(big.Int).SetString(r.M, 10)
	if !ok {
		return nil, fmt.Errorf("oracle: bad mantissa %q", r.M)
	}
	if r.S != 0 {
		m.Neg(m)
	}
	return FromMantExp(m, r.E), nil
}

// Close ends the subprocess.
func (o *Mpmath) Close() error {
	o.in.Close()
	return o.cmd.Wait()
}

// ToMantExp returns m and e with x = m·2^e exactly (x finite; m = 0 for a
// zero, whose sign is lost).
func ToMantExp(x *big.Float) (*big.Int, int64) {
	if x.Sign() == 0 {
		return new(big.Int), 0
	}
	e := x.MantExp(nil)
	p := int(x.Prec())
	t := new(big.Float).SetMantExp(x, p-e)
	m, _ := t.Int(nil)
	return m, int64(e - p)
}

// FromMantExp returns m·2^e exactly. Exponents beyond big.Float's range
// give ±Inf or ±0, as big.Float's own arithmetic does.
func FromMantExp(m *big.Int, e int64) *big.Float {
	v := new(big.Float).SetInt(m)
	if v.Sign() == 0 {
		return v
	}
	switch {
	case e+int64(m.BitLen()) > big.MaxExp:
		return v.SetInf(m.Sign() < 0)
	case e+int64(m.BitLen()) < big.MinExp:
		return v.SetInt64(0)
	}
	// v = m·2^-bitlen (in [1/2, 1)), then the remaining exponent, which
	// fits an int on 32-bit targets.
	n := m.BitLen()
	v.SetMantExp(v, -n)
	return v.SetMantExp(v, int(e+int64(n)))
}
