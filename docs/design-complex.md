# Design: arbitrary-precision complex numbers

Status: design, no code yet. This page decides the complex type of
`go-bignum/bigfloat`: where it lives, its API on `big.Float`, its rounding
contract, the algorithms of its elementary functions, the oracles that check
them, and the order of the pull requests. Almost every missing mpmath row of
[coverage.md](coverage.md) waits on it. Its double-precision counterpart is
[go-numerics/special docs/design-complex.md](https://github.com/go-numerics/special/blob/main/docs/design-complex.md);
the two share Annex G special values and branch cuts, nothing else (different
modules, different contracts).

## Decisions

1. **A type `Complex` in package `bigfloat` itself**, with the
   `math/big` method style (`z.Exp(x)` computes into the receiver and returns
   it), not a new package and not a new module.
2. **Correct rounding per component, MPC semantics**: each part of the result
   is the correctly rounded value of the exact part, at that part's own
   precision and in that part's own rounding mode, in all six
   `big.RoundingMode`s; `z.Re.Acc()` and `z.Im.Acc()` report the direction of
   each rounding.
3. **Special values and branch cuts from C11 Annex G**, as MPC does
   ("our goal is to compute the same result as specified for the
   corresponding function in the C standard", MPC manual §4.7). A result that
   Annex G gives with a NaN part panics with `ErrNaN`, since `big.Float` has
   no NaN (the existing convention of this package).
4. **Rigour from complex balls** (a pair of the existing real `ball`s, Arb's
   rectangular `acb` form) and **one Ziv loop deciding both parts**, each part
   accepted as soon as its own enclosure rounds to one value.
5. **Two independent oracles**: mpmath (the existing harness, extended to
   complex arguments) through the `Rounder`, in all six modes; and MPC run as
   a black box through gmpy2 for bit equality in the five modes it has.

### Why in package `bigfloat`, as a type with methods

- The complex functions are compositions of **enclosures** of real functions
  (e^x, sin, cos, sinh, log1p, atan2 at a working precision, as balls), not of
  correctly rounded real results. Those enclosures, `ball`, `fx`, `mag` and
  `ziv` are unexported. A separate package (`bigfloat/bigcomplex`) would need
  them exported — a speculative public API of ~40 internal functions — or the
  whole package moved under `internal/` first, which is a rewrite of every
  file for no user-visible gain. A new module adds versioning on top.
- Methods avoid the name clash with the real package-level functions
  (`bigfloat.Exp(z, x *big.Float)` stays; the complex one is
  `(*Complex).Exp`) without `C` prefixes, and match `big.Float` itself
  (`z.Sqrt(x)`, `z.Quo(x, y)`) and MPC's "output first" convention (MPC §4.3).
- Rejected: functions on `[2]*big.Float` pairs (no place for per-part
  precision defaults, aliasing rules unclear); generics (one `big.Float`
  element type, nothing to parametrise).

## Bibliography

Opened for this design:

- GNU MPC manual, edition 1.4.1, April 2026
  (<https://www.multiprecision.org/downloads/mpc-1.4.1.pdf>, GFDL): §4.1
  (precision of each part independent, `mpc_init3`), §4.4 (a complex rounding
  mode is a pair of real modes `MPC_RNDxy`), §4.5 (semantics: inputs are exact
  dyadic numbers; each part rounded according to its mode, which "uniquely
  defines the expected output of any function"), §4.6 (return value:
  `MPC_INEX_RE`, `MPC_INEX_IM`), §4.7 (special values, signed zeros and branch
  cuts as in the C standard; (+0)+(-0) = +0; a real operand's zero imaginary
  part is an exact unsigned zero, so (1 - 0i) + (+0) = 1 - 0i), §4.8 (no
  exceptions), §5.6-5.9 (the function list; `mpc_pow` zero-sign rule;
  `mpc_sqrt` "non-negative real part, and if its real part is zero, a
  non-negative imaginary part"; `mpc_log` imaginary part in (-π, π];
  `mpc_acosh` cut (-∞, 1); `mpc_agm` normalisation; `mpc_rootofunity`).
  The MPC sources (LGPL) were not read.
- ISO/IEC 9899:201x N1570 (C11) Annex G: G.3, G.5.1 (products and quotients
  with infinities), G.6 (branch cuts continuous on both sides with the sign of
  zero; `casin(z) = -i casinh(iz)` etc.; symmetries) and the special-value
  lists G.6.1.1-G.6.4.2. HTML rendering at
  <https://port70.net/~nsz/c/c11/n1570.html#G>.
- This repository: `docs/design.md` (contract, balls, fixed point, Ziv, exact
  cases), `docs/coverage.md`, `ziv.go` (`evalFn`, `zivLimit`, `decide`),
  `ball.go`, `oracle/` (`MpmathScript`, `Rounder`, `Exact`, `Self`).
- NumPy's `npy_math_complex.c.src` (msun, BSD-2) and Go's `math/cmplx`, read
  for the double-precision design; their cancellation-free formulas for
  `clog` near |z| = 1 and the Hull–Fairgrieve–Tang decomposition of the
  inverse functions are reused below in ball form.

Not obtained: Kahan (1987) and Hull, Fairgrieve, Tang (ACM TOMS 23(3), 1997)
themselves (cited through msun); the Arb/FLINT `acb` documentation (the
rectangular-ball form is taken from Johansson 2017, already cited in
design.md, not from the acb manual); Enge, Théveny, Zimmermann's papers on
MPC's error analysis (the manual was enough for the semantics).

## Contract

```go
// A Complex is x + iy with x = Re and y = Im, each an arbitrary-precision
// big.Float carrying its own precision and rounding mode. The zero value is
// 0 + 0i with precision 0 in both parts, ready to use.
type Complex struct {
	Re, Im big.Float
}
```

- **Result**: for a function f, `z.F(x)` sets `z.Re` to Re f(x) rounded to
  `z.Re.Prec()` bits in `z.Re.Mode()`, and `z.Im` likewise with its own
  precision and mode; operands are exact (MPC §4.5). All six modes
  (`ToNearestEven`, `ToNearestAway`, `ToZero`, `AwayFromZero`,
  `ToNegativeInf`, `ToPositiveInf`); MPC has five (no `ToNearestAway`).
- **Precision 0** (as `math/big`): a part of the receiver with precision 0
  takes the largest precision among all operand parts (64 for constants).
- **Accuracy**: `z.Re.Acc()`, `z.Im.Acc()` are `Below`, `Exact` or `Above`,
  truthfully (the `MPC_INEX_RE/IM` pair).
- **Aliasing**: the receiver may be any operand.
- **Over/underflow**: per part, with `big.Float`'s own semantics (±Inf when
  the rounded exponent exceeds `MaxExp`, ±0 below `MinExp`), decided once at
  the final scaling of each part, as the real functions do.
- **Special values**: inputs with ±Inf or ±0 parts follow Annex G G.6 to the
  sign of each zero; entries the standard leaves unspecified follow MPC's
  choice where its manual states one (the sign of an undetermined zero is +,
  §4.7), otherwise +. A result with a NaN part panics `ErrNaN`.
- **Branch cuts** (Annex G, MPC §5.8-5.9): `Sqrt`, `Log*`, `Pow` on
  (-∞, 0]; `Asin`, `Acos`, `Atanh` on (-∞, -1] ∪ [1, ∞); `Asinh`, `Atan` on
  i(-∞, -1] ∪ i[1, ∞); `Acosh` on (-∞, 1). An input exactly on a cut takes
  the side given by the sign of its zero part (`Sqrt(-4 - 0i) = -2i`).
- **Resource bounds**: `MaxReductionExponent` bounds the real parts given to
  sin/cos reductions as it does for the real functions; every function panics
  `ErrLimit` rather than exceed it; the Ziv loop's `zivLimit` applies per part
  and its panic is an internal error (a missed exact case), never a user
  error.

## API

```go
// Setting and inspecting
func NewComplex(re, im float64) *Complex                 // exact, precision 53
func (z *Complex) Set(x *Complex) *Complex                // rounds to z's precisions
func (z *Complex) SetParts(re, im *big.Float) *Complex
func (z *Complex) SetFloat64s(re, im float64) *Complex
func (z *Complex) SetComplex128(x complex128) *Complex
func (z *Complex) SetPrec(prec uint) *Complex             // both parts
func (z *Complex) SetMode(mode big.RoundingMode) *Complex // both parts
func (x *Complex) Complex128() (complex128, big.Accuracy, big.Accuracy)
func (x *Complex) IsInf() bool                            // Annex G G.3: one infinite part
func (x *Complex) String() string                         // "(re im)", MPC §5.4 format
func (x *Complex) Text(format byte, prec int) string
func (z *Complex) SetString(s string) (*Complex, bool)    // "(re im)" or a real

// Exact or correctly rounded arithmetic
func (z *Complex) Neg(x *Complex) *Complex
func (z *Complex) Conj(x *Complex) *Complex
func (z *Complex) Proj(x *Complex) *Complex
func (z *Complex) MulI(x *Complex, sign int) *Complex     // ±i·x, exact
func (z *Complex) Add(x, y *Complex) *Complex
func (z *Complex) Sub(x, y *Complex) *Complex
func (z *Complex) Mul(x, y *Complex) *Complex
func (z *Complex) Sqr(x *Complex) *Complex
func (z *Complex) Quo(x, y *Complex) *Complex
func (z *Complex) AddFloat(x *Complex, y *big.Float) *Complex // y's imaginary part is an exact unsigned 0
func (z *Complex) MulFloat(x *Complex, y *big.Float) *Complex
func (z *Complex) QuoFloat(x *Complex, y *big.Float) *Complex
func (z *Complex) Mul2Exp(x *Complex, k int) *Complex
func Norm(z *big.Float, x *Complex) *big.Float            // x.Re² + x.Im², correctly rounded
// |x| and arg x need no new function: Hypot(z, &x.Re, &x.Im) and
// Atan2(z, &x.Im, &x.Re) are already correctly rounded with Annex F zeros.

// Elementary functions (MPC §5.8-5.9; mpmath's complex forms)
func (z *Complex) Sqrt(x *Complex) *Complex
func (z *Complex) Root(x *Complex, n uint, k uint) *Complex  // k-th n-th root, mpmath root(x, n, k)
func (z *Complex) Exp(x *Complex) *Complex
func (z *Complex) Expm1(x *Complex) *Complex
func (z *Complex) Log(x *Complex) *Complex
func (z *Complex) Log1p(x *Complex) *Complex
func (z *Complex) Log2(x *Complex) *Complex
func (z *Complex) Log10(x *Complex) *Complex
func (z *Complex) Pow(x, y *Complex) *Complex
func (z *Complex) PowFloat(x *Complex, y *big.Float) *Complex
func (z *Complex) PowInt(x *Complex, n int64) *Complex
func (z *Complex) Sin(x *Complex) *Complex
func (z *Complex) Cos(x *Complex) *Complex
func (z *Complex) Tan(x *Complex) *Complex
func SinCos(s, c *Complex, x *Complex)                      // one reduction, two results
func (z *Complex) Sinh(x *Complex) *Complex
func (z *Complex) Cosh(x *Complex) *Complex
func (z *Complex) Tanh(x *Complex) *Complex
func (z *Complex) Asin(x *Complex) *Complex
func (z *Complex) Acos(x *Complex) *Complex
func (z *Complex) Atan(x *Complex) *Complex
func (z *Complex) Asinh(x *Complex) *Complex
func (z *Complex) Acosh(x *Complex) *Complex
func (z *Complex) Atanh(x *Complex) *Complex
func (z *Complex) Expj(x *big.Float) *Complex                // e^(ix), mpmath expj
func (z *Complex) ExpjPi(x *big.Float) *Complex              // e^(iπx), exact at half-integers
func (z *Complex) RootOfUnity(n, k uint64) *Complex          // e^(2πik/n), mpmath unitroots, MPC rootofunity
```

Later slices add, on the same type: `Agm` (MPC normalisation), the gamma
family, zeta, `hyper` and everything coverage.md lists as "complex type".

## Algorithms

Internal layer (unexported, `cball.go`): `type cball struct{ re, im ball }`
with `add`, `sub`, `mul` (four ball products; the radius bound of a product
of rectangles), `sqr` (x² - y², 2xy: one product fewer), `quo`, `conj`,
`mulI`, `norm` (exact in fixed point when the two exponents are close; the
smaller square folded into the radius beyond wp+64 binades, as the real
`addSub` does). A Ziv driver `ziv2(z *Complex, inPrec uint, eval evalFn2)`
with

```go
type evalFn2 func(wp uint) (b *cball, Kre, Kim int64, ok bool)
```

calls `eval`, decides each part with the existing `decide` at that part's
precision and mode, keeps a part once decided, and stops when both are;
each part has its own scaling exponent so that a tiny part next to a huge one
never leaves `big.Float`'s range. Parts known exact (below) are set before
the loop and never enter it.

| function | enclosure | exact parts detected before the loop |
|---|---|---|
| `Add`, `Sub`, `Neg`, `Conj`, `MulI`, `Mul2Exp`, `Proj` | none: `big.Float` operations per part (one rounding each) | all |
| `Mul`, `Sqr` | ac - bd and ad + bc computed **exactly** (products exact in `big.Float` at the sum of the precisions, then one correctly rounded sum each, with the exponent-gap fold of the real `addSub`) | all |
| `Quo` | (ac + bd)/(c² + d²) and (bc - ad)/(c² + d²): exact numerators and denominator, one correctly rounded division per part | a part is exact iff the exact quotient is dyadic of the target length: tested by the integer division remainder |
| `Norm` | exact sum of exact squares, one rounding | all |
| `Exp` | e^x as the real enclosure with its 2^K scaling, cos y and sin y from one reduction | Im = ±0 iff y = ±0 (sign of y), since sin y ≠ 0 for dyadic y ≠ 0 (π is irrational); Re never zero (cos y = 0 needs y ∈ π/2 + πℤ); exp(0) = 1 |
| `Expm1` | Re = expm1(x) cos y - 2 sin²(y/2) evaluated as a ball expression (the cancellation is paid for by Ziv's precision increase, bounded since the result is nonzero), Im = e^x sin y | Im = ±0 iff y = ±0; x = y = 0 |
| `Log`, `Log1p`, `Log2`, `Log10` | Re: ½ log1p(t) with t = x² + y² - 1 formed **exactly** (dyadic arithmetic) when \|z\| is in [1/2, 2], ½ log(x² + y²) otherwise (exact norm, real log enclosure); Im: atan2 enclosure; `Log1p` forms (1+x)² + y² - 1 = 2x + x² + y² exactly | Re = 0 iff x² + y² = 1, which for dyadic x, y happens only at (±1, 0) and (0, ±1) (a² + b² = 4^k has only the trivial solutions); Im = ±0 iff y = ±0 and x > 0; Log2/Log10 of powers of 2 / 10 on the positive axis as in the real functions |
| `Sqrt` | Algorithm 312 form: t = √((\|x\| + \|z\|)/2), other part y/(2t), as balls | exact roots: a candidate w at the target precision is squared exactly and compared with x (the `Root` approach of design.md); Re = +0 on the negative axis, Im = ±0 on the positive axis |
| `Root(x, n, k)` | e^((Log x + 2πik)/n) | exact candidates checked by exact n-th power |
| `Pow`, `PowFloat`, `PowInt` | `PowInt`: binary powering in complex balls, error tracked; general: Exp(y·Log x) with the product in balls | `PowInt` exact when the exact Gaussian-dyadic power fits (computed exactly with `big.Int` while the size stays below the bound of the real `Pow`); `PowFloat` with y = k/2^j through exact square roots, as the real `Pow`; x or y zero per MPC §5.8 (result 1 for 0^0, zero signs by the tiny-value rule); remaining cases: no exactness (a part that is exactly zero there is a missed case, caught by `zivLimit` in tests, not by users) |
| `Sin`, `Cos`, `SinCos`, `Sinh`, `Cosh` | sin x cosh y + i cos x sinh y etc., one reduction of x mod π/2 and one exponential for cosh/sinh; `Sin(z) = -i Sinh(iz)` (G.6p7) | Sin: Re = ±0 iff x = ±0, Im = ±0 iff y = ±0; Cos: Im = ±0 iff x = 0 or y = 0, Re never zero; same by symmetry for Sinh, Cosh |
| `Tan`, `Tanh` | Re = sin 2x / (cos 2x + cosh 2y), Im = sinh 2y / (cos 2x + cosh 2y), the denominator as one ball (no cancellation: cosh 2y >= 1 >= -cos 2x, equality only at y = 0) | Re = ±0 iff x = ±0; Im = ±0 iff y = ±0 |
| `Asinh`, `Asin`, `Acos`, `Acosh`, `Atanh`, `Atan` | the Hull–Fairgrieve–Tang decomposition in balls: A = (\|z+i\| + \|z-i\|)/2, B = y/A, A - 1 = f(x, 1+y) + f(x, 1-y) with f(a, b) = a²/(hypot(a, b) + b)/2, log1p and atan2 enclosures; `Atanh` as ¼ log1p(4x/((1-x)² + y²)) + ½ i atan2(2y, (1-x)(1+x) - y²) with the norms exact | the real-axis segment where the result is real (e.g. Asin(x ± 0i), \|x\| <= 1: Im = ±0) and the imaginary-axis segment where it is imaginary; Acos(1) = 0, Acosh(1) = 0, Atanh(0) = 0; the poles of Atanh at ±1 (±Inf, divide-by-zero in G.6.2.3) |
| `Expj`, `ExpjPi`, `RootOfUnity` | cos and sin enclosures (of πx for `ExpjPi`, of 2πk/n) | `ExpjPi`: x a multiple of 1/2 gives exact ±1, ±i; `RootOfUnity`: k/n in (1/4)ℤ |

Cancellation is handled by the enclosures, not by case analysis: a ball whose
radius exceeds its midpoint makes `decide` fail and the loop raises the
working precision. The case analysis above exists only where the exact
result is representable (Ziv can never decide those) and where the formula
would overflow at any precision (the real functions' 2^K scaling).

Performance choices, all algorithmic: one argument reduction shared by the
sine and cosine of a part (`SinCos` exposes it); cosh and sinh from one
exponential; exact norms in fixed point (`fx`) rather than `big.Float`; the
real enclosures evaluated once per working precision for both parts.

## Validation

- **Oracle 1, mpmath** (`oracle/`): `MpmathScript` gains complex requests,
  `{"f": name, "p": bits, "c": [[[mr, er], [mi, ei]], …]}` answered by
  `{"re": {…}, "im": {…}}` (each part as today: sign, decimal integer mantissa,
  binary exponent, or `{"i": ±1}`), with `mp.mpc` built from exact parts and
  the function list mapped to `mpmath`'s complex forms. The `Rounder` applies
  per part: far more bits than the target, widened by the declared slack,
  Ziv until both ends round alike in **all six modes**, confirmed at twice the
  precision; `Exact` (integer arithmetic, independent of the package) gains
  the complex exact cases of the table above.
- **Oracle 2, MPC as a black box**: the fixture generator also evaluates each
  point with gmpy2's `mpc` (which, from knowledge, wraps MPC; to be confirmed
  when the generator is written) in its five modes, and the test requires
  bit equality with it — a judge independent of mpmath. Running MPC through a
  Python wheel uses it, it does not read or copy its sources. `ToNearestAway`
  is checked against oracle 1 only.
- **Fixtures**: `testdata/c<name>.json.gz`, as the real ones: precisions 2,
  24, 53, 64, 113, 200, 1000 and random ones up to 4096, all six modes,
  points drawn per region:
  - near 0, near the branch points ±1, ±i (distances 2^-p … 2^-8);
  - **each branch cut on both sides**: Im = +0 and -0, Im = ±2^(-p-10),
    ±2^-8, for every function of the cut list above;
  - \|z\| near 1 (`Log`), parts of very different magnitudes (2^-1000 next
    to 2^1000), exponents near `MaxExp`/`MinExp`;
  - the exact cases of the table (each must come back `Exact` in both parts
    where both are exact);
  - the special-value grid {±0, ±tiny, ±1, ±huge, ±Inf}² against the G.6
    tables, with the NaN-part entries asserting the `ErrNaN` panic.
- **Bit-exact symmetries** on every point: f(conj z) = conj f(z) with the
  imaginary part's mode mirrored (`ToNegativeInf` ↔ `ToPositiveInf`); odd and
  even functions with both parts' directed modes mirrored.
- **Coverage**: 100 % of statements, as the real package.

## Performance

- Measured against mpmath's complex functions at 53, 113, 256, 1024, 4096 and
  16384 bits through the oracle's timing request (`"n"`, timed inside Python),
  and against MPC through gmpy2 for information; single core, interleaved
  rounds, never at a 1-minute load above 8.
- Target: at or below mpmath at every precision for every function (a slower
  row is an open defect); the MPC ratio is published, not gated (C against
  `math/big`).
- Allocation counts per call (exact under load) are part of the table; the
  `cball` scratch values are pooled as `zivWS` is.

## Fuzzing and security

- Fuzz targets per function over (re mantissa bytes, re exponent, im
  mantissa bytes, im exponent, precision ≤ 4096, mode), corpus committed under
  `testdata/fuzz/`: no panic other than the documented `ErrNaN`/`ErrLimit`;
  bit-exact symmetries; the rounded result inside an enclosure computed
  independently at 2p + 64 bits; `Acc()` consistent with that enclosure.
- **Bounds**: `MaxReductionExponent` applies to every real part reduced mod
  π/2 (Sin, Cos, Tan and their hyperbolic forms' imaginary parts, `Exp`'s
  imaginary part, `Expj`); `PowInt` and exact `Pow` candidates stop at the
  size bound of the real `Pow`; parsing (`SetString`) bounds the digit count
  and exponent as `big.Float.SetString` does. Bomb tests assert the refusal
  (`ErrLimit`) under `GOMEMLIMIT` and `ulimit -v`, never against the host's
  memory.
- **32-bit** (386, arm in CI): exponent arithmetic on `int64` (scaling
  exponents Kre, Kim, gaps between parts), never `int`.
- No dependency is added (`math/big` only); `govulncheck` already runs.

## File layout

```text
complex.go          type Complex, setters, conversion, String/Text/SetString
complex_arith.go    Neg, Conj, Proj, MulI, Add, Sub, Mul, Sqr, Quo, *Float, Mul2Exp, Norm
complex_exp.go      Exp, Expm1, Log, Log1p, Log2, Log10, Sqrt, Root
complex_pow.go      Pow, PowFloat, PowInt
complex_trig.go     Sin, Cos, SinCos, Tan, Sinh, Cosh, Tanh, Expj, ExpjPi, RootOfUnity
complex_inv.go      Asin, Acos, Atan, Asinh, Acosh, Atanh
complex_annexg.go   the special-value dispatch shared by all of the above
cball.go            complex balls, ziv2
complex_*_test.go, complex_fuzz_test.go, complex_bench_test.go
oracle/             complex requests, Rounder per part, Exact complex cases
testdata/c*.json.gz
docs/design-complex.md   this page
```

## Tasks (one pull request each, in order)

1. **Complex oracle**: complex requests in `MpmathScript`, the per-part
   `Rounder`, complex `Exact`, the generator's gmpy2/MPC column, fixture
   format.
2. **`Complex` type, setters and exact or once-rounded arithmetic**: `Neg`,
   `Conj`, `Proj`, `MulI`, `Mul2Exp`, `Add`, `Sub`, `Mul`, `Sqr`, `Quo`,
   `*Float` forms, `Norm`, `String`/`SetString`; the Annex G dispatch;
   `cball` and `ziv2`.
3. **`Exp`, `Expm1`, `Log`, `Log1p`, `Log2`, `Log10`, `Sqrt`, `Root`.**
4. **`Sin`, `Cos`, `SinCos`, `Tan`, `Sinh`, `Cosh`, `Tanh`, `Expj`,
   `ExpjPi`, `RootOfUnity`.**
5. **`Asin`, `Acos`, `Atan`, `Asinh`, `Acosh`, `Atanh`.**
6. **`Pow`, `PowFloat`, `PowInt`.**
7. **Perf report against mpmath and MPC, fuzz corpus**; coverage.md
   ("Complex numbers" rows done, the "complex type" dependencies resolved);
   README, landing and docs site; release (minor).

Then the slices of coverage.md that wait on it, each on `Complex`:
the gamma family (Bernoulli numbers, `gamma`, `loggamma` with its branch cut,
`rgamma`, `psi`, `polygamma`, …), zeta and polylog, the `hyper` engine (and
through it Bessel, Airy, `erf`, exponential integrals, orthogonal
polynomials), elliptic functions and theta/modular/Weierstrass, `agm`,
`polyroots`, `fft` over `big.Float`.
