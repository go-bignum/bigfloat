# Design of go-bignum/bigfloat

## The contract

Every function returns the **correctly rounded** value of its exact
mathematical result: the value that rounding the infinitely precise result to
the receiver's precision, in the receiver's rounding mode, would give. That
holds for every precision and for all six `big.RoundingMode`s
(`ToNearestEven`, `ToNearestAway`, `ToZero`, `AwayFromZero`, `ToNegativeInf`,
`ToPositiveInf`). `z.Acc()` reports truthfully whether the returned value is
below, equal to, or above the exact result. This is the MPFR contract
(Fousse, Hanrot, Lefèvre, Pélissier, Zimmermann, "MPFR: A multiple-precision
binary floating-point library with correct rounding", ACM TOMS 33(2), 2007,
doi:10.1145/1236463.1236468). It makes results reproducible and testable bit
for bit, on every architecture.

Conventions follow `math/big`: a receiver of precision 0 takes its operand's
precision (the larger of two), constants default to 64 bits, the receiver may
alias an operand. Exponent overflow gives ±Inf and underflow ±0, exactly as
`big.Float`'s own arithmetic: the result is ±Inf when the *rounded* exponent
exceeds `MaxExp`, with `Acc()` set as `big.Float` sets it. `big.Float` has no
NaN; a function whose result would be NaN panics with `bigfloat.ErrNaN`, as
`big.Float.Sqrt` panics with `big.ErrNaN` (whose message field is
unexported, hence a type of the same shape here). Special values follow C99
Annex F and IEEE 754 (`Atan2(±0, -0) = ±π`, `Pow(-1, ±Inf) = 1`,
`Root(-0, 2) = +0` as C23 `rootn`, …).

## Error model: balls

All approximations are **enclosures**. A `ball` is a midpoint `big.Float` and
a radius `mag` (Johansson, "Arb: efficient arbitrary-precision midpoint-radius
interval arithmetic", IEEE TC 66(8), 2017, doi:10.1109/TC.2017.2690633): each
operation returns a ball guaranteed to contain the exact result of the
operation applied to any points of its operands, the rounding of its own
midpoint added to the radius. The radius type `mag` is an upper bound m·2^e,
float64 m rounded upwards on every operation, int64 e, so radii never lose
range. Rigour is therefore by construction, not by a per-function error
analysis: a function is correct as soon as each of its steps is a valid ball
operation and each truncated series adds a valid bound on its tail.

The series cores run in **fixed point** (`fx`): an integer v and an error
count e in units of 2^-w, again an upper bound kept by rounding upwards. A
big.Int multiply-and-shift costs much less than a normalising big.Float
operation (Brent & Zimmermann, *Modern Computer Arithmetic*, CUP 2010, §4.4;
Arb does the same). Every fixed-point operation bounds `|computed − exact|`
for any operands within their errors, so input radii propagate rigorously.

Two choices remove whole classes of slow paths:

- an operand more than wp+64 binades below the other is folded into the
  radius instead of added: `big.Float` aligns operands by shifting by the
  full exponent gap, which for e^x − 2^-k at huge x meant gigabits;
- long arguments (precision far above the target) are rounded into balls
  first, their rounding becoming radius.

## Ziv's strategy

`ziv` evaluates an enclosure at a working precision wp, rounds both of its
ends at the target precision, and accepts the result when they agree and the
rounded value lies outside the ball (so that its side, hence `Acc()`, is
known); otherwise it raises wp by half and tries again (Ziv, "Fast evaluation
of elementary mathematical functions with correctly rounded last bit", ACM
TOMS 17(3), 1991, doi:10.1145/114697.116813). Each function's enclosure is of
f(x)·2^-K for an exponent K chosen by the function (e^x = 2^k·e^r), so that
the midpoint never leaves `big.Float`'s range and over/underflow is decided
once, at the final scaling, with `big.Float`'s own semantics.

An enclosure of a representable value, or of a rounding midpoint, never
rounds to a single value: **exact results never reach the loop**. Each
function detects them first:

- the trivial ones (exp(0), log(1), sin(0), acos(1), …);
- integer logarithms: log2 of powers of two, log10 of 10^k for k >= 0 (the
  only powers of ten a binary number can hold);
- roots and hypot are computed exactly with integers (below);
- powers: x^y with y = k/2^j is dyadic exactly when x is a perfect 2^j-th
  power and the k-th power is dyadic (a·k + b·2^j = 1 shows rationality of
  x^(k/2^j) implies that of x^(1/2^j)); found by at most about 32 + log2
  prec(x) exact square roots, then a bounded integer power;
- |x| = 1.

Everything else these functions return at a dyadic argument is irrational
(Lindemann–Weierstrass for exp, log, the trigonometric and hyperbolic
functions and their inverses; Gelfond–Schneider and the dyadic analysis above
for pow), so the loop terminates.

### Proxies: tiny arguments and saturation

Near-identity results (sin x ≈ x for tiny x, cos x ≈ 1, tanh x ≈ ±1 for large
x, expm1 x ≈ −1 for very negative x, x^y ≈ 1) sit so close to a representable
number that Ziv's loop would need about twice the argument's exponent in bits
— impossible for x = 2^-10^9. Instead the result is set from a **proxy**: an
exact number base ± 2^(exp(base)−P−1), P = max(p+2, prec(base)) + 2, which
lies strictly between base and every rounding breakpoint of precision p near
it, on the same side as the true value. Rounding the proxy gives the correct
value and the correct `Acc()`. The proxy applies when the distance to base is
provably below 2^(exp(base)−P+1), from a bound on the series remainder
(|sin x − x| <= |x|³/6, |e^x − 1| <= 2|x|, …). The same trick appends a
sticky bit to exact integer roots.

## Algorithms

| function | method |
|---|---|
| Sqrt, Cbrt, Root(n) | integer n-th root of the mantissa scaled to at least p+2 result bits, with a sticky bit when inexact: exact in every mode (math/big's own `Sqrt` is not correctly rounded in every mode). Newton from above, started from the exact root of the top half of the integer (recursively, down to a float64 estimate: the precision doubles from one level to the next, as mpmath's `nthroot_fixed`), so that one step reaches the root; x^n <= N then proves it. Beyond 2^20 bits of work: exp(log(x)/n), with an exactness test (odd part a perfect n-th power, n dividing the exponent). |
| Hypot | exact sum of squares after scaling; a proxy when the smaller argument is below half the precision |
| Exp | x = k·ln2 + r in fixed point against a cached fixed-point ln2, r reduced by 2^-s, Taylor series on m ≈ sqrt(n) accumulators (one division per term, one multiplication every m terms: mpmath's `exp_basecase` with m sums), s squarings (MCA §4.3.1); s ≈ sqrt(w/2) |
| Expm1 | the same with expm1 doublings E <- E(E+2) for |x| < 1/2, relative accuracy kept by scaling the fixed point to x's exponent |
| Log, Log2, Log10 | x = m·2^e, m in [1/sqrt2, sqrt2). Up to 8192 bits, mpmath's table reduction (`log_taylor_cached`): log m = log a + 2·atanh((m−a)/(m+a)) with a = n/2^9 the grid point below m, |t| < 2^-9 (18 bits per term), log a cached; within 2^-8 of 1 directly 2·atanh((m−1)/(m+1)), which keeps the relative accuracy of log m near 0. Above 8192 bits the direct quotient, |t| <= 0.172. Quotients formed exactly from the integer mantissa; e·ln2 added in fixed point; extra bits as m−1 has leading zeros |
| Log1p | 2·atanh(x/(2+x)) for |x| < 1/2 (formed exactly), log(1+x) otherwise |
| Pow | exact cases above; else exp(y·log|x|) with log|x| to wp + log2|y·log x| bits; overflow, underflow and near-1 proxies from a float64 bound of y·log|x| with a rigorous margin (|y|·2^-50 + |t|·2^-50), or from a 64-bit enclosure when that bound is near a threshold |
| Sin, Cos, Tan | x = q·π/2 + r with π to exp(x) + w bits (exact reduction, as MPFR's `mpfr_sin`; equivalent in precision to Payne–Hanek), q from one integer division against the cached fixed-point π; the reduction is redone with more bits when r has leading zeros. Up to scale 640, mpmath's `cos_sin_basecase`: sin and cos of d = r − n/2^8 (|d| < 2^-8) by their series, rotated by cached sin/cos(n/2^8); above, both series at r/2^j and j double-angle steps. The two series share the powers of r² on m accumulators |
| Atan, Atan2, Asin, Acos | up to 8192 bits, mpmath's `atan_taylor`: atan t = atan a + atan((t−a)/(1+a·t)), a = n/2^7 the grid point below |t|, atan a cached, 14 bits per term; above, j argument halvings t/(1+sqrt(1+t²)); the alternating series on m accumulators; |t| > 1 through π/2 − atan(1/t). Asin for |x| <= 0.86 by mpmath's 2·atan(x/(1+sqrt(1−x²))) in fixed point, acos = π/2 − asin there; beyond, asin = atan(x/sqrt((1−x)(1+x))) and acos likewise, 1 − x² formed exactly |
| Sinh, Cosh, Tanh | e^|x| = b·2^k and 2^-2k/b, combined without cancellation; expm1 for small arguments |
| Asinh, Acosh, Atanh | log1p forms free of cancellation (acosh(1+d) = log1p(d + sqrt(d(2+d))), d exact); log 2x plus a radius for huge arguments |
| Pi | Chudnovsky series by binary splitting (Haible & Papanikolaou, ANTS-III, LNCS 1423, 1998; MCA §4.9); each term ratio below 2^-45, so the tail after N terms is below 2^(2−45N) |
| E | Σ 1/k! by binary splitting, tail below 1/N! |
| Ln2 | 18 atanh(1/26) − 2 atanh(1/4801) + 8 atanh(1/8749), each by binary splitting |
| Ln10 | 3 ln2 + 2 atanh(1/9) |
| EulerGamma | Brent–McMillan algorithm B1 (Math. Comp. 34, 1980; MCA §4.12): U/V − γ = K0(2n)/I0(2n), below π e^-4n, with tails of U and V bounded from the computed balls |

Constants are cached at the highest precision asked so far and rounded down
from the cache; π and ln2 also as fixed-point balls per scale (up to 64
scales of at most 8192 bits).

### Speed without losing rigour

The table reductions keep the error model: each table entry is a
fixed-point ball computed by the table-free path at a power-of-two class
scale (next power of two above the scale asked, plus 32 bits, as mpmath's
`cache_prec_steps`), and truncating it to a smaller scale adds one ulp.
Entries are computed on first use and shared read-only; truncated copies
are memoised, at most 1024 per table. The series kernels write into
caller-provided fixed-point values drawn from size-classed `sync.Pool`s,
the Ziv decision reuses pooled `big.Float`s, the integer square root is an
allocation-free Newton iteration, and the float64 bounds of the error
counts move one ulp outwards with two flops instead of `math.Nextafter`
(x·(1 ± 2^-52) ± 2^-1074, at least one ulp for every non-negative x). A
53-bit call allocates a handful of objects.

Measured against mpmath in [perf.md](perf.md). Two kernels could serve
other fleet packages and belong lower down if a second user appears: the
fixed-point ball (`fx`) with its m-accumulator series, and the allocation
free integer square root and n-th root.

AGM-based logarithms and π (MCA §4.8) pay only from tens of thousands of
bits, beyond this package's present targets (up to a few thousand bits,
where the series with binary splitting or fixed point win); they are left for
when a user needs them, with the same ball framework.

Prior art: `ALTree/bigfloat` (MIT) offers Exp, Log and Pow on `big.Float`
without a correct-rounding guarantee; it was not used. MPFR (LGPL) is an
oracle reference only; no code was ported.

## Resource bounds

Inputs to these functions may be untrusted, so every cost is bounded and the
bound documented:

| what | bound | beyond it |
|---|---|---|
| Ziv working precision | 64·(p + input precision) + 65536 bits; worst cases of the table maker's dilemma need about 2p | panic (internal error): a missed exact case, never silently wrong |
| trigonometric argument exponent | `MaxReductionExponent` = 2^22 (π to four million bits) | panic `ErrLimit` |
| exact integer roots | 2^20 bits of work, then exp(log x / n) | — |
| exact integer powers | 2^20 bits | treated as inexact (still correctly rounded by Ziv) |
| exponent of exp/sinh/cosh arguments | |x| > 2^34 overflows or underflows at once | — |
| operands far apart | an operand wp+64 binades below the other becomes radius | — |
| long arguments | rounded to the working precision with their rounding in the radius | — |
| value tables (log, atan, sin, cos) | entries only up to scale 8192 (sin/cos 640), one per grid point and power-of-two class; at most 1024 truncated copies per table (about 1 MiB) | table-free path |
| fixed-point constant copies | 64 scales of at most 8192 bits | computed, not kept |

Other costs are proportional to the requested precision or to the
arguments' own size. The fuzz test found two unbounded paths, both fixed:
`Pow(-1, even integer)` looped to the Ziv limit, and the oracle's exactness
check built x^98304 as a five-megabit rational. It later found a defect of
the fixed-point reduction: the error of k·π/2 went through float64(k),
+Inf for |k| >= 2^1024 (Tan(1.5·2^1130) panicked on an internal invariant);
that error is now bounded after rescaling, never formed in float64.

## Verification

The oracle harness (`oracle`, `cmd/oracle`) drives mpmath 1.4.1 from Go over a
JSON-lines pipe, the Python side being one minimal script (`MpmathScript`)
that ends with its standard input; values travel as exact (mantissa,
exponent) pairs. A `Rounder` asks for 2p + 64 bits, widens the answer by the
oracle's declared slack (8 bits for mpmath), doubles the precision until both
ends of the interval round alike in all six modes, and **confirms** the
decision at twice that precision (an oracle less accurate than it claims at
one precision — mpmath rounding an argument longer than its working
precision, found during this work — is caught, not trusted). Exact results
are recognised by `oracle.Exact`, written with plain `big.Int`/`big.Rat`
independently of the package under test.

`cmd/oracle` generates the fixtures: for each function, random arguments per
region — tiny (on both sides of the proxy thresholds), small, the bulk of the
domain, huge (sin of 2^1000), near singularities (log near 1, acosh near 1,
atanh and asin near ±1, sin/cos/tan at the closest p-bit numbers to k·π/2),
the grid points of the table reductions and the fixed-point asin bound,
exact cases (perfect powers, powers of two and ten, Pythagorean pairs), and
the same regions with arguments of 4p+64 bits — at 24, 53, 64, 113, 200 and
1000 bits; constants at every precision from 2 to 300 bits and at 500, 1000
and 2000. 28 840 cases × 6 modes = 173 040 results, compared bit for bit
(value, sign of zero and accuracy) by `TestFixtures`, all exact.
`TestHostileMpmath` (when uv is installed; `TestHostileSelf` always) adds
exact grid points, ulp-wide perturbations of them, table/series edges and
huge or tiny exponents at eleven precisions straddling the table limits.

Also: special values and domain errors by table; unit tests of the error
model's refusals; `FuzzFunctions`, which compares every function, mode and
precision with the Rounder over the package itself at more than twice the
precision and admits only `ErrNaN` and `ErrLimit` panics (about 19 million
executions over two campaigns, three defects found and fixed, the crashers
kept in the corpus); 100% statement coverage with
`-race`; `GOARCH=386`.

## The oracle for the fleet

`oracle.Oracle` is the interface a validation suite needs: `Eval(fn, prec,
args...)` to about prec bits, a declared `Slack`, a `Name`. `Mpmath` and
`Self` implement it; `Rounder` turns either into correctly rounded results in
every mode; `Fixture`/`Case` store them, `Case.Check` compares an
implementation. A package such as go-numerics/special can replace its mpmath
calls by `oracle.Self{}` (no Python in CI) once its functions are expressed
through these elementary ones, keeping mpmath for the functions bigfloat does
not provide yet.
