# Coverage against mpmath and Mathematica

Inventory of what `go-bignum/bigfloat` (main, after v0.1.0) implements
against the public function list of mpmath, with Mathematica's
arbitrary-precision and number-theory functions as the delta and MPFR / Arb
capabilities as names only. `bigfloat` is the only repository of the
go-bignum organisation that holds code (the others are `brand`, `docs`, the
landing site and `.github`).

Status: **done** (same function, real arguments, and at least mpmath's
accuracy: here correctly rounded in all six modes), **partial** (what is
missing is named), **missing**. "Dep" is what a missing item needs beyond
this package. "Done (math/big)" means the standard library already provides
it and `bigfloat` adds nothing.

Sources. mpmath: `mpmath/__init__.py` at tag 1.4.1 (the version the oracle
uses) read line by line, plus the names that exist only on master
(`weierp` family, `rtheta`, `rtheta_jet`, `kleinjinv`, `g2g3from`,
`omega1omega2from`, `pinv`, `lll_gram`, `ulp`, `spherical_in`,
`spherical_kn`, `fft`, `invfft`), marked "(master)". Go side: README,
`go doc -all`, `docs/design.md` and the sources on main. Mathematica,
MPFR and Arb: from knowledge of their documented function families; the
Wolfram, MPFR and Arb reference pages were not opened for this inventory, so
every name in those columns must be re-checked against the manuals before a
slice is started.

How to read the real-only scope. `bigfloat` works on `*big.Float` and is
real-valued: where mpmath returns a complex number for a real argument
(`sqrt(-1)`, `log(-2)`, `acos(2)`, `asin(2)`, `acosh(0)`, `atanh(2)`,
`power(-8, 1/3)`) the Go function panics with `ErrNaN`. That is one row
("complex arguments") instead of a caveat on every function. Where Go has
no counterpart of a function, a missing row is not necessarily hard: the
dependency column says what it needs.

## Totals

One count per backticked mpmath name in the first column of a table,
computed by a throwaway Go program from this file (see "Counting").

| Section | done | partial | missing | names |
|---|---|---|---|---|
| Constants | 5 | 0 | 8 | 13 |
| Elementary functions | 26 | 1 | 22 | 49 |
| Arithmetic and utilities | 11 | 0 | 19 | 30 |
| Complex numbers | 0 | 0 | 11 | 11 |
| Gamma, factorial, beta, polygamma | 0 | 0 | 22 | 22 |
| Zeta, polylog and relatives | 0 | 0 | 19 | 19 |
| Combinatorial and number-theoretic | 0 | 0 | 13 | 13 |
| Hypergeometric and q-series | 0 | 0 | 23 | 23 |
| Error function, exponential integrals, Fresnel | 0 | 0 | 16 | 16 |
| Airy and Scorer | 0 | 0 | 6 | 6 |
| Elliptic integrals and AGM | 0 | 0 | 10 | 10 |
| Elliptic functions, theta, modular, Weierstrass | 0 | 0 | 19 | 19 |
| Orthogonal polynomials and other | 0 | 0 | 14 | 14 |
| Bessel and related | 0 | 0 | 29 | 29 |
| Waveforms, sigmoid, fft | 0 | 0 | 7 | 7 |
| **Functions and constants, subtotal** | **42** | **1** | **238** | **281** |
| Tools (precision, intervals, calculus, solvers, matrices, I/O) | 0 | 0 | 112 | 112 |

The Mathematica, MPFR and Arb tables are not counted.

## Constants

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `pi` | `Pi` | `Pi` | done | |
| `e` | `E` | `E` | done | |
| `ln2` | `Log[2]` | `Ln2` | done | |
| `ln10` | `Log[10]` | `Ln10` | done | |
| `euler` | `EulerGamma` | `EulerGamma` | done | |
| `phi` | `GoldenRatio` | | missing | none (`(1+sqrt 5)/2`, correctly rounded via `Sqrt`) |
| `catalan` | `Catalan` | | missing | none (series; see Brent-McMillan style or Zeta-based formula) |
| `khinchin` | `Khinchin` | | missing | `zeta` derivative series; slow convergence tables |
| `glaisher` | `Glaisher` | | missing | `zeta'(-1)` or the Kinkelin series |
| `apery` | `Zeta[3]` | | missing | none (Zeta(3) fast series) |
| `degree` | `Degree` | | missing | none (`Pi/180`) |
| `twinprime` | | | missing | none (mpmath tabulates it; needs a series for arbitrary digits) |
| `mertens` | | | missing | none (same) |

## Elementary functions

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `sqrt` | `Sqrt` | `Sqrt` | done | complex result for negative x: see complex row |
| `cbrt` | `CubeRoot` | `Cbrt` | done | |
| `root`, `nthroot` | `Surd`, `Power` | `Root` | done | k-th complex root of `root(x, n, k)` |
| `exp` | `Exp` | `Exp` | done | |
| `expm1` | `Exp[x]-1` | `Expm1` | done | |
| `ln`, `log` | `Log` | `Log` | done | `log(x, b)` two-argument form is `Log(x)/Log(b)`, not provided |
| `log2` | `Log2` is `Log[2,x]` | `Log2` | done | |
| `log10` | `Log10` | `Log10` | done | |
| `log1p` | `Log[1+x]` | `Log1p` | done | |
| `power` | `Power` | `Pow` | done | complex base or exponent, negative base with integer exponent is covered |
| `exp2` | `2^x` | | partial | `Pow(2, x)` works; no named function |
| `powm1` | | | missing | none (`x^y - 1` without cancellation) |
| `cos` | `Cos` | `Cos` | done | |
| `sin` | `Sin` | `Sin` | done | |
| `tan` | `Tan` | `Tan` | done | |
| `cosh` | `Cosh` | `Cosh` | done | |
| `sinh` | `Sinh` | `Sinh` | done | |
| `tanh` | `Tanh` | `Tanh` | done | |
| `acos` | `ArcCos` | `Acos` | done | |
| `asin` | `ArcSin` | `Asin` | done | |
| `atan` | `ArcTan` | `Atan` | done | |
| `asinh` | `ArcSinh` | `Asinh` | done | |
| `acosh` | `ArcCosh` | `Acosh` | done | |
| `atanh` | `ArcTanh` | `Atanh` | done | |
| `atan2` | `ArcTan[x,y]` | `Atan2` | done | |
| `hypot` | `Norm[{x,y}]` | `Hypot` | done | |
| `sec` | `Sec` | | missing | none (correctly rounded `1/cos` needs its own Ziv loop) |
| `csc` | `Csc` | | missing | none |
| `cot` | `Cot` | | missing | none |
| `sech` | `Sech` | | missing | none |
| `csch` | `Csch` | | missing | none |
| `coth` | `Coth` | | missing | none |
| `asec` | `ArcSec` | | missing | none |
| `acsc` | `ArcCsc` | | missing | none |
| `acot` | `ArcCot` | | missing | none |
| `asech` | `ArcSech` | | missing | none |
| `acsch` | `ArcCsch` | | missing | none |
| `acoth` | `ArcCoth` | | missing | none |
| `cospi` | `Cos[Pi x]` | | missing | none (exact at half-integers) |
| `sinpi` | `Sin[Pi x]` | | missing | none |
| `sinc` | `Sinc` | | missing | none |
| `sincpi` | `Sinc[Pi x]` | | missing | none |
| `cos_sin` | | | missing | none (one reduction, two results) |
| `cospi_sinpi` | | | missing | none |
| `unitroots` | `RootOfUnity`, `Exp[2 Pi I k/n]` | | missing | complex type |
| `cyclotomic` | `Cyclotomic` | | missing | polynomial evaluation at a point |
| `lambertw` | `ProductLog` | | missing | none for real branches; complex branches need the complex type |

## Arithmetic and utilities

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `fneg` | `Minus` | `big.Float.Neg` | done | done (math/big) |
| `fadd` | `Plus` | `big.Float.Add` | done | done (math/big) |
| `fsub` | `Subtract` | `big.Float.Sub` | done | done (math/big) |
| `fmul` | `Times` | `big.Float.Mul` | done | done (math/big) |
| `fdiv` | `Divide` | `big.Float.Quo` | done | done (math/big) |
| `fabs` | `Abs` | `big.Float.Abs` | done | done (math/big) |
| `sign` | `Sign` | `big.Float.Sign` | done | done (math/big) |
| `ldexp` | | `big.Float.SetMantExp` | done | done (math/big) |
| `frexp` | | `big.Float.MantExp` | done | done (math/big) |
| `isinf` | | `big.Float.IsInf` | done | done (math/big) |
| `isint` | `IntegerQ` | `big.Float.IsInt` | done | done (math/big) |
| `floor` | `Floor` | | missing | none (math/big has `Int`, which truncates) |
| `ceil` | `Ceiling` | | missing | none |
| `nint` | `Round` | | missing | none (ties-to-even or away) |
| `frac` | `FractionalPart` | | missing | none |
| `fmod` | `Mod` | | missing | none (exact for any exponent gap) |
| `fsum` | `Total` | | missing | none (exactly rounded sum) |
| `fdot` | `Dot` | | missing | none |
| `fprod` | `Times @@` | | missing | none |
| `ulp` (master) | | | missing | none |
| `degrees` | `x/Degree` | | missing | `Pi` (have) |
| `radians` | `x Degree` | | missing | `Pi` (have) |
| `absmin` | | | missing | none |
| `absmax` | | | missing | none |
| `almosteq` | | | missing | none |
| `isnormal` | | | missing | none |
| `isfinite` | | | missing | none |
| `isnan` | | | missing | not applicable: `big.Float` has no NaN |
| `nint_distance` | | | missing | none |
| `mag` | `Exponent`-like (`MantissaExponent`) | | missing | none |

## Complex numbers

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `mpc`, `make_mpc` | `Complex` | | missing | a `Complex` type on `*big.Float` pairs (new package or here); everything below depends on it |
| `re` | `Re` | | missing | complex type |
| `im` | `Im` | | missing | complex type |
| `conj` | `Conjugate` | | missing | complex type |
| `arg`, `phase` | `Arg` | | missing | complex type, `Atan2` (have) |
| `polar` | `AbsArg` | | missing | complex type |
| `rect` | | | missing | complex type, `Sin`, `Cos` (have) |
| `expj` | `Exp[I x]` | | missing | complex type |
| `expjpi` | `Exp[I Pi x]` | | missing | complex type |

Every function above and below whose mpmath form accepts a complex argument
(almost all of them) is real-only in Go until this type exists. Complex
`exp`, `log`, `sqrt`, `power`, the trigonometric and hyperbolic families and
their inverses are the first slice after the type, because the special
functions are built on them.

## Gamma, factorial, beta, polygamma

All depend on correctly rounded `Log`, `Exp`, `Sin`, `Pi` (have) and, for
the complex plane, the complex type. A shared kernel (Stirling series with
Bernoulli numbers, Spouge / Lanczos not acceptable for certified rounding)
serves most of them. `go-numerics/special` has float64 forms only.

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `gamma` | `Gamma` | | missing | Stirling with Bernoulli numbers, exact reflection; `Log`, `Exp` |
| `rgamma` | `1/Gamma` | | missing | `gamma` |
| `loggamma` | `LogGamma` | | missing | `gamma` kernel; complex branch cut |
| `factorial`, `fac` | `Factorial` | | missing | `gamma` (integer fast path via `big.Int`) |
| `fac2` | `Factorial2` | | missing | `gamma` |
| `beta` | `Beta` | | missing | `gamma` |
| `betainc` | `Beta[z1,z2,a,b]`, `BetaRegularized` | | missing | `hyp2f1` |
| `psi` | `PolyGamma` | | missing | `gamma` kernel |
| `polygamma` | `PolyGamma[n,x]` | | missing | `psi`, `zeta` |
| `digamma` | `PolyGamma[x]` | | missing | `psi` |
| `gammainc` | `Gamma[a,z0,z1]` | | missing | `hyp1f1`, continued fraction |
| `lower_gamma` | `Gamma[a,0,z]` | | missing | `gammainc` |
| `upper_gamma` | `Gamma[a,z]` | | missing | `gammainc` |
| `gammaprod` | | | missing | `gamma` |
| `binomial` | `Binomial` | | missing | `gamma` (`big.Int.Binomial` covers integer arguments) |
| `rf` | `Pochhammer` | | missing | `gamma` |
| `ff` | `FactorialPower` | | missing | `gamma` |
| `barnesg` | `BarnesG` | | missing | `loggamma`, `zeta` |
| `superfac` | `Superfactorial`-like (`BarnesG[n+2]`) | | missing | `barnesg` |
| `hyperfac` | `Hyperfactorial` | | missing | `barnesg` or `zeta'` |
| `harmonic` | `HarmonicNumber` | | missing | `psi` |

## Zeta, polylog and relatives

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `zeta` | `Zeta` | | missing | Euler-Maclaurin, Borwein, Riemann-Siegel; Bernoulli numbers |
| `altzeta` | `DirichletEta` | | missing | `zeta` |
| `stieltjes` | `StieltjesGamma` | | missing | quadrature or Euler-Maclaurin on derivatives |
| `hurwitz` | `Zeta[s,a]` | | missing | Euler-Maclaurin |
| `dirichlet` | `DirichletL` | | missing | `hurwitz` |
| `polylog` | `PolyLog` | | missing | series, `zeta`, `hyp` |
| `lerchphi` | `LerchPhi` | | missing | `hurwitz` family |
| `primezeta` | `PrimeZetaP` | | missing | `zeta`, Moebius sums |
| `secondzeta` | | | missing | `zeta`, quadrature |
| `riemannr` | `RiemannR` | | missing | `zeta`, `li` |
| `siegeltheta` | `RiemannSiegelTheta` | | missing | `loggamma` |
| `siegelz` | `RiemannSiegelZ` | | missing | `zeta` |
| `grampoint` | `GramPoint`-like | | missing | `siegeltheta` |
| `zetazero` | `ZetaZero` | | missing | `siegelz`; tables |
| `nzeros` | | | missing | `siegeltheta` |
| `backlunds` | | | missing | `nzeros` |
| `clsin` | `ClausenS`-like (`Im PolyLog`) | | missing | `polylog` |
| `clcos` | `ClausenC`-like | | missing | `polylog` |
| `polyexp` | | | missing | `polylog` |

## Combinatorial and number-theoretic

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `bernoulli` | `BernoulliB` | | missing | fast Bernoulli (zeta-based or recursion on `big.Int`) |
| `bernfrac` | `BernoulliB` (exact) | | missing | `big.Rat` |
| `bernpoly` | `BernoulliB[n,x]` | | missing | `bernoulli` |
| `eulerpoly` | `EulerE[n,x]` | | missing | Euler numbers |
| `eulernum` | `EulerE` | | missing | none (`big.Int`) |
| `bell` | `BellB` | | missing | none (`big.Int`) |
| `stirling1` | `StirlingS1` | | missing | none |
| `stirling2` | `StirlingS2` | | missing | none |
| `fib`, `fibonacci` | `Fibonacci` | | missing | none (closed form for non-integer index needs `Pow`, `Cos`) |
| `primepi` | `PrimePi` | | missing | sieve or Meissel-Lehmer |
| `primepi2` | | | missing | `primepi` |
| `mangoldt` | `MangoldtLambda` (via `FactorInteger`) | | missing | factorisation |

`unitroots`, `cyclotomic` and `lambertw` are listed under Elementary;
`identify`, `findpoly` and `pslq` under Tools.

## Hypergeometric and q-series

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `hyper` | `HypergeometricPFQ` | | missing | series with certified remainder, asymptotic and transformation regimes |
| `hyp0f1` | `Hypergeometric0F1` | | missing | `hyper` |
| `hyp1f1` | `Hypergeometric1F1` | | missing | `hyper` |
| `hyp1f2` | | | missing | `hyper` |
| `hyp2f1` | `Hypergeometric2F1` | | missing | `hyper`, transformations near 1 |
| `hyp2f2` | | | missing | `hyper` |
| `hyp2f0` | `HypergeometricPFQ[{a,b},{},z]` | | missing | divergent series, `hyperu` |
| `hyp2f3` | | | missing | `hyper` |
| `hyp3f2` | | | missing | `hyper` |
| `hyperu` | `HypergeometricU` | | missing | `hyper` |
| `hypercomb` | | | missing | `hyper` and limit handling |
| `meijerg` | `MeijerG` | | missing | `hyper` |
| `foxh` | `FoxH` | | missing | `meijerg` |
| `appellf1` | `AppellF1` | | missing | `hyper` |
| `appellf2` | | | missing | `hyper2d` |
| `appellf3` | | | missing | `hyper2d` |
| `appellf4` | | | missing | `hyper2d` |
| `hyper2d` | | | missing | `hyper` |
| `bihyper` | `HypergeometricPFQ` (bilateral) | | missing | `hyper` |
| `qp` | `QPochhammer` | | missing | none beyond `Pow` |
| `qhyper` | `QHypergeometricPFQ` | | missing | `qp` |
| `qgamma` | `QGamma` | | missing | `qp` |
| `qfac` | `QFactorial` | | missing | `qgamma` |

## Error function, exponential integrals, Fresnel

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `erf` | `Erf` | | missing | `hyp1f1` or direct series; asymptotic for large x |
| `erfc` | `Erfc` | | missing | `erf`, continued fraction |
| `erfi` | `Erfi` | | missing | `erf` |
| `erfinv` | `InverseErf` | | missing | `erf`, Newton with certified step |
| `npdf` | `PDF[NormalDistribution]` | | missing | `Exp` (have) |
| `ncdf` | `CDF[NormalDistribution]` | | missing | `erfc` |
| `expint` | `ExpIntegralE` | | missing | `gammainc` |
| `e1` | `ExpIntegralE[1,x]` | | missing | `expint` |
| `ei` | `ExpIntegralEi` | | missing | series, asymptotic |
| `li` | `LogIntegral` | | missing | `ei` |
| `ci` | `CosIntegral` | | missing | `ei` or direct series |
| `si` | `SinIntegral` | | missing | same |
| `chi` | `CoshIntegral` | | missing | same |
| `shi` | `SinhIntegral` | | missing | same |
| `fresnels` | `FresnelS` | | missing | `erf` or series |
| `fresnelc` | `FresnelC` | | missing | same |

## Airy and Scorer

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `airyai` | `AiryAi` | | missing | `hyp0f1`, Bessel |
| `airybi` | `AiryBi` | | missing | same |
| `airyaizero` | `AiryAiZero` | | missing | `airyai`, root finding |
| `airybizero` | `AiryBiZero` | | missing | `airybi`, root finding |
| `scorergi` | | | missing | `airybi`, quadrature |
| `scorerhi` | | | missing | same |

## Elliptic integrals and AGM

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `ellipk` | `EllipticK` | | missing | `agm` |
| `ellipe` | `EllipticE` | | missing | `agm`, or `elliprg` |
| `ellipf` | `EllipticF` | | missing | `elliprf` |
| `ellippi` | `EllipticPi` | | missing | `elliprj` |
| `elliprc` | `CarlsonRC` | | missing | duplication algorithm, certified error |
| `elliprj` | `CarlsonRJ` | | missing | same |
| `elliprf` | `CarlsonRF` | | missing | same |
| `elliprd` | `CarlsonRD` | | missing | same |
| `elliprg` | `CarlsonRG` | | missing | same |
| `agm` | `ArithmeticGeometricMean` | | missing | `Sqrt` (have), complex for general arguments |

## Elliptic functions, theta, modular, Weierstrass

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `jtheta` | `EllipticTheta` | | missing | complex type for general `z` |
| `ellipfun` | `JacobiSN`, `JacobiCN`, `JacobiDN`, ... | | missing | `jtheta` |
| `qfrom` | `EllipticNomeQ` | | missing | `ellipk` |
| `mfrom` | `InverseEllipticNomeQ` | | missing | `jtheta` |
| `kfrom` | | | missing | `mfrom` |
| `taufrom` | | | missing | complex type |
| `qbarfrom` | | | missing | complex type |
| `g2g3from` (master) | `WeierstrassInvariants` | | missing | complex type |
| `omega1omega2from` (master) | `WeierstrassHalfPeriods` | | missing | complex type |
| `rtheta` (master) | | | missing | `jtheta` |
| `rtheta_jet` (master) | | | missing | `rtheta` |
| `kleinj` | `KleinInvariantJ` | | missing | `eta`, complex type |
| `kleinjinv` (master) | `ModularLambda`-related | | missing | `kleinj` |
| `eta` | `DedekindEta` | | missing | complex type |
| `weierp` (master) | `WeierstrassP` | | missing | complex type |
| `weierpprime` (master) | `WeierstrassPPrime` | | missing | same |
| `weiersigma` (master) | `WeierstrassSigma` | | missing | same |
| `weierzeta` (master) | `WeierstrassZeta` | | missing | same |
| `weierpinv` (master) | `InverseWeierstrassP` | | missing | same |

## Orthogonal polynomials and other

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `jacobi` | `JacobiP` | | missing | `hyp2f1` |
| `chebyt` | `ChebyshevT` | | missing | `hyp2f1` or recurrence |
| `chebyu` | `ChebyshevU` | | missing | same |
| `legendre` | `LegendreP` | | missing | same |
| `legenp` | `LegendreP[n,m,x]` | | missing | `hyp2f1` |
| `legenq` | `LegendreQ` | | missing | same |
| `hermite` | `HermiteH` | | missing | `hyperu` / recurrence |
| `pcfd` | `ParabolicCylinderD` | | missing | `hyperu` |
| `pcfu` | | | missing | `pcfd` |
| `pcfv` | | | missing | `pcfd` |
| `pcfw` | | | missing | `pcfd` |
| `gegenbauer` | `GegenbauerC` | | missing | `hyp2f1` |
| `laguerre` | `LaguerreL` | | missing | `hyp1f1` |
| `spherharm` | `SphericalHarmonicY` | | missing | `legenp` |

## Bessel and related

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `besselj` | `BesselJ` | | missing | `hyp0f1`, asymptotic for large arguments |
| `j0`, `j1` | `BesselJ[0,x]`, `BesselJ[1,x]` | | missing | `besselj` |
| `besseli` | `BesselI` | | missing | `hyp0f1` |
| `bessely` | `BesselY` | | missing | `besselj`, limit for integer order |
| `besselk` | `BesselK` | | missing | same |
| `besseljzero` | `BesselJZero` | | missing | `besselj`, McMahon start |
| `besselyzero` | `BesselYZero` | | missing | `bessely` |
| `spherical_jn` (name from mpmath 1.4) | `SphericalBesselJ` | | missing | `besselj` |
| `spherical_yn` | `SphericalBesselY` | | missing | `bessely` |
| `spherical_in` (master) | `SphericalBesselI` | | missing | `besseli` |
| `spherical_kn` (master) | | | missing | `besselk` |
| `hankel1` | `HankelH1` | | missing | complex type |
| `hankel2` | `HankelH2` | | missing | complex type |
| `struveh` | `StruveH` | | missing | `hyper` |
| `struvel` | `StruveL` | | missing | `hyper` |
| `angerj` | `AngerJ` | | missing | `hyper` |
| `webere` | `WeberE` | | missing | `hyper` |
| `lommels1` | `LommelS1`-like | | missing | `hyper` |
| `lommels2` | | | missing | `hyper` |
| `whitm` | `WhittakerM` | | missing | `hyp1f1` |
| `whitw` | `WhittakerW` | | missing | `hyperu` |
| `ber` | `KelvinBer` | | missing | `besselj`, complex type |
| `bei` | `KelvinBei` | | missing | same |
| `ker` | `KelvinKer` | | missing | same |
| `kei` | `KelvinKei` | | missing | same |
| `coulombc` | | | missing | `loggamma` |
| `coulombf` | `CoulombF` | | missing | `hyp1f1`, complex type |
| `coulombg` | `CoulombG` | | missing | same |

## Waveforms, sigmoid, fft

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `squarew` | `SquareWave` | | missing | none |
| `trianglew` | `TriangleWave` | | missing | none |
| `sawtoothw` | `SawtoothWave` | | missing | none |
| `unit_triangle` | | | missing | none |
| `sigmoid` | `LogisticSigmoid` | | missing | `Exp` (have) |
| `fft` (master) | `Fourier` | | missing | complex type, go-fft over `big.Float` (a generic FFT) |
| `invfft` (master) | `InverseFourier` | | missing | same |

## Tools (not functions of one argument)

These are mpmath's calculus, solver, linear algebra, precision and I/O
tools. They are outside what `bigfloat` is (elementary functions on
`big.Float`); they are listed so the inventory is of the whole
`mpmath/__init__.py`. Each row is a family, one backticked name per mpmath
name. All missing; the dependency is the package where it would live.

| mpmath | Mathematica | Go | Status | Missing / dep |
|---|---|---|---|---|
| `mp`, `fp`, `iv`, `extraprec`, `extradps`, `workprec`, `workdps`, `autoprec`, `maxcalls`, `memoize` | `$MinPrecision`, `WorkingPrecision`, `Block` | `big.Float.SetPrec` per value, `Ziv` loop internal | missing | no global context by design; `autoprec` / `workprec` equivalents are a helper that re-evaluates until stable |
| `mpi` | `Interval`, `Around` | midpoint-radius ball in `ball.go` (unexported) | missing | export a ball type (Arb style); `iv` arithmetic and function enclosures |
| `quad`, `quadgl`, `quadts`, `quadosc`, `quadsubdiv`, `gauss_quadrature` | `NIntegrate` | | missing | `go-numerics/integrate` is float64 only; needs a big.Float quadrature (tanh-sinh, Gauss-Legendre) |
| `diff`, `diffs`, `diffs_prod`, `diffs_exp`, `diffun`, `differint`, `difference`, `taylor`, `pade`, `limit` | `ND`, `D` (numeric), `Series`, `PadeApproximant`, `NLimit` | | missing | none (finite differences at raised precision) |
| `nsum`, `nprod`, `sumem`, `sumap`, `richardson`, `shanks`, `levin`, `cohen_alt` | `NSum`, `NProduct`, `SequenceLimit` | | missing | none |
| `odefun`, `jacobian` | `NDSolve` | | missing | Taylor-series ODE; `jacobian` needs matrices |
| `invertlaplace`, `invlaptalbot`, `invlapstehfest`, `invlapdehoog` | `InverseLaplaceTransform` | | missing | complex type |
| `findroot`, `multiplicity`, `polyroots`, `polyval` | `FindRoot`, `NSolve`, `Roots` | | missing | complex type for `polyroots` |
| `fourier`, `fourierval`, `chebyfit` | `FourierCoefficient`, `MiniMaxApproximation` | | missing | `quad` |
| `pslq`, `identify`, `findpoly` | `FindIntegerNullVector`, `Recognize`, `RootApproximant` | | missing | `big.Int` lattice code; `lll_gram` below |
| `matrix`, `eye`, `diag`, `zeros`, `ones`, `hilbert`, `randmatrix`, `swap_row`, `extend`, `unitvector` | `Table`, `IdentityMatrix`, `DiagonalMatrix`, `HilbertMatrix` | | missing | a matrix over `big.Float`; go-ndarray has no big dtype |
| `norm`, `mnorm`, `lu`, `lu_solve`, `qr`, `qr_solve`, `inverse`, `pinv` (master), `residual`, `cholesky`, `cholesky_solve`, `det`, `cond`, `rank` | `Norm`, `LUDecomposition`, `QRDecomposition`, `Inverse`, `PseudoInverse`, `CholeskyDecomposition`, `Det`, `MatrixRank` | | missing | the matrix type above |
| `hessenberg`, `schur`, `eig`, `eig_sort`, `eigsy`, `eighe`, `eigh`, `svd_r`, `svd_c`, `svd` | `HessenbergDecomposition`, `SchurDecomposition`, `Eigensystem`, `SingularValueDecomposition` | | missing | matrix type, complex type |
| `expm`, `sqrtm`, `powm`, `logm`, `sinm`, `cosm` | `MatrixExp`, `MatrixPower`, `MatrixLog`, `MatrixFunction` | | missing | matrix type, complex type |
| `lll_gram` (master) | `LatticeReduce` | | missing | `big.Int` LLL |
| `nstr`, `nprint`, `chop`, `mpf`, `mpmathify`, `convert`, `rand`, `nan`, `inf`, `ninf`, `eps`, `j`, `fraction`, `linspace`, `arange`, `plot`, `cplot`, `splot`, `monitor`, `timing` | `NumberForm`, `Chop`, `RandomReal`, `Infinity`, `$MachineEpsilon`, `Range`, `Plot` | `big.Float.Text`, `SetString`, `SetPrec` | missing | formatting with a digit count (`Text('g', n)` covers it); randomness; plots out of scope |

## Mathematica: delta over mpmath

Arbitrary-precision and number-theory functions Mathematica has and mpmath
does not (or not as a named function). From knowledge; re-check against the
Wolfram Language reference before use.

| Mathematica | Go | Status | Dep |
|---|---|---|---|
| `N[expr, p]` (precision tracking, significance arithmetic) | | missing | a ball type carrying the radius (as Arb) |
| `Precision`, `Accuracy`, `SetPrecision` | `big.Float.Prec`, `SetPrec` | partial | per-value precision exists; no accuracy tracking |
| `RealDigits`, `IntegerDigits`, `FromDigits` | `big.Float.Text`, `big.Int.Text` | partial | digit-list API missing |
| `ContinuedFraction`, `FromContinuedFraction`, `Rationalize` | | missing | `big.Rat` |
| `PowerMod` | `big.Int.Exp` | done (math/big) | |
| `ModularInverse` | `big.Int.ModInverse` | done (math/big) | |
| `GCD`, `LCM`, `ExtendedGCD` | `big.Int.GCD` (with Bezout) | done (math/big) | LCM trivial |
| `JacobiSymbol`, `KroneckerSymbol` | `big.Jacobi` | done (math/big) | Kronecker missing |
| `PrimeQ` | `big.Int.ProbablyPrime` | done (math/big) | proof (`ProvablePrimeQ`) missing |
| `NextPrime`, `Prime`, `PrimePi` | | missing | sieve, Miller-Rabin loop |
| `FactorInteger`, `Divisors`, `DivisorSigma`, `EulerPhi`, `MoebiusMu`, `CarmichaelLambda`, `MultiplicativeOrder`, `PrimitiveRoot` | | missing | factorisation (Pollard rho, ECM, SIQS) |
| `ChineseRemainder` | | missing | none |
| `SquareFreeQ`, `PerfectNumberQ` | | missing | factorisation |
| `Multinomial`, `CatalanNumber`, `PartitionsP`, `PartitionsQ`, `LucasL` | | missing | none, `big.Int` |
| `Zeta[s]` derivative, `ZetaZero` | | missing | `zeta` |
| `InverseGammaRegularized`, `InverseBetaRegularized`, `InverseErfc` | | missing | `gammainc`, `betainc`, `erfc` |
| `GammaRegularized`, `BetaRegularized` | | missing | same |
| `Gudermannian`, `InverseGudermannian`, `Haversine`, `InverseHaversine` | | missing | none |
| `Dilog`-style (`PolyLog[2,x]`) | | missing | `polylog` |
| `InverseJacobiSN` and the other inverse Jacobi functions | | missing | `ellipfun` |
| `MathieuC`, `MathieuS`, `SpheroidalPS`, `HeunG`, `ZernikeR` | | missing | new algorithms |
| `WignerD`, `ClebschGordan`, `ThreeJSymbol` | | missing | exact rational arithmetic |
| `FindMinimum`, `NMinimize`, `NSolve` at high precision | | missing | go-numerics/optimize over `big.Float` |
| `Round`, `Floor`, `Ceiling` to a multiple | | missing | utilities section |
| `Interval` arithmetic, `IntervalIntersection` | | missing | ball type |

## MPFR and Arb capabilities (names only)

Not counted. From knowledge; the MPFR manual and the Arb documentation were
not opened.

| Library | Capability | Go | Status |
|---|---|---|---|
| MPFR | Correct rounding, six modes plus faithful, `inexact` flag | all six modes, `Acc()` | done |
| MPFR | `mpfr_sin_cos`, `sinh_cosh`, `sec`, `csc`, `cot`, `exp2`, `exp10`, `expm1`, `log1p`, `cbrt`, `rootn_ui`, `hypot`, `atan2` | most elementary | partial (missing: `sin_cos`, `sinh_cosh`, `sec` family, `exp10`, `exp2`, `rec_sqrt`) |
| MPFR | `sinpi`, `cospi`, `tanpi`, `asinpi`, `acospi`, `atanpi` (4.2+) | | missing |
| MPFR | `compound`, `powr`, `pown` (4.2+) | | missing |
| MPFR | `gamma`, `lngamma`, `lgamma`, `digamma`, `zeta`, `zeta_ui`, `erf`, `erfc`, `eint`, `li2`, `ai`, `j0`, `j1`, `jn`, `y0`, `y1`, `yn`, `agm`, `fac_ui` | | missing |
| MPFR | `fma`, `fms`, `fmma`, `sum`, `dot` (exactly rounded) | | missing |
| MPFR | `fmod`, `remainder`, `remquo`, `rint`, `round`, `trunc`, `ceil`, `floor`, `frac`, `modf`, `nextabove`, `nextbelow`, `nexttoward` | | missing |
| MPFR | Exponent range `emin`/`emax`, `subnormalize`, exception flags (underflow, overflow, nan, inexact, erange) | big.Float exponent range fixed; overflow gives Inf | missing |
| MPFR | `const_pi`, `const_log2`, `const_euler`, `const_catalan` | `Pi`, `Ln2`, `EulerGamma` | partial (no `Catalan`) |
| MPFR | `mpfr_get_str`, `mpfr_set_str`, `printf` | `big.Float.Text`, `SetString` | done (math/big) |
| MPFR | `mpfr_ubf` (unbounded floats) | | missing |
| MPFR | `mpfr_mul_2si`, `mpfr_div_2si`, `get_d_2exp` | `SetMantExp` | done (math/big) |
| MPFR / MPC | Complex numbers (`mpc_t`): all elementary functions, correctly rounded | | missing |
| Arb | Real ball `arb_t` (midpoint-radius) with certified error | `ball.go`, unexported | partial |
| Arb | Complex ball `acb_t` | | missing |
| Arb | `arb_poly`, `acb_poly` (series, composition, roots) | | missing |
| Arb | `arb_mat`, `acb_mat` (LU, solve, det, eigenvalues, exp) | | missing |
| Arb | `acb_hypgeom` (all hypergeometric, Bessel, Airy, gamma, erf, Ei, Ci, ...) | | missing |
| Arb | `acb_calc` rigorous integration, root isolation | | missing |
| Arb | `acb_elliptic`, `acb_modular`, `acb_dirichlet` (zeta, L-functions, characters, Riemann-Siegel) | | missing |
| Arb | `arb_fmpz_poly`, `acb_theta`, `arb_fpwrap` (double-precision wrappers) | | missing |
| Arb | `bernoulli`, `partitions` (Hardy-Ramanujan-Rademacher), `arith` module (number-theoretic functions) | | missing |
| FLINT (Arb's base) | `fmpz`, `fmpq`, `nmod`, `fmpz_poly`, `fmpz_mat`, factoring, primality proving | `math/big` partial | partial |
| Calcium | exact real and complex numbers (symbolic) | | missing |

## Proposed slices

1. **Elementary completions** (no dependency): `sec`/`csc`/`cot` and the
   hyperbolic and inverse families, `cospi`/`sinpi`/`sinc`/`sincpi`,
   `cos_sin`, `exp2`, `powm1`, `phi`, `degree`, `catalan`, `apery`,
   `log(x, b)`. Each needs its own Ziv loop to stay correctly rounded.
2. **Utilities** (no dependency): `floor`, `ceil`, `nint`, `frac`, `fmod`,
   `fsum`, `fdot`, `fprod`, `ulp`, `degrees`, `radians`, the predicates.
3. **Complex type** (new package `bigcomplex` or here): `mpc`, then complex
   `exp`, `log`, `sqrt`, `power`, trigonometric and hyperbolic functions.
   Everything below that returns complex values waits on it.
4. **Gamma family kernel**: Bernoulli numbers on `big.Int`, then `gamma`,
   `loggamma`, `rgamma`, `factorial`, `beta`, `psi`, `polygamma`,
   `harmonic`, `binomial`, `rf`/`ff`. Dependency: Bernoulli numbers.
5. **Zeta family**: `zeta`, `altzeta`, `hurwitz`, `polylog`, `lerchphi`,
   `stieltjes`, `bernpoly`, `euler*`. Dependency: slice 4.
6. **Hypergeometric engine**: `hyper` (series, asymptotic, transformations)
   then `hyp0f1`..`hyp3f2`, `hyperu`, orthogonal polynomials, Bessel, Airy,
   `erf` family, `gammainc`, exponential integrals. Dependency: slice 3 for
   complex arguments, slice 4 for parameters.
7. **Elliptic**: `agm`, `ellipk`/`ellipe`, Carlson `rf/rc/rd/rj/rg`, then
   theta, Jacobi, Weierstrass, modular. Dependency: slice 3.
8. **Ball type exported** (Arb style): `mpi` / `iv` equivalent. Dependency:
   `ball.go`, already inside.
9. **Number theory on `big.Int`**: `bell`, `stirling1/2`, `fib`, primes,
   factorisation, `PowerMod` family wrappers. No dependency.
10. **Tools** (`quad`, `findroot`, `odefun`, `nsum`, matrices): belong in
    `go-numerics/integrate`, `optimize` and `go-ndarray` once those accept a
    `big.Float` scalar; not in this package.

## Counting

The totals were produced by a small Go program that reads this file and
counts, for each `##` section, the backticked names in the first column by
the Status column. The program is not part of the repository (this change is
docs only); recount after editing a row.
