// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package main

import (
	"math/big"
	"math/rand/v2"

	"github.com/go-bignum/bigfloat"
)

// gen draws the arguments of one case at precision p.
type gen func(r *rand.Rand, p uint) []*big.Float

// region is a named argument distribution of a function.
type region struct {
	name string
	gen  gen
}

// randInt returns a uniform integer in [0, 2^bits).
func randInt(r *rand.Rand, bits uint) *big.Int {
	m := new(big.Int)
	for i := uint(0); i < bits; i += 64 {
		w := new(big.Int).SetUint64(r.Uint64())
		m.Or(m, w.Lsh(w, i))
	}
	return m.Rsh(m, (bits+63)/64*64-bits)
}

// mant returns a random p-bit number in [1/2, 1).
func mant(r *rand.Rand, p uint) *big.Float {
	m := randInt(r, p)
	m.SetBit(m, int(p-1), 1)
	x := new(big.Float).SetPrec(p).SetInt(m)
	return x.SetMantExp(x, -int(p))
}

// logu returns ±m·2^e, e uniform in [e1, e2]; sign 0 picks either.
func logu(r *rand.Rand, p uint, e1, e2 int, sign int) *big.Float {
	x := mant(r, p)
	x.SetMantExp(x, e1+r.IntN(e2-e1+1))
	if sign < 0 || sign == 0 && r.IntN(2) == 0 {
		x.Neg(x)
	}
	return x
}

// uni returns a p-bit number uniform in [lo, hi).
func uni(r *rand.Rand, p uint, lo, hi float64) *big.Float {
	u := mant(r, p)             // [1/2, 1)
	u.Sub(u, big.NewFloat(0.5)) // [0, 1/2)
	u.SetMantExp(u, 1)          // [0, 1)
	w := new(big.Float).SetPrec(p)
	w.Mul(u, new(big.Float).SetFloat64(hi-lo))
	return w.Add(w, new(big.Float).SetFloat64(lo))
}

// near returns c + d, d = ±m·2^e with e in [e1, e2], rounded to p bits
// and different from c.
func near(r *rand.Rand, p uint, c float64, e1, e2, sign int) *big.Float {
	cf := big.NewFloat(c)
	for {
		d := logu(r, p, e1, e2, sign)
		if x := new(big.Float).SetPrec(p).Add(d, cf); x.Cmp(cf) != 0 {
			return x
		}
	}
}

// ep scales an exponent bound with the precision.
func ep(p uint, f float64) int { return int(f * float64(p)) }

func one(f func(r *rand.Rand, p uint) *big.Float) gen {
	return func(r *rand.Rand, p uint) []*big.Float { return []*big.Float{f(r, p)} }
}

// tiny spans exponents from well inside the near-identity shortcut
// (exponent below about -p) to well outside it (-p/4).
func tiny(sign int) region {
	return region{"tiny", one(func(r *rand.Rand, p uint) *big.Float { return logu(r, p, -ep(p, 2.2)-8, -ep(p, 0.25), sign) })}
}

func small(sign int) region {
	return region{"small", one(func(r *rand.Rand, p uint) *big.Float { return logu(r, p, -ep(p, 0.25), 0, sign) })}
}

func uniform(name string, lo, hi float64) region {
	return region{name, one(func(r *rand.Rand, p uint) *big.Float { return uni(r, p, lo, hi) })}
}

func wide(name string, e1, e2, sign int) region {
	return region{name, one(func(r *rand.Rand, p uint) *big.Float { return logu(r, p, e1, e2, sign) })}
}

// nearC: within 2^-(p+4) .. 1/4 of c, on the side given by sign.
func nearC(name string, c float64, sign int) region {
	return region{name, one(func(r *rand.Rand, p uint) *big.Float { return near(r, p, c, -int(p)-4, -2, sign) })}
}

// nearPM1: ±(1 - d), d in [2^-p, 1/4].
func nearPM1() region {
	return region{"near±1", one(func(r *rand.Rand, p uint) *big.Float {
		x := near(r, p, 1, -int(p), -2, -1)
		if r.IntN(2) == 0 {
			x.Neg(x)
		}
		return x
	})}
}

func pow2(r *rand.Rand, p uint) *big.Float {
	x := new(big.Float).SetPrec(p).SetInt64(1)
	return x.SetMantExp(x, r.IntN(4001)-2000)
}

// pow10 returns 10^k representable at p bits (5^k needs 2.33k bits).
func pow10(r *rand.Rand, p uint) *big.Float {
	k := r.IntN(int(float64(p)/2.33) + 1)
	n := new(big.Int).Exp(big.NewInt(10), big.NewInt(int64(k)), nil)
	return new(big.Float).SetPrec(p).SetInt(n)
}

// perfectPower returns ±b^n·2^(n·s), b odd of about p/n bits: an exact
// n-th power representable at precision p.
func perfectPower(r *rand.Rand, p uint, n int, neg bool) *big.Float {
	b := randInt(r, uint(max(int(p)/n, 1)))
	b.SetBit(b, 0, 1)
	v := new(big.Int).Exp(b, big.NewInt(int64(n)), nil)
	x := new(big.Float).SetPrec(p).SetInt(v)
	x.SetMantExp(x, n*(r.IntN(41)-20))
	if neg {
		x.Neg(x)
	}
	return x
}

// nearHalfPiMultiple returns ±k·π/2 rounded to p bits, k random: the
// arguments where sin, cos or tan is closest to a zero or a pole.
func nearHalfPiMultiple(r *rand.Rand, p uint, kmax int64) *big.Float {
	w := p + 128
	h := bigfloat.Pi(new(big.Float).SetPrec(w))
	h.Mul(h, new(big.Float).SetPrec(w).SetInt64(r.Int64N(kmax)+1))
	h.SetMantExp(h, -1)
	if r.IntN(2) == 0 {
		h.Neg(h)
	}
	return new(big.Float).SetPrec(p).Set(h)
}

func intF(p uint, v int64) *big.Float { return new(big.Float).SetPrec(p).SetInt64(v) }

// regionsOf returns the argument regions of a function: tiny and small
// arguments, the bulk of the domain, huge arguments, the neighbourhood of
// singularities and special points, and the exact cases.
func regionsOf(fn string) []region {
	zero := region{"exact", one(func(_ *rand.Rand, p uint) *big.Float { return new(big.Float).SetPrec(p) })}
	oneR := region{"exact", one(func(_ *rand.Rand, p uint) *big.Float { return intF(p, 1) })}
	switch fn {
	case "exp":
		return []region{tiny(0), small(0), uniform("moderate", -20, 20), uniform("large", -1e5, 1e5), uniform("huge", -1e9, 1e9), zero}
	case "expm1":
		return []region{tiny(0), small(0), uniform("moderate", -40, 40), uniform("negative", -1e4, -40), uniform("large", 40, 1e5), zero}
	case "log", "log2", "log10":
		rs := []region{nearC("near1", 1, 0), uniform("moderate", 0, 100), wide("wide", -100000, 100000, 1), oneR}
		switch fn {
		case "log2":
			rs = append(rs, region{"exact", one(pow2)})
		case "log10":
			rs = append(rs, region{"exact", one(pow10)})
		}
		return rs
	case "log1p":
		return []region{tiny(0), small(0), nearC("near-1", -1, 1), uniform("moderate", 0, 1000), wide("wide", 1, 100000, 1), zero}
	case "sqrt", "cbrt":
		n, sign := 2, 1
		if fn == "cbrt" {
			n, sign = 3, 0
		}
		return []region{wide("wide", -100000, 100000, sign), uniform("moderate", 0, 10),
			{"exact", one(func(r *rand.Rand, p uint) *big.Float { return perfectPower(r, p, n, sign == 0 && r.IntN(2) == 0) })}}
	case "root":
		ns := []int{2, 3, 4, 5, 7, 10, 17, 100, 1001}
		return []region{
			{"wide", func(r *rand.Rand, p uint) []*big.Float {
				n := ns[r.IntN(len(ns))]
				sign := 1
				if n%2 == 1 {
					sign = 0
				}
				return []*big.Float{logu(r, p, -100000, 100000, sign), intF(p, int64(n))}
			}},
			{"exact", func(r *rand.Rand, p uint) []*big.Float {
				n := ns[r.IntN(6)]
				return []*big.Float{perfectPower(r, p, n, n%2 == 1 && r.IntN(2) == 0), intF(p, int64(n))}
			}},
			// n·(p+2) beyond the direct integer root: through exp(log(x)/n)
			{"large-n", func(r *rand.Rand, p uint) []*big.Float {
				n := int64(1<<20 + r.IntN(2))
				sign := 1
				if n%2 == 1 && r.IntN(2) == 0 {
					sign = -1 // odd n: negative arguments too
				}
				x := logu(r, p, -100000, 100000, sign)
				if r.IntN(4) == 0 { // exact: 2^(n·s)
					x = new(big.Float).SetPrec(p).SetMantExp(big.NewFloat(1), int(n)*(r.IntN(41)-20))
				}
				return []*big.Float{x, intF(p, n)}
			}},
		}
	case "hypot":
		return []region{
			{"wide", func(r *rand.Rand, p uint) []*big.Float {
				return []*big.Float{logu(r, p, -1000, 1000, 0), logu(r, p, -1000, 1000, 0)}
			}},
			{"gap", func(r *rand.Rand, p uint) []*big.Float {
				e := r.IntN(200) - 100
				return []*big.Float{logu(r, p, e, e+2, 0), logu(r, p, e-ep(p, 1.2), e, 0)}
			}},
			{"exact", func(r *rand.Rand, p uint) []*big.Float {
				bits := max(int(p)/2-1, 2)
				a := randInt(r, uint(bits)).Int64()%(1<<20) + 2
				b := randInt(r, uint(bits)).Int64()%(a-1) + 1
				return []*big.Float{intF(p, a*a-b*b), intF(p, 2*a*b)}
			}},
		}
	case "pow":
		return []region{
			{"moderate", func(r *rand.Rand, p uint) []*big.Float { return []*big.Float{uni(r, p, 0, 10), uni(r, p, -50, 50)} }},
			{"near1", func(r *rand.Rand, p uint) []*big.Float {
				// x = 1 ± 2^e·m, |y·log x| below 2^20: no overflow
				e := -4 - r.IntN(int(p)-3)
				return []*big.Float{near(r, p, 1, e, e, 0), logu(r, p, 0, 20-e, 0)}
			}},
			{"wide", func(r *rand.Rand, p uint) []*big.Float {
				return []*big.Float{logu(r, p, -1000, 1000, 1), uni(r, p, -3, 3)}
			}},
			{"negative", func(r *rand.Rand, p uint) []*big.Float {
				return []*big.Float{uni(r, p, -10, 0), intF(p, r.Int64N(81)-40)}
			}},
			{"exact", func(r *rand.Rand, p uint) []*big.Float {
				k := r.Int64N(8) + 1
				switch r.IntN(3) {
				case 0: // b^k
					b := intF(p, randInt(r, uint(max(int(p)/int(k), 2))).Int64()%(1<<40)|1)
					if r.IntN(2) == 0 {
						b.Neg(b)
					}
					return []*big.Float{b, intF(p, k)}
				case 1: // (s²)^((2k-1)/2)
					s := perfectPower(r, p, 2, false)
					y := intF(p, 2*k-1)
					return []*big.Float{s, y.SetMantExp(y, -1)}
				default: // 2^e to a negative power
					return []*big.Float{pow2(r, p), intF(p, -k)}
				}
			}},
		}
	case "sin", "cos", "tan":
		return []region{tiny(0), small(0), uniform("moderate", -100, 100), wide("large", 7, 60, 0), wide("huge", 60, 1000, 0),
			{"near-kpi/2", one(func(r *rand.Rand, p uint) *big.Float { return nearHalfPiMultiple(r, p, 1<<20) })}, zero}
	case "asin", "acos":
		return []region{tiny(0), uniform("moderate", -1, 1), nearPM1(), zero}
	case "atan":
		return []region{tiny(0), uniform("moderate", -10, 10), wide("wide", 0, 100000, 0), zero}
	case "atan2":
		return []region{
			{"wide", func(r *rand.Rand, p uint) []*big.Float {
				return []*big.Float{logu(r, p, -100, 100, 0), logu(r, p, -100, 100, 0)}
			}},
			{"ratio-tiny", func(r *rand.Rand, p uint) []*big.Float {
				return []*big.Float{logu(r, p, -ep(p, 2)-8, -ep(p, 0.3), 0), logu(r, p, 0, 2, 0)}
			}},
			{"ratio-huge", func(r *rand.Rand, p uint) []*big.Float {
				return []*big.Float{logu(r, p, 0, 2, 0), logu(r, p, -ep(p, 2)-8, -ep(p, 0.3), 0)}
			}},
		}
	case "sinh", "cosh":
		return []region{tiny(0), small(0), uniform("moderate", -50, 50), uniform("large", -1e5, 1e5), uniform("huge", -1e9, 1e9), zero}
	case "tanh":
		return []region{tiny(0), small(0), uniform("moderate", -20, 20),
			{"saturated", one(func(r *rand.Rand, p uint) *big.Float { return uni(r, p, float64(p)*0.3, float64(p)*0.45+2) })}, zero}
	case "asinh":
		return []region{tiny(0), small(0), uniform("moderate", -100, 100), wide("wide", 0, 100000, 0), zero}
	case "acosh":
		return []region{nearC("near1", 1, 1), uniform("moderate", 1, 100), wide("wide", 1, 100000, 1), oneR}
	case "atanh":
		return []region{tiny(0), small(0), nearPM1(), zero}
	}
	return []region{{"const", func(*rand.Rand, uint) []*big.Float { return nil }}}
}
