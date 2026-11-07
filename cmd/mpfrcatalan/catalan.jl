# Prints MPFR's mpfr_const_catalan (through Base.MathConstants.catalan, which
# ccalls it) at each precision of ARGS and in each rounding mode, one line per
# result:  <prec> <mode> <mantissa, hex integer of prec bits> <exponent>
# meaning value = mantissa · 2^(exponent - prec), mantissa/2^prec in [1/2, 1).
# (Everything runs at precision p: frexp and ldexp round to the default
# precision, 256 bits, otherwise.)
modes = [("N", RoundNearest), ("Z", RoundToZero), ("A", RoundFromZero),
         ("D", RoundDown), ("U", RoundUp)]
for a in ARGS
    p = parse(Int, a)
    setprecision(BigFloat, p) do
        for (name, r) in modes
            x = BigFloat(Base.MathConstants.catalan, r; precision = p)
            @assert precision(x) == p
            e = exponent(x) + 1
            mant = BigInt(ldexp(x, p - e))
            @assert ldexp(BigFloat(mant), e - p) == x
            println(p, " ", name, " ", string(mant, base = 16), " ", e)
        end
    end
end
println("# ", Base.MPFR.version())
