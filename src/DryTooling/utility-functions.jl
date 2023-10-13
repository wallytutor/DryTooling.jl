# -*- coding: utf-8 -*-

""" Automatic differentiable Heaviside function. """
heaviside(t) = 0.5 * (sign(t) + 1.0)

""" Returns 1 if ``x ∈ (a, b)``, 1/2 for `` x = a || x = b``, or 0 . """
interval(x; a = -Inf, b = Inf) = heaviside(x-a) - heaviside(x-b)

"""
    makestepwise1d(lo, hi, xc)

Creates an univariate function that is composed of two parts, the first
evaluated before a critical domain point `xc`, and the seconda above that
value. This is often required, for instance, for the evaluation of NASA
polynomials for thermodynamic properties. If `differentiable`, then the
returned function is compatible with symbolic argument as required when
using package `ModelingToolkit`, etc.
"""
function makestepwise1d(lo, hi, xc; differentiable = true)
    if differentiable
        f = x -> lo(x)*interval(x, b=xc)+hi(x)*interval(x, a=xc)
    else
        f = x -> (x < xc) ? lo(x) : hi(x)
    end
    return f
end

"Compute the power of `x` closest to `v`."
function closestpowerofx(v::Number; x::Number = 10)::Number
    rounder = x^floor(log(x, v))
    return convert(Int64, rounder * ceil(v/rounder))
end
