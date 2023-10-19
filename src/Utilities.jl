# -*- coding: utf-8 -*-
module Utilities

export closestpowerofx, axesunitscaler
export head, tail, body
export heaviside, interval
export makestepwise1d

function closestpowerofx(
        v::Number;
        x::Number = 10,
        roundf::Function = ceil
    )::Int64
    """ Compute the power of `x` closest to `v`. """
    rounder = x^floor(log(x, v))
    return convert(Int64, rounder * roundf(v / rounder))
end

function axesunitscaler(x::Number)::Tuple{String,Int64}
    """ Find scaling factor for multiples of 1000 units. """
	# Find the floor of log10 of number.
	m = convert(Int64, x |> log10 |> floor)

	# Get the order of magnitude number.
	n = div(m, 3)

	# Find scaling factor.
	p = 1000^(n)

	return (n == 0) ? ("", 1) : ("[×$(1000^n)]", p)
end

function head(z)
    """ Access view of array head. """
    return @view z[1:end-1]
end

function tail(z)
    """ Access view of array tail. """
    return @view z[2:end-0]
end

function body(z)
    """ Access view of array body. """
    return @view z[2:end-1]
end

function heaviside(t)
    """ Automatic differentiable Heaviside function. """
    return 0.5 * (sign(t) + 1.0)
end

function interval(x; a=-Inf, b=Inf)
    """ Returns 1 if ``x ∈ (a, b)``, 1/2 for `` x = a || x = b``, or 0 . """
    return heaviside(x - a) - heaviside(x - b)
end

function makestepwise1d(lo, hi, xc; differentiable=true)
    """
        makestepwise1d(lo, hi, xc)
    
    Creates an univariate function that is composed of two parts, the first
    evaluated before a critical domain point `xc`, and the seconda above that
    value. This is often required, for instance, for the evaluation of NASA
    polynomials for thermodynamic properties. If `differentiable`, then the
    returned function is compatible with symbolic argument as required when
    using package `ModelingToolkit`, etc.
    """
    if differentiable
        f = x -> lo(x) * interval(x, b=xc) + hi(x) * interval(x, a=xc)
    else
        f = x -> (x < xc) ? lo(x) : hi(x)
    end
    return f
end

end # module Utilities