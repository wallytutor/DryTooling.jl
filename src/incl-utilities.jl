# -*- coding: utf-8 -*-
export closestpowerofx, axesunitscaler
export head, tail, body
export heaviside, interval
export makestepwise1d

"""
    closestpowerofx(
        v::Number;
        x::Number = 10,
        roundf::Function = ceil
    )::Int64

Compute the integer power of `x` closest to `v` using `roundf` as
rouding method. This might be useful for automatic setting more
reasonable limits to plot axis or similar applications. Changing
the rouding method through `roundf` is also possible.

```jldoctest
julia> closestpowerofx(12.0; x = 10)
20
julia> closestpowerofx(12.0; x = 10, roundf = floor)
10
julia> closestpowerofx(12.0; x = 10, roundf = round)
10
```
"""
function closestpowerofx(
        v::Number;
        x::Number = 10,
        roundf::Function = ceil
    )::Int64
    rounder = x^floor(log(x, v))
    return convert(Int64, rounder * roundf(v / rounder))
end

"""
    axesunitscaler(x::Number)::Tuple{String, Int64}

Find scaling factor for multiples of 1000 units. Together with
`closestpowerofx` this can be used to produce better automatic
plot axes limits. The returned values provide the string for
modifying the axis label and the associated scaling factor.

**NOTE:** this function is not yet stable. In the future it will
instead return labels using symbols like `k`, `M`, `G` for the
units through a flag provided by the user.

```jldoctest
julia> axesunitscaler(1)
("", 1)
julia> axesunitscaler(1000)
("[×1000]", 1000)
julia> axesunitscaler(1000000)
("[×1000000]", 1000000)
```
"""
function axesunitscaler(x::Number)::Tuple{String, Int64}
	# Find the floor of log10 of number.
	m = convert(Int64, x |> log10 |> floor)

	# Get the order of magnitude number.
	n = div(m, 3)

	# Find scaling factor.
	p = 1000^(n)

	return (n == 0) ? ("", 1) : ("[×$(1000^n)]", p)
end

"""
    head(z)

Access view of array head. See also ``tail`` and ``body``.

```jldoctest
julia> head(1:4)
1:3
julia> head([1, 2, 3, 4])
3-element view(::Vector{Int64}, 1:3) with eltype Int64:
 1
 2
 3
```
"""
head(z) = @view z[1:end-1]

"""
    tail(z)

Access view of array tail. See also `head` and `body`.

```jldoctest
julia> tail([1, 2, 3, 4])
3-element view(::Vector{Int64}, 2:4) with eltype Int64:
 2
 3
 4
julia> tail(1:4)
2:4
```
"""
tail(z) = @view z[2:end-0]

"""
    body(z)

Access view of array body. See also `head` and `tail`.

```jldoctest
julia> body([1, 2, 3, 4])
2-element view(::Vector{Int64}, 2:3) with eltype Int64:
 2
 3
julia> body(1:4)
2:3
```
"""
body(z) = @view z[2:end-1]

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
