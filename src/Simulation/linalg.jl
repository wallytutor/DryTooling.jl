# -*- coding: utf-8 -*-
export TridiagonalProblem
export solve!
export change
export residual

""" Memory for a tridiagonal problem of rank `N`.

All tensors are filled with zeros upon creation. This is simply
a utility for memory allocation, no other operations are made.

$(TYPEDFIELDS)
"""
struct TridiagonalProblem <: AbstractMatrixProblem

    "Main problem matrix."
    A::Tridiagonal{Float64, Vector{Float64}}

    "Right-hand side vector."
    b::Vector{Float64}

    "Solution variable vector."
    x::Vector{Float64}

    "Auxiliary vector, *e.g.* for iterative problems."
    a::Vector{Float64}

    function TridiagonalProblem(N)
        A = Tridiagonal(zeros(N-1), zeros(N), zeros(N-1))
        return new(A, zeros(N), zeros(N), zeros(N))
    end
end

"""
    solve!(p::TridiagonalProblem)::Nothing

Solve problem ``x=A^{-1}b` updating internal memory.
"""
function solve!(p::TridiagonalProblem)::Nothing
    p.x[:] = p.A \ p.b
    return nothing
end

"""
    change(p::TridiagonalProblem)::Vector{Float64}

Change in solution ``A^{-1}b-x`` of problem without update ``x``.
"""
function change(p::TridiagonalProblem)::Vector{Float64}
    return p.A \ p.b - p.x
end

"""
    residual(p::TridiagonalProblem)::Vector{Float64}

Solution residual ``b-Ax`` of problem without update of ``x``.
"""
function residual(p::TridiagonalProblem)::Vector{Float64}
    return p.b - p.A * p.x
end

function Base.length(p::AbstractMatrixProblem)
    return length(p.x)
end
