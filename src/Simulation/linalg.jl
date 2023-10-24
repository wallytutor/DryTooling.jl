# -*- coding: utf-8 -*-
export TridiagonalProblem
export solve!
export change

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

Solve a tridiagonal problem and update internal memory.
"""
function solve!(p::TridiagonalProblem)::Nothing
    p.x[:] = p.A \ p.b
    return nothing
end

"""
    change(p::TridiagonalProblem)::Vector{Float64}

Compute change in solution of a tridiagonal problem without update.
"""
function change(p::TridiagonalProblem)::Vector{Float64}
    return p.A \ p.b - p.x
end

function Base.length(p::AbstractMatrixProblem)
    return length(p.x)
end
