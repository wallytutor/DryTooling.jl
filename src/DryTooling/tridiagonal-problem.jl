# -*- coding: utf-8 -*-

"Memory for a tridiagonal problem of size `N`."
struct TridiagonalProblem <: AbstractMatrixProblem
    A::Tridiagonal{Float64, Vector{Float64}}
    b::Vector{Float64}
    x::Vector{Float64}

    function TridiagonalProblem(N)
        A = Tridiagonal(zeros(N), zeros(N+1), zeros(N))
        return new(A, zeros(N+1), zeros(N+1))
    end
end
