# -*- coding: utf-8 -*-

struct RadialUserDefinedGrid1D <: AbstractGrid1D
    "Number of cells in domain."
    N::Int64

    "Radial coordinates of cells centers."
    r::Vector{Float64}

    "Radial coordinates of cells walls."
    w::Vector{Float64}

    function RadialUserDefinedGrid1D(r::Vector{Float64})
        w = vcat(0.0, (r[2:end] + r[1:end-1])/2, last(r))
        return new(length(r)-1, r, w)
    end
end

struct RadialEquispacedGrid1D <: AbstractGrid1D
    "Number of cells in domain."
    N::Int64

    "Radial coordinates of cells centers."
    r::Vector{Float64}

    "Radial coordinates of cells walls."
    w::Vector{Float64}

    function RadialEquispacedGrid1D(R::Float64, N::Int64)
        δr = R / N
        r = collect(0.0:δr:R)
        w = vcat(0.0, collect(0.5δr:δr:R-0.5δr), R)
        return new(N, r, w)
    end
end

struct RadialGeometricGrid1D <: AbstractGrid1D
    "Number of cells in domain."
    N::Int64

    "Radial coordinates of cells centers."
    r::Vector{Float64}

    "Radial coordinates of cells walls."
    w::Vector{Float64}

    function RadialGeometricGrid1D(R::Float64, N::Int64)
        α = (R + 1)^(1 / N)
        r = @. α^(0:N) - 1.0
        w = vcat(0.0, (r[2:end] + r[1:end-1])/2, R)
        return new(N, r, w)
    end
end
