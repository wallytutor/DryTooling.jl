# -*- coding: utf-8 -*-

struct CylinderGrid1DEquispaced <: AbstractGrid1D
    "Number of cells in domain."
    N::Int64

    "Radial coordinates of cells centers."
    r::Vector{Float64}

    "Radial coordinates of cells walls."
    w::Vector{Float64}

    function CylinderGrid1DEquispaced(R::Float64, N::Int64)
        δr = R / N
        r = collect(0.0:δr:R)
        w = vcat(0.0, collect(0.5δr:δr:R-0.5δr), R)
        return new(N, r, w)
    end
end

"Alias type for sphere equispaced discretization."
const SphereGrid1DEquispaced = CylinderGrid1DEquispaced
