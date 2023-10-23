# -*- coding: utf-8 -*-
module FiniteVolumes

export UserDefinedGrid1D
export equidistantcellsgrid1D
export geometriccellsgrid1D

using DryTooling

struct UserDefinedGrid1D <: AbstractGrid1D
    """

    Provides a very simple interface for the 1D grids used in standard finite volume
    solvers provided in the package. Constructor accepts a vector of coordinates and
    compute walls at mid-points between cells centers. First and last cells are over
    the boundaries and the models must consider this in the implementation.
    
    $(TYPEDFIELDS)
    """

    "Number of cells in domain, size of solution memory per variable."
    N::Int64

    "Radial coordinates of cells centers."
    r::Vector{Float64}

    "Radial coordinates of cells walls."
    w::Vector{Float64}

    function UserDefinedGrid1D(r::Vector{Float64})
        w = vcat(first(r), 0.5*(tail(r) + head(r)), last(r))
        return new(length(r), r, w)
    end
end

function equidistantcellsgrid1D(R::Float64, N::Int64)::UserDefinedGrid1D
    """
        equidistantcellsgrid(R::Float64, N::Int64)::UserDefinedGrid1D
    
    Helper for creating a `UserDefinedGrid1D` with all nodes equidistant.
    """
    return UserDefinedGrid1D(collect(0.0:R/N:R))
end

function geometriccellsgrid1D(R::Float64, N::Int64)::UserDefinedGrid1D
    """
        geometriccellsgrid(R::Float64, N::Int64)::UserDefinedGrid1D
    
    Helper for creating a `UserDefinedGrid1D` with all nodes equidistant.
    """
    return UserDefinedGrid1D(@. ((R+1)^(1/N))^(0:N)-1.0)
end

end # module FiniteVolumes