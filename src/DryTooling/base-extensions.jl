# -*- coding: utf-8 -*-

function Base.length(p::AbstractMatrixProblem)
    return length(p.x)
end

function Base.length(p::PorosityDescriptor)
    return length(p.ϕ)
end
