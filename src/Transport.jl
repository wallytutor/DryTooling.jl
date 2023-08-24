# -*- coding: utf-8 -*-
module Transport

export AbstractTransportModel
export LennardJonesTransport

""" Named access to transport models. """
const MODELS = Dict("gas" => :LennardJonesTransport)

""" Base type for transport models. """
abstract type AbstractTransportModel end

""" Lennard-Jones ideal gas transport model. """
struct LennardJonesTransport <: AbstractTransportModel
    geometry::String
    welldepth::Float64
    diameter::Float64
    polarizability::Float64
    rotationalrelaxation::Float64

    function LennardJonesTransport(transport)
        return new(
            transport["geometry"],
            get(transport, "well-depth", 0.0),
            get(transport, "diameter", 0.0),
            get(transport, "polarizability", 0.0),
            get(transport, "rotational-relaxation", 0.0)
        )
    end
end

end #(module Transport)