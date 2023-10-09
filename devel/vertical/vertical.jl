# -*- coding: utf-8 -*-
using CairoMakie
using Polynomials
using YAML

struct Substance
    "Mean molecular mass [g/mol]"
    W::Float64

    "Viscosity polynomial [Pa.s]"
    μ::Polynomial{Float64, :T}

    "Thermal conductivity polynomial [W/(m.K)]"
    k::Polynomial{Float64, :T}

    "Specific heat polynomial [J/(kg.K)]"
    c::Polynomial{Float64, :T}

    function Substance(;
            W::Float64,
            μ::Vector{Float64},
            k::Vector{Float64},
            c::Vector{Float64}
        )
        return new(
            W,
            Polynomial(μ, :T),
            Polynomial(k, :T),
            Polynomial(c, :T)
        )
    end
end

function Substance(d::Dict{Any, Any})
    return Substance(;
        W = d["W"],
        μ = d["mu"],
        k = d["kg"],
        c = d["cp"]
    )
end

struct Mixture
    "Number of components in system."
    K::Int64

    "Storage of substance objects."
    s::Vector{Substance}

    "Mixture temperature [K]"
    T::Float64

    "Mixture pressure [Pa]"
    P::Float64

    "Mixture mass fractions [-]"
    Y::Vector{Float64}

    function Mixture(; substances::Vector{Substance})
        return new(substances)
    end
end

function Mixture(d::Dict{Any, Any})
    return Mixture(; substances = map(Substance, values(d)))
end

# function settemperature(m::Mixture, T::Float64)
#     m.T = T
# end

data = YAML.load_file("mixtures.yaml")

mix = Mixture(data)