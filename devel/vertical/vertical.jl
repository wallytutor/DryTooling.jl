# -*- coding: utf-8 -*-
using CairoMakie
using Polynomials
using YAML

##############################################################################
# CONSTANTS
##############################################################################

const ZEROCELSIUS::Float64 = 273.15

const ONEATM::Float64 = 101_325.0

const GASCONSTANT::Float64 = 8314.462_618_153_24

##############################################################################
# ABSTRACT TYPES
##############################################################################

abstract type Substance end
abstract type Mixture end

##############################################################################
# STATELESS MODELS
##############################################################################

struct GasComponent <: Substance
    "Mean molecular mass [g/mol]"
    W::Float64

    "Viscosity polynomial [Pa.s]"
    μ::Polynomial{Float64, :T}

    "Thermal conductivity polynomial [W/(m.K)]"
    k::Polynomial{Float64, :T}

    "Specific heat polynomial [J/(kg.K)]"
    c::Polynomial{Float64, :T}

    function GasComponent(;
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

struct GasMixture <: Mixture
    "Number of components in system."
    K::Int64

    "Storage of gas component objects."
    s::Vector{GasComponent}

    function GasMixture(; components::Vector{GasComponent})
        return new(length(components), components)
    end
end

function GasComponent(d::Dict{Any, Any})
    return GasComponent(; W = d["W"], μ = d["mu"], k = d["kg"], c = d["cp"])
end

function GasMixture(d::Dict{Any, Any}; order::Any = nothing)
    order = isnothing(order) ? keys(d) : order
    components = map(k->GasComponent(d[k]), order)
    return GasMixture(; components = components)
end

##############################################################################
# METHODS
##############################################################################

function molecularmasses(m::GasMixture)::Vector{Float64}
    return map(x->x.W, m.s)
end

function meanmolecularmass(
        Y::Union{Vector{Float64},SubArray},
        W::Vector{Float64}
    )::Float64
    return 1.0 / sum(y / w for (y, w) in zip(Y, W))
end

function idealgasdensity(T::Float64, P::Float64, W::Float64)::Float64
    return P * W / (GASCONSTANT * T)
end

function idealgasdensity(
        T::Float64,
        P::Float64,
        Y::Union{Vector{Float64},SubArray};
        W::Vector{Float64}
    )::Float64
    return idealgasdensity(T, P, meanmolecularmass(Y, W))
end

function thermophysicalproperties(s::GasComponent, T::Float64)::Matrix{Float64}
    return [s.μ(T) s.k(T) s.c(T)]
end

function thermophysicalproperties(
        m::GasMixture,
        T::Float64,
        Y::Union{Vector{Float64},SubArray}
    )::Matrix{Float64}
    return sum(y*thermophysicalproperties(s, T) for (s, y) in zip(m.s, Y))
end

function mixtureproperties(
        T::Float64,
        P::Float64,
        Y::Union{Vector{Float64},SubArray};
        m::GasMixture,
        W::Vector{Float64}
    )::Matrix{Float64}
    ρ = idealgasdensity(T, P, Y; W=W)
    μ, k, c = thermophysicalproperties(m, T, Y)
    return [ρ μ k c]
end

function mixtureproperties(
        T::Vector{Float64},
        P::Vector{Float64},
        Y::Matrix{Float64};
        m::GasMixture,
        W::Vector{Float64}
    )::Vector{Matrix{Float64}}
    return mixtureproperties.(T, P, eachrow(Y); m = m, W = W)
end

##############################################################################
# DEVEL
##############################################################################

data = YAML.load_file("mixtures.yaml")
mix = GasMixture(data, order = ["fumes", "co2"])

W = molecularmasses(mix)

T = ZEROCELSIUS
P = ONEATM
Y = [1.0, 0.0]
ρ = idealgasdensity(T, P, Y; W = W)

K = 20
T = collect(range(300.0, 3000.0, K))
P = collect(range(300.0, 3000.0, K)) .+ ONEATM
Y = [range(0.0, 1.0, K) range(1.0, 0.0, K)]
ρ = idealgasdensity.(T, P, eachrow(Y); W = W)

K = 20
T = (ZEROCELSIUS + 25.0) * ones(K)
P = ONEATM * ones(K)
Y = [range(0.0, 1.0, K) range(1.0, 0.0, K)]
ρ, μ, k, c = zip(mixtureproperties(T, P, Y; m = mix, W = W)...)

K = 20
T = collect(range(300.0, 2000.0, K))
P = ONEATM * ones(K)
Y = [ones(Float64, K) zeros(Float64, K)]
ρ, μ, k, c = zip(mixtureproperties(T, P, Y; m = mix, W = W)...)
