# -*- coding: utf-8 -*-
using CairoMakie
using Polynomials
using YAML

const ZEROCELSIUS::Float64 = 273.15

const ONEATM::Float64 = 101_325.0

const GASCONSTANT::Float64 = 8314.462_618_153_24

abstract type Substance end
abstract type Mixture end

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

    "Gas mixture mass fractions [-]"
    Y::Vector{Float64}

    "Gas mixture temperature [K]"
    T::Base.RefValue{Float64}

    "Gas mixture pressure [Pa]"
    P::Base.RefValue{Float64}

    "Gas mixture pressure [Pa]"
    W::Base.RefValue{Float64}

    function GasMixture(;
            components::Vector{GasComponent},
            T::Float64 = ZEROCELSIUS,
            P::Float64 = ONEATM,
            Y::Vector{Float64} = Float64[]
        )
        K = length(components)

        if length(Y) != K
            @warn "Bad size of mass fractions array... resetting!"
            Y = zeros(Float64, K)
            Y[end] = 1.0
        end

        W = meanmolecularmass(Y, map(x->getproperty(x, :W), components))

        return new(K, components, Y, Ref(T), Ref(P), Ref(W))
    end
end

function GasComponent(d::Dict{Any, Any})
    return GasComponent(; W = d["W"], μ = d["mu"], k = d["kg"], c = d["cp"])
end

function GasMixture(d::Dict{Any, Any})
    return GasMixture(; components = map(GasComponent, values(d)))
end

function meanmolecularmass(Y::Vector{Float64}, W::Vector{Float64})::Float64
    return 1.0 / sum(y / w for (y, w) in zip(Y, W))
end

function meanmolecularmass(m::GasMixture)::Float64
    return 1.0 / sum(y / c.W for (y, c) in zip(m.Y, m.s))
end

function setstate(m::GasMixture, T::Float64, P::Float64, Y::Vector{Float64})
    @assert m.K == length(Y) "Mass fractions shape is wrong!"
    @assert sum(Y) ≈ 1.0     "Mass fractions do not add up to unit!"
    m.Y[1:end] = Y
    m.T[] = T
    m.P[] = P
    m.W[] = meanmolecularmass(m)
    return nothing
end

function densitymass(m::GasMixture)
    return m.P[] * m.W[] / (GASCONSTANT * m.T[])
end

mix = GasMixture(YAML.load_file("mixtures.yaml"))

setstate(mix, ZEROCELSIUS, ONEATM, [1.0, 0.0, 0.0])
densitymass(mix)

setstate(mix, 1000.0+ZEROCELSIUS, ONEATM, [1.0, 0.0, 0.0])
densitymass(mix)
