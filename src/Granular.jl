# -*- coding: utf-8 -*-
module Granular

using CairoMakie
using DifferentialEquations: ODEProblem, Tsit5
using DifferentialEquations: solve
using Distributions
using ModelingToolkit
using Printf
using Random
using Trapz: trapz

export PackedBedPorosityDescriptor
export SymbolicLinearKramersModel
export RotaryKilnBedSolution
export solvelinearkramersmodel
export plotlinearkramersmodel

struct PackedBedPorosityDescriptor
    """
        PackedBedPorosityDescriptor

    Provides description of porosity parameters with stochastic behavior.

    $(TYPEDFIELDS)
    """

    "Porosity volume fraction in medium [-]."
    ϕ::Union{Float64, Vector{Float64}}

    "Characteristic particle size in medium [m]."
    l::Union{Float64, Vector{Float64}}

    "Optional standard deviation of porosity volume fraction  [-]."
    σϕ::Union{Float64, Nothing}

    "Optional standard deviation of characteristic particle size [m]."
    σl::Union{Float64, Nothing}

    "Perimeter in reactor cross-section [m]."
    P::Union{Float64, Vector{Float64}}

    "Characteristic diameter of porosity channels [m]."
    D::Union{Float64, Vector{Float64}}

    "Reactor area used for scaling perimeter [m²]."
    A::Float64

    function PackedBedPorosityDescriptor(;
            ϕ::Float64,
            l::Float64,
            σϕ::Union{Float64, Nothing} = nothing,
            σl::Union{Float64, Nothing} = nothing,
            N::Union{Int64, Nothing} = nothing,
            ϕlims::Tuple{Float64, Float64} = (0.4, 0.8),
            llims::Tuple{Float64, Float64} = (0.0, 0.3),
            seed::Int64 = 42,
            area::Float64 = 1.0
        )
        if !any(isnothing, [σϕ, σl, N])
            Random.seed!(seed)
            ϕᵤ = rand(truncated(Normal(ϕ, σϕ), ϕlims...), N)
            lᵤ = rand(truncated(Normal(l, σl), llims...), N)
        else
            ϕᵤ = ϕ
            lᵤ = l
        end

        P = @. 6.0 * (1.0 - ϕᵤ) / lᵤ
        D = @. 4.0 * ϕᵤ / P

        return new(ϕᵤ, lᵤ, σϕ, σl, area * P, D, area)
    end
end

struct SymbolicLinearKramersModel
    """
        SymbolicLinearKramersModel

    Creates a reusable linear Kramers model for rotary kiln simulation.

    Implements the ordinary differential equation for prediction of bed
    height profile in a rotary kiln as proposed by Kramers and Croockewite
    (1952) [^1]. Its goal is to be used as a process support tool or to
    integrate more complex models requiring integration of the bed profile.

    [^1]: [Kramers et al., 1952](https://doi.org/10.1016/0009-2509(52)87019-8)

    $(TYPEDFIELDS)
    """

    "Symbolic kiln internal radius"
    R::Num

    "Symbolic kiln feed rate"
    Φ::Num

    "Symbolic kiln rotation rate"
    ω::Num

    "Symbolic kiln slope"
    β::Num

    "Symbolic solids repose angle"
    γ::Num

    "Symbolic kiln axial coordinates"
    z::Num

    "Symbolic bed height profile"
    h::Num

    "Problem ordinary differential equation"
    sys::ODESystem

    """ Symbolic model constructor. """
    function SymbolicLinearKramersModel()
        # Declare symbols and unknowns.
        @parameters z
        @parameters R Φ ω β γ
        @variables h(z)

        # Declare a derivative.
        Dz = Differential(z)

        # Compose problem right-hand side.
        C = (3//4) * tan(γ) * Φ / (π * R^3 * ω)
        f = C * ((h / R) * (2 - h / R))^(-3//2)

        # *Stack* equation.
        eqs = Dz(h) ~ f - tan(β) / cos(γ)

        # Assembly system for solution.
        @named sys = ODESystem(eqs)
        sys = structural_simplify(sys)

        return new(R, Φ, ω, β, γ, z, h, sys)
    end
end

struct RotaryKilnBedSolution
    """
        RotaryKilnBedSolution

    Description of a rotary kiln bed geometry computed from the solution
    of bed height along the kiln length. The main goal of the quantities
    computed here is their use with heat and mass transfer models for the
    simulation of rotary kiln process.

    Internal elements are initialized through the following constructor:

        RotaryKilnBedSolution(
            z::Vector{Float64},
            h::Vector{Float64},
            R::Float64,
            Φ::Float64
        )

    Where parameters are given as:

        - `z`: solution coordinates over length, [m].
        - `h`: bed profile solution over length, [m].
        - `R`: kiln internal radius, [m].
        - `Φ`: kiln feed rate, [m³/s].

    $(TYPEDFIELDS)
    """

    "Solution coordinates [m]"
    z::Vector{Float64}

    "Solution bed height [m]"
    h::Vector{Float64}

    "View angle from kiln center [rad]"
    θ::Vector{Float64}

    "Bed-freeboard cord length [m]"
    l::Vector{Float64}

    "Local bed cross section area [m²]"
    A::Vector{Float64}

    "Local loading based on height [-]"
    η::Vector{Float64}

    "Mean loading of kiln [%]"
    ηₘ::Float64

    "Bed integral volume [m³]"
    V::Float64

    "Residence time of particles"
    τ::Float64

    "Kiln slope [rad]"
    β::Float64

    function RotaryKilnBedSolution(z, h, β, R, Φ)
        L = z[end]
        θ = @. 2acos(1 - h / R)
        l = @. 2R * sin(θ / 2)
        A = @. (θ * R^2 - l * (R - h)) / 2
        η = @. (θ - sin(θ)) / 2π
        ηₘ = 100trapz(z, η) / L

        # Integrate mid-point volume approximation.
        Aₘ = (1//2) * (A[1:end-1] + A[2:end])
        δz = z[2:end] - z[1:end-1]
        V = sum(@. Aₘ * δz)

        # Residence time is bed volume divided by flow rate.
        τ = V  / Φ

        # Construct solution object.
        return new(z, h, θ, l, A, η, ηₘ, V, τ, β)
    end
end

function Base.length(p::PackedBedPorosityDescriptor)
    return length(p.ϕ)
end

function solvelinearkramersmodel(;
        model::SymbolicLinearKramersModel,
        L::Float64,
        R::Float64,
        Φ::Float64,
        ω::Float64,
        β::Float64,
        γ::Float64,
        d::Float64,
        solver::Any = Tsit5(),
        rtol::Float64 = 1.0e-08,
        atol::Float64 = 1.0e-08
    )
    """
        solvelinearkramersmodel(;
            model::SymbolicLinearKramersModel,
            L::Float64,
            R::Float64,
            Φ::Float64,
            ω::Float64,
            β::Float64,
            γ::Float64,
            d::Float64,
            solver::Any = Tsit5(),
            rtol::Float64 = 1.0e-08,
            atol::Float64 = 1.0e-08
        )

    Integrates an instance of `SymbolicLinearKramersModel`.

    **Important:** inputs must be provided in international system (SI) units
    as a better physical practice. The only exception is the rotation rate `ω`
    provided in revolution multiples. If the discharge end is held by a dam,
    its height must be provided instead of the particle size, as it is used
    as the ODE initial condition.

    - `model`: a symbolic kiln model.
    - `L`: kiln length, [m].
    - `R`: kiln internal radius, [m].
    - `Φ`: kiln feed rate, [m³/s].
    - `ω`: kiln rotation rate, [rev/s].
    - `β`: kiln slope, [rad].
    - `γ`: solids repose angle, [rad].
    - `d`: particle size or dam height, [m].
    """
    h = [model.h => d]
    p = [model.R => R,
         model.Φ => Φ,
         model.ω => ω,
         model.β => β,
         model.γ => γ]

    prob = ODEProblem(model.sys, h, (0.0, L), p, jac = true)
    sol = solve(prob, solver, reltol = rtol, abstol = atol)
    bed = RotaryKilnBedSolution(sol.t, sol[1, :], β, R, Φ)
    return bed
end

function plotlinearkramersmodel(
        model::RotaryKilnBedSolution;
        normz::Bool = false,
        normh::Bool = false
    )::Figure
    """
        plotlinearkramersmodel(
            model::SolutionLinearKramersModel;
            normz::Bool = false,
            normh::Bool = false
        )::Any

    Display plot of model solution for rotary kiln bed profile. Arguments
    `normz` and `normh` control whether z-coordinate and bed height must
    be normalized, respectively.
    """
    z = model.z
    h = tan(model.β) * z + model.h

    z = normz ? (100z / maximum(z[end])) : z
    h = normh ? (100h / maximum(h[end])) : 100h

    unitz = normz ? "%" : "m"
    unith = normh ? "%" : "cm"

    η = @sprintf("%.1f", model.ηₘ)
    τ = @sprintf("%.0f", model.τ / 60)

    title  = "Loading $(η)% | Residence $(τ) min"
    xlab  = "Coordinate [$(unitz)]"
    ylab  = "Bed height [$(unith)]"

    xlims = (normz) ? (0.0, 100.0) : (0.0, model.z[end])
    ylims = (normh) ? (0.0, 100.0) : (0.0, round(maximum(h)+1))

    fig = Figure()
    ax = Axis(fig[1, 1], title = title, xlabel = xlab, ylabel = ylab,
                    xticks = range(xlims..., 6), yticks = range(ylims..., 6))
    lines!(ax, z, h, color = :red, label = "Profile")
    limits!(ax, xlims, ylims)
    axislegend(position = :lt)

    return fig
end

end # module Granular