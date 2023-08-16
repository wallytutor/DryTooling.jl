# -*- coding: utf-8 -*-
"""
    module Kramers

Implements the differential equation for prediction of bed height profile
in a rotary kiln as proposed by Kramers and Croockewite (1952) [^1].

[^1]: [Kramers et al., 1952](https://doi.org/10.1016/0009-2509(52)87019-8)
"""
module Kramers

using ModelingToolkit
using Plots
using Printf

using DifferentialEquations: ODEProblem
using DifferentialEquations: Tsit5
using DifferentialEquations: solve
using Trapz: trapz

export RotaryKilnBedGeometry
export SymbolicLinearKramersModel
export SolutionLinearKramersModel
export solvelinearkramersmodel
export plotlinearkramersmodel

"""
    RotaryKilnBedGeometry

Description of a rotary kiln bed geometry computed from the solution
of bed height along the kiln length. The main goal of the quantities
computed here is their use with heat and mass transfer models for the
simulation of rotary kiln process.
"""
struct RotaryKilnBedGeometry
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

    "Mean loading of kiln [-]"
    ηₘ::Float64

    "Bed integral volume [m³]"
    V::Float64

    """
        RotaryKilnBedGeometry(
            z::Vector{Float64},
            h::Vector{Float64},
            R::Float64,
            L::Float64
        )

    - `z`: solution coordinates over length, [m].
    - `h`: bed profile solution over length, [m].
    - `R`: kiln internal radius, [m].
    - `L`: kiln length, [m].
    """
    function RotaryKilnBedGeometry(
            z::Vector{Float64},
            h::Vector{Float64},
            R::Float64,
            L::Float64
        )
        θ = @. 2acos(1 - h / R)
        l = @. 2R * sin(θ / 2)
        A = @. (θ * R^2 - l * (R - h)) / 2
        η = @. (θ - sin(θ)) / 2π
        ηₘ = 100trapz(z, η) / L

        # Integrate mid-point volume approximation.
        Aₘ = (1//2) * (A[1:end-1] + A[2:end])
        δz = z[2:end] - z[1:end-1]
        V = sum(@. Aₘ * δz)

        return new(z, h, θ, l, A, η, ηₘ, V)
    end
end

"""
        SymbolicLinearKramersModel

    Creates a reusable linear Kramers model for rotary kiln simulation.
"""
struct SymbolicLinearKramersModel
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

    """
        SymbolicLinearKramersModel()

    Symbolic model constructor.
    """
    function SymbolicLinearKramersModel()
        # Declare symbols and unknowns.
        @parameters z
        @parameters R Φ ω β γ
        @variables h(z)

        # Declare a derivative.
        D = Differential(z)

        # Compose problem right-hand side.
        C = (3//4) * tan(γ) * Φ / (π * R^3 * ω)
        f = C * ((h / R) * (2 - h / R))^(-3//2)

        # *Stack* equation.
        eqs = D(h) ~ f - tan(β) / cos(γ)

        # Assembly system for solution.
        @named sys = ODESystem(eqs)
        sys = structural_simplify(sys)

        return new(R, Φ, ω, β, γ, z, h, sys)
    end
end

"""
    SolutionLinearKramersModel

Solve and process results from a symbolic linear Kramers model.
"""
struct SolutionLinearKramersModel
    "Object containing profile results"
    bed::RotaryKilnBedGeometry

    "Residence time of particles"
    τ::Float64

    """
        SolutionLinearKramersModel(;
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

    Integrates an instace of `SymbolicLinearKramersModel`.

    **Note:** if the discharge end is hold by a dam, its height
    must be provided instead of the particle size, as it is used
    as the ODE initial condition.

    - `L`: kiln length, [m].
    - `R`: kiln internal radius, [m].
    - `Φ`: kiln feed rate, [m³/s].
    - `ω`: kiln rotation rate, [rev/s].
    - `β`: kiln slope, [rad].
    - `γ`: solids repose angle, [rad].
    - `d`: particle size or dam height, [m].
    """
    function SolutionLinearKramersModel(;
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
        # Map initial condition (dam/particle size).
        h₀ = [model.h => d]

        # Map model parameters.
        p = [model.R => R,
             model.Φ => Φ,
             model.ω => ω,
             model.β => β,
             model.γ => γ]

        # Create and solve problem.
        prob = ODEProblem(model.sys, h₀, (0.0, L), p, jac = true)
        sol = solve(prob, solver, reltol = rtol, abstol = atol)

        # Bed geometry processing.
        bed = RotaryKilnBedGeometry(sol.t, sol[1, :], R, L)

        # Residence time is bed volume divided by flow rate.
        τ = bed.V  / Φ

        return new(bed, τ)
    end
end

"""
    solvelinearkramersmodel(;
        L::Float64,
        D::Float64,
        Φ::Float64,
        ω::Float64,
        β::Float64,
        γ::Float64,
        d::Float64,
        model::Union{SymbolicLinearKramersModel,Nothing} = nothing
    )::SolutionLinearKramersModel

- `L`: kiln length, [m].
- `D`: kiln internal diameter, [m].
- `Φ`: kiln feed rate, [m³/h].
- `ω`: kiln rotation rate, [rev/min].
- `β`: kiln slope, [°].
- `γ`: solids repose angle, [°].
- `d`: particle size or dam height, [mm].

**Important:** inputs are in customary units based in international
system (SI) - no imperial units supported - while all the computed
values are provided in SI units as a better physical practice. If `d`
represents a particle size, take care that it is provided in millimeters
for ease of exchange with common dam sizes.
"""
function solvelinearkramersmodel(;
        L::Float64,
        D::Float64,
        Φ::Float64,
        ω::Float64,
        β::Float64,
        γ::Float64,
        d::Float64,
        model::Union{SymbolicLinearKramersModel,Nothing} = nothing
    )::SolutionLinearKramersModel
    if isnothing(model)
        model = SymbolicLinearKramersModel()
    end

    return SolutionLinearKramersModel(
            model = model,
            L = L,
            R = D / 2.0,
            Φ = Φ / 3600.0,
            ω = ω / 60.0,
            β = deg2rad(β),
            γ = deg2rad(γ),
            d = d / 1000.0
    )
end

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
function plotlinearkramersmodel(
        model::SolutionLinearKramersModel;
        normz::Bool = false,
        normh::Bool = false
    )::Plots.Plot
    z = model.bed.z
    h = model.bed.h
    τ = model

    z = normz ? (100z / maximum(z[end])) : z
    h = normh ? (100h / maximum(h[end])) : 100h

    unitz = normz ? "%" : "m"
    unith = normh ? "%" : "cm"

    xlims  = normz ? (0.0, 100.0) : (0.0, model.bed.z[end])
    xticks = normz ? (0.0:20.0:100.0) : nothing

    ylims  = normz ? (0.0, 100.0) : nothing
    yticks = normh ? (0.0:20.0:100.0) : nothing

    p = plot()

    plot!(p, z, h, linewidth = 3, label = nothing)

    η = @sprintf("%.1f", model.bed.ηₘ)
    τ = @sprintf("%.0f", model.τ / 60)

    plot!(p,
          title  = "Loading $(η)% | Residence $(τ) min",
          xaxis  = "Coordinate [$(unitz)]",
          yaxis  = "Bed height [$(unith)]",
          xlims  = xlims,
          xticks = xticks,
          ylims  = ylims,
          yticks = yticks
    )

    return p
end

end # (module Kramers)