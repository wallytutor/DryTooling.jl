# -*- coding: utf-8 -*-
module Kramers

using ModelingToolkit
using Plots
using Printf

using DocStringExtensions: TYPEDFIELDS
using DifferentialEquations: ODEProblem
using DifferentialEquations: Tsit5
using DifferentialEquations: solve
using Trapz: trapz

export SymbolicLinearKramersModel
export RotaryKilnBedSolution
export solvelinearkramersmodel
export plotlinearkramersmodel

#############################################################################
# SymbolicLinearKramersModel
#############################################################################

"""
    SymbolicLinearKramersModel

Creates a reusable linear Kramers model for rotary kiln simulation.

$(TYPEDFIELDS)
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
end

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
    Dz = Differential(z)

    # Compose problem right-hand side.
    C = (3//4) * tan(γ) * Φ / (π * R^3 * ω)
    f = C * ((h / R) * (2 - h / R))^(-3//2)

    # *Stack* equation.
    eqs = Dz(h) ~ f - tan(β) / cos(γ)

    # Assembly system for solution.
    @named sys = ODESystem(eqs)
    sys = structural_simplify(sys)

    return SymbolicLinearKramersModel(R, Φ, ω, β, γ, z, h, sys)
end

#############################################################################
# RotaryKilnBedSolution
#############################################################################

"""
    RotaryKilnBedSolution

Description of a rotary kiln bed geometry computed from the solution
of bed height along the kiln length. The main goal of the quantities
computed here is their use with heat and mass transfer models for the
simulation of rotary kiln process.

$(TYPEDFIELDS)
"""
struct RotaryKilnBedSolution
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

    "Residence time of particles"
    τ::Float64
end

"""
    RotaryKilnBedSolution(
        z::Vector{Float64},
        h::Vector{Float64},
        R::Float64,
        Φ::Float64
    )

- `z`: solution coordinates over length, [m].
- `h`: bed profile solution over length, [m].
- `R`: kiln internal radius, [m].
- `Φ`: kiln feed rate, [m³/s].
"""
function RotaryKilnBedSolution(z, h, R, Φ)
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
    return RotaryKilnBedSolution(z, h, θ, l, A, η, ηₘ, V, τ)
end

#############################################################################
# Functions
#############################################################################

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
    h = [model.h => d]
    p = [model.R => R,
         model.Φ => Φ,
         model.ω => ω,
         model.β => β,
         model.γ => γ]

    prob = ODEProblem(model.sys, h, (0.0, L), p, jac = true)
    sol = solve(prob, solver, reltol = rtol, abstol = atol)
    bed = RotaryKilnBedSolution(sol.t, sol[1, :], R, Φ)
    return bed
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
        model::RotaryKilnBedSolution;
        normz::Bool = false,
        normh::Bool = false
    )::Plots.Plot
    z = model.z
    h = model.h

    z = normz ? (100z / maximum(z[end])) : z
    h = normh ? (100h / maximum(h[end])) : 100h

    unitz = normz ? "%" : "m"
    unith = normh ? "%" : "cm"

    xlims  = normz ? (0.0, 100.0) : (0.0, model.z[end])
    xticks = normz ? (0.0:20.0:100.0) : nothing

    ylims  = normz ? (0.0, 100.0) : nothing
    yticks = normh ? (0.0:20.0:100.0) : nothing

    η = @sprintf("%.1f", model.ηₘ)
    τ = @sprintf("%.0f", model.τ / 60)

    p = plot()
    plot!(p, z, h, linewidth = 3, label = nothing)
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