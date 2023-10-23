# -*- coding: utf-8 -*-
module Granular

using CairoMakie
using DifferentialEquations: ODEProblem, Tsit5
using DifferentialEquations: solve
using Distributions
using DocStringExtensions: TYPEDFIELDS
using ModelingToolkit
using Printf
using Random
using Trapz: trapz

export PackedBedPorosityDescriptor
export SymbolicLinearKramersModel
export RotaryKilnBedSolution
export solvelinearkramersmodel
export plotlinearkramersmodel

"""
Provides description of porosity parameters with stochastic behavior.

$(TYPEDFIELDS)

# Usage

`PackedBedPorosityDescriptor` can be used to describe the geometry
of exchange section of a packed bed for a single set of arguments.

```jldoctest
julia> using DryTooling.Granular

julia> PackedBedPorosityDescriptor(; ϕ = 0.65, l = 0.10, area = 1.0)
PackedBedPorosityDescriptor(P = 21.000000 m, D = 0.123810 m)
```

It can also be used to describe randomly varying reactors, what is
a more realistic thing to do when using this structure to simulate
real world systems.

```jldoctest
julia> using DryTooling.Granular

julia> PackedBedPorosityDescriptor(;
            ϕ  = 0.65, l  = 0.10,
            σϕ = 0.03, σl = 0.01,
            N = 2,
            ϕlims = (0.4, 0.8),
            llims = (0.0, 0.3),
            seed = 42,
            area = 1.0
        )
PackedBedPorosityDescriptor(
    P from  21.455749 m to  24.370742 m
    D from   0.125589 m to   0.102353 m
)
```
"""
struct PackedBedPorosityDescriptor
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

"""
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

"""
General geometric description of a bed from Kramers equation solution.

$(TYPEDFIELDS)

# Arguments

Internal elements are initialized through the following constructor:

```julia
RotaryKilnBedSolution(z::Vector{Float64}, h::Vector{Float64}, R::Float64, Φ::Float64)
```

Where parameters are given as:

- `z`: solution coordinates over length, [m].
- `h`: bed profile solution over length, [m].
- `R`: kiln internal radius, [m].
- `Φ`: kiln feed rate, [m³/s].

An outer constructor is also provided for managing the integration of an
instance of `SymbolicLinearKramersModel`. This is the recommended usage
that is illustrated below.

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
- `solver`: Solver for `DifferentialEquations`. Defaults to `Tsit5`.
- `rtol`: Relative integration tolerance. Defaults to 1.0e-08.
- `atol`: Absolute integration tolerance. Defaults to 1.0e-08.

# Usage

Data in next example is an SI conversion of an example from Kramers (1952).

```jldoctest
julia> using DryTooling.Granular

julia> L = 13.715999999999998;  # Kiln length [m]

julia> D = 1.8897599999999999;  # Kiln diameter [m]

julia> β = 2.3859440303888126;  # Kiln slope [°]

julia> γ = 45.0;                # Repose angle [°]

julia> d = 1.0;                 # Particle/dam size [mm]

julia> Φ = 10.363965852671996;  # Feed rate [m³/h]

julia> ω = 3.0300000000000002;  # Rotation rate [rev/min]

julia> bed = RotaryKilnBedSolution(;
            model = SymbolicLinearKramersModel(),
            L     = L,
            R     = D / 2.0,
            Φ     = Φ / 3600.0,
            ω     = ω / 60.0,
            β     = deg2rad(β),
            γ     = deg2rad(γ),
            d     = d / 1000.0
        );

julia> bed
RotaryKilnBedSolution(τ = 13.169938 min, ηₘ = 5.913271 %)

julia> bed.τ
790.1963002204092
```
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

function RotaryKilnBedSolution(;
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
    return RotaryKilnBedSolution(sol.t, sol[1, :], β, R, Φ)
end

function Base.length(p::PackedBedPorosityDescriptor)
    return length(p.ϕ)
end

function Base.show(io::IO, obj::PackedBedPorosityDescriptor)
    if any(isnothing, [obj.σϕ, obj.σl])
        P = @sprintf("%.6f m", obj.P)
        D = @sprintf("%.6f m", obj.D)
        print(io, "PackedBedPorosityDescriptor(P = $(P), D = $(D))")
    else
        P = map(x->@sprintf("%10.6f m", x), obj.P)
        D = map(x->@sprintf("%10.6f m", x), obj.D)
        print(io, """\
        PackedBedPorosityDescriptor(
            P from $(P[1]) to $(P[2])
            D from $(D[1]) to $(D[2])
        )\
        """)
    end
end

function Base.show(io::IO, obj::RotaryKilnBedSolution)
    τ = @sprintf("%.6f min", obj.τ/60)
    ηₘ = @sprintf("%.6f", obj.ηₘ)
    print(io, "RotaryKilnBedSolution(τ = $(τ), ηₘ = $(ηₘ) %)")
end

"""
    plotlinearkramersmodel(
        model::RotaryKilnBedSolution;
        normz::Bool = false,
        normh::Bool = false
    )::Figure

Standardized plotting of `RotaryKilnBedSolution` bed profile. It
supports normalization of axes throught keywords `normz` for axial
coordinate and `normh` for bed depth.
"""
function plotlinearkramersmodel(
        model::RotaryKilnBedSolution;
        normz::Bool = false,
        normh::Bool = false
    )::Figure
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