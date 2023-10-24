# -*- coding: utf-8 -*-
export SymbolicLinearKramersModel
export RotaryKilnBedSolution
export solvelinearkramersmodel
export plotlinearkramersmodel
export dimlessNΦ
export dimlessNₖ
export sullivansηₘ
export perrayresidence
export kramersnlapprox

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

"Kramers (1952) dimensionless group NΦ."
function dimlessNΦ(R, β, ω, Φ, γ)
    return Φ * sin(γ) / (ω * R^3 * tan(β))
end

"Kramers (1952) dimensionless group Nₖ."
function dimlessNₖ(L, R, β, γ)
    return R * cos(γ) / (L * tan(β))
end

"Sullivans approximation to kiln filling."
function sullivansηₘ(R, β, ω, Φ, γ)
    return 3.8 * dimlessNΦ(R, β, ω, Φ, γ) * sqrt(γ) / sin(γ)
end

"Compute residence time from Peray's equation."
function perrayresidence(L, ω, D, β)
    return 0.19 * L / (ω * D * tan(β))
end

"Nonlinear formulation of Kramers model approximate solution."
function kramersnlapprox(; z, R, Φ, ω, β, γ, d)
    L = z[end]
    N = length(z)

    NΦ = dimlessNΦ(R, β, ω, Φ, γ)
    Nₖ = dimlessNₖ(L, R, β, γ)

    C₁ = R * NΦ
    C₂ = 3C₁ / (4π * 1.24)
    C₃ = C₁ / (L * NΦ * Nₖ)

    optim = JuMP.Model(Ipopt.Optimizer)
    JuMP.set_silent(optim)

    @JuMP.variable(optim, h[1:N])
    @JuMP.NLconstraint(
        optim,
        [i = 1:N],
        C₂ * log((d - C₂) / (h[i] - C₂)) - C₃ * z[i] - h[i] + d == 0,
    )
    @JuMP.NLconstraint(optim, [i = 1:N], h[i] >= 0.0)
    JuMP.optimize!(optim)

    return RotaryKilnBedSolution(z, JuMP.value.(h), β, R, Φ)
end
