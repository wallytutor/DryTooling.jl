# -*- coding: utf-8 -*-
using CairoMakie
using Distributions
using LinearAlgebra
using Polynomials
using Random
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
abstract type MatrixProblem end
abstract type AbstractModel end

##############################################################################
# UTILITIES
##############################################################################

"Compute the power of `x` closest to `v`."
function closestpowerofx(v::Number; x::Number = 10)::Number
    rounder = x^floor(log(x, v))
    return convert(Int64, rounder * ceil(v/rounder))
end

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
# REACTOR DATA
##############################################################################

"Geometric description of vertical reactor bounding volume."
struct VerticalReactorGeometry
    "Reactor total height [m]"
    H::Float64

    "Reactor cross-section depth [m]"
    D::Float64

    "Reactor cross-section width [m]"
    W::Float64

    "Reactor cross-section perimeter [m]"
    P::Float64

    "Reactor cross-section area [m²]"
    A::Float64

    "Reactor total volume [m³]"
    V::Float64

    function VerticalReactorGeometry(; H, D, W)
        P = 2 * (D + W)
        A = D * W
        V = A * H
        return new(H, D, W, P, A, V)
    end
end

##############################################################################
# POROSITY EVALUATION
##############################################################################

"Provides description of porosity parameters with stochastic behavior."
struct PorosityDescriptor
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

    function PorosityDescriptor(;
            μϕ::Float64,
            μl::Float64,
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
            ϕᵤ = rand(truncated(Normal(μϕ, σϕ), ϕlims...), N)
            lᵤ = rand(truncated(Normal(μl, σl), llims...), N)
        else
            ϕᵤ = μϕ
            lᵤ = μl
        end

        P = @. 6.0 * (1.0 - ϕᵤ) / lᵤ
        D = @. 4.0 * ϕᵤ / P

        return new(ϕᵤ, lᵤ, σϕ, σl, area * P, D, area)
    end
end

function Base.length(p::PorosityDescriptor)
    return length(p.ϕ)
end

##############################################################################
# RESIDUALS CONTROL
##############################################################################

"Manage residuals storage during a simulation."
mutable struct ResidualsRaw
    inner::Int64
    outer::Int64
    counter::Int64
    innersteps::Vector{Int64}
    residuals::Vector{Float64}
    function ResidualsRaw(inner::Int64, outer::Int64)
        innersteps = -ones(Int64, outer)
        residuals = -ones(Float64, outer * inner)
        return new(inner, outer, 0, innersteps, residuals)
    end
end

"Post-processed simulation residuals."
struct ResidualsProcessed
    counter::Int64
    innersteps::Vector{Int64}
    residuals::Vector{Float64}
    finalsteps::Vector{Int64}
    finalresiduals::Vector{Float64}

    function ResidualsProcessed(r::ResidualsRaw)
        innersteps = r.innersteps[r.innersteps .> 0.0]
        residuals = r.residuals[r.residuals .> 0.0]

        finalsteps = cumsum(innersteps)
        finalresiduals = residuals[finalsteps]

        return new(r.counter, innersteps, residuals,
                   finalsteps, finalresiduals)
    end
end

"Feed residuals in an inner relaxation loop."
function feedinnerresidual(r::ResidualsRaw, ε::Float64)
    # TODO: add resizing test here!
    r.counter += 1
    r.residuals[r.counter] = ε
end

##############################################################################
# RESIDUALS METRICS
##############################################################################

"Maximum relative change in a solution array."
function maxrelativevariation(
        x::Vector{Float64},
        Δx::Vector{Float64}
    )::Float64
    return maximum(abs.(Δx ./ x))
end

##############################################################################
# GENERIC MEMORY ALLOCATION
##############################################################################

"Memory for a tridiagonal problem of size `N`."
struct TridiagonalProblem <: MatrixProblem
    A::Tridiagonal{Float64, Vector{Float64}}
    b::Vector{Float64}
    x::Vector{Float64}

    function TridiagonalProblem(N)
        A = Tridiagonal(zeros(N), zeros(N+1), zeros(N))
        return new(A, zeros(N+1), zeros(N+1))
    end
end

function Base.length(p::MatrixProblem)
    return length(p.x)
end

"Integrates and arbitrary problem stepping over non-linear relaxations."
function relaxationouterloop(;
        model::AbstractModel,
        updaterouter::Function,
        updaterinner::Function,
        tend::Float64,
        tau::Float64,
        iters::Int64 = 10,
        relax::Float64 = 0.5,
        tol::Float64 = 1.0e-08,
        metric::Function = maxrelativevariation
    )::ResidualsProcessed
    times = 0.0:tau:tend
    residual = ResidualsRaw(iters, length(times))

    # TODO use `ts` for time dependent *things*!
    for (nouter, ts) in enumerate(times)
        updaterouter(model, nouter, ts)

        residual.innersteps[nouter] = relaxationinnerloop(;
            model    = model,
            updater  = updaterinner,
            residual = residual,
            iters    = iters,
            relax    = relax,
            tol      = tol,
            metric   = metric,
        )

        # store solution here
    end

    return ResidualsProcessed(residual)
end

"Solve an arbitrary problem through successive relaxations."
function relaxationinnerloop(;
        model::AbstractModel,
        updater::Function,
        residual::ResidualsRaw,
        iters::Int64 = 10,
        relax::Float64 = 0.5,
        tol::Float64 = 1.0e-08,
        metric::Function = maxrelativevariation
    )::Int64
    for niter in 1:iters
        updater(model)
        ε = relaxationstep(model.problem, relax, metric)
        feedinnerresidual(residual, ε)
        if ε <= tol
            return niter
        end
    end

    @warn "Did not converge after $(iters) steps"
    return iters
end

"Applies relaxation to solution and compute residual."
function relaxationstep(
        p::MatrixProblem,
        relax::Float64,
        metric::Function
    )::Float64
    Δx = (1.0 - relax) * (p.A \ p.b - p.x)
    p.x[:] += Δx
    return metric(p.x, Δx)
end

##############################################################################
# DIFFUSION MODELS
##############################################################################

"Thermal diffusion in a sphere represented in temperature space."
struct SphereTemperatureModel <: AbstractModel
    problem::TridiagonalProblem
    z::Vector{Float64}
    α::Vector{Float64}
    β::Vector{Float64}
    U::Float64
    T∞::Float64
    k::Function
    t::Float64
    τ::Float64
    Q::Vector{Float64}

    function SphereTemperatureModel(;
            N::Int64,
            R::Float64,
            ρ::Float64,
            c::Float64,
            h::Float64,
            T∞::Float64,
            T₀::Float64,
            k::Function,
            t::Union{Float64, Nothing} = nothing,
            M::Union{Int64, Nothing} = nothing
        )
        # Compute final integration time if required.
        tend = isnothing(t) ? 2*ρ*c*R^2/k(0.5(T∞+T₀)) : t
        steps = isnothing(M) ? convert(Int64, round(tend/10)) : M

        # Space discretization.
        δr = R / N
        z = collect(0.0:δr:R)
        w = collect(0.5δr:δr:R-0.5δr)
        r = vcat(0.0, w, R)

        # Increments.
        δ = z[2:end-0] - z[1:end-1]
        τ = tend / steps

        # To handle boundaries use ``r`` for computing α.
        α = @. (r[2:end-0]^3 - r[1:end-1]^3)*(ρ*c)/(3τ)

        # For β use only internal walls ``w``.
        β = @. w^2 / δ

        # Heat transfer coefficient multiplied by area.
        U = 4π * R^2 * h

        # Create linear problem memory.
        problem = TridiagonalProblem(N)
        problem.x[:] .= T₀

        Q = zeros(steps+2)
        Q[1] = U * (T∞ - T₀)

        return new(problem, z, α, β, U, T∞, k, tend, τ, Q)
    end
end

##############################################################################
# PLOTTING
##############################################################################

"Plot problem residuals over iterations or steps."
function plotresiduals(
        r::ResidualsProcessed;
        ε::Union{Float64, Nothing} = nothing,
        showinner::Bool = false,
        yticks::Any = nothing,
        xbase::Number = 20
    )::Figure
    function getxticks(xv)
        δi = closestpowerofx(xv/10; x = xbase)
        imax = closestpowerofx(xv; x = xbase)
        return 0:δi:imax
    end

    xs = r.finalsteps
    ys = log10.(r.finalresiduals)

    fig = Figure(resolution = (720, 500))
    ax = Axis(fig[1, 1], yscale = identity)
    ax.ylabel = "log10(Residual)"
    ax.title = "Maximum of $(maximum(r.innersteps)) iterations per step"

    if !isnothing(yticks)
        ax.yticks = yticks
        ylims!(ax, extrema(yticks))
    end

    if showinner
        xg = 1:r.counter
        yg = log10.(r.residuals)
        xticks = getxticks(xg[end])
        ax.xlabel = "Global iteration"
        lines!(ax, xg, yg, color = :black, linewidth = 0.5)
        scatter!(ax, xs, ys, color = :red)
    else
        xticks = getxticks(length(xs))
        ax.xlabel = "Outer iteration"
        lines!(ax, xs, ys, color = :red)
    end

    ax.xticks = xticks
    xlims!(ax, extrema(xticks))

    if !isnothing(ε)
        hlines!(ax, log10(ε), color = :blue, linewidth = 2)
    end

    return fig
end

##############################################################################
# DEVEL
##############################################################################

# data = YAML.load_file("mixtures.yaml")
# mix = GasMixture(data, order = ["fumes", "co2"])

# W = molecularmasses(mix)

# T = ZEROCELSIUS
# P = ONEATM
# Y = [1.0, 0.0]
# ρ = idealgasdensity(T, P, Y; W = W)

# K = 20
# T = collect(range(300.0, 3000.0, K))
# P = collect(range(300.0, 3000.0, K)) .+ ONEATM
# Y = [range(0.0, 1.0, K) range(1.0, 0.0, K)]
# ρ = idealgasdensity.(T, P, eachrow(Y); W = W)

# K = 20
# T = (ZEROCELSIUS + 25.0) * ones(K)
# P = ONEATM * ones(K)
# Y = [range(0.0, 1.0, K) range(1.0, 0.0, K)]
# ρ, μ, k, c = zip(mixtureproperties(T, P, Y; m = mix, W = W)...)

# K = 20
# T = collect(range(300.0, 2000.0, K))
# P = ONEATM * ones(K)
# Y = [ones(Float64, K) zeros(Float64, K)]
# ρ, μ, k, c = zip(mixtureproperties(T, P, Y; m = mix, W = W)...)
