# -*- coding: utf-8 -*-
module PlugFlow

using DryTooling: AbstractPhysicalModel
using DryTooling: TridiagonalProblem

##############################################################################
# STATELESS MODELS
##############################################################################

# struct GasComponent <: Substance
#     "Mean molecular mass [g/mol]"
#     W::Float64

#     "Viscosity polynomial [Pa.s]"
#     μ::Polynomial{Float64, :T}

#     "Thermal conductivity polynomial [W/(m.K)]"
#     k::Polynomial{Float64, :T}

#     "Specific heat polynomial [J/(kg.K)]"
#     c::Polynomial{Float64, :T}

#     function GasComponent(;
#             W::Float64,
#             μ::Vector{Float64},
#             k::Vector{Float64},
#             c::Vector{Float64}
#         )
#         return new(
#             W,
#             Polynomial(μ, :T),
#             Polynomial(k, :T),
#             Polynomial(c, :T)
#         )
#     end
# end

# struct GasMixture <: Mixture
#     "Number of components in system."
#     K::Int64

#     "Storage of gas component objects."
#     s::Vector{GasComponent}

#     function GasMixture(; components::Vector{GasComponent})
#         return new(length(components), components)
#     end
# end

# function GasComponent(d::Dict{Any, Any})
#     return GasComponent(; W = d["W"], μ = d["mu"], k = d["kg"], c = d["cp"])
# end

# function GasMixture(d::Dict{Any, Any}; order::Any = nothing)
#     order = isnothing(order) ? keys(d) : order
#     components = map(k->GasComponent(d[k]), order)
#     return GasMixture(; components = components)
# end

# ##############################################################################
# # METHODS
# ##############################################################################

# function molecularmasses(m::GasMixture)::Vector{Float64}
#     return map(x->x.W, m.s)
# end

# function meanmolecularmass(
#         Y::Union{Vector{Float64},SubArray},
#         W::Vector{Float64}
#     )::Float64
#     return 1.0 / sum(y / w for (y, w) in zip(Y, W))
# end

# function idealgasdensity(T::Float64, P::Float64, W::Float64)::Float64
#     return P * W / (1000GAS_CONSTANT * T)
# end

# function idealgasdensity(
#         T::Float64,
#         P::Float64,
#         Y::Union{Vector{Float64},SubArray};
#         W::Vector{Float64}
#     )::Float64
#     return idealgasdensity(T, P, meanmolecularmass(Y, W))
# end

# function thermophysicalproperties(s::GasComponent, T::Float64)::Matrix{Float64}
#     return [s.μ(T) s.k(T) s.c(T)]
# end

# function thermophysicalproperties(
#         m::GasMixture,
#         T::Float64,
#         Y::Union{Vector{Float64},SubArray}
#     )::Matrix{Float64}
#     return sum(y*thermophysicalproperties(s, T) for (s, y) in zip(m.s, Y))
# end

# function mixtureproperties(
#         T::Float64,
#         P::Float64,
#         Y::Union{Vector{Float64},SubArray};
#         m::GasMixture,
#         W::Vector{Float64}
#     )::Matrix{Float64}
#     ρ = idealgasdensity(T, P, Y; W=W)
#     μ, k, c = thermophysicalproperties(m, T, Y)
#     return [ρ μ k c]
# end

# function mixtureproperties(
#         T::Vector{Float64},
#         P::Vector{Float64},
#         Y::Matrix{Float64};
#         m::GasMixture,
#         W::Vector{Float64}
#     )::Vector{Matrix{Float64}}
#     return mixtureproperties.(T, P, eachrow(Y); m = m, W = W)
# end

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
# DIFFUSION MODELS
##############################################################################

"Thermal diffusion in a sphere represented in temperature space."
struct SphereTemperatureModel <: AbstractPhysicalModel
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
    T::Matrix{Float64}

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
        α = @. (r[2:end-0]^3 - r[1:end-1]^3) * ρ * c / (3τ)

        # For β use only internal walls ``w``.
        β = @. w^2 / δ

        # Heat transfer coefficient multiplied by area.
        U = 4π * R^2 * h

        # Create linear problem memory.
        problem = TridiagonalProblem(N)
        problem.x[:] .= T₀

        # Memory for storing surface fluxes.
        Q = zeros(steps+2)

        # Memory for intermediate solution.
        T = zeros(steps+2, N+1)

        return new(problem, z, α, β, U, T∞, k, tend, τ, Q, T)
    end
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
end # (module PlugFlow)