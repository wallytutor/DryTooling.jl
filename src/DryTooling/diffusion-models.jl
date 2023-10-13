# -*- coding: utf-8 -*-
export SphereTemperatureModel
export solve

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

"Interface for solving a `SphereTemperatureModel` instance."
function CommonSolve.solve(model::SphereTemperatureModel; kwargs...)
	return relaxationouterloop(; model = model, kwargs...)
end
