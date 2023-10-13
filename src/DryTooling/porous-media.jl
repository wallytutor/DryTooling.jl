# -*- coding: utf-8 -*-

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
