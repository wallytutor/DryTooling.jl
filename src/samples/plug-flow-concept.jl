# -*- coding: utf-8 -*-
"""
Thermal plug-flow model concept
===============================

- [x] Finite volume matrix formulation.
- [ ] Box-model with ModelingToolkit.
- [ ] Solution in enthalpy space.
"""
using Plots
using Roots
using SparseArrays: spdiagm

struct FvmLinearSpace1D
    """ Cell length [m]. """
    δ::Float64

    """ Centers coordinates [m]. """
    zc::Vector{Float64}

    """ Boundaries coordinates [m]. """
    zw::Vector{Float64}

    function FvmLinearSpace1D(L::Float64, N::Int64)
        δ = L / N
        zc₀ = δ / 2
        zc₁ = L - δ / 2
        zc = collect(zc₀:δ:zc₁)
        zw = collect(0.0:δ:L)
        return new(δ, zc, zw)
    end
end

struct FvmThermalModelAllConstant
    T::Vector{Float64}

    function FvmThermalModelAllConstant(
            space::FvmLinearSpace1D,
            R::Float64,
            ρ::Float64,
            u::Float64,
            cₚ::Float64,
            h::Float64,
            Tₚ::Float64,
            Tₛ::Float64
        )
        Aᵥ = π * R^2
        Aₛ = 2 * π * R * space.δ
        
        a₁ = ρ * u * Aᵥ * cₚ
        a₂ = h * Aₛ
        
        A⁺ = a₁ + a₂ / 2
        A⁻ = a₁ - a₂ / 2
        
        C₁ = a₂ * Tₛ
        C₂ = C₁ + 2 * A⁻ * Tₚ
        
        d₀ = A⁺ * ones(N)
        d₀[1] = 2 * a₁
        
        d₁ = -A⁻ * ones(N - 1)
        
        C = C₁ * ones(N)
        C[1] = C₂
        
        M = spdiagm(0 => d₀, -1 => d₁)

        return new(M \ C)
    end
end

N = 1000
L = 10.0

Tₚ = 300.0
Tₛ = 1000.0

R = 0.015
ρ = 1.0
u = 1.0
cₚ = 1000.0
h = 10.0

# space = FvmLinearSpace1D(L, N)
# solution = FvmThermalModelAllConstant(space, R, ρ, u, cₚ, h, Tₚ, Tₛ)
# plot(space.zc, solution.T)


# - Guess solution in temperature.
# - Compute RHS and solve for enthalpy.
# - Find temperature root of enthalpy.
# - Repeat.

ṁ = 1.0
ĥ = 10.0
Aₛ = 2 * π * R * space.δ

enthalpy(T) = 1000 * T + 1000
h₀ = enthalpy(Tₚ)

space = FvmLinearSpace1D(L, N)
T = Tₛ * ones(N)

hₛ = enthalpy.(T)

find_zero((T->enthalpy(T)-hₛ[1]), (0.9T[1]))