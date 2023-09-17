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
            ĥ::Float64,
            Tₚ::Float64,
            Tₛ::Float64
        )
        Aᵥ = π * R^2
        Aₛ = 2 * π * R * space.δ
        
        a₁ = ρ * u * Aᵥ * cₚ
        a₂ = ĥ * Aₛ
        
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

N = 100
L = 10.0

Tₚ = 300.0
Tₛ = 1000.0

R = 0.015
ρ = 1.0
u = 1.0
ĥ = 10.0
cₚ = 1000.0

space = FvmLinearSpace1D(L, N)
solution = FvmThermalModelAllConstant(space, R, ρ, u, cₚ, ĥ, Tₚ, Tₛ)
# plot(space.zc, solution.T)

# - Guess solution in temperature.
# - Compute RHS and solve for enthalpy.
# - Find temperature root of enthalpy.
# - Repeat.

h(T) = cₚ * T
find_temperature(Tₖ, hₖ) = find_zero(T->h(T)-hₖ, Tₖ)
underrelax(aold, anew, α) = (1 - α) * anew + α * aold

Aᵥ = π * R^2
Aₛ = 2 * π * R * space.δ
ṁ = ρ * u * Aᵥ

α = 0.8
M = spdiagm(0 => ones(N), -1 => -ones(N-1))
C₁ = ĥ * Aₛ / ṁ

# Random initialization of expected solution.
# Force boundary condition (fix this later with computed value!)
T_old = Tₛ * ones(N+1)
T_old[1] = Tₚ

atol = 1.0e-08
maxiter = 100

niter = 0
residuals = []
T_new = similar(T_old)

# XXX: otherwise it always have a huge residual.
T_new[1] = Tₚ

while niter < maxiter
    T_mid = (T_old[1:end-1] + T_old[2:end])/2
    x = C₁ * (Tₛ .- T_mid)

    # Apply boundary condition (h1 = (C*(Ts - T*) + 2h0)/2).
    x[1] = (x[1] + 2 * h(Tₚ)) / 2

    T_new[2:end] = broadcast(find_temperature, T_old[2:end], M \ x)
    T_old[2:end] = underrelax(T_old[2:end], T_new[2:end], α)

    ε = maximum(abs2.(T_old - T_new))
    push!(residuals, ε)

    if ε <= atol
        break
    end

    niter += 1
end

p = plot()
plot!(space.zc, solution.T)
plot!(space.zc, T_old[2:end])
