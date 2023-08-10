# -*- coding: utf-8 -*-
module Kramers

using ModelingToolkit

using DifferentialEquations: ODEProblem
using DifferentialEquations: solve
using Trapz: trapz

export RotaryKilnBedGeometry
export SymbolicLinearKramersModel
export SolutionLinearKramersModel
export solvelinearkramersmodel

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

		return new(z, h, θ, l, A, η, ηₘ)
	end
end

struct SymbolicLinearKramersModel
	R::Num
	Φ::Num
	ω::Num
	β::Num
	γ::Num
	z::Num
	h::Num
	sys::ODESystem
	
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

struct SolutionLinearKramersModel
	bed::RotaryKilnBedGeometry
    τ::Float64
	
	function SolutionLinearKramersModel(;
		model::SymbolicLinearKramersModel,
		L::Float64,
		R::Float64,
		Φ::Float64,
		ω::Float64,
		β::Float64,
		γ::Float64,
		d::Float64,
        solver::Any = ODE.Tsit5(),
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

function solvelinearkramersmodel(;
        L::Float64,
        D::Float64,
        Φ::Float64,
        ω::Float64,
        β::Float64,
        γ::Float64,
        d::Float64,
        model::SymbolicLinearKramersModel = nothing
    )::RotaryKilnBedGeometry
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
            d = d
    )
end

end # (module Kramers)