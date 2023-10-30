# -*- coding: utf-8 -*-
import Pkg
import Revise

if Base.current_project() != Base.active_project()
    Pkg.activate(Base.current_project())
    Pkg.resolve()
    Pkg.instantiate()

    using CairoMakie
    using CommonSolve
    using CommonSolve: solve
    using DocStringExtensions: TYPEDFIELDS
    using Trapz: trapz
    using DryTooling
    using DryTooling.Simulation
    using DryTooling: Temperature1DModelStorage
    using DryTooling: interfaceconductivity1D
end

struct Sphere1DEnthalpyModel <: AbstractDiffusionModel1D
    "Energy equation in an 1-D sphere formulated in enthalpy."

    "Grid over which problem will be solved."
    grid::AbstractGrid1D

    "Memory for model linear algebra problem."
    problem::TridiagonalProblem

    "Model coefficient α."
    α::Vector{Float64}

    "Constant part of model coefficient β."
    β′::Vector{Float64}

    "Thermal conductivity in terms of temperature."
    κ::Function

    "Global heat transfer coefficient ``U=hR``."
    U::Function

    "Surface environment temperature."
    B::Function

    "Enthalpy function of temperature."
    H::Function

    "Time-step used in integration."
    τ::Base.RefValue{Float64}

    "Memory storage for solution retrieval."
    mem::Base.RefValue{Temperature1DModelStorage}

    "Residuals tracking during solution."
    res::Base.RefValue{TimeSteppingSimulationResiduals}

    "Surface area scaling factor."
    scale::Float64

    function Sphere1DEnthalpyModel(;
            grid::AbstractGrid1D,
            h::Union{Function,Float64},
            B::Union{Function,Float64},
            κ::Union{Function,Float64},
            H::Function,
            ρ::Float64
        )
        hu = (typeof(h) <: Function) ? h : (t) -> h
        Bu = (typeof(B) <: Function) ? B : (t) -> B
        κu = (typeof(κ) <: Function) ? κ : (T) -> κ

        problem = TridiagonalProblem(grid.N)

        rₙ = tail(grid.w)
        rₛ = head(grid.w)
        α = @. ρ * (rₙ^3 - rₛ^3) / 3.0

        rₙ = tail(grid.r)
        rₛ = head(grid.r)
        wⱼ = body(grid.w)
        β′ = @. wⱼ^2 / (rₙ - rₛ)

        R = last(grid.r)^2
        U = (t) -> hu(t) * R
        τ = Ref(-Inf)

        mem = Ref(Temperature1DModelStorage(0, 0))
        res = Ref(TimeSteppingSimulationResiduals(1, 0, 0))

        return new(grid, problem, α, β′, κu, U, Bu, H, τ, mem, res, 4π)
    end
end

function initialize!(
        m::Sphere1DEnthalpyModel,
        t::Float64,
        τ::Float64;
        T::Union{Float64,Nothing} = nothing,
        M::Int64 = 50
    )::Nothing
    "Set initial condition of thermal diffusion model."
    # XXX: this is the same as for LocalAbstractTemperature1DModel!
    # Try to instantiate that function instead of the next calls!
    # ---
    nsteps = convert(Int64, round(t / τ))
    m.τ[] = Base.step(range(0.0, t, nsteps))
    m.res[] = TimeSteppingSimulationResiduals(1, M, nsteps)
    m.mem[] = Temperature1DModelStorage(m.grid.N, nsteps)
    
    # Do not reinitialize a problem.
    if !isnothing(T)
        m.problem.x[:] .= T
    end
    # ---

    # Update constant coefficient with time-step.
    m.α[:] = @. m.α / m.τ[] 

    return nothing
end

function DryTooling.Simulation.fouter!(
        m::Sphere1DEnthalpyModel, t::Float64, n::Int64
    )::Nothing
    "Time-step dependent updater for model."
    # XXX: for now evaluating B.C. at mid-step, fix when going full
    # semi-implicit generalization!
    U = m.U(t + m.τ[]/2)
    B = m.B(t + m.τ[]/2)

    # Follow surface heat flux and store partial solutions.
    # XXX: note the factor 4π because U = r²h only and A = 4πr²!!!!
    # XXX: note the factor 2π because U = rh only and A = 2πrl!!!!
    m.mem[].t[n] = t
    m.mem[].Q[n] = m.scale * U * (B - last(m.problem.x))
    m.mem[].T[n, 1:end] = m.problem.x

    @. m.problem.b[1:end] = map(m.H, m.problem.x)
    m.problem.b[end] += U * B / m.α[end]
    return nothing
end

function DryTooling.Simulation.finner!(
        m::Sphere1DEnthalpyModel, t::Float64, n::Int64
    )::Nothing
    "Non-linear iteration updater for model."
    κ = interfaceconductivity1D(m.κ.(m.problem.x))
    β = κ .* m.β′

    # XXX: for now evaluating B.C. at mid-step, fix when going full
    # semi-implicit generalization!
    U = m.U(t + m.τ[]/2)

    m.problem.A.dl[1:end] = -β
    m.problem.A.du[1:end] = -β
    m.problem.A.d[1:end]  = m.α

    m.problem.A.d[2:end-1] += tail(β) + head(β)
    m.problem.A.d[1]       += first(β)
    m.problem.A.d[end]     += last(β) + U / m.α[end]
    return nothing
end

function DryTooling.Simulation.fsolve!(
        m::Sphere1DEnthalpyModel, t::Float64, n::Int64, α::Float64
    )::Float64
    "Solve problem for one non-linear step."
    # ε = relaxationstep!(m.problem, α, maxabsolutechange)
    # addresidual!(m.res[], [ε])
    # return ε
end
