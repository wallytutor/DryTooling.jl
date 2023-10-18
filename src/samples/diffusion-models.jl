# -*- coding: utf-8 -*-
import Pkg
Pkg.activate(Base.current_project())
Pkg.instantiate()

using Revise
using CairoMakie
using ExtendableGrids: geomspace
import DryTooling as dry

function workingsofsimulationresiduals()
    """
        workingsofsimulationresiduals()

    Function to illustrate how to use `SimulationResiduals` in a loop.
    This might be a tutorial for new solver conception, for instance.
    """
    r = dry.SimulationResiduals(2, 5, 10)

    for kouter in 1:7
        for kinner in 1:rand(1:5)
            r.innersteps[kouter] = kinner
            dry.step(r, rand(r.N))
        end
    end

    s = dry.SimulationResiduals(r)
    return dry.plotsimulationresiduals(s; showinner = true)
end

function mwediffusionmodel(mtype)
    """
        mwediffusionmodel()

    Minimum example for solving heat equation with constant properties.
    """
    model = mtype(;
        grid = dry.equidistantcellsgrid(0.05, 20),
        h    = 20.0,
        B    = 1500.0,
        κ    = 2.0,
        ρ    = 3000.0,
        c    = 900.0
    )
    dry.solve(model; t = 200.0, τ = 2.0, T = 300.0, M = 50)

    s = dry.SimulationResiduals(model.res[])
    fig, ax, p = dry.plotsimulationresiduals(s, showinner = false)

    return model, fig, ax, p
end

function compareradialdiffusion()
    """
        compareradialdiffusion()
    
    Vomparison of a radial and spherical coodinates solution.
    """
    retsc = mwediffusionmodel(dry.Cylinder1DTemperatureModel)
    retss = mwediffusionmodel(dry.Sphere1DTemperatureModel)
    c = retsc[1]
    s = retss[1]
    
    fig = Figure(resolution = (720, 500))
    ax = Axis(fig[1, 1], yscale = identity)
    lines!(ax, 100c.grid.r, c.problem.x, label = "Cylinder")
    lines!(ax, 100s.grid.r, s.problem.x, label = "Sphere")
    ax.xlabel = "Radial coordinate [cm]"
    ax.ylabel = "Temperature [K]"
    ax.xticks = 0.0:1.0:100last(c.grid.r)
    axislegend(ax; position = :lt)

    return fig
end

function illustrategridconvergence()
    """
        illustrategridconvergence()

    This function aims at serving as a tutorial for initialization and
    use of standard diffusion models. Here a set of grids is used to
    solve the same problem, illustrating the advantage of selecting a
    grid with some *enconded physics*.
    """
    function creatediffusionsample(grid)
        c = 900.0
        ρ = 3000.0
        h = (t) -> (t < 2000) ? 20000.0 : 0.0
        B = (t) -> 1500.0 + 200sin(2π * t / 240)
        κ = (T) -> 2.0 + 0.0001T
        return dry.Cylinder1DTemperatureModel(; grid, h, B, κ, ρ, c)
    end

    N = 100
    R = 0.05

    grids = [
        dry.equidistantcellsgrid(R, N),
        dry.geometriccellsgrid(R, N),
        dry.UserDefinedGrid1D(geomspace(0, R, 0.003, 0.0005))
    ]

    fig = Figure(resolution = (720, 500))
    ax = Axis(fig[1, 1], yscale = identity)

    for grid in grids
        model = creatediffusionsample(grid)
        @time dry.solve(model; t = 2400.0, τ = 2.0, T = 300.0, M = 50)
        lines!(ax, 100grid.r, model.problem.x, label = "N = $(grid.N)")
    end

    ax.xlabel = "Radial coordinate [cm]"
    ax.ylabel = "Temperature [K]"
    ax.xticks = 0.0:1.0:100R
    xlims!(ax, (-0.2, 100R+0.2))
    axislegend(ax; position = :lt)
    return fig, ax
end

rets = workingsofsimulationresiduals()

fig = compareradialdiffusion()

rets = illustrategridconvergence()


using CommonSolve
using DryTooling
using DryTooling: GAS_CONSTANT
using DryTooling: AbstractDiffusionModel1D, AbstractGrid1D
using DryTooling: TridiagonalProblem
using DryTooling: Temperature1DModelStorage, SimulationResiduals
using DryTooling: head, tail, step, maxabsolutechange
using DryTooling: relaxationstep!, interfaceconductivity1D
using DryTooling: fouter!, finner!, fsolve!, advance!
using Trapz: trapz

struct AusteniteCarburizing1DModel <: AbstractDiffusionModel1D
    "Carbon diffusion in a plate represented in temperature space."

    "Grid over which problem will be solved."
    grid::AbstractGrid1D
    
    "Memory for model linear algebra problem."
    problem::TridiagonalProblem
    
    "Constant part of model coefficient α."
    α′::Vector{Float64}
    
    "Constant part of model coefficient β."
    β′::Vector{Float64}
    
    "Diffusivity terms of composition."
    D::Function
    
    "Global mass transfer coefficient."
    h::Function
    
    "Environment equivalent concentration."
    C::Function
    
    "Time-step used in integration."
    τ::Base.RefValue{Float64}
    
    "Memory storage for solution retrieval."
    mem::Base.RefValue{Temperature1DModelStorage}

    "Residuals tracking during solution."
    res::Base.RefValue{SimulationResiduals}

    function AusteniteCarburizing1DModel(;
            grid::AbstractGrid1D,
            h::Union{Function,Float64},
            C::Union{Function,Float64},
            T::Float64
        )
        A(xc) = 4.84e-05exp(-38.0xc) / (1.0 - 5.0xc)
        E(xc) = 155_000.0 - 570_000.0xc
        D(xc) = A(xc) * exp(-E(xc) / (GAS_CONSTANT * T))

        hu = (typeof(h) <: Function) ? h : (t) -> h
        Cu = (typeof(C) <: Function) ? C : (t) -> C

        problem = TridiagonalProblem(grid.N)
        
        α′ = tail(grid.w) - head(grid.w)
        β′ = 1 ./ (tail(grid.r) - head(grid.r))

        τ = Ref(-Inf)
        mem = Ref(Temperature1DModelStorage(0, 0))
        res = Ref(SimulationResiduals(1, 0, 0))

        return new(grid, problem, α′, β′, D, hu, Cu, τ, mem, res)
    end
end

function initialize!(
        m::AusteniteCarburizing1DModel,
        t::Float64,
        τ::Float64,
        x::Float64;
        M::Int64 = 50
    )::Nothing
    "Set initial condition of thermal diffusion model."
    nsteps = convert(Int64, round(t / τ))
    m.τ[] = Base.step(range(0.0, t, nsteps))
    m.res[] = SimulationResiduals(1, M, nsteps)
    m.mem[] = Temperature1DModelStorage(m.grid.N, nsteps)
    m.problem.x[:] .= x
    return nothing
end

function DryTooling.fouter!(
        m::AusteniteCarburizing1DModel, t::Float64, n::Int64
    )::Nothing
    "Time-step dependent updater for model."
    # XXX: for now evaluating B.C. at mid-step, fix when going full
    # semi-implicit generalization!
    h = m.h(t + m.τ[]/2)
    C = m.C(t + m.τ[]/2)

    # Follow surface mass flux and store partial solutions.
    m.mem[].Q[n] = h * (C - last(m.problem.x))
    m.mem[].T[n, 1:end] = m.problem.x

    @. m.problem.b[1:end] = (m.α′ / m.τ[]) * m.problem.x
    m.problem.b[end] += h * C
    return nothing
end

function DryTooling.finner!(
        m::AusteniteCarburizing1DModel, t::Float64, n::Int64
    )::Nothing
    "Non-linear iteration updater for model."
    D = interfaceconductivity1D(m.D.(m.problem.x))
    β = D .* m.β′
    α = m.α′./ m.τ[]

    # XXX: for now evaluating B.C. at mid-step, fix when going full
    # semi-implicit generalization!
    h = m.h(t + m.τ[]/2)

    m.problem.A.dl[1:end] = -β
    m.problem.A.du[1:end] = -β
    m.problem.A.d[1:end]  = α

    m.problem.A.d[2:end-1] += tail(β) + head(β)
    m.problem.A.d[1]       += first(β)
    m.problem.A.d[end]     += last(β) + h
    return nothing
end

function DryTooling.fsolve!(
        m::AusteniteCarburizing1DModel, t::Float64, n::Int64, α::Float64
    )::Float64
    "Solve problem for one non-linear step."
    ε = relaxationstep!(m.problem, α, maxabsolutechange)
    step(m.res[], [ε])
    return ε
end

function solve(
        m::AusteniteCarburizing1DModel;
        t::Float64,
        τ::Float64,
        x::Float64,
        M::Int64 = 50,
        α::Float64 = 0.1,
        ε::Float64 = 1.0e-10
    )::Nothing
    "Interface for solving a `Cylinder1DTemperatureModel` instance."
    initialize!(m, t, τ, x, M = M)
    advance!(m; α, ε, M)
    return nothing
end

function masstomolefraction(w)
    return w * (w / 0.012 + (1 - w) / 0.055)^(-1) / 0.012
end

function moletomassfraction(x)
    return 0.012 * x / (0.012*x + (1 - x) * 0.055)
end

function massintake(z, y0, yf; ρ = 7890.0)
	σ₀ = y0 * last(z) / (1.0 - y0)
	σ₁ = trapz(z, @. yf / (1.0 - yf))
	return 1000ρ  * (σ₁ - σ₀)
end

function carburize(L, N, t, τ, T, h, y0, ys, ; M = 50)
    x = masstomolefraction(y0)
    C = masstomolefraction(ys)
    grid = dry.equidistantcellsgrid(L, N)
    model = AusteniteCarburizing1DModel(; grid, h, C, T)
    @time solve(model; t, τ, x, M = M, α = 0.05)
    return model
end

N = 100
L = 0.002
T = 1173.15
τ = 10.0

# t = 6*3600.0
# y0 = 0.0016
# ys = 0.0100
# h = (t) -> (t <= 7200.0) ? 6.0e-04 : 0.0

t = 5*3600.0
y0 = 0.0023
ys = 0.0095
h = (t) -> (t <= 7200.0) ? 1.0e-03 : 0.0

let 
    model = carburize(L, N, t, τ, T, h, y0, ys, ; M = 50)
    
    z = model.grid.r
    y = moletomassfraction.(model.problem.x)
    m = massintake(model.grid.r, y0, y)

    fig = Figure(resolution = (720, 500))
    ax = Axis(fig[1, 1], yscale = identity)
    lines!(ax, 1000z, 100reverse(y))
    
    ax.title  = "Mass intake $(round(m, digits = 2)) g/m²"
    ax.xlabel = "Coordinate [cm]"
    ax.ylabel = "Mass percentage [%]"
    ax.xticks = 0.0:0.2:1.2
    ax.yticks = 0.0:0.1:0.6

    xlims!(ax, extrema(ax.xticks.val))
    ylims!(ax, extrema(ax.yticks.val))
    fig
end