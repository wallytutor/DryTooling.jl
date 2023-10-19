# -*- coding: utf-8 -*-
import Pkg
Pkg.activate(Base.current_project())
Pkg.instantiate()

using Revise
using CairoMakie
using ExtendableGrids: geomspace
using DryTooling
using DryTooling.Residuals
using DryTooling.Grids
using DryTooling.Utilities
using DryTooling.HeatConduction
using DryTooling.DiffusionInSolids
using DryTooling.DiffusionInSolids: carburizemoletomassfraction
import DryTooling as dry

function workingsofsimulationresiduals()
    """
        workingsofsimulationresiduals()

    Function to illustrate how to use `SimulationResiduals` in a loop.
    This might be a tutorial for new solver conception, for instance.
    """
    r = SimulationResiduals(2, 5, 10)

    for kouter in 1:7
        for kinner in 1:rand(1:5)
            r.innersteps[kouter] = kinner
            addresidual!(r, rand(r.N))
        end
    end

    s = SimulationResiduals(r)
    return plotsimulationresiduals(s; showinner = true)
end

function mwediffusionmodel(mtype)
    """
        mwediffusionmodel()

    Minimum example for solving heat equation with constant properties.
    """
    model = mtype(;
        grid = equidistantcellsgrid1D(0.05, 20),
        h    = 20.0,
        B    = 1500.0,
        κ    = 2.0,
        ρ    = 3000.0,
        c    = 900.0
    )
    dry.solve(model; t = 200.0, τ = 2.0, T = 300.0, M = 50)

    s = SimulationResiduals(model.res[])
    fig, ax, p = plotsimulationresiduals(s, showinner = false)

    return model, fig, ax, p
end

function compareradialdiffusion()
    """
        compareradialdiffusion()
    
    Vomparison of a radial and spherical coodinates solution.
    """
    retsc = mwediffusionmodel(Cylinder1DTemperatureModel)
    retss = mwediffusionmodel(Sphere1DTemperatureModel)
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
        return Cylinder1DTemperatureModel(; grid, h, B, κ, ρ, c)
    end

    N = 100
    R = 0.05

    grids = [
        equidistantcellsgrid1D(R, N),
        geometriccellsgrid1D(R, N),
        UserDefinedGrid1D(geomspace(0, R, 0.003, 0.0005))
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
    grid = equidistantcellsgrid1D(L, N)
    model = carburize(grid, t, τ, T, h, y0, ys, ; M = 50)
    
    z = model.grid.r
    y = carburizemoletomassfraction.(model.problem.x)
    m = interstitialmassintake(model.grid.r, y0, y)

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