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

function mwediffusionmodel()
    """
        mwediffusionmodel()

    Minimum example for solving heat equation with constant properties.
    """
    model = dry.Cylinder1DTemperatureModel(;
        grid = dry.equidistantcellsgrid(0.05, 20),
        h    = 20.0,
        B    = 1500.0,
        κ    = 2.0,
        ρ    = 3000.0,
        c    = 900.0
    )
    dry.solve(model; t = 2400.0, τ = 2.0, T = 300.0, M = 50)

    s = dry.SimulationResiduals(model.res[])
    fig, ax, p = dry.plotsimulationresiduals(s, showinner = false)

    return model, fig, ax, p
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

rets = mwediffusionmodel()

rets = illustrategridconvergence()
