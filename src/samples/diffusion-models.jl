# -*- coding: utf-8 -*-
import Pkg
Pkg.activate(Base.current_project())
Pkg.instantiate()

using Revise
using CairoMakie
using DryTooling: RadialUserDefinedGrid1D
using DryTooling: RadialEquispacedGrid1D
using DryTooling: RadialGeometricGrid1D
using DryTooling: Cylinder1DTemperatureModel
using DryTooling: solve

function illustrategridconvergence()
    N = 100
    R = 0.05
    c = 900.0
    ρ = 3000.0
    h = (t) -> (t < 2000) ? 20000.0 : 0.0
    B = (t) -> 1500.0 + 200sin(2π * t / 240)
    κ = (T) -> 2.0 + 0.0001T

    grids = [
        RadialEquispacedGrid1D(R, N),
        RadialGeometricGrid1D(R, N),
        RadialUserDefinedGrid1D(geomspace(0, R, 0.003, 0.0005))
    ]

    fig = Figure(resolution = (720, 500))
    ax = Axis(fig[1, 1], yscale = identity)

    for grid in grids
        f = typeof(grid)
        model = Cylinder1DTemperatureModel(grid, h, B, κ, ρ, c)

        @time r = solve(model;
            T     = 300.0,
            τ     = 2.0,
            t     = 2400.0,
            iters = 20,
            relax = 0.001,
            tol   = 1.0e-10
        )

        lines!(ax, 100grid.r, model.problem.x,
            label = "$(f) | N = $(grid.N)")
    end

    ax.xlabel = "Radial coordinate [cm]"
    ax.ylabel = "Temperature [K]"
    ax.xticks = 0.0:1.0:100R
    xlims!(ax, (-0.2, 100R+0.2))
    axislegend(ax; position = :lt)
    fig
end

fig = illustrategridconvergence()

