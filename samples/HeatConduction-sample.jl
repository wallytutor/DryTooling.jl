# -*- coding: utf-8 -*-
"""
DryTooling.HeatConduction sample
================================

In this script we illustrate how to use the main functionalities of module
DryTooling.HeatConduction. Sample 1 compares solutions with different 1D
solvers available in the module. Sample 2 shows the support to arbitrary
space discretizations.
"""

using GLMakie
using ExtendableGrids: geomspace
using DryTooling.Grids
using DryTooling.HeatConduction
using DryTooling.HeatConduction: solve
using DryTooling.Simulation: plotsimulationresiduals

begin # (sample 1)
    # These are the common parameters used in all simulation cases.
    # It is possible to use numerical or functional arguments.
    args_num = (
        h = 20.0,
        B = 1500.0,
        κ = 2.0,
        ρ = 3000.0,
        c = 900.0
    )

    args_fun = (
        h = (t) -> 20.0,
        B = (t) -> 1500.0 + 20 * rand(),
        κ = (T) -> 2.0 + 0.01 * T,
        ρ = 3000.0,
        c = 900.0
    )

    # Comment/uncomment line to select example arguments.
    # args = args_num
    args = args_fun

    # Common parameters used when calling the solver.
    solve_pars = (
        t = 600.0,
        τ = 2.0,
        T = 300.0,
        α = 0.01,
        ε = 1.0e-10,
        M = 20
    )

    # A simple 1D grid for both the cylinder or sphere.
    grid = equidistantcellsgrid1D(0.05, 20)

    # Create models.
    model_cyl = Cylinder1DTemperatureModel(; grid, args...)
    model_sph = Sphere1DTemperatureModel(; grid, args...)

    # Solve models, *i.e.* integrate over time.
    solve(model_cyl; solve_pars...)
    solve(model_sph; solve_pars...)

    # Visualize residuals if using under-relaxation.
    if solve_pars.α > 0.0
        fig_cyl = plotsimulationresiduals(model_cyl.res[]; ε = solve_pars.ε)[1]
        fig_sph = plotsimulationresiduals(model_sph.res[]; ε = solve_pars.ε)[1]
    end

    # Compare solutions: the sphere should heat faster!
    fig = let
        fig = Figure(resolution = (720, 500))
        ax = Axis(fig[1, 1], yscale = identity)
        lines!(ax, 100model_cyl.grid.r, model_cyl.problem.x, label = "Cylinder")
        lines!(ax, 100model_sph.grid.r, model_sph.problem.x, label = "Sphere")
        ax.xlabel = "Radial coordinate [cm]"
        ax.ylabel = "Temperature [K]"
        ax.xticks = 0.0:1.0:100last(model_cyl.grid.r)
        axislegend(ax; position = :lt)
        fig
    end
end

begin # (sample 2)
    N = 100
    R = 0.05

    grids = [
        equidistantcellsgrid1D(R, N),
        geometriccellsgrid1D(R, N),
        UserDefinedGrid1D(geomspace(0, R, 0.003, 0.0005))
    ]

    fig = let
        fig = Figure(resolution = (720, 500))
        ax = Axis(fig[1, 1], yscale = identity)

        for grid in grids
            c = 900.0
            ρ = 3000.0
            h = (t) -> (t < 2000) ? 20000.0 : 0.0
            B = (t) -> 1500.0 + 200sin(2π * t / 240)
            κ = (T) -> 2.0 + 0.0001T
            model = Cylinder1DTemperatureModel(; grid, h, B, κ, ρ, c)
            @time solve(model; t = 2400.0, τ = 2.0, T = 300.0, M = 50)
            lines!(ax, 100grid.r, model.problem.x, label = "N = $(grid.N)")
        end

        ax.xlabel = "Radial coordinate [cm]"
        ax.ylabel = "Temperature [K]"
        ax.xticks = 0.0:1.0:100R
        xlims!(ax, (-0.2, 100R+0.2))
        axislegend(ax; position = :lt)

        fig
    end
end
