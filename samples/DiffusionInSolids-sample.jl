# -*- coding: utf-8 -*-
"""
DryTooling.DiffusionInSolids sample
===================================

Module under development... documentation comming soon!
"""

import Pkg
Pkg.activate(Base.current_project())
Pkg.instantiate()

using Revise
using CairoMakie
using DryTooling.Residuals
using DryTooling.Grids
using DryTooling.DiffusionInSolids
using DryTooling.DiffusionInSolids: carburizemoletomassfraction

begin
    N = 100
    L = 0.002
    T = 1173.15
    τ = 10.0

    t = 5*3600.0
    y0 = 0.0023
    ys = 0.0095
    h = (t) -> (t <= 7200.0) ? 1.0e-03 : 0.0

    grid = equidistantcellsgrid1D(L, N)
    model = carburize(grid, t, τ, T, h, y0, ys, ; M = 50)

    z = model.grid.r
    y = carburizemoletomassfraction.(model.problem.x)
    m = interstitialmassintake(model.grid.r, y0, y)

    fig = let
        fig = Figure(resolution = (720, 500))
        ax = Axis(fig[1, 1], yscale = identity)
        lines!(ax, 1000z, 100reverse(y))
        ax.title  = "Mass intake $(round(m, digits = 2)) g/m²"
        ax.xlabel = "Coordinate [cm]"
        ax.ylabel = "Mass percentage [%]"
        ax.xticks = 0.0:0.2:1.2
        ax.yticks = 0.2:0.1:0.6
        xlims!(ax, extrema(ax.xticks.val))
        ylims!(ax, extrema(ax.yticks.val))
        fig
    end

    fig = plotsimulationresiduals(model.res[]; showinner = false)[1]
end