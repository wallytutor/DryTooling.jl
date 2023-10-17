# -*- coding: utf-8 -*-
import Pkg
Pkg.activate(Base.current_project())
Pkg.instantiate()

using Revise
import DryTooling as dry

R = 0.05
tol = 1.0e-10

grid = dry.CylinderGrid1DEquispaced(R, 100)

h = (t) -> (t < 2300) ? 20.0 : 0.0
B = 1500.0
κ = 2.0
ρ = 3000.0
c = 900.0

model = dry.Cylinder1DTemperatureModel(grid, h, B, κ, ρ, c)

@time residuals = dry.solve(model;
    T     = 300.0,
    τ     = 2.0,
    t     = 2400.0,
    iters = 20,
    relax = 0.0001,
    tol   = tol,
    metric = dry.maxrelativevariation
)
