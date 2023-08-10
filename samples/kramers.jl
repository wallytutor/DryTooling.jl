# -*- coding: utf-8 -*-
using DryTooling.Kramers
using Plots

# Conversions for dealing with paper's imperial units,
in1_to_m1(v) = 0.0254 * v
ft1_to_m1(v) = in1_to_m1(12.0) * v
ft3_to_m3(v) = ft1_to_m1(1.0)^3 * v

# From last example of reference paper.
sol = solvelinearkramersmodel(
    L = ft1_to_m1(45.0),
    D = 2 * ft1_to_m1(3.1),
    Φ = ft3_to_m3(6.1) * 60,
    ω =  0.0505 * 60.0,
    β = rad2deg(atan(0.5 / 12)),
    γ = 45.0,
    d = 0.001,
    model = nothing
)

