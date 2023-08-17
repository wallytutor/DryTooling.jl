# -*- coding: utf-8 -*-
using Plots
using DryTooling.Kramers
import DryTooling.Kramers as kramers

function gettestbed()
    R = 1.0e+00
    Φ = 1.0e-02
    z = collect(0.0:0.1:10.0)
    h = (0.5R) * ones(size(z))
    return RotaryKilnBedSolution(z, h, R, Φ)
end

function testsolverSI(L, D, Φ, ω, β, γ, d)
    bed = kramers.solvelinearkramersmodelSI(;
        model = SymbolicLinearKramersModel(),
        L = L,
        R = D / 2.0,
        Φ = Φ / 3600.0,
        ω = ω / 60.0,
        β = deg2rad(β),
        γ = deg2rad(γ),
        d = d / 1000.0
    )

    p = plotlinearkramersmodel(bed, normz=true, normh=true)
    display(p)

    return bed
end

function testsolver(L, D, Φ, ω, β, γ, d)
    bed = solvelinearkramersmodel(
        L = L,
        D = D,
        Φ = Φ,
        ω = ω,
        β = β,
        γ = γ,
        d = d,
        model = nothing
    )

    p = plotlinearkramersmodel(bed, normz=true, normh=true)
    display(p)

    return bed
end

# Data from last example of reference paper.

# Kiln length [m].
L = 13.715999999999998

# Kiln diameter [m].
D = 1.8897599999999999

# Feed rate [m³/h].
Φ = 10.363965852671996

# Rotation rate [rev/min].
ω = 3.0300000000000002

# Kiln slope [°].
β = 2.3859440303888126

# Bed repose angle [°].
γ = 45.0

# Particle size [mm].
d = 1.0

bed1 = testsolverSI(L, D, Φ, ω, β, γ, d)
bed2 = testsolver(L, D, Φ, ω, β, γ, d)