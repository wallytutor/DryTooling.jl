# -*- coding: utf-8 -*-
using Plots
using DryTooling.Kramers
import DryTooling.Kramers as kramers

function solvekramers1952()
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

    return solvelinearkramersmodel(;
        model = SymbolicLinearKramersModel(),
        L     = L,
        R     = D / 2.0,
        Φ     = Φ / 3600.0,
        ω     = ω / 60.0,
        β     = deg2rad(β),
        γ     = deg2rad(γ),
        d     = d / 1000.0
    )
end

bed = solvekramers1952()
p = plotlinearkramersmodel(bed, normz=true, normh=true)
display(p)
