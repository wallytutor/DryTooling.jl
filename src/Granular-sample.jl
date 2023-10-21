# -*- coding: utf-8 -*-
"""
DryTooling.Granular sample
==========================

"""

import Pkg
Pkg.activate(Base.current_project())
Pkg.instantiate()

using Revise
using GLMakie
using DryTooling.Granular

begin # (sample 1)
    @info """\
    PackedBedPorosityDescriptor can be used to describe the geometry \
    of exchange section of a packed bed for a single set of arguments.
    """
    p1 = PackedBedPorosityDescriptor(;
        ϕ = 0.65,
        l = 0.10, 
        area = 1.0
    )
end

begin # (sample 2)
    @info """\
    PackedBedPorosityDescriptor can also be used to describe randomly \
    varying reactors, what is a more realistic thing to do when using \
    this structure to simulate real world systems.
    """
    p2 = PackedBedPorosityDescriptor(;
        ϕ = 0.65,
        l = 0.10, 
        σϕ = 0.03,
        σl = 0.01,
        N = 2,
        ϕlims = (0.4, 0.8),
        llims = (0.0, 0.3),
        seed = 42,
        area = 1.0
    )
end

begin # (sample 3)
    # Geometry
    Lᵣ = 13.715999999999998  # Kiln length [m]
    Dᵣ = 1.8897599999999999  # Kiln diameter [m]
    βᵣ = 2.3859440303888126  # Kiln slope [°]

    # Material
    γᵣ = 45.0                # Repose angle [°]
    dᵣ = 1.0                 # Particle/dam size [mm]

    # Process
    Φᵣ = 10.363965852671996  # Feed rate [m³/h]
    ωᵣ = 3.0300000000000002  # Rotation rate [rev/min]

    bed = solvelinearkramersmodel(;
        model = SymbolicLinearKramersModel(),
        L     = L,
        R     = D / 2.0,
        Φ     = Φ / 3600.0,
        ω     = ω / 60.0,
        β     = deg2rad(β),
        γ     = deg2rad(γ),
        d     = d / 1000.0
    )
    fig = plotlinearkramersmodel(bed, normz=false, normh=false)
end