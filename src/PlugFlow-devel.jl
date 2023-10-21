# -*- coding: utf-8 -*-
"""
DryTooling.PlugFlow sample
==========================

Module under development... documentation comming soon!
"""

import Pkg
Pkg.activate(Base.current_project())
Pkg.instantiate()

using Revise
using CairoMakie
using DifferentialEquations: solve
using DryTooling.Constants
using DryTooling.Grids
using LaTeXStrings
using ModelingToolkit
using Symbolics
using Symbolics: scalarize
using YAML

# XXX: these will be placed in FluidModels later
using DryTooling: concentration
using DryTooling: densitymass
using DryTooling: idealgasdensity
using DryTooling: meanmolecularmass
using DryTooling: molefraction2massfraction
using DryTooling: massfraction2molefraction
using DryTooling: IdealGasMixture

struct Graf2007AcetyleneKinetics
    r::Vector{Num}
    ω::Vector{Num}

    function Graf2007AcetyleneKinetics(T::Num, X::Symbolics.Arr{Num, 1})
        RT = GAS_CONSTANT * T

        # Reactions coefficients.
        ν = [
            -1  1 -1  1 -1 -2  2 -1  0
            -1  1 -3  3  1  0  0  0  3
             1 -1  0  0  0  0  0  0  0
             0  0  2 -2  0  0  0  0  0
             0  0  0  0  0  1 -1 -1  0
             0  0  0  0  0  0  0  1 -1
             0  0  0  0  2  0  0  0  6
             0  0  0  0  0  0  0  0  0
        ]

        # Rate constants.
        k = [
            4.4e+03 * exp(-1.0300e+05 / RT)
            3.8e+07 * exp(-2.0000e+05 / RT)
            1.4e+05 * exp(-1.5000e+05 / RT)
            8.6e+06 * exp(-1.9500e+05 / RT)
            5.5e+06 * exp(-1.6500e+05 / RT)
            1.2e+05 * exp(-1.2070e+05 / RT)
            1.0e+15 * exp(-3.3520e+05 / RT)
            1.8e+03 * exp(-6.4500e+04 / RT)
            1.0e+03 * exp(-7.5000e+04 / RT)
        ]

        # Reaction rates in molar units.
        # TODO compute using mole fractions and prepend with (P/RT)^k!
        r = k .* [
            abs(X[1]) * abs(X[2])^0.36
            abs(X[3])^0.50
            abs(X[1])^0.35 * abs(X[2])^0.22
            abs(X[4])^0.21
            abs(X[1])^1.90 / (1.0 + 18.0*X[2])
            abs(X[1])^1.60
            abs(X[5])^0.75
            abs(X[1])^1.30 * X[5]^0.60
            abs(X[6])^0.75 / (1.0 + 22.0*X[2])
        ]

        # Species production rates.
        ω = scalarize(ν * r)

        return new(r, ω)
    end
end

loaddatabase(fpath) = YAML.load_file(realpath(fpath))
standarddensity(M) = idealgasdensity(ZERO_CELSIUS, ONE_ATM, M)
sccmtomdot(q, M) = standarddensity(M) * q / 6.0e+07

# Dimensions of reactor [m].
R = 0.014
L = 0.5

# Load database for retrieving thermodynamic data.
mech = "data/CT-hydrocarbon-dalmazsi-2017-mech.yaml"
data = loaddatabase(joinpath(@__DIR__, mech))

# Create gas phase only once.
selected = ["C2H2", "H2", "C2H4", "CH4", "C4H4", "C6H6", "Cs", "N2"]
mix = IdealGasMixture(data, selected)

# Create a grid for the reactor.
grid = equidistantcellsgrid1D(L, 20)

# Cross-section of reactor [m²].
A = (π * R^2) * ones(Float64, grid.N)

# Interpolated area at cell centers.
Ac = (1/2)*(A[1:end-1] + A[2:end])

# Create symbolic parameters (user-provided).
pars = @parameters T P ṁ

# Create symbolic model independent and solution variables.
@variables z (Y(z))[1:mix.nspecies]

# Also create the observables (rates) as variables.
@variables (RHS(z))[1:mix.nspecies]

# Create differential.
D = Differential(z)

# Compute intermediates.
W = mix.molecularmasses
X = concentration(mix, T, P, Y)

# Get mechanism with rates.
kin = Graf2007AcetyleneKinetics(T, X)

# Product ρu is constant in this case.
# TODO: generalize to A(z).
ρu = ṁ / A[1]
scaler = @. W / ρu

# Assembly problem equations.
rhs = @. kin.ω * scaler
eqs = [scalarize(D.(Y) .~ RHS)...,
       scalarize(RHS .~ rhs)...]

@named model = ODESystem(eqs, z, [Y..., RHS...], pars)
system = structural_simplify(model)

Y0 = let
    # Mole fraction of acetylene in system.
    x1 = 0.36

    # Add acetylene impurities to initialization.
    X0 = zeros(mix.nspecies)
    X0[1] = 0.980 * x1
    X0[4] = 0.002 * x1
    X0[end] = 1 - sum(X0[1:end-1])

    # Convert to mass fractions for the model.
    molefraction2massfraction(X0, W)
end


M = meanmolecularmass(Y0, W)

p = [
    ṁ => sccmtomdot(222.0, 1000M),
    P => 5000.0,
    T => 1173.15
]

prob = ODEProblem(system, Y0, (0.0, L), p)
sol = solve(prob; reltol=1.0e-03, abstol=1.0e-08, saveat = grid.r)

density(Y) = 
X = map((y)->massfraction2molefraction(y, W), sol.u)
ρ = map((y)->densitymass(mix, p[3].second, p[2].second, y), sol.u)
u = @. (p[1].second / (ρ * A[1]))

X = mapreduce(permutedims, vcat, X)

fig =  let 
    f = Figure(resolution = (800, 600))

    f[1, 1] = GridLayout()
    f[1, 2] = GridLayout()
    f[2, 1] = GridLayout()
    f[2, 2] = GridLayout()

    ax1 = Axis(f[1, 1])
    ax2 = Axis(f[1, 2])
    ax3 = Axis(f[2, 1])
    ax4 = Axis(f[2, 2])

    linkxaxes!(ax1, ax2, ax3, ax4)

    lines!(ax1, 100grid.r, 100X[:, 1], label = L"\mathrm{C_2H_2}")
    lines!(ax2, 100grid.r, 100X[:, 2], label = L"\mathrm{H_2}")
    lines!(ax2, 100grid.r, 100X[:, 3], label = L"\mathrm{C_2H_4}")
    lines!(ax2, 100grid.r, 100X[:, 4], label = L"\mathrm{CH_4}")
    lines!(ax3, 100grid.r, 100X[:, 5], label = L"\mathrm{C_4H_4}")
    lines!(ax3, 100grid.r, 100X[:, 6], label = L"\mathrm{C_6H_6}")
    lines!(ax3, 100grid.r, 100X[:, 7], label = L"\mathrm{C(s)}")
    lines!(ax4, 100grid.r, u)

    ax1.ylabel = "Mole percentage [%]"
    ax2.ylabel = "Mole percentage [%]"
    ax3.ylabel = "Mole percentage [%]"
    ax4.ylabel = "Velocity [m/s]"

    ax3.xlabel = "Coordinate [cm]"
    ax4.xlabel = "Coordinate [cm]"

    ax1.xticks = 0.0:10:50.0
    ax2.xticks = 0.0:10:50.0
    ax3.xticks = 0.0:10:50.0
    ax4.xticks = 0.0:10:50.0

    xlims!(ax1, extrema(ax1.xticks.val))
    xlims!(ax2, extrema(ax1.xticks.val))
    xlims!(ax3, extrema(ax1.xticks.val))
    xlims!(ax4, extrema(ax1.xticks.val))

    axislegend(ax1; position = :rt)
    axislegend(ax2; position = :lt)
    axislegend(ax3; position = :lt)

    f
end
