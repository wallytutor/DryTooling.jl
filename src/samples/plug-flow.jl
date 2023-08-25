# -*- coding: utf-8 -*-
import DryTooling as dry
using DifferentialEquations: solve
using Symbolics: scalarize
using ModelingToolkit
using Plots
using YAML

loaddatabase(fpath) = YAML.load_file(realpath(fpath))
data = loaddatabase(joinpath(pwd(), "src/samples/plug-flow.yaml"))

# Create gas phase only once.
selected = ["C2H2", "H2", "C2H4", "CH4", "C4H4", "C6H6", "N2"]
mix = dry.IdealGasMixture(data, selected)

# Radius of reactor [m].
R = 0.014

# Reactor length [m].
L = 1.0

# Number of discrete cells.
n = 100

# Cell length [m].
δ = L / n

# Space coordinates of reactor [m].
zn = collect(0:δ:L)

# Cross-section of reactor [m²].
A = (π * R^2) * ones(length(zn))

# Interpolated area at cell centers.
Ac = (1/2)*(A[1:end-1] + A[2:end])


# function pfrfactory()
# Create symbolics.
pars = @parameters T P ṁ
@variables z (Y(z))[1:mix.nspecies]
#(RHS(z))[1:mix.nspecies]

# Create differential.
D = Differential(z)

# Compute intermediates.
W = mix.molecularmasses
X = dry.concentration(mix, T, P, Y)
RT = dry.GAS_CONSTANT * T

# Reaction rates in molar units.
r = [
    4.4e+03 * exp(-1.0300e+05 / RT) * X[1] * X[2]^0.36
    3.8e+07 * exp(-2.0000e+05 / RT) * X[3]^0.50
    1.4e+05 * exp(-1.5000e+05 / RT) * X[1]^0.35 * X[2]^0.22
    8.6e+06 * exp(-1.9500e+05 / RT) * X[4]^0.21
    1.2e+05 * exp(-1.2070e+05 / RT) * X[1]^1.60
    1.0e+15 * exp(-3.3520e+05 / RT) * X[5]^0.75
    1.8e+03 * exp(-6.4500e+04 / RT) * X[1]^1.30 * X[5]^0.60
]

# Reactions coefficients.
ν = [
    -1   1  -1   1  -2   2  -1
    -1   1  -3   3   0   0   0
     1  -1   0   0   0   0   0
     0   0   2  -2   0   0   0
     0   0   0   0   1  -1  -1
     0   0   0   0   0   0   1
     0   0   0   0   0   0   0
]

# Species production rates.
ω = scalarize(ν * r)

# Product ρu is constant in this case.
# TODO: generalize to A(z).
ρu = ṁ / A[1]

rhs = @. ω * W / ρu
eqs = scalarize(D.(Y) .~ rhs)

    # return eqs
# end

# eqs = pfrfactory()
@named model = ODESystem(eqs, z, [Y...], pars)

Yz = zeros(mix.nspecies)
Yz[1] = 0.35
Yz[7] = 0.65

p = [ṁ => 1.0e-05, P => dry.ONE_ATM, T => 1073.15]
prob = ODEProblem(model, Yz, (0.0, 0.5), p)
sol = solve(prob)
plot(sol)

# @variables z (A(z))[1:2] (RHS(z)

# D = Differential(z)

# eq1 = 