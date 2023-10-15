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
selected = ["C2H2", "H2", "C2H4", "CH4", "C4H4", "C6H6", "Cs", "N2"]
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

# Create symbolics.
pars = @parameters T P ṁ
@variables z (Y(z))[1:mix.nspecies] (RHS(z))[1:mix.nspecies]

# Create differential.
D = Differential(z)

# Compute intermediates.
W = mix.molecularmasses
X = dry.concentration(mix, T, P, Y)
RT = dry.GAS_CONSTANT * T

# Reaction rates in molar units.
r = [
    4.4e+03 * exp(-1.0300e+05 / RT) * abs(X[1]) * abs(X[2])^0.36
    3.8e+07 * exp(-2.0000e+05 / RT) * abs(X[3])^0.50
    1.4e+05 * exp(-1.5000e+05 / RT) * abs(X[1])^0.35 * abs(X[2])^0.22
    8.6e+06 * exp(-1.9500e+05 / RT) * abs(X[4])^0.21
    5.5e+06 * exp(-1.6500e+05 / RT) * abs(X[1])^1.90 / (1.0 + 18.0*X[2])
    1.2e+05 * exp(-1.2070e+05 / RT) * abs(X[1])^1.60
    1.0e+15 * exp(-3.3520e+05 / RT) * abs(X[5])^0.75
    1.8e+03 * exp(-6.4500e+04 / RT) * abs(X[1])^1.30 * X[5]^0.60
    1.0e+03 * exp(-7.5000e+04 / RT) * abs(X[6])^0.75 / (1.0 + 22.0*X[2])
]

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

# Species production rates.
ω = scalarize(ν * r)

# # Product ρu is constant in this case.
# # TODO: generalize to A(z).
ρu = ṁ / A[1]
scaler = @. W / ρu

rhs = @. ω * scaler
eqs = [scalarize(D.(Y) .~ RHS)..., scalarize(RHS .~ rhs)...]

@named model = ODESystem(eqs, z, [Y..., RHS...], pars)

system = structural_simplify(model)

y1 = 0.35
Y0 = zeros(mix.nspecies)
Y0[1] = 0.980 * y1
Y0[4] = 0.002 * y1
Y0[end] = 1 - sum(Y0[1:end-1])

p = [ṁ => 1.0e-05, P => 5000.0, T => 1173.15]
prob = ODEProblem(system, Y0, (0.0, L), p)
sol = solve(prob; reltol=1.0e-03, abstol=1.0e-08)
# plot(sol[RHS])

molefractions(Y) = dry.massfraction2molefraction(Y, W)
density(Y) = dry.densitymass(mix, p[3].second, p[2].second, Y)

X = [molefractions(y) for y in sol.u]
ρ = [density(y) for y in sol.u]
u = (p[1].second / (ρ * A[1]))'

x1 = [x[2] for x in X]
plot(sol.t, x1)
# dry.massfraction2molefraction(W, sol.u[2])
