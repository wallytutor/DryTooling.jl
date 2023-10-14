# -*- coding: utf-8 -*-
module DryTooling

using CommonSolve
using Distributions
using LinearAlgebra
using ModelingToolkit
using Polynomials
using Plots
using Printf
using Random
using YAML

using DocStringExtensions: TYPEDFIELDS
using DifferentialEquations: ODEProblem, Tsit5
using DifferentialEquations: solve
using Trapz: trapz

# TODO: remove dependence from Plots!
import CairoMakie as CM

include("DryTooling/abstract-types.jl")
include("DryTooling/constants.jl")
include("DryTooling/utility-functions.jl")

include("DryTooling/simulation-residuals.jl")
include("DryTooling/tridiagonal-problem.jl")
include("DryTooling/iterative-solver.jl")

include("DryTooling/chemical-elements.jl")
include("DryTooling/chemical-conversions.jl")
include("DryTooling/gas-phase-models.jl")
include("DryTooling/solid-phase-models.jl")

include("DryTooling/kramers-model.jl")
include("DryTooling/porous-media.jl")

include("DryTooling/grid-generation.jl")
include("DryTooling/diffusion-models.jl")
include("DryTooling/plug-flow-reactors.jl")

include("DryTooling/base-extensions.jl")

# include("CanteraAPI.jl")

end
