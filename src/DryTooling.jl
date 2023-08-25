# -*- coding: utf-8 -*-
module DryTooling

using ModelingToolkit
using Plots
using Printf

using DocStringExtensions: TYPEDFIELDS
using DifferentialEquations: ODEProblem, Tsit5
using DifferentialEquations: solve
using Trapz: trapz

include("DryTooling/constants.jl")
include("DryTooling/utility-functions.jl")
include("DryTooling/chemical-elements.jl")
include("DryTooling/gas-phase-models.jl")
include("DryTooling/kramers-model.jl")
# include("CanteraAPI.jl")

end
