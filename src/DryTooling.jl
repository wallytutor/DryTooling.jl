# -*- coding: utf-8 -*-
module DryTooling

begin # module DryTooling
    using CairoMakie
    using CommonSolve
    using DocStringExtensions: TYPEDFIELDS
    using LinearAlgebra
    using Polynomials
    using Roots
    using YAML

    include("DryTooling/abstract.jl")
    include("DryTooling/constants.jl")
    include("DryTooling/utilities.jl")
end

module Cantera
    using Libdl
    using Logging
    using Printf

    include("Cantera/core.jl")
    include("Cantera/pointers.jl")
    include("Cantera/interfaces.jl")
end

# XXX
include("tmp.jl")

module Simulation
    using CairoMakie
    using DocStringExtensions: TYPEDFIELDS
    using LinearAlgebra
    using DryTooling

    include("Simulation/residuals.jl")
    include("Simulation/linalg.jl")
    include("Simulation/nlstepping.jl")
end # module Simulation

module FiniteVolumes
    using CommonSolve
    using CommonSolve: solve
    using DocStringExtensions: TYPEDFIELDS
    using Trapz: trapz
    using DryTooling
    using DryTooling.Simulation
    # using DryTooling.Simulation: fouter!, finner!, fsolve!, timepoints
    # TODO this will migrate!
    using DryTooling: Temperature1DModelStorage
    using DryTooling: interfaceconductivity1D

    include("FiniteVolumes/grid-generation.jl")
    include("FiniteVolumes/heat-conduction.jl")
    include("FiniteVolumes/diffusion-in-solids.jl")
end # module FiniteVolumes

module FluidModels  
end # module FluidModels

module Granular
    using CairoMakie
    using DifferentialEquations: ODEProblem, Tsit5
    using DifferentialEquations: solve
    using Distributions
    using DocStringExtensions: TYPEDFIELDS
    using ModelingToolkit
    using Printf
    using Random
    using Trapz: trapz

    include("Granular/porous-media.jl")
    include("Granular/rotary-kiln.jl")
end # module Granular

module Kinetics
    using ModelingToolkit
    using Symbolics
    using Symbolics: scalarize
    using DryTooling
    # TODO this will migrate!
    using DryTooling: meanmolecularmass

    include("Kinetics/core.jl")
    include("Kinetics/hard-coded-mechs.jl")
end # module Kinetics

module PlugFlow
    using CommonSolve
    using CommonSolve: solve
    using DifferentialEquations: solve
    using ModelingToolkit
    using Symbolics
    using Symbolics: scalarize
    using DryTooling
    using DryTooling.Kinetics

    include("PlugFlow/core.jl")
end # module PlugFlow

using DryTooling.Simulation
using DryTooling.FiniteVolumes
using DryTooling.FluidModels
using DryTooling.Granular
using DryTooling.Kinetics
using DryTooling.PlugFlow

end
