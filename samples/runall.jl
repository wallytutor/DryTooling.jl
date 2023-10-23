# -*- coding: utf-8 -*-
import Pkg
using Revise

if Base.current_project() != Base.active_project()
    Pkg.activate(Base.current_project())
    Pkg.resolve()
    Pkg.instantiate()
end

# include("CanteraAPI-sample.jl")
include("DiffusionInSolids-sample.jl")
include("HeatConduction-sample.jl")
include("Kinetics-sample.jl")
include("PlugFlow-sample.jl")

