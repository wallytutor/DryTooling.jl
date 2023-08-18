var documenterSearchIndex = {"docs":
[{"location":"kramers/#DryTooling.Kramers","page":"Kramers","title":"DryTooling.Kramers","text":"","category":"section"},{"location":"kramers/","page":"Kramers","title":"Kramers","text":"Implements the ordinary differential equation for prediction of bed height profile in a rotary kiln as proposed by Kramers and Croockewite (1952) [1]. Its goal is to be used as a process support tool or to integrate more complex models requiring integration of the bed profile.","category":"page"},{"location":"kramers/","page":"Kramers","title":"Kramers","text":"","category":"page"},{"location":"kramers/","page":"Kramers","title":"Kramers","text":"Modules = [\n    DryTooling.Kramers,\n]","category":"page"},{"location":"kramers/#DryTooling.Kramers.RotaryKilnBedSolution","page":"Kramers","title":"DryTooling.Kramers.RotaryKilnBedSolution","text":"RotaryKilnBedSolution\n\nDescription of a rotary kiln bed geometry computed from the solution of bed height along the kiln length. The main goal of the quantities computed here is their use with heat and mass transfer models for the simulation of rotary kiln process.\n\nz::Vector{Float64}: Solution coordinates [m]\nh::Vector{Float64}: Solution bed height [m]\nθ::Vector{Float64}: View angle from kiln center [rad]\nl::Vector{Float64}: Bed-freeboard cord length [m]\nA::Vector{Float64}: Local bed cross section area [m²]\nη::Vector{Float64}: Local loading based on height [-]\nηₘ::Float64: Mean loading of kiln [%]\nV::Float64: Bed integral volume [m³]\nτ::Float64: Residence time of particles\n\n\n\n\n\n","category":"type"},{"location":"kramers/#DryTooling.Kramers.RotaryKilnBedSolution-NTuple{4, Any}","page":"Kramers","title":"DryTooling.Kramers.RotaryKilnBedSolution","text":"RotaryKilnBedSolution(\n    z::Vector{Float64},\n    h::Vector{Float64},\n    R::Float64,\n    Φ::Float64\n)\n\nz: solution coordinates over length, [m].\nh: bed profile solution over length, [m].\nR: kiln internal radius, [m].\nΦ: kiln feed rate, [m³/s].\n\n\n\n\n\n","category":"method"},{"location":"kramers/#DryTooling.Kramers.SymbolicLinearKramersModel","page":"Kramers","title":"DryTooling.Kramers.SymbolicLinearKramersModel","text":"SymbolicLinearKramersModel\n\nCreates a reusable linear Kramers model for rotary kiln simulation.\n\nR::Symbolics.Num: Symbolic kiln internal radius\nΦ::Symbolics.Num: Symbolic kiln feed rate\nω::Symbolics.Num: Symbolic kiln rotation rate\nβ::Symbolics.Num: Symbolic kiln slope\nγ::Symbolics.Num: Symbolic solids repose angle\nz::Symbolics.Num: Symbolic kiln axial coordinates\nh::Symbolics.Num: Symbolic bed height profile\nsys::ModelingToolkit.ODESystem: Problem ordinary differential equation\n\n\n\n\n\n","category":"type"},{"location":"kramers/#DryTooling.Kramers.SymbolicLinearKramersModel-Tuple{}","page":"Kramers","title":"DryTooling.Kramers.SymbolicLinearKramersModel","text":"SymbolicLinearKramersModel()\n\nSymbolic model constructor.\n\n\n\n\n\n","category":"method"},{"location":"kramers/#DryTooling.Kramers.plotlinearkramersmodel-Tuple{DryTooling.Kramers.RotaryKilnBedSolution}","page":"Kramers","title":"DryTooling.Kramers.plotlinearkramersmodel","text":"plotlinearkramersmodel(\n    model::SolutionLinearKramersModel;\n    normz::Bool = false,\n    normh::Bool = false\n)::Any\n\nDisplay plot of model solution for rotary kiln bed profile. Arguments normz and normh control whether z-coordinate and bed height must be normalized, respectively.\n\n\n\n\n\n","category":"method"},{"location":"kramers/#DryTooling.Kramers.solvelinearkramersmodel-Tuple{}","page":"Kramers","title":"DryTooling.Kramers.solvelinearkramersmodel","text":"solvelinearkramersmodel(;\n    model::SymbolicLinearKramersModel,\n    L::Float64,\n    R::Float64,\n    Φ::Float64,\n    ω::Float64,\n    β::Float64,\n    γ::Float64,\n    d::Float64,\n    solver::Any = Tsit5(),\n    rtol::Float64 = 1.0e-08,\n    atol::Float64 = 1.0e-08\n)\n\nIntegrates an instance of SymbolicLinearKramersModel.\n\nImportant: inputs must be provided in international system (SI) units  as a better physical practice. The only exception is the rotation rate ω provided in revolution multiples. If the discharge end is held by a dam, its height must be provided instead of the particle size, as it is used as the ODE initial condition.\n\nmodel: a symbolic kiln model.\nL: kiln length, [m].\nR: kiln internal radius, [m].\nΦ: kiln feed rate, [m³/s].\nω: kiln rotation rate, [rev/s].\nβ: kiln slope, [rad].\nγ: solids repose angle, [rad].\nd: particle size or dam height, [m].\n\n\n\n\n\n","category":"method"},{"location":"kramers/","page":"Kramers","title":"Kramers","text":"[1]: Kramers et al., 1952","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = DryTooling","category":"page"},{"location":"#DryTooling","page":"Home","title":"DryTooling","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for DryTooling.","category":"page"}]
}