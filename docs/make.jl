# -*- coding: utf-8 -*-
using DryTooling
using Documenter

DocMeta.setdocmeta!(DryTooling, :DocTestSetup, :(using DryTooling); recursive=true)

format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical  = "https://wallytutor.github.io/DryTooling.jl",
    repolink   = "https://github.com/wallytutor/DryTooling.jl",
    edit_link  = "main",
    assets     = String[],
)

makedocs(;
    modules  = [DryTooling],
    authors  = "Walter Dal'Maz Silva <walter.dalmazsilva.manager@gmail.com> & contributors",
    repo     = "https://github.com/wallytutor/DryTooling.jl/blob/{commit}{path}#{line}",
    sitename = "DryTooling.jl",
    format   = format,
    pages    = [
        "Home"                  => "index.md",
        "Table of contents"     => "toc.md",

        "Module Finite Volumes" => "FiniteVolumes/index.md",
        "Module Granular"       => "Granular/index.md",

        "Module Cantera"        => "Cantera/index.md",
        "Module Thermodynamics" => "Thermodynamics/index.md",
        "Module Kinetics"       => "Kinetics/index.md",
        "Module Plug Flow"      => "PlugFlow/index.md",
        
        "Module Simulation"     => "Simulation/index.md",
        "DryTooling Core"       => "DryTooling/index.md",

        "Reference API"         => "api.md",
    ]
)

deploydocs(;
    repo="github.com/wallytutor/DryTooling.jl",
    devbranch="main",
)
