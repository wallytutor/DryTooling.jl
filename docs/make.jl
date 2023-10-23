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
        "Module Granular"       => "granular.md",
        "Module Simulation"     => "simulation.md",
        "Module Finite Volumes" => "finitevolumes.md",
        "Cantera Wrapper"       => "canterawrapper.md",
        "DryTooling Core"       => "drytooling.md",
        "Reference API"         => "apireference.md",
    ]
)

deploydocs(;
    repo="github.com/wallytutor/DryTooling.jl",
    devbranch="main",
)
