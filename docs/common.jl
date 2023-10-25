# -*- coding: utf-8 -*-
using Documenter
using DocumenterCitations
using DryTooling

bib_filepath = joinpath(@__DIR__, "src/references.bib")
bib = CitationBibliography(bib_filepath)

DocMeta.setdocmeta!(DryTooling, :DocTestSetup, :(using DryTooling); recursive=true)

# format = Documenter.LaTeX()

format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical  = "https://wallytutor.github.io/DryTooling.jl",
    repolink   = "https://github.com/wallytutor/DryTooling.jl",
    edit_link  = "main",
    assets     = String[],
)

pages  = [
    "Home"                  => "index.md",
    
    "Module Finite Volumes" => [
        "Documentation" => "FiniteVolumes/index.md",
        "Examples"      => "FiniteVolumes/samples.md",
        "Theory guide"  => "FiniteVolumes/theory.md",
    ],

    "Module Granular"       => [
        "Documentation" => "Granular/index.md",
        "Examples"      => "Granular/samples.md",
        "Theory guide"  => "Granular/theory.md",
    ],
    
    "Module Cantera"        => "Cantera/index.md",
    "Module Thermodynamics" => "Thermodynamics/index.md",
    "Module Kinetics"       => "Kinetics/index.md",
    "Module Plug Flow"      => "PlugFlow/index.md",
    
    "Module Simulation"     => "Simulation/index.md",
    "DryTooling Core"       => "DryTooling/index.md",
    
    "Reference API"         => "api.md",
    "Table of contents"     => "toc.md",
    "References"         => "references.md",
]

makedocs(;
    modules  = [DryTooling],
    authors  = "Walter Dal'Maz Silva <walter.dalmazsilva.manager@gmail.com> & contributors",
    repo     = "https://github.com/wallytutor/DryTooling.jl/blob/{commit}{path}#{line}",
    sitename = "DryTooling.jl",
    format   = format,
    plugins  = [bib],
    pages    = pages
)
