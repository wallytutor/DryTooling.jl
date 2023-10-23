# -*- coding: utf-8 -*-
import Pkg
using Revise

if Base.current_project() != Base.active_project()
    Pkg.activate(Base.current_project())
    Pkg.instantiate()
end

using DryTooling
using Documenter

# DocMeta.setdocmeta!(DryTooling, :DocTestSetup, :(using DryTooling); recursive=false)

makedocs(;
    modules=[DryTooling],
    authors="Walter Dal'Maz Silva <walter.dalmazsilva.manager@gmail.com> & contributors",
    repo="https://github.com/wallytutor/DryTooling.jl/blob/{commit}{path}#{line}",
    sitename="DryTooling.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://wallytutor.github.io/DryTooling.jl",
        repolink="https://github.com/wallytutor/DryTooling.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home"            => "index.md",
        "Module Granular" => "granular.md",
        "Cantera API"     => "DocsCanteraAPI.md",
        "Reference API"   => "references.md",
    ],
)
