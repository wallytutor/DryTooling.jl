using DryTooling
using Documenter

DocMeta.setdocmeta!(DryTooling, :DocTestSetup, :(using DryTooling); recursive=true)

makedocs(;
    modules=[DryTooling],
    authors="Walter Dal'Maz Silva <walter.dalmazsilva.manager@gmail.com> & contributors",
    repo="https://github.com/wallytutor/DryTooling.jl/blob/{commit}{path}#{line}",
    sitename="DryTooling.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://wallytutor.github.io/DryTooling.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home"      => "index.md",
        "General"   => "DocsGeneral.md",
        "IdealGas"  => "DocsIdealGas.md",
        "Transport" => "DocsTransport.md",
        "Kramers"   => "DocsKramers.md",
    ],
)

deploydocs(;
    repo="github.com/wallytutor/DryTooling.jl",
    devbranch="main",
)
