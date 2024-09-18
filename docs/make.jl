using CornerPlotting
using Documenter

DocMeta.setdocmeta!(CornerPlotting, :DocTestSetup, :(using CornerPlotting); recursive=true)

makedocs(;
    modules=[CornerPlotting],
    authors="Pablo Marchant, Alina Istrate, Reinhold Willcox",
    sitename="CornerPlotting.jl",
    format=Documenter.HTML(;
        canonical="https://orlox.github.io/CornerPlotting.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/orlox/CornerPlotting.jl",
    devbranch="main",
)
