using CornerPlotting
using Documenter
using Literate

# Parse examples using Literate
pkg_path = pkgdir(CornerPlotting)

function ignore_code_blocks(content)
    content = replace(content, "##\n" => "\n")  # remove code blocks
    content = replace(content, "###" => "##")  # make level 3 headers level 2
end

Literate.markdown(pkg_path * "/examples/tutorial.jl", pkg_path * "/docs/src/", preprocess=ignore_code_blocks)

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
        "Tutorial" => "tutorial.md",
    ],
)

deploydocs(;
    repo="github.com/orlox/CornerPlotting.jl",
    devbranch="main",
)
