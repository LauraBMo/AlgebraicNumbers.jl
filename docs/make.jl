using AlgebraicNumbers
using Documenter

DocMeta.setdocmeta!(AlgebraicNumbers, :DocTestSetup, :(using AlgebraicNumbers); recursive=true)

makedocs(;
    modules=[AlgebraicNumbers],
    authors="Laura Brustenga i Moncusí <brust@math.ku.dk> and contributors",
    repo="https://github.com/LauraBMo/AlgebraicNumbers.jl/blob/{commit}{path}#{line}",
    sitename="AlgebraicNumbers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://LauraBMo.github.io/AlgebraicNumbers.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LauraBMo/AlgebraicNumbers.jl",
)
