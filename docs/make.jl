using HydroChem
using Documenter

DocMeta.setdocmeta!(HydroChem, :DocTestSetup, :(using HydroChem); recursive=true)

makedocs(;
    modules=[HydroChem],
    authors="Siddhant A. Deshmukh <siddhant593@gmail.com>",
    repo="https://github.com/chongchonghe/Hydro.jl/blob/{commit}{path}#{line}",
    sitename="Hydro.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://chongchonghe.github.io/Hydro.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(;
    repo="github.com/chongchonghe/Hydro.jl",
    devbranch="main"
)
