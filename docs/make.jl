using Documenter
using GlobalHeightFunctions

makedocs(
    modules = [GlobalHeightFunctions],
    authors = "PenguinxCutCell contributors",
    sitename = "GlobalHeightFunctions.jl",
    format = Documenter.HTML(
        canonical = "https://PenguinxCutCell.github.io/GlobalHeightFunctions.jl",
        repolink = "https://github.com/PenguinxCutCell/GlobalHeightFunctions.jl",
        collapselevel = 2,
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "reference.md",
    ],
    pagesonly = true,
    warnonly = true,
    remotes = nothing,
)

if get(ENV, "CI", "") == "true"
    deploydocs(
        repo = "github.com/PenguinxCutCell/GlobalHeightFunctions.jl",
        push_preview = true,
    )
end
