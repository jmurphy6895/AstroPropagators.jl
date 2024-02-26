using Documenter
using AstroForceModels

makedocs(
    modules = [AstroForceModels],
    format=Documenter.HTML(; prettyurls=!("local" in ARGS), highlights=["yaml"], ansicolor=true), 
    sitename = "AstroForceModels.jl",
    authors = "Jordan Murphy",
    pages = [
        "Home" => "index.md",
        "Usage" => "man/usage.md",
        "API" => "man/api.md",
        "Library" => "lib/library.md",
    ],
)

deploydocs(
    repo = "github.com/jmurphy6895/AstroForceModels.jl.git",
    target = "build",
)
