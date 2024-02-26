using Documenter
using HAMMERHEAD

makedocs(
    modules = [HAMMERHEAD],
    format=Documenter.HTML(; prettyurls=!("local" in ARGS), highlights=["yaml"], ansicolor=true), 
    sitename = "HAMMERHEAD.jl",
    authors = "Jordan Murphy",
    pages = [
        "Home" => "index.md",
        "Usage" => "man/usage.md",
        "API" => "man/api.md",
        "Library" => "lib/library.md",
    ],
)

deploydocs(
    repo = "github.com/jmurphy6895/HAMMERHEAD.jl.git",
    target = "build",
)
