push!(LOAD_PATH,"../src/")

using Documenter, LatticeUtilities

DocMeta.setdocmeta!(LatticeUtilities, :DocTestSetup, :(using LatticeUtilities); recursive=true)

makedocs(
    sitename = "LatticeUtilities.jl",
    modules  = [LatticeUtilities],
    doctest  = true,
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages    = [
        "index.md",
        "Getting Started" => "getting_started.md",
        "Examples" => "examples.md",
        "API" => "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/cohensbw/LatticeUtilities.jl.git"
)
