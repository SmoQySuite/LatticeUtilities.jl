push!(LOAD_PATH,"../src/")

using Documenter, LatticeUtilities

DocMeta.setdocmeta!(LatticeUtilities, :DocTestSetup, :(using LatticeUtilities); recursive=true)

makedocs(
    sitename = "LatticeUtilities.jl",
    modules  = [LatticeUtilities],
    pages    = [
        "Home"            => "index.md",
        "Getting Started" => "getting_started.md",
        "Examples"        => "examples.md",
        "API"             => "api.md"
    ]
)

deploydocs(
    repo = "github.com/cohensbw/LatticeUtilities.jl.git",
)
