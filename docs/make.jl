using Documenter
using LatticeUtilities

push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "LatticeUtilities",
    pages = [
        "index.md" 
    ],
    format = Documenter.HTML(),
    modules = [LatticeUtilities]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
