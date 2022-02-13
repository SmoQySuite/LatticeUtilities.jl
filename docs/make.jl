push!(LOAD_PATH,"../src/")

using Documenter, LatticeUtilities

DocMeta.setdocmeta!(LatticeUtilities, :DocTestSetup, :(using LatticeUtilities); recursive=true)

makedocs(
    sitename = "LatticeUtilities",
    format = Documenter.HTML(),
    modules = [LatticeUtilities],
    doctest = true
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
