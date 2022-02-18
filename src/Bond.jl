"""
    Bond

Defines a bond defintion in lattice.

# Fields
$(TYPEDFIELDS)
"""
struct Bond

    "Dimension of system bond is in."
    D::Int

    "Initial and final orbital species respectively."
    o::Vector{Int}

    "Displacement in unit cells."
    Δl::Vector{Int}
end

"""
    Bond(Δl::AbstractVector{Int},o::AbstractVector{Int})

Constrcut a [`Bond`](@ref)
"""
function Bond(o::AbstractVector{Int},Δl::AbstractVector{Int})

    D = length(Δl)
    return  Bond(D,o,Δl)
end