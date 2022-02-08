"""
    Lattice

A type defining a finite lattice in arbitary dimensions.

# Fields
- `D::Int`: number of spatial dimensions.
- `N::Int`: number of unit cells.
- `L::Vector{Int}`: linear extent of lattice in the direction of each lattice vector.
- `periodic::Vector{Bool}`: whether the lattice is periodic in the direction of each lattice vector.
"""
struct Lattice

    "number of spatial dimensions."
    D::Int

    "number of unit cells."
    N::Int

    "linear extent of lattice in the direction of each lattice vector."
    L::Vector{Int}

    "whether the lattice is periodic in the direction of each lattice vector."
    periodic::Vector{Bool}
end

"""
    Lattice(L::Vector{Int},periodic::Vector{Bool})

Constructs a `Lattice`.
"""
function Lattice(L::Vector{Int},periodic::Vector{Bool})

    @assert 1 <= length(L) <= 3
    @assert length(L) == length(periodic)
    @assert all(l -> l > 0, L)
     
    # dimension of lattice
    D = length(L)
    # number of unit cells in lattice
    N = prod(L)

    return Lattice(D,N,L,periodic)
end