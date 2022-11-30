"""
    Bond{D}

Defines a bond in a `D` dimensional lattice.

# Fields

- `orbitals::NTuple{2,Int}`: Orbital pair associated with bond. Bond goes from `orbitals[1]` to `orbitals[2]`.
- `displacement::SVector{D,Int}`: Displacement in unit cells in the direction of each lattice vector associated with bond.
"""
struct Bond{D}
    
    # initial and final orbitals
    orbitals::NTuple{2,Int}

    # displacement in Unit Cells
    displacement::SVector{D,Int}
end

"""
    Bond(orbitals, displacement)

    Bond(; orbitals, displacement)

Construct a [`Bond`](@ref)
"""
Bond(orbitals, displacement) = Bond(SVector{2,Int}(orbitals), SVector{length(displacement),Int}(displacement))
Bond(; orbitals, displacement) = Bond(orbitals, displacement)


"""
    Base.show(io::IO, bond::Bond{D}) where {D}
    
    Base.show(io::IO, ::MIME"text/plain", bond::Bond{D}) where {D}

Show lattice.
"""
Base.show(io::IO, bond::Bond{D}) where {D} = print(io, "Bond{$(bond.D)}(orbitals=$(bond.orbitals), displacement=$(bond.displacement))")
function Base.show(io::IO, ::MIME"text/plain", bond::Bond{D}) where {D}

    (; orbitals, displacement) = bond
    @printf io "[Bond]\n\n"
    @printf io "dimensions   = %d\n" D
    @printf io "orbitals     = %s\n" string(orbitals)
    @printf io "displacement = %s\n" string(displacement)

    return nothing
end


"""
    Base.:(==)(b1::Bond{D1}, b2::Bond{D2}) where {D1, D2}

Tests if two bonds `b1` and `b2` are equivalent.
"""
function Base.:(==)(b1::Bond{D1}, b2::Bond{D2}) where {D1, D2}

    return (D1==D2) && (b1.orbitals==b2.orbitals) && (b1.displacement==b2.displacement)
end