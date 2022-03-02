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
    orbitals::Vector{Int}

    "Displacement in unit cells."
    displacement::Vector{Int}
end

"""
    Bond(orbitals::AbstractVector{Int}, displacement::AbstractVector{Int})

Constrcut a [`Bond`](@ref)
"""
function Bond(orbitals::AbstractVector{Int},displacement::AbstractVector{Int})

    D = length(displacement)
    return  Bond(D,orbitals,displacement)
end

Bond(; orbitals, displacement) = Bond(orbitals,displacement)


"""
    Base.show(io::IO, bond::Bond)
    Base.show(io::IO, ::MIME"text/plain", bond::Bond)

Show lattice.
"""
Base.show(io::IO, bond::Bond) = print(io,"Bond(", (bond.D),", ",bond.o,", ",bond.Δl)
function Base.show(io::IO, ::MIME"text/plain", bond::Bond)

    (; D, orbitals, displacement) = bond
    print(io, "Bond:\n")
    print(io, " • D = ", D, "\n")
    print(io, " • orbitals = ", orbitals, "\n")
    print(io, " • displacement = ", displacement)
    return nothing
end