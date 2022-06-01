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
    orbitals::Tuple{Int,Int}

    "Displacement in unit cells."
    displacement::Vector{Int}
end

"""
    Bond(orbitals::Tuple{Int,Int}, displacement::AbstractVector{Int})

Constrcut a [`Bond`](@ref)
"""
function Bond(orbitals::Tuple{Int,Int}, displacement::AbstractVector{Int})

    @assert all(i -> i > 0, orbitals)
    D = length(displacement)
    return  Bond(D, orbitals, displacement)
end

"""
    Bond(orbitals::AbstractVector{Int}, displacement::AbstractVector{Int})

Constrcut a [`Bond`](@ref)
"""
function Bond(orbitals::AbstractVector{Int}, displacement::AbstractVector{Int})

    @assert length(orbitals) == 2
    return  Bond(Tuple(i for i in orbitals), displacement)
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


"""
    Base.:==(b₁::Bond,b₂::Bond)

Tests if two bonds `b₁` and `b₂` are equivalent.
"""
function Base.:(==)(b₁::Bond,b₂::Bond)

    D₁  = b₁.D
    Δl₁ = b₁.displacement
    o₁  = b₁.orbitals
    D₂  = b₂.D
    Δl₂ = b₂.displacement
    o₂  = b₂.orbitals

    equal = false
    if D₁==D₂
        if Δl₁==Δl₂ && o₁==o₂
            equal = true
        elseif all(Δl₁[i]==-Δl₂[i] for i in 1:D₁) && o₁[1]==o₂[2] && o₁[2]==o₂[1]
            equal = true
        end
    end

    return equal
end