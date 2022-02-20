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


"""
    Base.show(io::IO, bond::Bond)
    Base.show(io::IO, ::MIME"text/plain", bond::Bond)

Show lattice.
"""
Base.show(io::IO, bond::Bond) = print(io,"Bond(D=$(bond.D), ",bond.o,", ",bond.Δl)
function Base.show(io::IO, ::MIME"text/plain", bond::Bond)

    (; D, o, Δl) = bond
    println(io, "Bond:")
    println(io, " - D  = ", D)
    println(io, " - o  = ", o)
    println(io, " - Δl = ", Δl)
    return nothing
end