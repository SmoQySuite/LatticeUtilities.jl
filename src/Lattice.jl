"""
    Lattice

A type defining a finite lattice in arbitary dimensions.

# Fields
$(TYPEDFIELDS)
"""
struct Lattice

    "Number of spatial dimensions."
    D::Int

    "Number of unit cells."
    N::Int

    "Linear extent of lattice in the direction of each lattice vector."
    L::Vector{Int}

    "Whether the lattice is periodic in the direction of each lattice vector."
    periodic::Vector{Bool}

    "Storage space for representing a location or displacement in the lattice."
    lvec::Vector{Int}
end

"""
    Lattice(L::Vector{Int}, periodic::Vector{Bool})

Constructs a [`Lattice`](@ref).
"""
function Lattice(L::Vector{Int},periodic::Vector{Bool})

    @assert 1 <= length(L) <= 3
    @assert length(L) == length(periodic)
    @assert all(l -> l > 0, L)
     
    # dimension of lattice
    D = length(L)
    # number of unit cells in lattice
    N = prod(L)
    # location/displacment array
    lvec = zeros(Int,D)

    return Lattice(D,N,L,periodic,lvec)
end

Lattice(; L, periodic) = Lattice(L,periodic)


"""
    Base.show(io::IO, lattice::Lattice)
    Base.show(io::IO, ::MIME"text/plain", lattice::Lattice)

Show lattice.
"""
Base.show(io::IO, lattice::Lattice) = print(io,"Lattice(D=$(lattice.D), N=$(lattice.N), ",
                                            lattice.L, ", ", lattice.periodic)
function Base.show(io::IO, ::MIME"text/plain", lattice::Lattice)

    (; D, N, L, periodic) = lattice
    println(io, "Lattice:")
    println(io, " • D = $D")
    println(io, " • N = $N")
    println(io, " • L = ", L)
    println(io, " • periodic = ", periodic)
    return nothing
end


"""
    valid_location(l::AbstractVector{Int}, lattice::Lattice)

Determine if `l` describes a valid location in the lattice.
"""
function valid_loc(l::AbstractVector{Int}, lattice::Lattice)

    (; D, N, L) = lattice

    isvalid = true
    # check if location valid in each direction
    for d in 1:D
        if !(0<=l[d]<L[d])
            isvalid = false
            break
        end
    end

    return isvalid
end

valid_loc(; l, lattice) = valid_loc(l, lattice)


"""
    pbc!(l::AbstractVector{Int}, lattice::Lattice)

Apply periodic boundary to unit cell location `l`.
"""
function pbc!(l::AbstractVector{Int}, lattice::Lattice)

    (; D, N, L, periodic) = lattice
    @assert length(l) == D
    for d in 1:D
        # check if given direction in lattice is periodic
        if periodic[d]
            # make sure each value is positive
            if l[d] < 0
                l[d] = l[d] + L[d] * (abs(l[d])÷L[d] + 1)
            end
            # apply periodic boundary conditions
            l[d] = mod(l[d], L[d])
        end
    end
    return nothing
end

pbc!(; l, lattice) = pbc!(l, lattice)


"""
    unitcell_to_loc!(l::AbstractVector{Int}, u::Int, lattice::Lattice)

Calculate the location `l` of a unit cell `u`.
"""
function unitcell_to_loc!(l::AbstractVector{Int}, u::Int, lattice::Lattice)

    (; D, N, L) = lattice
    @assert length(l) == D

    for d in D:-1:1
        N    = N ÷ L[d]
        l[d] = (u-1) ÷ N
        u    = mod1(u,N)
    end

    return nothing
end

unitcell_to_loc!(; l, u, lattice) = unitcell_to_loc!(l, u, lattice)

"""
    unitcell_to_loc(u::Int, lattice::Lattice)

Return the location of unit cell `u`.
"""
function unitcell_to_loc(u::Int, lattice::Lattice)

    l = zeros(Int,lattice.D)
    unitcell_to_loc!(l,u,lattice)
    return l
end

unitcell_to_loc(; u, lattice) = unitcell_to_loc(u, lattice)


"""
    loca_to_unitcell(l::AbstractVector{Int}, lattice::Lattice)

Return the unit cell found at location `l` in the lattice.
"""
function loc_to_unitcell(l::AbstractVector{Int}, lattice::Lattice)

    (; D, N, L) = lattice
    @assert length(l) == D
    @assert valid_loc(l, lattice)

    u = 1
    for d in D:-1:1
        N = N ÷ L[d]
        u = u + N * l[d]
    end

    return u
end

loc_to_unitcell(; l, lattice) = loc_to_unitcell(l, lattice)