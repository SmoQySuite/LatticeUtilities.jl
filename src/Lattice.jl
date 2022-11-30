"""
    Lattice{D}

A type defining a finite lattice in `D` dimensions.

# Fields

- `N::Int`: Number of unit cells in finite lattice.
- `L::SVector{D,Int}`: Size of finite lattice in the direction of each lattice vector in unit cells.
- `periodic::SVector{D,Bool}`: Specifies whether each lattice vector direction hosts periodic or open boundary conditions.
- `lvec::MVector{D,Int}`: Private temporary storage vector to contain intermediate location and displacement vectors.
"""
struct Lattice{D}

    N::Int
    L::SVector{D,Int}
    periodic::SVector{D,Bool}
    lvec::MVector{D,Int}
end

"""
    Lattice(L, periodic)

    Lattice(; L, periodic)

Constructs a [`Lattice`](@ref).
"""
Lattice(L, periodic) = Lattice(prod(L), SVector{length(L),Int}(L), SVector{length(periodic),Bool}(periodic), MVector{length(L),Int}(zeros(Int,length(L))))
Lattice(; L, periodic) = Lattice(L ,periodic)


"""
    Base.show(io::IO, lattice::Lattice{D}) where {D}
    
    Base.show(io::IO, ::MIME"text/plain", lattice::Lattice{D}) where {D}

Show lattice.
"""
Base.show(io::IO, lattice::Lattice{D}) where {D} = print(io,"Lattice{$(D)}(N=$(lattice.N), L=$(lattice.L), periodic=$(lattice.periodic))")
function Base.show(io::IO, ::MIME"text/plain", lattice::Lattice{D}) where {D}

    (; N, L, periodic) = lattice
    @printf io "[Lattice]\n\n"
    @printf io "dimensions   = %d\n" D
    @printf io "n_unit_cells = %d\n" N
    @printf io "size         = %s\n" string(L)
    @printf io "periodic     = [%s]\n" join(periodic, ", ")
    return nothing
end


"""
    valid_loc(l, lattice::Lattice{D})::Bool where {D}

Determine if `l` describes a valid location in the lattice.
"""
function valid_loc(l, lattice::Lattice{D})::Bool where {D}

    return all(i -> 0<=l[i]<lattice.L[i], 1:D) && length(l)==D
end

valid_loc(; l, lattice) = valid_loc(l, lattice)


"""
    pbc!(l::AbstractVector{Int}, lattice::Lattice{D}) where {D}

Apply periodic boundary to unit cell location `l`.
"""
function pbc!(l::AbstractVector{Int}, lattice::Lattice{D}) where {D}

    (; L, periodic) = lattice
    @fastmath @inbounds for d in eachindex(l)
        # check if given direction in lattice is periodic
        if periodic[d]
            # apply periodic boundary conditions
            l[d] = mod(l[d], L[d])
        end
    end
    return nothing
end

pbc!(; l, lattice) = pbc!(l, lattice)


"""
    unitcell_to_loc!(l::AbstractVector{Int}, u::Int, lattice::Lattice{D}) where {D}

Calculate the location `l` of a unit cell `u`.
"""
function unitcell_to_loc!(l::AbstractVector{Int}, u::Int, lattice::Lattice{D}) where {D}

    (; N, L) = lattice

    @fastmath @inbounds for d in D:-1:1
        N    = N ÷ L[d]
        l[d] = (u-1) ÷ N
        u    = mod1(u,N)
    end

    return nothing
end

unitcell_to_loc!(; l, u, lattice) = unitcell_to_loc!(l, u, lattice)

"""
    unitcell_to_loc(u::Int, lattice::Lattice{D})::SVector{D,Int} where {D}

Return the location of unit cell `u` as an instance of type `SVector{D,Int}`.
"""
function unitcell_to_loc(u::Int, lattice::Lattice{D})::SVector{D,Int} where {D}

    l = zeros(Int,lattice.D)
    unitcell_to_loc!(l,u,lattice)
    return SVector{D,Int}(l)
end

unitcell_to_loc(; u, lattice) = unitcell_to_loc(u, lattice)


"""
    loc_to_unitcell(l, lattice::Lattice{D}) where {D}

Return the unit cell found at location `l` in the lattice.
"""
function loc_to_unitcell(l, lattice::Lattice{D}) where {D}

    (; N, L) = lattice

    u = 1
    @fastmath @inbounds for d in D:-1:1
        N = N ÷ L[d]
        u = u + N * l[d]
    end

    return u
end

loc_to_unitcell(; l, lattice) = loc_to_unitcell(l, lattice)