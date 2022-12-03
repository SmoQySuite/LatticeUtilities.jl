@doc raw"""
    UnitCell{D,T<:AbstractFloat,N}

A type defining a unit cell in `D` spatial dimensions.

# Fields

- `n::Int`: Number of orbitals in the unit cell.
- `lattice_vecs::SMatrix{D,D,T,N}`: Matrix where columns give the lattice vectors.
- `reciprocal_vecs::SMatrix{D,D,T,N}`: Matrix where columns give the reciprocal lattice vectors.
- `basis_vecs::Vector{SVector{D,T}}`: Vector of basis vectors specifying the location of each orbital in unit cell.
"""
struct UnitCell{D,T<:AbstractFloat,N}
    
    # number of orbitals per unit cell
    n::Int

    # columns of matrix are lattice vectors
    lattice_vecs::SMatrix{D,D,T,N}

    # columns of matrix are reciprocal lattice vectors
    reciprocal_vecs::SMatrix{D,D,T,N}

    # list (vector) of basis vectors
    basis_vecs::Vector{SVector{D,T}}
end


@doc raw"""
    UnitCell(lattice_vecs::AbstractMatrix{T}, basis_vecs::AbstractMatrix{T}) where {T<:AbstractFloat}

    UnitCell(lattice_vecs::AbstractVector{Vector{T}}, basis_vecs::AbstractVector{Vector{T}}) where {T<:AbstractFloat}

    UnitCell(; lattice_vecs, basis_vecs)

Returns an instance of the type [`UnitCell`](@ref).
"""
function UnitCell(lattice_vecs::AbstractMatrix{T}, basis_vecs::AbstractMatrix{T}) where {T<:AbstractFloat}

    D = size(lattice_vecs, 1)
    n = size(basis_vecs, 2)
    reciprocal_vecs = 2π * inv(lattice_vecs)'

    static_basis_vecs = SVector{D,T}[]
    for i in 1:n
        push!(static_basis_vecs, SVector{D}(@view basis_vecs[:,i]))
    end

    return UnitCell(n, SMatrix{D,D}(lattice_vecs), SMatrix{D,D}(reciprocal_vecs), static_basis_vecs)
end

function UnitCell(lattice_vecs::AbstractVector{Vector{T}}, basis_vecs::AbstractVector{Vector{T}}) where {T<:AbstractFloat}

    return UnitCell(hcat(lattice_vecs...), hcat(basis_vecs...))
end

UnitCell(; lattice_vecs, basis_vecs) = UnitCell(lattice_vecs, basis_vecs)


"""
    Base.show(io::IO, uc::UnitCell{D,T}) where {D,T}
    
    Base.show(io::IO, ::MIME"text/plain", uc::UnitCell{D,T}) where {D,T}

Show unit cell.
"""
Base.show(io::IO, uc::UnitCell{D,T}) where {D,T} = print(io,"UnitCell{$D,$T}(n=$(uc.n))")
function Base.show(io::IO, ::MIME"text/plain", uc::UnitCell{D,T}) where {D,T}

    (; n, lattice_vecs, reciprocal_vecs, basis_vecs) = uc
    @printf io "[UnitCell]\n\n"
    @printf io "dimensions = %d\n" D
    @printf io "orbitals = %d\n\n" n
    @printf io "[UnitCell.lattice_vecs]\n\n"
    for d in axes(lattice_vecs,2)
        a = @view lattice_vecs[:,d]
        @printf io "a_%d = %s\n" d string(a)
    end
    @printf io "\n"
    @printf io "[UnitCell.reciprocal_vecs]\n\n"
    for d in axes(reciprocal_vecs,2)
        a = @view reciprocal_vecs[:,d]
        @printf io "b_%d = %s\n" d string(a)
    end
    @printf io "\n"
    @printf io "[UnitCell.basis_vecs]\n"
    for d in eachindex(basis_vecs)
        @printf io "\nb_%d = %s" d string(basis_vecs[d])
    end
    return nothing
end


"""
    loc_to_pos!(r::AbstractVector{T}, l, unit_cell::UnitCell{D,T}) where {D,T}

Calculate the position `r` of a unit cell at location `l`.
"""
function loc_to_pos!(r::AbstractVector{T}, l, unit_cell::UnitCell{D,T}) where {D,T}

    (; lattice_vecs) = unit_cell

    fill!(r,0.0)
    @fastmath @inbounds for d in eachindex(r)
        a_d = @view lattice_vecs[:,d]
        @. r += l[d] * a_d
    end

    return nothing
end

"""
    loc_to_pos!(r::AbstractVector{T}, l, o::Int, unit_cell::UnitCell{D,T}) where {D,T}

Calculate the position `r` of an orbital `o` at location `l`.
"""
function loc_to_pos!(r::AbstractVector{T}, l, o::Int, unit_cell::UnitCell{D,T}) where {D,T}

    (; basis_vecs) = unit_cell
    loc_to_pos!(r, l, unit_cell)
    @. r += basis_vecs[o]

    return nothing
end

"""
    loc_to_pos(l, unit_cell::UnitCell{T}) where {T}

Return the position `r` of a unit cell at location `l`.
"""
function loc_to_pos(l, unit_cell::UnitCell{D,T}) where {D,T}

    r = zeros(T,D)
    loc_to_pos!(r, l, unit_cell)

    return SVector{D,T}(r)
end

"""
    loc_to_pos(l, s::Int, unit_cell::UnitCell{D,T}) where {D,T}

Return the position `r` of a orbital `o` at location `l` as a vector or type `SVector{D,T}`.
"""
function loc_to_pos(l, o::Int, unit_cell::UnitCell{D,T}) where {D,T}

    r = zeros(T,D)
    loc_to_pos!(r,l,o,unit_cell)

    return SVector{D,T}(r)
end


"""
    displacement_to_vec!(Δr::AbstractVector{T}, Δl, o_init::Int, o_final::Int, unit_cell::UnitCell{D,T}) where {D,T}

Computes the position space displacement vector `Δr` corresponding to a displacement definition given by initial and final
orbitals `o₁` and `o₂` in the unit cell respectively, along with a displacement in unit cells `Δl`.

# Arguments
- `Δr::AbstractVector{T}`: displacement vector in position space.
- `Δl`: displacement in unit cells.
- `o_init::Int`: initial orbital in unit cell.
- `o_final::Int`: final orbital in unit cell.
- `unit_cell::UnitCell{D,T}`: unit cell.
"""
function displacement_to_vec!(Δr::AbstractVector{T}, Δl, o_init::Int, o_final::Int, unit_cell::UnitCell{D,T}) where {D,T}

    (; lattice_vecs, basis_vecs) = unit_cell

    fill!(Δr,0.0)
    @fastmath @inbounds for d in eachindex(Δr)
        a_d = @view lattice_vecs[:,d]
        @. Δr += Δl[d] * a_d
    end
    @. Δr += basis_vecs[o_final] - basis_vecs[o_init]

    return nothing
end

displacement_to_vec!(; Δr, Δl, o₁, o₂, unit_cell) = displacement_to_vec!(Δr, Δl, o₁, o₂, unit_cell)

"""
    displacement_to_vec(Δl, o_init::Int, o_final::Int, unit_cell::UnitCell{D,T})::SVector{D,T} where {D,T}

Returns the position space displacement vector `Δr` corresponding to a displacement definition given by initial and final
orbitals `o_init` and `o_final` in the unit cell respectively, along with a displacement in unit cells `Δl`.
"""
function displacement_to_vec(Δl, o_init::Int, o_final::Int, unit_cell::UnitCell{D,T})::SVector{D,T} where {D,T}
    
    Δr = zeros(T, D)
    displacement_to_vec!(Δr, Δl, o_init, o_final, unit_cell)

    return SVector{D,T}(Δr)
end

displacement_to_vec(; Δl, o_init, o_final, unit_cell) = displacement_to_vec(Δl, o_init, o_final, unit_cell)