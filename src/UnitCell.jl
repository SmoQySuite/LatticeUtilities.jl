"""
    UnitCell{T<:AbstractFloat}

A type defining a unit cell.

# Fields
$(TYPEDFIELDS)
"""
struct UnitCell{T<:AbstractFloat}
    
    "Number of spatial dimensions."
    D::Int

    "Orbitals per unit cell."
    n::Int

    "Matrix where the columns are the lattice vectors."
    lattice_vecs::Matrix{T}

    "Matrix where the columns are the reciprocal lattice vectors."
    reciprocal_vecs::Matrix{T}

    "Matrix where the columns are then basis vectors."
    basis_vecs::Matrix{T}
end

"""
    UnitCell(lattice_vecs::Matrix{T}, basis_vecs::Matrix{T}) where {T<:AbstractFloat}

Constrcuts a [`UnitCell`](@ref).
"""
function UnitCell(lattice_vecs::Matrix{T}, basis_vecs::Matrix{T}) where {T<:AbstractFloat}

    @assert 1 <= size(lattice_vecs,1) <= 3
    @assert size(basis_vecs,2) >= 1
    @assert size(lattice_vecs,1)==size(basis_vecs,1)
    @assert size(lattice_vecs,1)==size(lattice_vecs,2)

    # dimension of unit cell
    D = size(lattice_vecs,1)

    # number of orbitals in unit cell
    n = size(basis_vecs,2)

    # calculating reciprocal lattice vectors corresponding to lattice vectors
    reciprocal_vecs = 2π * inv(lattice_vecs)'

    return UnitCell{T}(D, n, lattice_vecs, reciprocal_vecs, basis_vecs)
end

"""
    UnitCell(lattice_vecs::AbstractVector{Vector{T}}, basis_vecs::AbstractVector{Vector{T}})
        where {T<:AbstractFloat}

Constrcuts a [`UnitCell`](@ref).
"""
function UnitCell(lattice_vecs::AbstractVector{Vector{T}}, basis_vecs::AbstractVector{Vector{T}}) where {T<:AbstractFloat}

    return UnitCell(hcat(lattice_vecs...), hcat(basis_vecs...))
end


"""
    Base.show(io::IO, uc::UnitCell{T}) where {T}
    Base.show(io::IO, ::MIME"text/plain", uc::UnitCell{T}) where {T}

Show unit cell.
"""
Base.show(io::IO, uc::UnitCell{T}) where {T} = print(io,"UnitCell{$T}(D=$(uc.D), n=$(uc.n))")
function Base.show(io::IO, ::MIME"text/plain", uc::UnitCell{T}) where {T}

    (; D, n, lattice_vecs, reciprocal_vecs, basis_vecs) = uc
    print(io, "UnitCell{$T}:\n")
    print(io," - D = $D\n")
    print(io," - n = $n\n")
    print(io," - lattice_vecs =\n")
    show(io,"text/plain",lattice_vecs)
    print(io,"\n")
    print(io," - reciprocal_vecs =\n")
    show(io,"text/plain",reciprocal_vecs)
    print(io,"\n")
    print(io," - basis_vecs =\n")
    show(io,"text/plain",basis_vecs)
    return nothing
end


"""
    loc_to_pos!(r::AbstractVector{T}, l::AbstractVector{Int}, unit_cell::UnitCell{T})
        where {T}

Calculate the position `r` of a unit cell at location `l`.
"""
function loc_to_pos!(r::AbstractVector{T}, l::AbstractVector{Int}, unit_cell::UnitCell{T}) where {T}

    (; D, lattice_vecs) = unit_cell
    @assert length(r) == length(l) == D

    fill!(r,0.0)
    for d in 1:D
        @views @. r += l[d] * lattice_vecs[:,d]
    end

    return nothing
end

"""
    loc_to_pos!(r::AbstractVector{T}, l::AbstractVector{Int}, s::Int, unit_cell::UnitCell{T})
        where {T}

Calculate the position `r` of a orbital `o` at location `l`.
"""
function loc_to_pos!(r::AbstractVector{T}, l::AbstractVector{Int}, o::Int, unit_cell::UnitCell{T}) where {T}

    loc_to_pos!(r,l,unit_cell)
    @views @. r += unit_cell.basis_vecs[:,o]

    return nothing
end

"""
    loc_to_pos(l::AbstractVector{Int}, unit_cell::UnitCell{T})
        where {T}

Return the position `r` of a unit cell at location `l`.
"""
function loc_to_pos(l::AbstractVector{Int}, unit_cell::UnitCell{T}) where {T}

    r = zeros(T,unit_cell.D)
    loc_to_pos!(r,l,unit_cell)

    return r
end

"""
    loc_to_pos(l::AbstractVector{Int}, s::Int, unit_cell::UnitCell{T})
        where {T}

Return the position `r` of a orbital `o` at location `l`.
"""
function loc_to_pos(l::AbstractVector{Int}, o::Int, unit_cell::UnitCell{T}) where {T}

    r = zeros(T,unit_cell.D)
    loc_to_pos!(r,l,o,unit_cell)

    return r
end


"""
    displacement_to_vec!(Δr::AbstractVector{T}, Δl::AbstractVector{Int},
        o₁::Int, o₂::Int, unit_cell::UnitCell{T}) where {T}

Computes the position space displacement vector `Δr` corresponding to a displacement definition given by initial and final
orbitals `o₁` and `o₂` in the unit cell respectively, along with a displacement in unit cells `Δl`.

# Arguments
- `Δr::AbstractVector{T}`: displacement vector in position space.
- `Δl::AbstractVector{Int}`: displacement in unit cells.
- `o₁::Int`: initial orbital in unit cell.
- `o₂::Int`: final orbital in unit cell.
- `unit_cell::UnitCell{T}`: unit cell.
"""
function displacement_to_vec!(Δr::AbstractVector{T}, Δl::AbstractVector{Int}, o₁::Int, o₂::Int, unit_cell::UnitCell{T}) where {T}
    
    (; D, n, lattice_vecs, basis_vecs) = unit_cell
    @assert length(Δr) == length(Δl) == D
    @assert 1 <= o₁ <= n
    @assert 1 <= o₂ <= n

    fill!(Δr,0.0)
    for d in 1:D
        @views @. Δr += Δl[d] * lattice_vecs[:,d]
    end
    @views @. Δr += basis_vecs[:,o₂] - basis_vecs[:,o₁]

    return nothing
end

"""
    displacement_to_vec(Δl::AbstractVector{Int}, o₁::Int, o₂::Int,
        unit_cell::UnitCell{T}) where {T}

Returns the position space displacement vector `Δr` corresponding to a displacement definition given by initial and final
orbitals `o₁` and `o₂` in the unit cell respectively, along with a displacement in unit cells `Δl`.
"""
function displacement_to_vec(Δl::AbstractVector{Int}, o₁::Int, o₂::Int, unit_cell::UnitCell{T}) where {T}
    
    Δr = zeros(T,unit_cell.D)
    displacement_to_vec!(Δr, Δl, o₁, o₂, unit_cell)

    return Δr
end