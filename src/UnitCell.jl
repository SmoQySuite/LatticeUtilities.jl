"""
    UnitCell{T<:AbstractFloat}

A type defining a unit cell.

# Fields
- `D::Int`: number of spatial dimensions.
- `n::Int`: sites per unit cell.
- `lattice_vecs::Matrix{T}`: matrix where the columns are the lattice vectors.
- `reciprocal_vecs::Matrix{T}`: matrix where the columns are the reciprocal lattice vectors.
- `basis_vecs::Matrix{T}`: matrix where the columns are the basis vectors.
"""
struct UnitCell{T<:AbstractFloat}
    
    "number of spatial dimensions."
    D::Int

    "sites per unit cell."
    n::Int

    "matrix where the columns are the lattice vectors."
    lattice_vecs::Matrix{T}

    "matrix where the columns are the reciprocal lattice vectors."
    reciprocal_vecs::Matrix{T}

    "matrix where the columns are then basis vectors."
    basis_vecs::Matrix{T}
end

"""
    UnitCell(lattice_vecs::Matrix{T}, basis_vecs::Matrix{T}) where {T<:AbstractFloat}

Constrcuts a `UnitCell`.

# Example
Define the unit cell for a Kagome lattice.
```jldoctest; output = false
lattice_vecs = [1.0  1/2;
                0.0 √3/2]
basis_vecs   = [0.0 1/2  1/4;
                0.0 0.0 √3/4]
unit_cell    = UnitCell(lattice_vecs, basis_vecs)

# output

UnitCell{Float64}:
- D = 2
- n = 3
- lattice_vecs =
2×2 Matrix{Float64}:
 6.28319  -3.6276
 0.0       7.2552
- reciprocal_vecs =
2×2 Matrix{Float64}:
 1.0  0.5
 0.0  0.866025
- basis_vecs =
2×3 Matrix{Float64}:
 0.0  0.5  0.25
 0.0  0.0  0.433013
```
"""
function UnitCell(lattice_vecs::Matrix{T}, basis_vecs::Matrix{T}) where {T<:AbstractFloat}

    @assert 1 <= size(lattice_vecs,1) <= 3
    @assert size(basis_vecs,2) >= 1
    @assert size(lattice_vecs,1)==size(basis_vecs,1)
    @assert size(lattice_vecs,1)==size(lattice_vecs,2)

    # dimension of unit cell
    D = size(lattice_vecs,1)

    # number of sites in unit cell
    n = size(basis_vecs,2)

    # calculating reciprocal lattice vectors corresponding to lattice vectors
    reciprocal_vecs = 2π*inv(lattice_vecs)

    return UnitCell{T}(D, n, lattice_vecs, reciprocal_vecs, basis_vecs)
end

"""
    UnitCell(lattice_vecs::AbstractVector{Vector{T}}, basis_vecs::AbstractVector{Vector{T}})
        where {T<:AbstractFloat}

Constrcuts a `UnitCell`.

# Example
Define the unit cell for a Kagome lattice.
```jldoctest; output = false
lattice_vecs = [ [1.0, 0.0], [1/2, √3/2] ]
basis_vecs   = [ [0.0, 0.0], [1/2, 0.0], [1/4, √3/4] ]
unit_cell    = UnitCell(lattice_vecs, basis_vecs)

# output

UnitCell{Float64}:
- D = 2
- n = 3
- lattice_vecs =
2×2 Matrix{Float64}:
 6.28319  -3.6276
 0.0       7.2552
- reciprocal_vecs =
2×2 Matrix{Float64}:
 1.0  0.5
 0.0  0.866025
- basis_vecs =
2×3 Matrix{Float64}:
 0.0  0.5  0.25
 0.0  0.0  0.433013
```
"""
function UnitCell(lattice_vecs::AbstractVector{Vector{T}}, basis_vecs::AbstractVector{Vector{T}}) where {T<:AbstractFloat}

    return UnitCell(hcat(lattice_vecs...), hcat(basis_vecs...))
end

"""
    Base.show(io::IO, uc::UnitCell{T}) where {T}
    Base.show(io::IO, ::MIME"text/plain", uc::UnitCell{T}) where{T}

Show unit cell.
"""
Base.show(io::IO, uc::UnitCell{T}) where {T} = print(io,"UnitCell{$T}(D=$(uc.D), n=$(uc.n))")
function Base.show(io::IO, ::MIME"text/plain", uc::UnitCell{T}) where{T}

    (; D, n, lattice_vecs, reciprocal_vecs, basis_vecs) = uc
    print(io, "UnitCell{$T}:\n")
    print(io,"- D = $D\n")
    print(io,"- n = $n\n")
    print(io,"- lattice_vecs =\n")
    show(io,"text/plain",reciprocal_vecs)
    print(io,"\n")
    print(io,"- reciprocal_vecs =\n")
    show(io,"text/plain",lattice_vecs)
    print(io,"\n")
    print(io,"- basis_vecs =\n")
    show(io,"text/plain",basis_vecs)
    return nothing
end

"""
    get_r!(r::AbstractVector{T}, l::AbstractVector{Int}, unit_cell::UnitCell{T})
        where {T}

Calculate the position `r` of a unit cell at location `l`.
"""
function get_r!(r::AbstractVector{T}, l::AbstractVector{Int}, unit_cell::UnitCell{T}) where {T}

    (; D, lattice_vecs) = unit_cell
    @assert length(r) == length(l) == D

    fill!(r,0.0)
    for d in 1:D
        @views @. r += l[d] * lattice_vecs[:,d]
    end

    return nothing
end

"""
    get_r!(r::AbstractVector{T}, l::AbstractVector{Int}, s::Int, unit_cell::UnitCell{T})
        where {T}

Calculate the position `r` of a site `s` at location `l`.
"""
function get_r!(r::AbstractVector{T}, l::AbstractVector{Int}, s::Int, unit_cell::UnitCell{T}) where {T}

    get_r!(r,l,unit_cell)
    @views @. r += unit_cell.basis_vecs[:,s]

    return nothing
end

"""
    get_r(l::AbstractVector{Int}, unit_cell::UnitCell{T})
        where {T}

Return the position `r` of a unit cell at location `l`.
"""
function get_r(l::AbstractVector{Int}, unit_cell::UnitCell{T}) where {T}

    r = zeros(T,unit_cell.D)
    get_r!(r,l,unit_cell)

    return r
end

"""
    get_r(l::AbstractVector{Int}, s::Int, unit_cell::UnitCell{T})
        where {T}

Return the position `r` of a site `s` at location `l`.
"""
function get_r(l::AbstractVector{Int}, s::Int, unit_cell::UnitCell{T}) where {T}

    r = zeros(T,unit_cell.D)
    get_r!(r,l,s,unit_cell)

    return r
end


"""
    Δl_to_Δr!(Δr::AbstractVector{T}, Δl::AbstractVector{Int}, s₁::Int, s₂::Int,
        unit_cell::UnitCell{T}) where {T}

Computes the position space displacement vector `Δr` corresponding to a displacement definition given by initial and final
sites `s₁` and `s₂` in the unit cell respectively, along with a displacement in unit cells `Δl`.

# Arguments
- `Δr::AbstractVector{T}`: displacement vector in position space.
- `Δl::AbstractVector{Int}`: displacement in unit cells.
- `s₁::Int`: initial site in unit cell.
- `s₂::Int`: final site in unit cell.
- `unit_cell::UnitCell{T}`: unit cell.
"""
function Δl_to_Δr!(Δr::AbstractVector{T}, Δl::AbstractVector{Int}, s₁::Int, s₂::Int, unit_cell::UnitCell{T}) where {T}
    
    (; D, n, lattice_vecs, basis_vecs) = unit_cell
    @assert length(Δr) == length(Δx) == D
    @assert 1 <= s₁ <= n
    @assert 1 <= s₂ <= n

    fill!(Δr,0.0)
    for d in in 1:D
        @views @. Δr += Δl[d] * lattice_vecs[:,d]
    end
    @views @. Δr += basis_vecs[:,s₂] - basis_vecs[:,s₁]

    return nothing
end

"""
    Δl_to_Δr(Δl::AbstractVector{Int}, s₁::Int, s₂::Int, unit_cell::UnitCell{T})
        where {T}

Returns the position space displacement vector `Δr` corresponding to a displacement definition given by initial and final
sites `s₁` and `s₂` in the unit cell respectively, along with a displacement in unit cells `Δl`.
"""
function Δl_to_Δr(Δl::AbstractVector{Int}, s₁::Int, s₂::Int, unit_cell::UnitCell{T}) where {T}
    
    Δr = zeros(T,unit_cell.D)
    Δl_to_Δr!(Δr, Δl, s₁, s₂, unit_cell)

    return Δr
end