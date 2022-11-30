@doc raw"""
    Base.ndims(unit_cell::UnitCell{D}) where {D}

    Base.ndims(lattice::Lattice{D}) where {D}

    Base.ndims(bond::Bond{D}) where {D}

Return the spatial dimensions of a [`UnitCell`](@ref), [`Lattice`](@ref) or [`Bond`](@ref).
"""
Base.ndims(unit_cell::UnitCell{D}) where {D} = D
Base.ndims(lattice::Lattice{D}) where {D} = D
Base.ndims(bond::Bond{D}) where {D} = D

@doc raw"""
    norbits(unit_cell::UnitCell)

Return the number of orbitals in a unit cell.
"""
norbits(unit_cell::UnitCell) = unit_cell.n

@doc raw"""
    nsites(unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

    nsites(; unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

Return the total number of sites/orbitals in a finite lattice.
"""
nsites(unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D} = unit_cell.n * lattice.N
nsites(; unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D} = nsites(unit_cell, lattice)


"""
    valid_site(s::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

    valid_site(; s, unit_cell, lattice) = valid_site(s, unit_cell, lattice)

Return whether `s` is a valid site index.
"""
function valid_site(s::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

    N = nsites(unit_cell, lattice)
    return 0 < s <= N
end

valid_site(; s, unit_cell, lattice) = valid_site(s, unit_cell, lattice)

"""
    site_to_unitcell(s::Int, unit_cell::UnitCell)

Return the unit cell `u` containing lattice site `s`.
"""
function site_to_unitcell(s::Int, unit_cell::UnitCell)

    return (s-1) ÷ unit_cell.n + 1 
end

site_to_unitcell(; s, unit_cell) = site_to_unitcell(s, unit_cell)


"""
    site_to_orbital(s::Int, unit_cell::UnitCell)

Return the orbtial species of site `s`.
"""
function site_to_orbital(s::Int, unit_cell::UnitCell)

    return mod1(s, unit_cell.n)
end


"""
    site_to_loc!(l::AbstractVector{Int}, s::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

    site_to_loc!(; l::AbstractVector{Int}, s::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

For a given site `s` in the lattice, calculate the location `l` of the unit cell it is in
and return the orbital species `o` of the site.
"""
function site_to_loc!(l::AbstractVector{Int}, s::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

    o = site_to_orbital(s, unit_cell)
    u = site_to_unitcell(s, unit_cell)
    unitcell_to_loc!(l, u, lattice)
    return o
end

site_to_loc!(; l::AbstractVector{Int}, s::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D} = site_to_loc!(l, s, unit_cell, lattice)

"""
    site_to_loc(s::Int, unit_cell::UnitCell{D,T}, lattice::Lattice{D}) where {D,T}

For a given site `s` in the lattice, return the location `l` of the unit cell it is in
and the orbital species `o`.
"""
function site_to_loc(s::Int, unit_cell::UnitCell{D,T}, lattice::Lattice{D}) where {D,T}

    @assert unit_cell.D == lattice.D
    l = zeros(Int,unit_cell.D)
    o = site_to_loc!(l, s, unit_cell, lattice)
    return (SVector{D,T}(l), o)
end

site_to_loc(; s, unit_cell, lattice) = site_to_loc(s, unit_cell, lattice)


"""
    loc_to_site(l, o::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

Given a unit cell location `l` and orbital species `o`, return the corresponding
site `s` in the lattice. If the location is not valid owing to open boundary conditions
then return `s = 0`.
"""
function loc_to_site(l, o::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}
    
    if valid_loc(l, lattice)
        u = loc_to_unitcell(l, lattice)
        s = loc_to_site(u, o, unit_cell)
    else
        s = 0
    end
    
    return s
end

loc_to_site(; l, o, unit_cell, lattice) = loc_to_site(l, o, unit_cell, lattice)

"""
    loc_to_site(u::Int, o::Int, unit_cell::UnitCell{D}) where {D}
    loc_to_site(; u::Int, o::Int, unit_cell::UnitCell{D}) where {D}

Given a unit cell index `u` and orbital `o`, return the correspond site `s`.
"""
function loc_to_site(u::Int, o::Int, unit_cell::UnitCell{D}) where {D}

    (; n) = unit_cell
    s = n * (u-1) + o
    return s
end

loc_to_site(; u::Int, o::Int, unit_cell::UnitCell{D}) where {D} = loc_to_site(u, o, unit_cell, lattice)


"""
    site_to_site(s::Int, Δl, o::Int, unit_cell::UnitCell, lattice::Lattice)

Given an initial site `s`, and a displacement in unit cells `Δl` and a terminating orbital
species `o`, return the resulting site `s′` in the lattice. If the displacement is
not allowed as a result of open boundary conditions, then  `s′=0` is returned.
"""
function site_to_site(s::Int, Δl, o::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

    l = lattice.lvec

    # get unit cell location containing s₁
    o′ = site_to_loc!(l, s, unit_cell, lattice)

    # displace unit cell location
    @. l += Δl

    # apply periodic boundary conditions
    pbc!(l, lattice)

    # get final site
    s′ = loc_to_site(l, o, unit_cell, lattice)

    return s′
end

site_to_site(; s₁, Δl, o₂, unit_cell, lattice) = site_to_site(s₁, Δl, o₂, unit_cell, lattice)


"""
    simplify!(Δl::AbstractVector{Int}, lattice::Lattice{D}) where {D}

Simplify displacement `Δl` so that it is as short as possible accounting
for periodic boundary conditions where necessary.
"""
function simplify!(Δl::AbstractVector{Int}, lattice::Lattice{D}) where {D}

    (; L, periodic) = lattice

    # simplify to shortest displacement accounting
    # for periodic boundary conditions
    @fastmath @inbounds for d in eachindex(Δl)
        if periodic[d] && abs(Δl[d]) > (L[d]÷2)
            Δl[d] -= sign(Δl[d]) * L[d]
        end
    end

    return nothing
end

simplify!(; Δl, lattice) = simplify!(Δl, lattice)

"""
    simplify(bond::Bond{D}, lattice::Lattice{D}) where {D}

Simplify a bond so that the displacement is the shortest possible accounting
for periodic boundary conditions where necessary, returning the new bond.
"""
function simplify(bond::Bond{D}, lattice::Lattice{D}) where {D}

    Δl = zeros(Int, D)
    copyto!(Δl, bond.displacement)
    simplify!(Δl, lattice)

    return Bond(bond.orbitals, Δl)
end

simplify(; bond, lattice) = simplify!(bond, lattice)


"""
    function sites_to_displacement!(Δl::AbstractVector{Int}, s::Int, s′::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

When getting displaced from site `s` to `s′`, calculate the corresponding displacement in unit cells `Δl`.
"""
function sites_to_displacement!(Δl::AbstractVector{Int}, s::Int, s′::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

    # get initial and final orbital
    o[1] = site_to_orbital(s₁, unit_cell)
    o[2] = site_to_orbital(s₂, unit_cell)

    # calculate displacement in unit cells
    l = lattice.lvec
    l′ = Δl
    site_to_loc!(l, s, unit_cell, lattice)
    site_to_loc!(l′, s′, unit_cell, lattice)
    @. Δl = l′ - l

    # simplify the displacement
    simplify!(Δl,lattice)

    return nothing
end

sites_to_displacement!(; Δl, s, s′, unit_cell, lattice) = sites_to_displacement!(Δl, s, s′, unit_cell, lattice)


"""
    function sites_to_displacement(s::Int, s′::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

When getting displaced from site `s` to `s′`, return the corresponding displacement in unit cells `Δl::MVector{D,Int}`.
"""
function sites_to_displacement(s::Int, s′::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

    Δl = MVector{D,Int}(undef)
    sites_to_displacement!(Δl, s, s′, unit_cell, lattice)

    return Δl
end

sites_to_displacement(; s, s′, unit_cell, lattice) = sites_to_displacement(s, s′, unit_cell, lattice)


"""
    sites_to_bond(s::Int, s′::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

Return the [`Bond`](@ref) associated with getting displaced from site `s` to `s′`. 
"""
function sites_to_bond(s::Int, s′::Int, unit_cell::UnitCell{D}, lattice::Lattice{D}) where {D}

    Δl = sites_to_displacement(s, s′, unit_cell, lattice)
    o  = site_to_orbital(s, unit_cell)
    o′ = site_to_orbital(s′, unit_cell)

    return Bond((o,o′), Δl)
end

sites_to_bond(; s, s′, unit_cell, lattice) = sites_to_bond(s, s′, unit_cell, lattice)


"""
    calc_k_point!(k_point::AbstractVector{T}, k_loc::AbstractVector{Int},
                  unit_cell::UnitCell{D,T}, lattice::Lattice{D}) where {D,T}

Calculate the k-point `k_point` corresponding to the k-point location `k_loc`.
"""
function calc_k_point!(k_point::AbstractVector{T}, k_loc::AbstractVector{Int},
                       unit_cell::UnitCell{D,T}, lattice::Lattice{D}) where {D,T}


    (; reciprocal_vecs) = unit_cell
    (; L, periodic) = lattice

    fill!(k_point,0.0)
    @fastmath @inbounds for d in eachindex(k_point)
        L_d = max( L[d]*periodic[d] , 1 )
        @views @. k_point += k_loc[d] * reciprocal_vecs[:,d] / L_d
    end

    return nothing
end

# calc_k_point!(; k_point, k_loc, unit_cell, lattice) = calc_k_point!(k_point, k_loc, unit_cell, lattice)

# """
#     calc_k_point(k_point::AbstractVector{T}, k_loc::AbstractVector{Int},
#         unit_cell::UnitCell{T}, lattice::Lattice) where {T}

# Return the k-point `k_point` corresponding to the k-point location `k_loc`.
# """
# function calc_k_point(k_loc::AbstractVector{Int}, unit_cell::UnitCell{T}, lattice::Lattice) where {T}
    
#     k_point = zeros(T,lattice.D)
#     calc_k_point!(k_point, k_loc, unit_cell, lattice)

#     return k_point
# end

# calc_k_point(; k_loc, unit_cell, lattice) = calc_k_point(k_loc, unit_cell, lattice)


# """
#     calc_k_points!(k_points::AbstractArray{T}, unit_cell::UnitCell{T},
#         lattice::Lattice) where {T}

# Calculate the k-point grid `k_points` assicated with a finite lattice.
# """
# function calc_k_points!(k_points::AbstractArray{T}, unit_cell::UnitCell{T}, lattice::Lattice) where {T}
    
#     @assert unit_cell.D == lattice.D
#     (; D, reciprocal_vecs) = unit_cell
#     (; L, N, periodic)     = lattice
#     k_loc                  = lattice.lvec
    
#     for ci in CartesianIndices( size(k_points)[2:D+1] )
#         for d in 1:D
#             k_loc[d] = ci[d] - 1
#         end
#         k_point = @view k_points[:,ci]
#         calc_k_point!(k_point, k_loc, unit_cell, lattice)
#     end

#     return nothing
# end

# calc_k_points!(; k_points, unit_cell, lattice) = calc_k_points!(k_points, unit_cell, lattice)

# """
#     calc_k_points!(unit_cell::UnitCell{T}, lattice::Lattice) where {T}

# Return the k-point grid assicated with a finite lattice.
# For a `D` dimensional lattice, a `D+1` dimensional array will be returned.
# If the system has open boundary conditions in a given direction, it will treat the linear
# extent of the system in that direction as equalling `L=1` for the purposes of calucting the k-points.
# """
# function calc_k_points(unit_cell::UnitCell{T}, lattice::Lattice) where {T}
    
#     (; D, L, periodic) = lattice
#     k_points = zeros(T, D, (max(L[d]*periodic[d] , 1) for d in 1:D)...)
#     calc_k_points!(k_points, unit_cell, lattice)

#     return k_points
# end

# calc_k_points(; unit_cell, lattice) = calc_k_points(unit_cell, lattice)


# """
#     bond_to_vec!(Δr::AbstractVector{T}, bond::Bond, unit_cell::UnitCell{T}) where {T}

# Calculate the displacement vector associated with a `bond`.
# """
# function bond_to_vec!(Δr::AbstractVector{T}, bond::Bond, unit_cell::UnitCell{T}) where {T}

#     (; displacement, orbitals) = bond
#     displacement_to_vec!(Δr,displacement,orbitals[1],orbitals[2],unit_cell)
#     return nothing
# end

# bond_to_vec!(; Δr, bond, unit_cell) = bond_to_vec!(Δr, bond, unit_cell)

# """
#     bond_to_vec(bond::Bond, unit_cell::UnitCell{T})

# Return the displacement vector associated with a `bond`.
# """
# function bond_to_vec(bond::Bond, unit_cell::UnitCell{T}) where {T}

#     Δr = zeros(T,bond.D)
#     bond_to_vec!(Δr, bond, unit_cell)
#     return Δr
# end

# bond_to_vec(; bond, unit_cell) = bond_to_vec(bond, unit_cell)


# """
#     build_neighbor_table(bond::Bond, unit_cell::UnitCell, lattice::Lattice)

# Construct the neighbor table corresponding to `bond`.
# """
# function build_neighbor_table(bond::Bond, unit_cell::UnitCell, lattice::Lattice)

#     @assert bond.D == unit_cell.D == lattice.D
#     @assert length(bond.orbitals) == 2
#     @assert 0 < bond.orbitals[1] <= unit_cell.n
#     @assert 0 < bond.orbitals[2] <= unit_cell.n

#     (; D, L, N, periodic)      = lattice
#     (; n)                      = unit_cell
#     (; displacement, orbitals) = bond

#     # initialize empty neighbor table
#     neighbor_table = Vector{Int}[]

#     # iterate over all unit cells
#     for u in 1:N
#         # get initial site
#         s₁ = loc_to_site(u, orbitals[1], unit_cell, lattice)
#         # get final site
#         s₂ = site_to_site(s₁, displacement, orbitals[2], unit_cell, lattice)
#         # check if final site was found
#         if s₂ != 0
#             # add to neighbor table
#             push!(neighbor_table, [s₁,s₂])
#         end
#     end

#     return hcat(neighbor_table...)
# end

# build_neighbor_table(; bond, unit_cell, lattice) = build_neighbor_table(bond, unit_cell, lattice)

# """
#     build_neighbor_table(bonds::AbstractVector{Bond}, unit_cell::UnitCell, lattice::Lattice)

# Construct the neighbor table corresponding to `bonds`.
# """
# function build_neighbor_table(bonds::AbstractVector{Bond}, unit_cell::UnitCell, lattice::Lattice)

#     neighbor_tables = Matrix{Int}[]
#     for i in eachindex(bonds)
#         neighbor_table = build_neighbor_table(bonds[i], unit_cell, lattice)
#         push!(neighbor_tables, neighbor_table)
#     end

#     return hcat(neighbor_tables...)
# end

# build_neighbor_table(; bonds, unit_cell, lattice) = build_neighbor_table(bonds, unit_cell, lattice)


# """
#     map_neighbor_table(neighbor_table::Matrix{Int})

# For a given neighbor table, return a dictionary that reports the bonds and neighbors
# associated with each site in the lattice. If `neighbor_table` is modified,
# then a new map must be constructed.
# """
# function map_neighbor_table(neighbor_table::Matrix{Int})

#     @assert size(neighbor_table,1) == 2
#     Nsites = maximum(neighbor_table)
#     Nbonds = size(neighbor_table,2)
#     nt_map = Dict( i => (bonds=Int[], neighbors=Int[]) for i in 1:Nsites)
#     for b in 1:Nbonds
#         i           = neighbor_table[1,b]
#         j           = neighbor_table[2,b]
#         info_i      = nt_map[i]
#         info_j      = nt_map[j]
#         bonds_i     = info_i.bonds
#         neighbors_i = info_i.neighbors
#         bonds_j     = info_j.bonds
#         neighbors_j = info_j.neighbors
#         push!(bonds_i, b)
#         push!(bonds_j, b)
#         push!(neighbors_i, j)
#         push!(neighbors_j, i)
#     end

#     return nt_map
# end


# """
#     function translational_avg!(fg::AbstractArray{Complex{T}}, f::AbstractArray{Complex{T}},
#         g::AbstractArray{Complex{T}}; restore::Bool=true) where {T<:AbstractFloat}

# Let `f[i]` and `g[j]` be two distinct multi-dimensional arrays, where `i` and `j` represent
# an index into them and also correspond to the position of a unit cell in a periodic finite lattice
# This method then computes in-place the product `(f⋅g)[i-j]` that is averaged over translation symmetry.
# If `restore = true` then `f` and `g` are left unchanged, otherwise they will be left modified.
# """
# function translational_avg!(fg::AbstractArray{Complex{T}}, f::AbstractArray{Complex{T}},
#     g::AbstractArray{Complex{T}}; restore::Bool=true) where {T<:AbstractFloat}
    
#     @assert size(fg) == size(f) "$(size(fg)) == $(size(f))"
#     @assert size(fg) == size(f) "$(size(fg)) == $(size(g))"

#     fft!(f)
#     fft!(g)
#     N  = length(f)
#     g′ = fg
#     circshift!(g′, g, size(fg,d)-1 for d in 1:ndims(fg))
#     reverse!(g′)
#     @. fg = f * g′ / N
#     ifft!(fg)
#     if restore
#         ifft!(f)
#         ifft!(g)
#     end
#     return nothing
# end

# translational_avg!(; fg, f, g, restore) = translational_avg!(fg, f, g, restore)