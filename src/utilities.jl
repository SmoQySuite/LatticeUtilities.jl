"""
    get_num_sites(unit_cell::UnitCell,lattice::Lattice)

Returns the number of sites `Nₛ` in a finite lattice.
"""
function get_num_sites(unit_cell::UnitCell,lattice::Lattice)

    @assert unit_cell.D == lattice.D
    return unit_cell.n * lattice.N
end


"""
    valid_site(s::Int,unit_cell::UnitCell,lattice::Lattice)

Return whether `s` is a valid site index.
"""
function valid_site(s::Int,unit_cell::UnitCell,lattice::Lattice)

    @assert unit_cell.D == lattice.D
    Nₛ = get_num_sites(unit_cell, lattice)
    return 0 < s <= Nₛ
end


"""
    site_to_unitcell(s::Int,unit_cell::UnitCell,lattice::Lattice)

Return the unit cell `u` containing lattice site `s`.
"""
function site_to_unitcell(s::Int,unit_cell::UnitCell,lattice::Lattice)

    @assert valid_site(s,unit_cell,lattice)

    return (s-1) ÷ unit_cell.n + 1 
end


"""
    site_to_orbital(s::Int,unit_cell::UnitCell,lattice::Lattice)

Return the orbtial species of site `s`.
"""
function site_to_orbital(s::Int,unit_cell::UnitCell,lattice::Lattice)

    return mod1(s,unit_cell.n)
end


"""
    site_to_loc!(loc::AbstractVector{Int},s::Int,unit_cell::UnitCell,lattice::Lattice)

For a given site `s` in the lattice, calculate the location `l` of the unit cell it is in
and return the orbital species `o` of the site.
"""
function site_to_loc!(l::AbstractVector{Int},s::Int,unit_cell::UnitCell,lattice::Lattice)

    o = site_to_orbital(s,unit_cell,lattice)
    u = site_to_unitcell(s,unit_cell,lattice)
    unitcell_to_loc!(l,u,lattice)
    return o
end

"""
    site_to_loc(s::Int,unit_cell::UnitCell,lattice::Lattice)

For a given site `s` in the lattice, return the location `l` of the unit cell it is in
and the orbital species `o`.
"""
function site_to_loc(s::Int,unit_cell::UnitCell,lattice::Lattice)

    @assert unit_cell.D == lattice.D
    l = zeros(Int,unit_cell.D)
    o   = site_to_loc!(l,s,unit_cell,lattice)
    return (l, o)
end


"""
    loc_to_site(l::AbstractVector{Int},o::Int,unit_cell::UnitCell,lattice::Lattice)

Given a unit cell location `l` and orbital species `o`, return the corresponding
site `s` in the lattice.
"""
function loc_to_site(l::AbstractVector{Int},o::Int,unit_cell::UnitCell,lattice::Lattice)

    @assert unit_cell.D == lattice.D
    @assert 0 < o <= unit_cell.n
    @assert valid_location(l,lattice)

    u = loc_to_unitcell(l,lattice)
    s = unit_cell.n * (u-1) + o
    return s
end


"""
    site_to_site(s₁::Int,Δl::AbstractVector{Int},o₂::Int,unit_cell::UnitCell,lattice::Lattice)

Given an initial site `s₁`, and a displacement in unit cells `Δl` and a terminating orbital
species `o₂`, return the resulting site `s₂` in the lattice.
"""
function site_to_site(s₁::Int,Δl::AbstractVector{Int},o₂::Int,unit_cell::UnitCell,lattice::Lattice)

    (; D, n) = unit_cell
    l = lattice.lvec

    # check that initial site index is valid
    @assert valid_site(s₁, unit_cell, lattice)

    # get unit cell location containing s₁
    o₁ = site_to_loc!(l, s₁, unit_cell, lattice)

    # displace unit cell location
    @. l += Δl

    # apply periodic boundary conditions
    pbc!(l, lattice)

    # check if valid unit cell location
    @assert valid_location(l, lattiice)

    # get final site
    s₂ = loc_to_site(l, o₂, unit_cell, lattice)

    # check that final site index is valid
    @assert valid_site(s₂, unit_cell, lattice)

    return s₂
end


"""
    calc_k_point!(k_point::AbstractVector{T}, k_loc::AbstractVector{Int}, unit_cell::UnitCell{T},
        lattice::Lattice) where {T}

Calculate the k-point `k_point` corresponding to the k-point index `k_loc`.
"""
function calc_k_point!(k_point::AbstractVector{T}, k_loc::AbstractVector{Int}, unit_cell::UnitCell{T}, lattice::Lattice) where {T}


    @assert length(k_point) == length(k_loc) == length(unit_cell.D) == length(lattice.D)
    (; reciprocal_vecs) = unit_cell
    (; D, L, periodic) = lattice

    for d in 1:D
        l = max( L[d]*periodic[d] , 1 )
        @assert 0 <= k_loc[d] < l
        @views @. k_point = k_loc[d] * reciprocal_vecs[:,d] / l
    end

    return nothing
end

"""
    calc_k_point(k_point::AbstractVector{T}, k_loc::AbstractVector{Int}, unit_cell::UnitCell{T},
        lattice::Lattice) where {T}

Return the k-point `k_point` corresponding to the k-point index `k_loc`.
"""
function calc_k_point(k_loc::AbstractVector{Int}, unit_cell::UnitCell{T}, lattice::Lattice) where {T}
    
    k_point = zeros(T,lattice.D)
    calc_k_point!(k_point,k_loc,unit_cell)

    return k_point
end


"""
    calc_k_points!(k_points::AbstractArray{T}, unit_cell::UnitCell{T},
        lattice::Lattice) where {T}

Calculate the k-point grid `k_points` assicated with a finite lattice.
"""
function calc_k_points!(k_points::AbstractArray{T}, unit_cell::UnitCell{T}, lattice::Lattice) where {T}
    
    @assert unit_cell.D == lattice.D
    (; D, reciprocal_vecs) = unit_cell
    (; L, N, periodic)     = lattice
    k_loc                  = lattice.lvec
    
    for ci in CartesianIndices( size(k_points)[2:D+1] )
        for d in 1:D
            k_loc[d] = ci[d]
        end
        k_point = @view k_points[:,ci]
        calc_k_point!(k_point, k_loc, unit_cell, lattice)
    end

    return nothing
end

"""
    calc_k_points!(unit_cell::UnitCell{T}, lattice::Lattice) where {T}

Return the k-point grid assicated with a finite lattice.
For a `D` dimensional lattice, a `D+1` dimensional array will be returned.
If the system has open boundary conditions in a given direction, it will treat the linear
extent of the system in that direction as equalling `L=1` for the purposes of calucting the k-points.
"""
function calc_k_points(unit_cell::UnitCell{T}, lattice::Lattice) where {T}
    
    k_points = zeros(T, D, (max(L[d]*periodic[d] , 1) for d in 1:D)...)
    calc_k_points!(k_points, unit_cell, lattice)

    return k_points
end