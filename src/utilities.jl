"""
    get_Nₛ(unit_cell::UnitCell,lattice::Lattice)

Returns the number of sites `Nₛ` in a finite lattice.
"""
function get_Nₛ(unit_cell::UnitCell,lattice::Lattice)

    @assert unit_cell.D == lattice.D
    return unit_cell.n * lattice.N
end


"""
    valid_site(s::Int,unit_cell::UnitCell,lattice::Lattice)

Return whether `s` is a valid site index.
"""
function valid_site(s::Int,unit_cell::UnitCell,lattice::Lattice)

    @assert unit_cell.D == lattice.D
    Nₛ = get_Nₛ(unit_cell, lattice)
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

For a given site `s` in the lattice, calculate the location `loc` of the unit cell it is in
and return the orbital species `o` of the site.
"""
function site_to_loc!(loc::AbstractVector{Int},s::Int,unit_cell::UnitCell,lattice::Lattice)

    o = site_to_orbital(s,unit_cell,lattice)
    u = site_to_unitcell(s,unit_cell,lattice)
    unitcell_to_loc!(loc,u,lattice)
    return o
end

"""
    site_to_loc(s::Int,unit_cell::UnitCell,lattice::Lattice)

For a given site `s` in the lattice, return the location `loc` of the unit cell it is in
and the orbital species `o`.
"""
function site_to_loc(s::Int,unit_cell::UnitCell,lattice::Lattice)

    @assert unit_cell.D == lattice.D
    loc = zeros(Int,unit_cell.D)
    o   = site_to_loc!(loc,s,unit_cell,lattice)
    return (loc, o)
end


"""
    loc_to_site(loc::AbstractVector{Int},o::Int,unit_cell::UnitCell,lattice::Lattice)

Given a unit cell location `loc` and orbital species `o`, return the corresponding
site `s` in the lattice.
"""
function loc_to_site(loc::AbstractVector{Int},o::Int,unit_cell::UnitCell,lattice::Lattice)

    @assert unit_cell.D == lattice.D
    @assert 0 < o <= unit_cell.n
    @assert valid_location(loc,lattice)

    u = loc_to_unitcell(loc,lattice)
    s = unit_cell.n * (u-1) + o
    return s
end


"""
    site_to_site(s::Int,Δl::AbstractVector{Int},o::Int,unit_cell::UnitCell,lattice::Lattice)

Given an initial site `s₁`, and a displacement in unit cell `dl` and a terminating orbital
species `o₂`, return the resulting site `s₂` in the lattice.
"""
function site_to_site(s₁::Int,dl::AbstractVector{Int},o₂::Int,unit_cell::UnitCell,lattice::Lattice)

    (; D, n) = unit_cell
    (; lvec) = lattice

    # check that initial site index is valid
    @assert valid_site(s₁, unit_cell, lattice)

    # get unit cell location containing intial site
    o₁ = site_to_loc!(lvec, s, unit_cell, lattice)

    # displace unit cell location
    @. lvec += dl

    # apply periodic boundary conditions
    pbc!(lvec, lattice)

    # check if valid unit cell location
    @assert valid_location(loc, lattiice)

    # get final site
    s₂ = loc_to_site(loc, o₂, unit_cell, lattice)

    # check that final site index is valid
    @assert valid_site(s₂, unit_cell, lattice)

    return s₂
end


"""
    calc_k_points(unit_cell::UnitCell{T}, lattice::Lattice) where {T}

Return the k-points grid assicated with a finite lattice.
For a `D` dimensional lattice, a `D+1` dimensional array will be returned.
If the system has open boundary conditions in a given direction, it will treat the linear
extent of the system in that direction as equalling `L=1` for the purposes of calucting the k-points.
"""
function calc_k_points(unit_cell::UnitCell{T}, lattice::Lattice) where {T}
    
    @assert unit_cell.D == lattice.D
    (; D, reciprocal_vecs) = unit_cell
    (; L, periodic)        = lattice

    k = zeros(T, D, (max(L[d]*periodic[d] , 1) for d in 1:D)...)
    
    for ci in CartesianIndices(k)
        for d in 1:D
            k[ci] += (ci[d+1]-1)/(size(k,d+1)) * reciprocal_vecs[ci[1],d]
        end
    end

    return k
end