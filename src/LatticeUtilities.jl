module LatticeUtilities

using LinearAlgebra
using StaticArrays
using Printf

# rexport used static arrays types from the StaticArrays package
export SMatrix, SVector, MVector

include("UnitCell.jl")
export UnitCell
export loc_to_pos!, loc_to_pos
export displacement_to_vec!, displacement_to_vec

include("Lattice.jl")
export Lattice
export valid_location, pbc!
export unitcell_to_loc!, unitcell_to_loc, loc_to_unitcell

include("Bond.jl")
export Bond

include("utilities.jl")
export norbits, nsites
export get_num_sites, valid_site
export site_to_unitcell, site_to_orbital
export site_to_loc!, site_to_loc
export loc_to_site, site_to_site
export calc_k_point!, calc_k_point
export calc_k_points!, calc_k_points
export bond_to_vec!, bond_to_vec
export build_neighbor_table
export simplify!
export sites_to_displacement!, sites_to_displacement
export sites_to_bond!, sites_to_bond
export map_neighbor_table

end # module
