var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [LatticeUtilities]","category":"page"},{"location":"api/#LatticeUtilities.Bond","page":"API","title":"LatticeUtilities.Bond","text":"Bond\n\nDefines a bond defintion in lattice.\n\nFields\n\nD::Int64\nDimension of system bond is in.\no::Vector{Int64}\nInitial and final orbital species respectively.\nΔl::Vector{Int64}\nDisplacement in unit cells.\n\n\n\n\n\n","category":"type"},{"location":"api/#LatticeUtilities.Bond-Tuple{AbstractVector{Int64}, AbstractVector{Int64}}","page":"API","title":"LatticeUtilities.Bond","text":"Bond(Δl::AbstractVector{Int},o::AbstractVector{Int})\n\nConstrcut a Bond\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.Lattice","page":"API","title":"LatticeUtilities.Lattice","text":"Lattice\n\nA type defining a finite lattice in arbitary dimensions.\n\nFields\n\nD::Int64\nNumber of spatial dimensions.\nN::Int64\nNumber of unit cells.\nL::Vector{Int64}\nLinear extent of lattice in the direction of each lattice vector.\nperiodic::Vector{Bool}\nWhether the lattice is periodic in the direction of each lattice vector.\nlvec::Vector{Int64}\nStorage space for representing a location or displacement in the lattice.\n\n\n\n\n\n","category":"type"},{"location":"api/#LatticeUtilities.Lattice-Tuple{Vector{Int64}, Vector{Bool}}","page":"API","title":"LatticeUtilities.Lattice","text":"Lattice(L::Vector{Int},periodic::Vector{Bool})\n\nConstructs a Lattice.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.UnitCell","page":"API","title":"LatticeUtilities.UnitCell","text":"UnitCell{T<:AbstractFloat}\n\nA type defining a unit cell.\n\nFields\n\nD::Int64\nNumber of spatial dimensions.\nn::Int64\nOrbitals per unit cell.\nlattice_vecs::Matrix{T} where T<:AbstractFloat\nMatrix where the columns are the lattice vectors.\nreciprocal_vecs::Matrix{T} where T<:AbstractFloat\nMatrix where the columns are the reciprocal lattice vectors.\nbasis_vecs::Matrix{T} where T<:AbstractFloat\nMatrix where the columns are then basis vectors.\n\n\n\n\n\n","category":"type"},{"location":"api/#LatticeUtilities.UnitCell-Union{Tuple{T}, Tuple{AbstractArray{Vector{T}, 1}, AbstractArray{Vector{T}, 1}}} where T<:AbstractFloat","page":"API","title":"LatticeUtilities.UnitCell","text":"UnitCell(lattice_vecs::AbstractVector{Vector{T}}, basis_vecs::AbstractVector{Vector{T}})\n    where {T<:AbstractFloat}\n\nConstrcuts a UnitCell.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.UnitCell-Union{Tuple{T}, Tuple{Matrix{T}, Matrix{T}}} where T<:AbstractFloat","page":"API","title":"LatticeUtilities.UnitCell","text":"UnitCell(lattice_vecs::Matrix{T}, basis_vecs::Matrix{T}) where {T<:AbstractFloat}\n\nConstrcuts a UnitCell.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.show-Tuple{IO, Lattice}","page":"API","title":"Base.show","text":"Base.show(io::IO, lattice::Lattice)\n\nShow lattice.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.show-Union{Tuple{T}, Tuple{IO, UnitCell{T}}} where T","page":"API","title":"Base.show","text":"Base.show(io::IO, uc::UnitCell{T}) where {T}\nBase.show(io::IO, ::MIME\"text/plain\", uc::UnitCell{T}) where {T}\n\nShow unit cell.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.bond_to_vec!-Union{Tuple{T}, Tuple{AbstractVector{T}, Bond, UnitCell{T}}} where T","page":"API","title":"LatticeUtilities.bond_to_vec!","text":"bond_to_vec!(Δr::AbstractVector{T},bond::Bond,unit_cell::UnitCell{T}) where {T}\n\nCalculate the displacement vector associated with a bond.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.bond_to_vec-Union{Tuple{T}, Tuple{Bond, UnitCell{T}}} where T","page":"API","title":"LatticeUtilities.bond_to_vec","text":"bond_to_vec(bond::Bond,unit_cell::UnitCell{T})\n\nReturn the displacement vector associated with a bond.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.build_neighbor_table-Tuple{AbstractVector{Bond}, UnitCell, Lattice}","page":"API","title":"LatticeUtilities.build_neighbor_table","text":"build_neighbor_table(bonds::AbstractVector{Bond}, unit_cell::UnitCell, lattice::Lattice)\n\nConstruct the neighbor table corresponding to bonds.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.build_neighbor_table-Tuple{Bond, UnitCell, Lattice}","page":"API","title":"LatticeUtilities.build_neighbor_table","text":"build_neighbor_table(bond::Bond, unit_cell::UnitCell, lattice::Lattice)\n\nConstruct the neighbor table corresponding to bond.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.calc_k_point!-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{Int64}, UnitCell{T}, Lattice}} where T","page":"API","title":"LatticeUtilities.calc_k_point!","text":"calc_k_point!(k_point::AbstractVector{T}, k_loc::AbstractVector{Int},\n    unit_cell::UnitCell{T}, lattice::Lattice) where {T}\n\nCalculate the k-point k_point corresponding to the k-point location k_loc.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.calc_k_point-Union{Tuple{T}, Tuple{AbstractVector{Int64}, UnitCell{T}, Lattice}} where T","page":"API","title":"LatticeUtilities.calc_k_point","text":"calc_k_point(k_point::AbstractVector{T}, k_loc::AbstractVector{Int},\n    unit_cell::UnitCell{T}, lattice::Lattice) where {T}\n\nReturn the k-point k_point corresponding to the k-point location k_loc.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.calc_k_points!-Union{Tuple{T}, Tuple{AbstractArray{T}, UnitCell{T}, Lattice}} where T","page":"API","title":"LatticeUtilities.calc_k_points!","text":"calc_k_points!(k_points::AbstractArray{T}, unit_cell::UnitCell{T},\n    lattice::Lattice) where {T}\n\nCalculate the k-point grid k_points assicated with a finite lattice.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.calc_k_points-Union{Tuple{T}, Tuple{UnitCell{T}, Lattice}} where T","page":"API","title":"LatticeUtilities.calc_k_points","text":"calc_k_points!(unit_cell::UnitCell{T}, lattice::Lattice) where {T}\n\nReturn the k-point grid assicated with a finite lattice. For a D dimensional lattice, a D+1 dimensional array will be returned. If the system has open boundary conditions in a given direction, it will treat the linear extent of the system in that direction as equalling L=1 for the purposes of calucting the k-points.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.displacement_to_vec!-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{Int64}, Int64, Int64, UnitCell{T}}} where T","page":"API","title":"LatticeUtilities.displacement_to_vec!","text":"displacement_to_vec!(Δr::AbstractVector{T}, Δl::AbstractVector{Int},\n    o₁::Int, o₂::Int, unit_cell::UnitCell{T}) where {T}\n\nComputes the position space displacement vector Δr corresponding to a displacement definition given by initial and final orbitals o₁ and o₂ in the unit cell respectively, along with a displacement in unit cells Δl.\n\nArguments\n\nΔr::AbstractVector{T}: displacement vector in position space.\nΔl::AbstractVector{Int}: displacement in unit cells.\no₁::Int: initial orbital in unit cell.\no₂::Int: final orbital in unit cell.\nunit_cell::UnitCell{T}: unit cell.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.displacement_to_vec-Union{Tuple{T}, Tuple{AbstractVector{Int64}, Int64, Int64, UnitCell{T}}} where T","page":"API","title":"LatticeUtilities.displacement_to_vec","text":"displacement_to_vec(Δl::AbstractVector{Int}, o₁::Int, o₂::Int,\n    unit_cell::UnitCell{T}) where {T}\n\nReturns the position space displacement vector Δr corresponding to a displacement definition given by initial and final orbitals o₁ and o₂ in the unit cell respectively, along with a displacement in unit cells Δl.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.get_num_sites-Tuple{UnitCell, Lattice}","page":"API","title":"LatticeUtilities.get_num_sites","text":"get_num_sites(unit_cell::UnitCell,lattice::Lattice)\n\nReturns the number of sites Nₛ in a finite lattice.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.loc_to_pos!-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{Int64}, Int64, UnitCell{T}}} where T","page":"API","title":"LatticeUtilities.loc_to_pos!","text":"loc_to_pos!(r::AbstractVector{T}, l::AbstractVector{Int}, s::Int, unit_cell::UnitCell{T})\n    where {T}\n\nCalculate the position r of a orbital o at location l.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.loc_to_pos!-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{Int64}, UnitCell{T}}} where T","page":"API","title":"LatticeUtilities.loc_to_pos!","text":"loc_to_pos!(r::AbstractVector{T}, l::AbstractVector{Int}, unit_cell::UnitCell{T})\n    where {T}\n\nCalculate the position r of a unit cell at location l.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.loc_to_pos-Union{Tuple{T}, Tuple{AbstractVector{Int64}, Int64, UnitCell{T}}} where T","page":"API","title":"LatticeUtilities.loc_to_pos","text":"loc_to_pos(l::AbstractVector{Int}, s::Int, unit_cell::UnitCell{T})\n    where {T}\n\nReturn the position r of a orbital o at location l.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.loc_to_pos-Union{Tuple{T}, Tuple{AbstractVector{Int64}, UnitCell{T}}} where T","page":"API","title":"LatticeUtilities.loc_to_pos","text":"loc_to_pos(l::AbstractVector{Int}, unit_cell::UnitCell{T})\n    where {T}\n\nReturn the position r of a unit cell at location l.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.loc_to_site-Tuple{AbstractVector{Int64}, Int64, UnitCell, Lattice}","page":"API","title":"LatticeUtilities.loc_to_site","text":"loc_to_site(l::AbstractVector{Int},o::Int,unit_cell::UnitCell,lattice::Lattice)\n\nGiven a unit cell location l and orbital species o, return the corresponding site s in the lattice. If the location is not valid owing to open boundary conditions then return s = 0.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.loc_to_site-Tuple{Int64, Int64, UnitCell, Lattice}","page":"API","title":"LatticeUtilities.loc_to_site","text":"loc_to_site(u::Int,o::Int,unit_cell::UnitCell,lattice::Lattice)\n\nGiven a unit cell index u and orbital o, return the correspond site s.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.loc_to_unitcell-Tuple{AbstractVector{Int64}, Lattice}","page":"API","title":"LatticeUtilities.loc_to_unitcell","text":"loc_to_unitcell(loc::AbstractVector{Int},lattice::Lattice)\n\nReturn the unit cell u found at location l in the lattice.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.pbc!-Tuple{AbstractVector{Int64}, Lattice}","page":"API","title":"LatticeUtilities.pbc!","text":"pbc!(loc::AbstractVector{Int}, lat::Lattice)\n\nApply periodic boundary to unit cell location loc.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.site_to_loc!-Tuple{AbstractVector{Int64}, Int64, UnitCell, Lattice}","page":"API","title":"LatticeUtilities.site_to_loc!","text":"site_to_loc!(loc::AbstractVector{Int},s::Int,unit_cell::UnitCell,lattice::Lattice)\n\nFor a given site s in the lattice, calculate the location l of the unit cell it is in and return the orbital species o of the site.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.site_to_loc-Tuple{Int64, UnitCell, Lattice}","page":"API","title":"LatticeUtilities.site_to_loc","text":"site_to_loc(s::Int,unit_cell::UnitCell,lattice::Lattice)\n\nFor a given site s in the lattice, return the location l of the unit cell it is in and the orbital species o.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.site_to_orbital-Tuple{Int64, UnitCell, Lattice}","page":"API","title":"LatticeUtilities.site_to_orbital","text":"site_to_orbital(s::Int,unit_cell::UnitCell,lattice::Lattice)\n\nReturn the orbtial species of site s.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.site_to_site-Tuple{Int64, AbstractVector{Int64}, Int64, UnitCell, Lattice}","page":"API","title":"LatticeUtilities.site_to_site","text":"site_to_site(s₁::Int,Δl::AbstractVector{Int},o₂::Int,unit_cell::UnitCell,lattice::Lattice)\n\nGiven an initial site s₁, and a displacement in unit cells Δl and a terminating orbital species o₂, return the resulting site s₂ in the lattice. If the displacement is not allowed as a result of open boundary conditions, then  s₂=0 is returned.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.site_to_unitcell-Tuple{Int64, UnitCell, Lattice}","page":"API","title":"LatticeUtilities.site_to_unitcell","text":"site_to_unitcell(s::Int,unit_cell::UnitCell,lattice::Lattice)\n\nReturn the unit cell u containing lattice site s.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.sort_neighbor_table!-Tuple{Matrix{Int64}}","page":"API","title":"LatticeUtilities.sort_neighbor_table!","text":"sort_neighbor_table!(neighbor_table::Matrix{Int})\n\nSorts neighbor_table so that the first row is in strictly ascending order, and for fixed values in the first row, the second row is also in ascending order. Also returns the inverse of the sorting perumtation, so original order of neighbors in neighbor_table can be easily recovered.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.sorted_neighbor_table_perm!-Tuple{Matrix{Int64}}","page":"API","title":"LatticeUtilities.sorted_neighbor_table_perm!","text":"sorted_neighbor_table_perm!(neighbor_table::Matrix{Int})\n\nReturns the permutation that sorts neighbor_table so that the first row is in strictly ascending order, and for fixed values in the first row, the second row is also in ascending order. This method also modifies the neighbor_table such that the smaller index in each column is always in the first row.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.unitcell_to_loc!-Tuple{AbstractVector{Int64}, Int64, Lattice}","page":"API","title":"LatticeUtilities.unitcell_to_loc!","text":"unitcell_to_loc!(loc::AbstractVector{Int},u::Int,lattice::Lattice)\n\nCalculate the location l of a unit cell u.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.unitcell_to_loc-Tuple{Int64, Lattice}","page":"API","title":"LatticeUtilities.unitcell_to_loc","text":"unitcell_to_loc(u::Int,lattice::Lattice)\n\nReturn the location l of a unit cell u.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.valid_location-Tuple{AbstractVector{Int64}, Lattice}","page":"API","title":"LatticeUtilities.valid_location","text":"valid_location(loc::AbstractVector{Int}, lat::Lattice)\n\nDetermine if loc describes a valid location in the lattice.\n\n\n\n\n\n","category":"method"},{"location":"api/#LatticeUtilities.valid_site-Tuple{Int64, UnitCell, Lattice}","page":"API","title":"LatticeUtilities.valid_site","text":"valid_site(s::Int,unit_cell::UnitCell,lattice::Lattice)\n\nReturn whether s is a valid site index.\n\n\n\n\n\n","category":"method"},{"location":"getting_started/#Getting-Started","page":"Getting Started","title":"Getting Started","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"First the LatticeUtilities.jl package needs to be imported.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"julia> using LatticeUtilities","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"As an initial example, we construct UnitCell type to represent the unit cell for a cubic lattice:","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"julia> cubic = UnitCell([[1.,0.,0.], [0.,1.,0.], [0.,0.,1.]],\n                        [[0.,0.,0.]])\nUnitCell{Float64}:\n- D = 3\n- n = 1\n- lattice_vecs =\n3×3 Matrix{Float64}:\n 1.0  0.0  0.0\n 0.0  1.0  0.0\n 0.0  0.0  1.0\n- reciprocal_vecs =\n3×3 Matrix{Float64}:\n 6.28319  0.0      0.0\n 0.0      6.28319  0.0\n 0.0      0.0      6.28319\n- basis_vecs =\n3×1 Matrix{Float64}:\n 0.0\n 0.0\n 0.0","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Next we construct an instance of the type Lattice that describes the size of a finite lattice for a given unit cell definition. In the example below we assume periodic boundary conditions in the direction of all three lattice vectors:","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"julia> lattice = Lattice([4,4,4],[true,true,true])\nLattice:\n- D = 3\n- N = 64\n- L = [4, 4, 4]\n- periodic = Bool[1, 1, 1]","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Bonds or edeges in a lattice are represented by the Bond type. Considering just nearest neighbors, there are three bonds that need to be defined:","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"julia> bond_x = Bond([1,1],[1,0,0])\nBond(3, [1, 1], [1, 0, 0])\n\njulia> bond_y = Bond([1,1],[0,1,0])\nBond(3, [1, 1], [0, 1, 0])\n\njulia> bond_z = Bond([1,1],[0,0,1])\nBond(3, [1, 1], [0, 0, 1])","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Using these three bond defintions, we can construct the corresponding neighbor table using the build_neighbor_table method:","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"julia> neighbor_table = build_neighbor_table([bond_x,bond_y,bond_z],cubic,lattice)\n2×192 Matrix{Int64}:\n 1  2  3  4  5  6  7  8   9  10  11  …  56  57  58  59  60  61  62  63  64\n 2  3  4  1  6  7  8  5  10  11  12      8   9  10  11  12  13  14  15  16","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"This page includes several examples on how to use the LatticeUtilities.jl package.","category":"page"},{"location":"examples/#Kagome-Lattice","page":"Examples","title":"Kagome Lattice","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"As a first example we will represent a Kagome lattice. The lattice vectors for the Kagome lattice are","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"beginalign*\nmathbfa_1 = left(10right)\nmathbfa_2 = left(frac12fracsqrt32right)\nendalign*","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"and with three orbitals per unit cell, the basis vectors are","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"beginalign*\nmathbfr_rm 1 = left(00right)\nmathbfr_rm 2 = left(frac120right)\nmathbfr_rm 3 = left(frac14fracsqrt34right)\nendalign*","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The reciprocal lattice vectors are then","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"beginalign*\nmathbfb_1 = 2pi left(1-frac1sqrt3right)\nmathbfb_2 = 2pi left(0frac2sqrt3right)\nendalign*","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We begin by constructing an instance of the type UnitCell to represent the Kagome lattice unit cell:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> kagome = UnitCell([[1.0,0.0], [1/2,√3/2]],\n                         [[0.0,0.0], [1/2,0.0], [1/4,√3/4]])\nUnitCell{Float64}:\n- D = 2\n- n = 3\n- lattice_vecs =\n2×2 Matrix{Float64}:\n 1.0  0.5\n 0.0  0.866025\n- reciprocal_vecs =\n2×2 Matrix{Float64}:\n  6.28319  0.0\n -3.6276   7.2552\n- basis_vecs =\n2×3 Matrix{Float64}:\n 0.0  0.5  0.25\n 0.0  0.0  0.433013","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"It is straightforward to calculate position space vectors of the form mathbfr = n_1 mathbfa_1 + n_2 mathbfa_2 + mathbfr_alpha using the method loc_to_pos, where mathbfr_alpha is one of the three possible basis vectors:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> n₁,n₂,α = 1,1,2\n(1, 1, 2)\n\njulia> loc_to_pos([n₁,n₂],α,kagome)\n2-element Vector{Float64}:\n 2.0\n 0.8660254037844386","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Displacement vectors of the form mathbfr = n_1 mathbfa_1 + n_2 mathbfa_2 + (mathbfr_beta - mathbfr_alpha) can also be caluculated using the displacement_to_vec method:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> n₁,n₂,α,β = 1,1,2,3\n(1, 1, 2, 3)\n\njulia> displacement_to_vec([n₁,n₂],α,β,kagome)\n2-element Vector{Float64}:\n 1.25\n 1.299038105676658","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Note that there are memory allocation free variants of many methods, such as loc_to_pos! and displacement_to_vec!, with names ending in ! that work by modifying the first array argument in place.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The next step is to construct an instance of the type Lattice, which can be used in conjunction with UnitCell instance of like dimenion D, to represent a finite lattice. In this example we will consider a 3 times 3 unit cell Kagome lattice with periodic boundary conditions in the direction of both lattice vectors:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> lattice = Lattice([3,3], [true,true])\nLattice:\n- D = 2\n- N = 9\n- L = [3, 3]\n- periodic = Bool[1, 1]","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Given an initial site, a displacement in unit cells, and a terminating orbital species, we can calculate the final site using the site_to_site method:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> site_to_site(17, [1,1], 2, kagome, lattice)\n20","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We can also calculate the k-points associated with out finite lattice using the calc_k_points method:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> calc_k_points(kagome,lattice)\n2×3×3 Array{Float64, 3}:\n[:, :, 1] =\n 0.0   2.0944   4.18879\n 0.0  -1.2092  -2.4184\n\n[:, :, 2] =\n 0.0     2.0944  4.18879\n 2.4184  1.2092  0.0\n\n[:, :, 3] =\n 0.0     2.0944  4.18879\n 4.8368  3.6276  2.4184","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Individual k-points can be caluclated using the calc_k_point and calc_k_point! methods.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"If we concern ourselves with just nearest neighbors, there are six bond defintions that need to be defined:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> bond_1 = Bond([1,2],[0,0])\nBond(2, [1, 2], [0, 0])\n\njulia> bond_2 = Bond([1,3],[0,0])\nBond(2, [1, 3], [0, 0])\n\njulia> bond_3 = Bond([2,3],[0,0])\nBond(2, [2, 3], [0, 0])\n\njulia> bond_4 = Bond([2,1],[1,0])\nBond(2, [2, 1], [1, 0])\n\njulia> bond_5 = Bond([3,1],[0,1])\nBond(2, [3, 1], [0, 1])\n\njulia> bond_6 = Bond([3,2],[-1,1])\nBond(2, [3, 2], [-1, 1])","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now we are ready to build the corresponding neighbor table:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> neighbor_table = build_neighbor_table([bond_1, bond_2, bond_3,\n                                              bond_4, bond_5, bond_6],\n                                              kagome, lattice)\n2×54 Matrix{Int64}:\n 1  4  7  10  13  16  19  22  25  1  4  …   3   6   9  12  15  18  21  24  27\n 2  5  8  11  14  17  20  23  26  3  6     17  11  14  26  20  23   8   2   5","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In certain situations it is useful to sort neighbor_table so that the first row is in strictly ascending order, and for constant value in the first row the second row is also ascending. This sorting can be performed using the sort_neighbor_table! method:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> inv_perm = sort_neighbor_table!(neighbor_table); neighbor_table\n2×54 Matrix{Int64}:\n 1  1  1   1  2  2   2   3   3  4  4  …  20  20  22  22  23  23  25  25  26\n 2  3  8  21  3  4  24  10  17  5  6     21  22  23  24  24  25  26  27  27","category":"page"},{"location":"#LatticeUtilities.jl","page":"LatticeUtilities.jl","title":"LatticeUtilities.jl","text":"","category":"section"},{"location":"","page":"LatticeUtilities.jl","title":"LatticeUtilities.jl","text":"Documentation for the LatticeUtilities.jl package. This package exports a suite of types and methods useful for defining arbitrary lattice geometries, and the construction of neighbor tables.","category":"page"}]
}
