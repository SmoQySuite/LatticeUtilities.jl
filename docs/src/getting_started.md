# Getting Started

First the [`LatticeUtilities.jl`](https://cohensbw.github.io/LatticeUtilities.jl/dev/)
package needs to be imported.

```jldoctest getting_started
julia> using LatticeUtilities
```

As an initial example, we construct an instance of the [`UnitCell`](@ref) type
to represent the unit cell for a cubic lattice:

```jldoctest getting_started
julia> cubic = UnitCell(lattice_vecs = [[1.,0.,0.], [0.,1.,0.], [0.,0.,1.]],
                        basis_vecs   = [[0.,0.,0.]])
UnitCell{Float64}:
 - D = 3
 - n = 1
 - lattice_vecs =
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0
 - reciprocal_vecs =
3×3 Matrix{Float64}:
 6.28319  0.0      0.0
 0.0      6.28319  0.0
 0.0      0.0      6.28319
 - basis_vecs =
3×1 Matrix{Float64}:
 0.0
 0.0
 0.0
```

Next we construct an instance of the type [`Lattice`](@ref) that describes the size
of a finite lattice for a given unit cell definition. In the example below we assume
periodic boundary conditions in the direction of all three lattice vectors:

```jldoctest getting_started
julia> lattice = Lattice(L = [4,4,4], periodic = [true,true,true])
Lattice:
 - D = 3
 - N = 64
 - L = [4, 4, 4]
 - periodic = Bool[1, 1, 1]
```

Bonds or edeges in a lattice are represented by the [`Bond`](@ref) type.
Considering just nearest neighbors, there are three bonds that need to be defined:

```jldoctest getting_started
julia> bond_x = Bond(orbitals = [1,1], displacement = [1,0,0])
Bond:
 - D = 3
 - orbitals = [1, 1]
 - displacement = [1, 0, 0]

julia> bond_y = Bond([1,1], [0,1,0])
Bond:
 - D = 3
 - orbitals = [1, 1]
 - displacement = [0, 1, 0]

julia> bond_z = Bond([1,1], [0,0,1])
Bond:
 - D = 3
 - orbitals = [1, 1]
 - displacement = [0, 0, 1]
```

Using these three bond defintions, we can construct the corresponding neighbor table
using the [`build_neighbor_table`](@ref) method:

```jldoctest getting_started
julia> neighbor_table = build_neighbor_table([bond_x,bond_y,bond_z], cubic, lattice)
2×192 Matrix{Int64}:
 1  2  3  4  5  6  7  8   9  10  11  …  56  57  58  59  60  61  62  63  64
 2  3  4  1  6  7  8  5  10  11  12      8   9  10  11  12  13  14  15  16
```