# Getting Started

First the [`LatticeUtilities.jl`](https://SmoQySuite.github.io/LatticeUtilities.jl/dev/)
package needs to be imported.

```jldoctest getting_started
julia> using LatticeUtilities
```

As an initial example, we construct an instance of the [`UnitCell`](@ref) type
to represent the unit cell for a square lattice:

```jldoctest getting_started
julia> square = UnitCell(lattice_vecs = [[1.,0.],[0.,1.]],
                         basis_vecs   = [[0.,0.]])
[UnitCell]

dimensions = 2
orbitals = 1

[UnitCell.lattice_vecs]

a_1 = [1.0, 0.0]
a_2 = [0.0, 1.0]

[UnitCell.reciprocal_vecs]

b_1 = [6.283185307179586, 0.0]
b_2 = [0.0, 6.283185307179586]

[UnitCell.basis_vecs]

b_1 = [0.0, 0.0]
```

Next we construct an instance of the type [`Lattice`](@ref) that describes the size
of a finite lattice for a given unit cell definition. In the example below we assume
periodic boundary conditions in the direction of both lattice vectors:

```jldoctest getting_started
julia> lattice = Lattice(L = [4,4], periodic = [true,true])
[Lattice]

dimensions   = 2
n_unit_cells = 16
size         = [4, 4]
periodic     = [true, true]
```

Bonds (or edges) in a lattice are represented by the [`Bond`](@ref) type.
Considering just nearest neighbors, there are two bonds that need to be defined:

```jldoctest getting_started
julia> bond_x = Bond(orbitals = (1,1), displacement = [1,0])
[[Bond]]

dimensions   = 2
orbitals     = [1, 1]
displacement = [1, 0]

julia> bond_y = Bond((1,1), [0,1])
[[Bond]]

dimensions   = 2
orbitals     = [1, 1]
displacement = [0, 1]
```

Using these two bond definitions, we can construct the corresponding neighbor table
using the [`build_neighbor_table`](@ref) method:

```jldoctest getting_started
julia> neighbor_table = build_neighbor_table([bond_x,bond_y], square, lattice)
2×32 Matrix{Int64}:
 1  2  3  4  5  6  7  8   9  10  11  …   8   9  10  11  12  13  14  15  16
 2  3  4  1  6  7  8  5  10  11  12     12  13  14  15  16   1   2   3   4
```