# Examples

This page includes several examples on how to use the
[`LatticeUtilities.jl`](https://cohensbw.github.io/LatticeUtilities.jl/dev/) package.

## Kagome Lattice

As a first example we will represent a Kagome lattice.
The lattice vectors for the Kagome lattice are

```math
\begin{align*}
\mathbf{a}_{1} &= \left(1,0\right)\\
\mathbf{a}_{2} &= \left(\frac{1}{2},\frac{\sqrt{3}}{2}\right),
\end{align*}
```

and with three orbitals per unit cell, the basis vectors are

```math
\begin{align*}
\mathbf{r}_{{\rm 1}} &= \left(0,0\right)\\
\mathbf{r}_{{\rm 2}} &= \left(\frac{1}{2},0\right)\\
\mathbf{r}_{{\rm 3}} &= \left(\frac{1}{4},\frac{\sqrt{3}}{4}\right).
\end{align*}
```

The reciprocal lattice vectors are then

```math
\begin{align*}
\mathbf{b}_{1} &= 2\pi \left(1,-\frac{1}{\sqrt{3}}\right)\\
\mathbf{b}_{2} &= 2\pi \left(0,\frac{2}{\sqrt{3}}\right).
\end{align*}
```

We begin by constructing an instance of the type [`UnitCell`](@ref) to represent
the Kagome lattice unit cell:

```jldoctest kagome; setup = :(using LatticeUtilities)
julia> kagome = UnitCell(lattice_vecs = [[1.0,0.0], [1/2,√3/2]],
                         basis_vecs   = [[0.0,0.0], [1/2,0.0], [1/4,√3/4]])
UnitCell{Float64}:
 • D = 2
 • n = 3
 • lattice_vecs =
2×2 Matrix{Float64}:
 1.0  0.5
 0.0  0.866025
 • reciprocal_vecs =
2×2 Matrix{Float64}:
  6.28319  0.0
 -3.6276   7.2552
 • basis_vecs =
2×3 Matrix{Float64}:
 0.0  0.5  0.25
 0.0  0.0  0.433013
```

It is straightforward to calculate position space vectors of the form
``\mathbf{r} = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2 + \mathbf{r}_\alpha``
using the method [`loc_to_pos`](@ref), where ``\mathbf{r}_\alpha`` is one of the
three possible basis vectors:

```jldoctest kagome
julia> n₁,n₂,α = 1,1,2
(1, 1, 2)

julia> loc_to_pos([n₁,n₂],α,kagome)
2-element Vector{Float64}:
 2.0
 0.8660254037844386
```

Displacement vectors of the form
``\mathbf{r} = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2 + (\mathbf{r}_\beta - \mathbf{r}_\alpha)``
can also be calculate using the [`displacement_to_vec`](@ref) method:

```jldoctest kagome
julia> n₁,n₂,α,β = 1,1,2,3
(1, 1, 2, 3)

julia> displacement_to_vec([n₁,n₂],α,β,kagome)
2-element Vector{Float64}:
 1.25
 1.299038105676658
```

Note that there are memory allocation free variants of many methods, such as [`loc_to_pos!`](@ref)
and [`displacement_to_vec!`](@ref), with names ending in `!` that work by modifying
the first array argument in place.

The next step is to construct an instance of the type [`Lattice`](@ref), which can be used in conjunction
with a [`UnitCell`](@ref) instance of like dimension `D`, to represent a finite lattice.
In this example we will consider a ``3 \times 3`` unit cell Kagome lattice with periodic boundary conditions
in the direction of both lattice vectors:

```jldoctest kagome
julia> lattice = Lattice(L = [3,3], periodic = [true,true])
Lattice:
 • D = 2
 • N = 9
 • L = [3, 3]
 • periodic = Bool[1, 1]
```

Given an initial site, a displacement in unit cells, and a terminating orbital species,
we can calculate the final site using the [`site_to_site`](@ref) method:

```jldoctest kagome
julia> site_to_site(17, [1,1], 2, kagome, lattice)
20
```

The set of relevant ``\mathbf{k}``-points associated with a finite lattice in
three dimensions is given by
```math
\mathbf{k} = \frac{n_1}{L_1}\mathbf{b}_1 + \frac{n_2}{L_2}\mathbf{b}_2 + \frac{n_3}{L_3}\mathbf{b}_3
```
where each index runs from ``n_i = 0, 1, \dots, L_i-1``.\
This set of ``\mathbf{k}``-points can be constructed using the [`calc_k_points`](@ref) method:

```jldoctest kagome
julia> calc_k_points(kagome,lattice)
2×3×3 Array{Float64, 3}:
[:, :, 1] =
 0.0   2.0944   4.18879
 0.0  -1.2092  -2.4184

[:, :, 2] =
 0.0     2.0944  4.18879
 2.4184  1.2092  0.0

[:, :, 3] =
 0.0     2.0944  4.18879
 4.8368  3.6276  2.4184
```

Individual ``\mathbf{k}``-points can be calculated using the [`calc_k_point`](@ref)
and [`calc_k_point!`](@ref) methods.

If we concern ourselves with just nearest neighbor bonds,
we need to define six instances of [`Bond`](@ref) type:

```jldoctest kagome
julia> bond_1 = Bond(orbitals = [1,2], displacement = [0,0])
Bond:
 • D = 2
 • orbitals = (1, 2)
 • displacement = [0, 0]

julia> bond_2 = Bond([1,3], [0,0])
Bond:
 • D = 2
 • orbitals = (1, 3)
 • displacement = [0, 0]

julia> bond_3 = Bond([2,3], [0,0])
Bond:
 • D = 2
 • orbitals = (2, 3)
 • displacement = [0, 0]

julia> bond_4 = Bond([2,1], [1,0])
Bond:
 • D = 2
 • orbitals = (2, 1)
 • displacement = [1, 0]

julia> bond_5 = Bond([3,1], [0,1])
Bond:
 • D = 2
 • orbitals = (3, 1)
 • displacement = [0, 1]

julia> bond_6 = Bond([3,2], [-1,1])
Bond:
 • D = 2
 • orbitals = (3, 2)
 • displacement = [-1, 1]
```

Now we are ready to build the corresponding neighbor table:
```jldoctest kagome
julia> neighbor_table = build_neighbor_table([bond_1, bond_2, bond_3,
                                              bond_4, bond_5, bond_6],
                                              kagome, lattice)
2×54 Matrix{Int64}:
 1  4  7  10  13  16  19  22  25  1  4  …   3   6   9  12  15  18  21  24  27
 2  5  8  11  14  17  20  23  26  3  6     17  11  14  26  20  23   8   2   5
```

It is frequently useful to have a reference that lists the neighbors
associated with each site in the lattice, as well as the bonds, which
are specified by a column index into `neighbor_table`. This information
is returned by the function [`map_neighbor_table`](@ref):

```jldoctest kagome
julia> neighbor_table_map = map_neighbor_table(neighbor_table)
Dict{Int64, NamedTuple{(:bonds, :neighbors), Tuple{Vector{Int64}, Vector{Int64}}}} with 27 entries:
  5  => (bonds = [2, 20, 29, 54], neighbors = [4, 6, 7, 27])
  16 => (bonds = [6, 15, 32, 39], neighbors = [17, 18, 14, 9])
  20 => (bonds = [7, 25, 34, 50], neighbors = [19, 21, 22, 15])
  12 => (bonds = [13, 22, 40, 49], neighbors = [10, 11, 19, 26])
  24 => (bonds = [17, 26, 44, 53], neighbors = [22, 23, 4, 2])
  8  => (bonds = [3, 21, 30, 52], neighbors = [7, 9, 1, 21])
  17 => (bonds = [6, 24, 33, 46], neighbors = [16, 18, 10, 3])
  1  => (bonds = [1, 10, 30, 43], neighbors = [2, 3, 8, 21])
  19 => (bonds = [7, 16, 36, 40], neighbors = [20, 21, 26, 12])
  22 => (bonds = [8, 17, 34, 41], neighbors = [23, 24, 20, 15])
  23 => (bonds = [8, 26, 35, 51], neighbors = [22, 24, 25, 18])
  6  => (bonds = [11, 20, 38, 47], neighbors = [4, 5, 13, 11])
  11 => (bonds = [4, 22, 31, 47], neighbors = [10, 12, 13, 6])
  9  => (bonds = [12, 21, 39, 48], neighbors = [7, 8, 16, 14])
  14 => (bonds = [5, 23, 32, 48], neighbors = [13, 15, 16, 9])
  3  => (bonds = [10, 19, 37, 46], neighbors = [1, 2, 10, 17])
  7  => (bonds = [3, 12, 29, 45], neighbors = [8, 9, 5, 27])
  25 => (bonds = [9, 18, 35, 42], neighbors = [26, 27, 23, 18])
  4  => (bonds = [2, 11, 28, 44], neighbors = [5, 6, 2, 24])
  ⋮  => ⋮
```

## Cubic Lattice and Boundary Conditions

In this section we consider a simple cubic lattice.
As usual, we begin by defining a [`UnitCell`](@ref):

```jldoctest cubic; setup = :(using LatticeUtilities)
julia> cubic = UnitCell(lattice_vecs = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],
                        basis_vecs   = [[0.,0.,0.]])
UnitCell{Float64}:
 • D = 3
 • n = 1
 • lattice_vecs =
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0
 • reciprocal_vecs =
3×3 Matrix{Float64}:
 6.28319  0.0      0.0
 0.0      6.28319  0.0
 0.0      0.0      6.28319
 • basis_vecs =
3×1 Matrix{Float64}:
 0.0
 0.0
 0.0
```

Next we define four instances of the type [`Lattice`](@ref).
The first one represents a lattice that is periodic in the direction
of all three lattice vectors:

```jldoctest cubic
julia> lattice_ppp = Lattice([4,4,4],[true,true,true])
Lattice:
 • D = 3
 • N = 64
 • L = [4, 4, 4]
 • periodic = Bool[1, 1, 1]
```

The second defines a lattice that is periodic only in the direction of the
first two lattice vectors, ``\hat{\mathbf{x}}`` and ``\hat{\mathbf{y}}``:

```jldoctest cubic
julia> lattice_ppo = Lattice([4,4,4],[true,true,false])
Lattice:
 • D = 3
 • N = 64
 • L = [4, 4, 4]
 • periodic = Bool[1, 1, 0]
```

The third defines a lattice that is periodic only in the direction
of the first lattice vector, ``\hat{\mathbf{x}}``:

```jldoctest cubic
julia> lattice_poo = Lattice([4,4,4],[true,false,false])
Lattice:
 • D = 3
 • N = 64
 • L = [4, 4, 4]
 • periodic = Bool[1, 0, 0]
```

Lastly, we define a lattice with strictly open boundary conditions:

```jldoctest cubic
julia> lattice_ooo = Lattice([4,4,4],[false,false,false])
Lattice:
 • D = 3
 • N = 64
 • L = [4, 4, 4]
 • periodic = Bool[0, 0, 0]
```

Note that the set of ``\mathbf{k}``-points generated by the function
[`calc_k_points`](@ref) changes with the boundary conditions.
In particular, if a lattice has open boundary conditions
in the direction of the lattice vector ``\mathbf{a}_i``, then it
treats ``L_i = 1`` when calculating the relevant set of
``\mathbf{k}``-points associated with a lattice:

```jldoctest cubic
julia> size( calc_k_points(cubic,lattice_ppp) )
(3, 4, 4, 4)

julia> size( calc_k_points(cubic,lattice_ppo) )
(3, 4, 4, 1)

julia> size( calc_k_points(cubic,lattice_poo) )
(3, 4, 1, 1)

julia> size( calc_k_points(cubic,lattice_ooo) )
(3, 1, 1, 1)

julia> calc_k_points(cubic,lattice_ooo)
3×1×1×1 Array{Float64, 4}:
[:, :, 1, 1] =
 0.0
 0.0
 0.0
```

As you can see, in the case that the lattice only has open boundary
conditions, then just the ``\mathbf{k} = 0`` point is returned.

Next we will consider how the boundary conditions effect how neighbor tables
are computed. Considering just nearest neighbors, we have three bonds to define:

```jldoctest cubic
julia> bond_x = Bond(orbitals = [1,1], displacement = [1,0,0])
Bond:
 • D = 3
 • orbitals = (1, 1)
 • displacement = [1, 0, 0]

julia> bond_y = Bond([1,1], [0,1,0])
Bond:
 • D = 3
 • orbitals = (1, 1)
 • displacement = [0, 1, 0]

julia> bond_z = Bond([1,1], [0,0,1])
Bond:
 • D = 3
 • orbitals = (1, 1)
 • displacement = [0, 0, 1]
```

Base on these three bonds, the neighbor table generated using [`build_neighbor_table`](@ref)
will depend on the boundary conditions used:

```jldoctest cubic
julia> nt_ppp = build_neighbor_table([bond_x,bond_y,bond_z],
                                      cubic, lattice_ppp)
2×192 Matrix{Int64}:
 1  2  3  4  5  6  7  8   9  10  11  …  56  57  58  59  60  61  62  63  64
 2  3  4  1  6  7  8  5  10  11  12      8   9  10  11  12  13  14  15  16

julia> nt_ppo = build_neighbor_table([bond_x,bond_y,bond_z],
                                      cubic, lattice_ppo)
2×176 Matrix{Int64}:
 1  2  3  4  5  6  7  8   9  10  11  …  40  41  42  43  44  45  46  47  48
 2  3  4  1  6  7  8  5  10  11  12     56  57  58  59  60  61  62  63  64

julia> nt_poo = build_neighbor_table([bond_x,bond_y,bond_z],
                                      cubic, lattice_poo)
2×160 Matrix{Int64}:
 1  2  3  4  5  6  7  8   9  10  11  …  40  41  42  43  44  45  46  47  48
 2  3  4  1  6  7  8  5  10  11  12     56  57  58  59  60  61  62  63  64

julia> nt_ooo = build_neighbor_table([bond_x,bond_y,bond_z],
                                      cubic, lattice_ooo)
2×144 Matrix{Int64}:
 1  2  3  5  6  7   9  10  11  13  14  …  40  41  42  43  44  45  46  47  48
 2  3  4  6  7  8  10  11  12  14  15     56  57  58  59  60  61  62  63  64
```

## Honeycomb Lattice and Bond Operations

In this section we consider the Honeycomb, or Hexagonal, lattice
with just nearest neighbor interactions. The lattice vectors for
the honeycomb lattice are

```math
\begin{align*}
\mathbf{a}_{1} =&	\left(\frac{3}{2},\frac{\sqrt{3}}{2}\right)\\
\mathbf{a}_{2} =&	\left(\frac{3}{2},-\frac{\sqrt{3}}{2}\right),
\end{align*}
```

and basis vectors are

```math
\begin{align*}
\mathbf{r}_{1} =&	\left(0,0\right)\\
\mathbf{r}_{2} =&	\left(1,0\right).
\end{align*}
```

The corresponding reciprocal lattice vectors are

```math
\begin{align*}
\mathbf{b}_{1} =&	2\pi\left(\frac{1}{3},\frac{\sqrt{3}}{3}\right)\\
\mathbf{b}_{2} =&	2\pi\left(\frac{1}{3},-\frac{\sqrt{3}}{3}\right).
\end{align*}
```

As usual, we begin by constructing an instance of the type [`UnitCell`](@ref)
to reflect this:

```jldoctest honeycomb; setup = :(using LatticeUtilities)
julia> honeycomb = UnitCell(lattice_vecs = [[3/2,√3/2],[3/2,-√3/2]],
                            basis_vecs   = [[0.,0.],[1.,0.]])
UnitCell{Float64}:
 • D = 2
 • n = 2
 • lattice_vecs =
2×2 Matrix{Float64}:
 1.5        1.5
 0.866025  -0.866025
 • reciprocal_vecs =
2×2 Matrix{Float64}:
 2.0944   2.0944
 3.6276  -3.6276
 • basis_vecs =
2×2 Matrix{Float64}:
 0.0  1.0
 0.0  0.0
```

And again, we define an instance of the type [`Lattice`](@ref) to
represent a finite lattice:

```jldoctest honeycomb
julia> lattice = Lattice([3,3],[true,true])
Lattice:
 • D = 2
 • N = 9
 • L = [3, 3]
 • periodic = Bool[1, 1]
```

If we consider just nearest neighbor relations, three types of bonds
need to be defined:

```jldoctest honeycomb
julia> bond_1 = Bond(orbitals = [1,2], displacement = [0,0])
Bond:
 • D = 2
 • orbitals = (1, 2)
 • displacement = [0, 0]

julia> bond_2 = Bond([1,2], [-1,0])
Bond:
 • D = 2
 • orbitals = (1, 2)
 • displacement = [-1, 0]

julia> bond_3 = Bond([1,2], [0,-1])
Bond:
 • D = 2
 • orbitals = (1, 2)
 • displacement = [0, -1]
```

However, these three bond definitions are not unique.
For instance, it would be similarly correct to define the nearest neighbor
relations using these three bond definitions:

```jldoctest honeycomb
julia> bond_1′ = Bond(orbitals = [2,1], displacement = [0,0])
Bond:
 • D = 2
 • orbitals = (2, 1)
 • displacement = [0, 0]

julia> bond_2′ = Bond([2,1], [1,0])
Bond:
 • D = 2
 • orbitals = (2, 1)
 • displacement = [1, 0]

julia> bond_3′ = Bond([2,1], [0,1])
Bond:
 • D = 2
 • orbitals = (2, 1)
 • displacement = [0, 1]
```

While the explicit definitions are different here, we can test for
equivalencies.

```jldoctest honeycomb
julia> bond_1 == bond_1′
true

julia> bond_2 == bond_2′
true

julia> bond_3 == bond_3′
true

julia> bond_1 == bond_2′
false
```

Note that while not explicitly enforced, for comparison operations like this
to be valid each bond definition need to be defined assuming shared lattice
vector definitions, which themselves are not unique.

Note that is also possible to [`simplify!`](@ref) a bond, accounting for
the size of a finite lattice and also periodic boundary conditions
where necessary:

```jldoctest honeycomb
julia> bond′ = Bond([1,2],[2,0])
Bond:
 • D = 2
 • orbitals = (1, 2)
 • displacement = [2, 0]

julia> simplify!(bond′,lattice)

julia> bond′
Bond:
 • D = 2
 • orbitals = (1, 2)
 • displacement = [-1, 0]
```

Lastly, given a pair of sites in a finite lattice, it is very easy
to determine the bond definition connecting the two sites using the
[`sites_to_bond`](@ref) function:

```jldoctest honeycomb
julia> bond_1_2 = sites_to_bond(1, 2, honeycomb, lattice)
Bond:
 • D = 2
 • orbitals = (1, 2)
 • displacement = [0, 0]

julia> bond_1_3 = sites_to_bond(1, 3, honeycomb, lattice)
Bond:
 • D = 2
 • orbitals = (1, 1)
 • displacement = [1, 0]
```