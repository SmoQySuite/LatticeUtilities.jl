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
julia> kagome = UnitCell([[1.0,0.0], [1/2,√3/2]],
                         [[0.0,0.0], [1/2,0.0], [1/4,√3/4]])
UnitCell{Float64}:
 - D = 2
 - n = 3
 - lattice_vecs =
2×2 Matrix{Float64}:
 1.0  0.5
 0.0  0.866025
 - reciprocal_vecs =
2×2 Matrix{Float64}:
  6.28319  0.0
 -3.6276   7.2552
 - basis_vecs =
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
can also be caluculated using the [`displacement_to_vec`](@ref) method:

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
with [`UnitCell`](@ref) instance of like dimenion `D`, to represent a finite lattice.
In this example we will consider a ``3 \times 3`` unit cell Kagome lattice with periodic boundary conditions
in the direction of both lattice vectors:

```jldoctest kagome
julia> lattice = Lattice([3,3], [true,true])
Lattice:
 - D = 2
 - N = 9
 - L = [3, 3]
 - periodic = Bool[1, 1]
```

Given an initial site, a displacement in unit cells, and a terminating orbital species,
we can calculate the final site using the [`site_to_site`](@ref) method:

```jldoctest kagome
julia> site_to_site(17, [1,1], 2, kagome, lattice)
20
```

We can also calculate the k-points associated with out finite lattice using the [`calc_k_points`](@ref) method:

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

Individual k-points can be caluclated using the [`calc_k_point`](@ref) and [`calc_k_point!`](@ref) methods.

If we concern ourselves with just nearest neighbors, there are
six bond defintions that need to be defined:

```jldoctest kagome
julia> bond_1 = Bond([1,2],[0,0])
Bond:
 - D  = 2
 - o  = [1, 2]
 - Δl = [0, 0]

julia> bond_2 = Bond([1,3],[0,0])
Bond:
 - D  = 2
 - o  = [1, 3]
 - Δl = [0, 0]

julia> bond_3 = Bond([2,3],[0,0])
Bond:
 - D  = 2
 - o  = [2, 3]
 - Δl = [0, 0]

julia> bond_4 = Bond([2,1],[1,0])
Bond:
 - D  = 2
 - o  = [2, 1]
 - Δl = [1, 0]

julia> bond_5 = Bond([3,1],[0,1])
Bond:
 - D  = 2
 - o  = [3, 1]
 - Δl = [0, 1]

julia> bond_6 = Bond([3,2],[-1,1])
Bond:
 - D  = 2
 - o  = [3, 2]
 - Δl = [-1, 1]
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

In certain situations it is useful to sort `neighbor_table` so that the first row
is in strictly ascending order, and for constant value in the first row the second
row is also ascending. This sorting can be performed using the [`sort_neighbor_table!`](@ref)
method:

```jldoctest kagome
julia> inv_perm = sort_neighbor_table!(neighbor_table); neighbor_table
2×54 Matrix{Int64}:
 1  1  1   1  2  2   2   3   3  4  4  …  20  20  22  22  23  23  25  25  26
 2  3  8  21  3  4  24  10  17  5  6     21  22  23  24  24  25  26  27  27
```