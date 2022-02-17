# Examples

This page includes several examples on how to use the `LatticeUtilities.jl` package.

## Kagome Lattice Example

As an first example we will represent a Kagome lattice.
The lattice vectors for the Kagome lattice are

```math
\begin{align*}
\mathbf{a}_{1} &= \left(1,0\right)\\
\mathbf{a}_{2} &= \left(\frac{1}{2},\frac{\sqrt{3}}{2}\right),
\end{align*}
```

and with three orbitals per unit cell the basis vectors are

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

The first thing to construct an instance of the type [`UnitCell`](@ref) to represent
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
``\mathbf{r} = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2 + \mathbf{r}_\beta - \mathbf{r}_\alpha``
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
and [`displacement_to_vec!`](@ref), with names ending with `!` that work by modifying
the first array argument in place.