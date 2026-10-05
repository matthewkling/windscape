# Self-contribution of each source cell in a stream-mode random walk

Computes `G_cc`, the stream-mode residence time in cell `c` per unit of
release from cell `c` (see
[`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md)).
This is each source's contribution to its own value, which can be
subtracted for a leave-one-out surface, e.g. to avoid a site predicting
itself.

## Usage

``` r
rw_self_retention(
  rose,
  half_life = Inf,
  timescale = 1,
  latitude_correction = TRUE,
  cells = NULL,
  exact = FALSE,
  chunk = 200
)
```

## Arguments

- rose:

  A `wind_rose`.

- half_life:

  Half-life of airborne mass, in hours. Must match the value used for
  the stream-mode walk; see
  [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md).

- timescale:

  Time step scaling factor; see
  [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md).
  Affects only the approximate values; exact values are independent of
  `timescale`.

- latitude_correction:

  Logical. Must match the value used for the stream-mode walk; see
  [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md).

- cells:

  Cell numbers, or a two-column matrix of coordinates, at which to
  compute `G_cc`. Optional for the approximation (which defaults to
  every cell); required if `exact = TRUE`.

- exact:

  Logical: compute exact values rather than the approximation? Exact
  values are a sparse solve per cell (sharing one factorization),
  practical for a set of occupied cells but not for every cell of a
  large grid. The approximation is a lower bound. Default `FALSE`.

- chunk:

  With `exact = TRUE`: number of cells to solve for at once. Larger
  values are faster but use more memory.

## Value

`G_cc` in hours per unit of release: a SpatRaster of approximate values
for every cell if `cells` is NULL, otherwise a numeric vector with one
value per cell. For a leave-one-out surface,
`residence_loo = residence - init * G_cc` and
`deposition_loo = log(2) / half_life * residence_loo`, using the
`residence` layer from
[`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md)
and `init` values at the source cells.

## Details

The approximation `G_cc ~ (1 - lambda) * t / (1 - (1 - lambda) * p_cc)`
counts only paths in which mass never leaves the cell (`p_cc` is the
retention probability; see
[`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md)).
It omits mass that leaves and returns, which is substantial in a
nearest-neighbor walk, especially in fast cells where `p_cc` is near
zero.
