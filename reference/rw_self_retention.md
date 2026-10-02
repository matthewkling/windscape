# Self-contribution of each source cell in a stream-mode random walk

Computes G_cc, the steady-state airborne mass in cell c per unit of
per-step release from cell c, i.e. the diagonal of G = (I - (1 - lambda)
P')^-1. Multiplying by a cell's release gives its contribution to its
own value, which can be subtracted for a leave-one-out surface:
`n_loo = n - n0 * G_cc`.

## Usage

``` r
rw_self_retention(
  rose,
  half_life = Inf,
  cells = NULL,
  exact = FALSE,
  timescale = 1,
  chunk = 200
)
```

## Arguments

- rose:

  A `wind_rose`.

- half_life:

  Half-life of airborne mass, in hours; see
  [random_walk](https://matthewkling.github.io/windscape/reference/random_walk.md).
  Must match the value used for the stream-mode walk.

- cells:

  Optional cell numbers, or a two-column coordinate matrix. Required if
  `exact = TRUE`.

- exact:

  Logical: compute exact values rather than the diagonal approximation?

- timescale:

  Time step scaling factor; see
  [random_walk](https://matthewkling.github.io/windscape/reference/random_walk.md).
  Must match the value used for the stream-mode walk.

- chunk:

  Number of right-hand sides per solve when `exact = TRUE`.

## Value

If `cells` is NULL, a SpatRaster of approximate G_cc. Otherwise a
numeric vector with one value per cell.

## Details

The approximation `G_cc ~ 1 / (1 - (1 - lambda) p_cc)` counts only paths
in which mass never leaves the cell, and is therefore a lower bound. It
omits mass that leaves and returns, which is substantial in a
nearest-neighbor walk, especially in fast cells where the retention
probability p_cc is near zero. The exact value requires one sparse solve
per cell (sharing a single factorization), so is practical for a set of
occupied cells but not for every cell of a large grid.
