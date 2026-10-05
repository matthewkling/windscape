# Pairwise random walk connectivity among sites

For every pair of sites, computes how strongly particles released at one
site reach the other under a stream-mode random walk (see
[`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md)):
by default, the probability that a particle released at the first site
is deposited at the second. This is the random walk counterpart to
[`pairwise_least_cost()`](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md).

## Usage

``` r
pairwise_random_walk(
  rose,
  sites,
  half_life = NULL,
  value = c("deposition", "residence"),
  density = TRUE,
  timescale = 1,
  latitude_correction = TRUE,
  chunk = 200
)
```

## Arguments

- rose:

  A `wind_rose`.

- sites:

  A two-column matrix (or data frame) of site coordinates, or a
  `SpatVector` of points.

- half_life:

  Half-life of airborne mass, in hours of transport time; see
  [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md).
  Required for `value = "deposition"`, which is zero without deposition.

- value:

  `"deposition"` (the default) or `"residence"`; see Value.

- density:

  Logical: divide by the area of the destination site's grid cell,
  giving values per km^2? Default `TRUE`. Without this, values depend on
  cell size: a larger cell receives more of the particles passing
  nearby. On a longitude/latitude grid, cell area shrinks toward the
  poles (by about 25 percent from 30 to 50 degrees latitude), so raw
  values are biased toward lower-latitude destinations even within a
  single analysis. Dividing by area removes that bias, and also makes
  values comparable across grid resolutions (though see Details).

- timescale, latitude_correction:

  See
  [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md).
  Results do not depend on `timescale`.

- chunk:

  Number of sites to solve for at once. Larger values are faster but use
  more memory.

## Value

A square matrix with one row and column per site, giving connectivity
from the row's site (the origin) to the column's site (the destination).
Large values mean strong connectivity, the opposite orientation from the
travel times of
[`pairwise_least_cost()`](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md).
With `value = "deposition"`, entries are the probability that a particle
released at the origin is deposited in the destination's grid cell (per
km^2 if `density = TRUE`). With `value = "residence"`, they are the
time, in hours, that a unit release at the origin spends airborne over
the destination's cell (per km^2 if `density = TRUE`). Because wind is
directional, the matrix is generally asymmetric. The diagonal is each
site's self-connectivity: release that is deposited (or stays airborne)
in its own cell, which can be large; see
[`rw_self_retention()`](https://matthewkling.github.io/windscape/reference/rw_self_retention.md).
Sites in the same grid cell get identical rows and columns.

## Details

Each origin's row is the destination values of a stream-mode random walk
released from that origin, so it equals the result of
`random_walk(rose, origin, mode = "stream", ...)` read at each
destination. All origins share one matrix factorization, so the cost
grows slowly with the number of sites.

`density = TRUE` removes the direct effect of cell size on how much a
destination receives, but the random walk's spread itself depends
somewhat on grid resolution (see
[`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md)),
so values are not fully resolution-independent.

## Examples

``` r
rose <- windscape_example("wind_rose")
sites <- cbind(c(-110, -105, -100, -95), c(40, 42, 38, 44))
pairwise_random_walk(rose, sites, half_life = 24)
#>              [,1]         [,2]         [,3]         [,4]
#> [1,] 1.249959e-04 9.325450e-07 7.104573e-10 1.513116e-10
#> [2,] 3.550546e-15 5.577307e-05 4.317434e-09 2.141014e-09
#> [3,] 4.053452e-17 6.038675e-11 7.086970e-05 1.498152e-08
#> [4,] 3.191795e-21 1.232927e-12 2.449646e-12 7.976103e-05
```
