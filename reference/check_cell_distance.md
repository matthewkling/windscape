# Check how grid cells distort distances among sites

Compares pairwise distances among sites with distances between the
centers of the grid cells the sites fall in, and prints a report.
Connectivity models that work from cell to cell treat each site as the
center of its cell: the random walk functions
([`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md),
[`pairwise_random_walk()`](https://matthewkling.github.io/windscape/reference/pairwise_random_walk.md)),
and the least-cost functions when sites are snapped to cell centers
([`least_cost_surface()`](https://matthewkling.github.io/windscape/reference/least_cost_surface.md),
[`least_cost_paths()`](https://matthewkling.github.io/windscape/reference/least_cost_paths.md),
and `pairwise_least_cost(snap = TRUE)`). For these, sites separated by
only a few cells have distorted distances and directions, and sites in
the same cell can't be distinguished at all.
[`pairwise_least_cost()`](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md)
with its default `snap = FALSE` places sites at their actual locations,
so this check does not apply to it.

## Usage

``` r
check_cell_distance(x, ll, return = FALSE)
```

## Arguments

- x:

  A SpatRaster on the grid used for connectivity modeling, e.g. a
  `wind_rose`.

- ll:

  A two-column matrix of site coordinates.

- return:

  Logical: return the matrix of ratios? Default `FALSE`, which only
  prints the report.

## Value

Prints the number of site pairs, the number (and percentage) in the same
grid cell, and the distribution of discrepancies between cell and point
distances, as percentages of point distance. If `return = TRUE`, also
returns a matrix of the ratios of cell distances to point distances
(`NaN` for a site with itself, and 0 for distinct sites in the same
cell).

## Details

Where many site pairs are affected, options are to use wind data on a
finer grid, or to
[`downscale()`](https://matthewkling.github.io/windscape/reference/downscale.md)
the wind rose (see its documentation for how downscaling changes random
walk results).
