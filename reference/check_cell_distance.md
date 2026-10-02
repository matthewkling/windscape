# Compare point distances to cell distances

Check for discrepancies between pairwise distances between a set of
points, versus pairwise distances between the centers of the raster grid
cells that the points fall within. There will always be some differences
due to the continuous versus discrete nature of the data types, but if
discrepancies are large as a percentage of the distances, then they will
contribute substantial noise that could invalidate wind connectivity
estimates. This function reports how many site pairs are in the same
grid cell (an infinite difference that makes wind calculations
impossible) as well as information about the distribution of
discrepancies as percentages of distance. The more site pairs there are
with nontrivial discrepancies, the less reliable a wind connectivity
model will be for these sites. Problems identified in this check can be
solved using the
[downscale](https://matthewkling.github.io/windscape/reference/downscale.md)
or [vrcd](https://matthewkling.github.io/windscape/reference/vrcd.md)
functions.

## Usage

``` r
check_cell_distance(x, ll, return = FALSE)
```

## Arguments

- x:

  SpatRaster (e.g. a \`wind_rose\`)

- ll:

  two-column matrix of site coordinates

- return:

  Logical: return the matrix of ratios? Default is FALSE, which only
  prints a report.

## Value

If `return = TRUE`, a matrix of the ratios of cell distances to point
distances
