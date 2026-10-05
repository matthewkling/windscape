# Pairwise least-cost travel times among sites

Calculate pairwise wind cost-distances (e.g. travel times) or flow rates
(the inverse of cost-distances) among a set of sites, using the least
cost path algorithm. This function is a wrapper around
[costDistance](https://AgrDataSci.github.io/gdistance/reference/costDistance-methods.html).
For the random walk counterpart, see
[`pairwise_random_walk()`](https://matthewkling.github.io/windscape/reference/pairwise_random_walk.md).

## Usage

``` r
pairwise_least_cost(graph, sites, adjust = TRUE, rate = FALSE)
```

## Arguments

- graph:

  A
  [wind_graph](https://matthewkling.github.io/windscape/reference/wind_graph.md).

- sites:

  A two-column matrix of point coordinates.

- adjust:

  Whether to scale results to correct for discrepancies between
  point-to-point distances and cell-to-cell distances. Default is TRUE.

- rate:

  Whether to return values as "rates" instead of the default "cost
  distances". Rates are the inverse of cost distances, representing flow
  rather than travel time.

## Value

A square matrix with one row and column per site: the least-cost travel
time (or, with `rate = TRUE`, its inverse) from the row's site to the
column's site. Small travel times mean strong connectivity. Because wind
graphs are directed, the matrix is generally asymmetric. The diagonal is
zero.

## Details

Paths are restricted to the eight neighbor directions, so cost distances
are overestimated for routes between neighbor bearings. For uniform wind
on square cells the maximum error is about 8 percent (routes 22.5
degrees off a grid axis). On a longitude/latitude grid, cells narrow
east-west toward the poles and neighbor bearings become uneven, so the
maximum error grows with latitude, to roughly 18 percent at 60 degrees.
The `adjust` option corrects for point-versus-cell distance
discrepancies, not for this directional bias.
