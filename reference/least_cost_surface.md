# Accumulated wind cost surface

Calculate the accumulated wind cost-distance (e.g. travel times) or flow
rate (the inverse of cost-distance) from one or more sites to every grid
cell across the domain, using the least cost path algorithm. This
function is a wrapper around
[accCost](https://AgrDataSci.github.io/gdistance/reference/accCost.html).

## Usage

``` r
least_cost_surface(graph, sites, rate = FALSE)
```

## Arguments

- graph:

  A
  [wind_graph](https://matthewkling.github.io/windscape/reference/wind_graph.md).

- sites:

  A two-column matrix of point coordinates.

- rate:

  Whether to return values as "rates" instead of the default "cost
  distances". Rates are the inverse of cost distances, representing flow
  rather than travel time.

## Value

A single-layer SpatRaster. With `rate = FALSE`, the layer is named
`hours` and gives the travel cost from the nearest site to each cell, in
hours if `trans = 1` in
[`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
and wind speeds are in m/s (otherwise in relative units). With
`rate = TRUE`, it is named `rate` and gives the inverse.

## Details

Paths are restricted to the eight neighbor directions, so cost distances
are overestimated for routes between neighbor bearings; see
[pairwise_least_cost](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md)
for magnitudes.
