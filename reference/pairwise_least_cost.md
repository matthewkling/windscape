# Pairwise least-cost travel times among sites

Calculate pairwise wind cost-distances (e.g. travel times) or flow rates
(the inverse of cost-distances) among a set of sites, using the least
cost path algorithm. For the random walk counterpart, see
[`pairwise_random_walk()`](https://matthewkling.github.io/windscape/reference/pairwise_random_walk.md).

## Usage

``` r
pairwise_least_cost(graph, sites, snap = FALSE, rate = FALSE)
```

## Arguments

- graph:

  A
  [wind_graph](https://matthewkling.github.io/windscape/reference/wind_graph.md).

- sites:

  A two-column matrix of point coordinates.

- snap:

  Logical: snap sites to the centers of their grid cells? Default
  `FALSE`, which places sites at their actual locations within cells;
  see details.

- rate:

  Whether to return values as "rates" instead of the default "cost
  distances". Rates are the inverse of cost distances, representing flow
  rather than travel time.

## Value

A square matrix with one row and column per site: the least-cost travel
time (or, with `rate = TRUE`, its inverse) from the row's site to the
column's site, in hours if `trans = 1` in
[`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
and wind speeds are in m/s. Small travel times mean strong connectivity.
Because wind graphs are directed, the matrix is generally asymmetric.
The diagonal is zero. Sites outside the graph's extent get `NA`, with a
warning.

## Details

Least-cost travel on a wind graph moves between neighboring cell
centers, so a grid alone places each site at the center of its cell. For
sites a few cells apart or less, that distorts both the distance and the
direction between them, and sites in the same cell would be zero hours
apart. By default (`snap = FALSE`), sites are instead added to the graph
at their actual locations. Each site is linked to the centers of its own
and the eight surrounding cells, and directly to any other site in those
cells, by edges whose travel time is computed exactly for the wind where
the edge starts, treated as uniform along the edge. That wind is the
cell's own at a cell center, and is interpolated between cell centers at
a site. Paths then run from site to site through these edges and the
grid. Results change continuously as sites move, rather than jumping at
cell boundaries. With `snap = TRUE`, sites are snapped to cell centers,
as in
[`gdistance::costDistance()`](https://AgrDataSci.github.io/gdistance/reference/costDistance-methods.html).

The exact travel time through uniform wind is the continuum limit of
least-cost travel on an eight-neighbor grid: the time to cover a
displacement using the cheapest combination of the cell's eight flow
vectors (conductance times the displacement to each neighbor; see
[`net_flow()`](https://matthewkling.github.io/windscape/reference/net_flow.md)),
which uses at most two of them. Equivalently, the region reachable in
one hour is the convex hull of the flow vectors. Site edges therefore
share the grid's metric: for sites at cell centers, results are close to
those from `snap = TRUE`, and slightly lower where a site edge combines
two directions more efficiently than the grid can.

Paths are restricted to the eight neighbor directions, so cost distances
are overestimated for routes between neighbor bearings. For uniform wind
on square cells the maximum error is about 8 percent (routes 22.5
degrees off a grid axis). On a longitude/latitude grid, cells narrow
east-west toward the poles and neighbor bearings become uneven, so the
maximum error grows with latitude, to roughly 18 percent at 60 degrees.
Site edges have the same directional bias, so that it is consistent
across distances.
