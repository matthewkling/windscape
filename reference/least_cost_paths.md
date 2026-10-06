# Least-cost paths

Traces the least-cost (fastest) wind paths between sites, as sequences
of grid cell centers. The result has the same structure as
[`wind_trails()`](https://matthewkling.github.io/windscape/reference/wind_trails.md)
output, so it can be drawn with
[`geom_wind_path()`](https://matthewkling.github.io/windscape/reference/geom_wind_path.md).
Paths run between `sites` and the points in `to`: downwind, from the
sites to the points, or upwind, from the points to the sites. Without
`to`, the points are a regular grid across the domain, so the paths show
the network of fastest routes downwind (or upwind) of the sites.

## Usage

``` r
least_cost_paths(
  rose,
  sites,
  to = NULL,
  direction = "downwind",
  pairs = c("all", "nearest", "matched"),
  n = 50,
  ...
)
```

## Arguments

- rose:

  A
  [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md),
  or a `wind_graph` built in advance whose direction matches
  `direction`; see
  [`pairwise_least_cost()`](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md).

- sites, to:

  Sites, and the points to trace paths to (downwind) or from (upwind):
  two-column matrices (or data frames) of longitude and latitude, or
  `SpatVector`s of points. If `to` is `NULL` (the default), the points
  are about `n` points on a regular grid across the domain, and points
  that can't be reached are dropped silently.

- direction:

  Either `"downwind"` (the default), for paths from `sites` to `to`, or
  `"upwind"`, for paths from `to` to `sites`. With explicit `to` points,
  the two give the same paths for swapped arguments; `direction` matters
  most with the grid of points and with `pairs = "nearest"`.

- pairs:

  Which pairs to trace: `"all"` (the default) traces a path between
  every site and every point in `to`; `"nearest"` traces one path for
  each site, to or from the point in `to` with the lowest travel time;
  `"matched"` pairs `sites` and `to` row by row, which requires them to
  have the same number of rows; it isn't available without `to`.

- n:

  Approximate number of grid points, when `to` is `NULL`.

- ...:

  Further arguments passed to
  [`wind_graph()`](https://matthewkling.github.io/windscape/reference/wind_graph.md),
  such as `wrap`. Not allowed when `rose` is already a `wind_graph`.

## Value

A data frame with one row per path vertex, ordered along each path in
the direction of travel (from the site to the point for downwind paths,
and from the point to the site for upwind paths):

- `trail`: line id, for drawing. Each path is one trail, unless it
  crosses the east-west seam of a wrapped global grid (see `wrap` in
  [`wind_graph()`](https://matthewkling.github.io/windscape/reference/wind_graph.md)),
  where a new trail starts so that lines don't cross the map.

- `site`, `to`: row numbers of the path's site in `sites` and point in
  `to`, which together identify the path.

- `step`: vertex number along the path, starting at 0 at its upwind end.

- `hours`: cumulative travel time along the path, in hours (if
  `trans = 1` in
  [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
  and wind speeds are in m/s). Its final value on each path equals the
  pair's travel time from
  [`pairwise_least_cost()`](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md)
  (with `snap = TRUE`).

- `x`, `y`: longitude and latitude of the grid cell center.

Pairs with no path (e.g. because wind never blows from one toward the
other), and sites and points outside the grid, are omitted with a
warning. Pairs whose site and point fall in the same grid cell are
omitted, since their path has no length.

## Details

Paths are computed with
[`gdistance::shortestPath()`](https://AgrDataSci.github.io/gdistance/reference/shortestPath.html).
Because wind is directional, the fastest path from A to B generally
differs from the fastest path from B to A.

Paths move between neighboring grid cells, in eight directions, so they
appear as segments at 45-degree angles, including occasional staircase
patterns where the optimal route lies between two neighbor directions.
This reflects the model's structure rather than the plot. Paths from one
site form a tree: once two paths meet, they share the rest of their
route back to the site, so shared segments show the main transport
corridors. Because diagonal steps between cell centers can cross without
passing through a common cell, two branches of the tree occasionally
appear to cross.

## Examples

``` r
rose <- windscape_example("wind_rose")
site <- cbind(-105, 40)
destinations <- cbind(c(-95, -115, -100, -110), c(45, 35, 32, 48))
paths <- least_cost_paths(rose, site, destinations)
head(paths)
#>   trail site to step    hours         x        y
#> 1     1    1  1    0  0.00000 -105.1579 39.84127
#> 2     1    1  1    1  6.39227 -104.8421 39.84127
#> 3     1    1  1    2 27.63580 -104.5263 40.15873
#> 4     1    1  1    3 35.55941 -104.2105 40.15873
#> 5     1    1  1    4 44.14945 -103.8947 40.15873
#> 6     1    1  1    5 53.24461 -103.5789 40.15873

library(ggplot2)
ggplot(paths, aes(x, y)) +
  geom_wind_path(aes(color = hours)) +
  coord_quickmap()


# the network of fastest routes from the site, and to it
down <- least_cost_paths(rose, site, n = 200)
up <- least_cost_paths(rose, site, n = 200, direction = "upwind")
ggplot(rbind(cbind(down, direction = "downwind"), cbind(up, direction = "upwind")),
       aes(x, y)) +
  geom_wind_path(aes(color = hours), arrow = NULL) +
  facet_wrap(~direction) +
  coord_quickmap()
```
