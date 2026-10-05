# Least-cost paths through a wind graph

Traces the least-cost (fastest) paths between sites through a wind
graph, as sequences of grid cell centers. The result has the same
structure as
[`wind_trails()`](https://matthewkling.github.io/windscape/reference/wind_trails.md)
output, so it can be drawn with
[`geom_wind_trail()`](https://matthewkling.github.io/windscape/reference/geom_wind_trail.md).

## Usage

``` r
least_cost_paths(graph, from, to, pairs = c("all", "nearest", "matched"))
```

## Arguments

- graph:

  A `wind_graph`, created with
  [`wind_graph()`](https://matthewkling.github.io/windscape/reference/wind_graph.md).

- from, to:

  Origin and destination sites: two-column matrices (or data frames) of
  longitude and latitude, or `SpatVector`s of points.

- pairs:

  Which origin-destination pairs to trace: `"all"` (the default) traces
  a path from every origin to every destination; `"nearest"` traces one
  path from each origin, to the destination it can reach at the lowest
  cost; `"matched"` pairs `from` and `to` row by row, which requires
  them to have the same number of sites.

## Value

A data frame with one row per path vertex, ordered from origin to
destination along each path:

- `trail`: path id.

- `from`, `to`: row numbers of the path's origin in `from` and
  destination in `to`.

- `step`: vertex number along the path, starting at 0 at the origin.

- `hours`: cumulative travel cost from the origin, in hours (if
  `trans = 1` in
  [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
  and wind speeds are in m/s). Its final value on each path equals the
  pair's
  [`pairwise_least_cost()`](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md)
  (with `snap = TRUE`).

- `x`, `y`: longitude and latitude of the grid cell center.

Pairs with no path (where the destination can't be reached, e.g. because
wind never blows toward it) are omitted with a warning. Pairs whose
origin and destination fall in the same grid cell are omitted, since
their path has no length.

## Details

Paths are computed with
[`gdistance::shortestPath()`](https://AgrDataSci.github.io/gdistance/reference/shortestPath.html).
They follow the graph's direction: with a downwind graph (the default in
[`wind_graph()`](https://matthewkling.github.io/windscape/reference/wind_graph.md)),
a path from `from` to `to` is the fastest downwind route. With an upwind
graph, it is the fastest downwind route from `to` to `from`, traced in
reverse. Because wind graphs are directed, the path from A to B
generally differs from the path from B to A.

Paths move between neighboring grid cells, in eight directions, so they
appear as segments at 45-degree angles, including occasional staircase
patterns where the optimal route lies between two neighbor directions.
This reflects the model's structure rather than the plot. Paths from one
origin often share segments, showing the main transport corridors.

## Examples

``` r
rose <- windscape_example("wind_rose")
graph <- wind_graph(rose)
site <- cbind(-105, 40)
destinations <- cbind(c(-95, -115, -100, -110), c(45, 35, 32, 48))
paths <- least_cost_paths(graph, site, destinations)
head(paths)
#>   trail from to step    hours         x        y
#> 1     1    1  1    0  0.00000 -105.1579 39.84127
#> 2     1    1  1    1  6.39227 -104.8421 39.84127
#> 3     1    1  1    2 27.63580 -104.5263 40.15873
#> 4     1    1  1    3 35.55941 -104.2105 40.15873
#> 5     1    1  1    4 44.14945 -103.8947 40.15873
#> 6     1    1  1    5 53.24461 -103.5789 40.15873

library(ggplot2)
ggplot(paths, aes(x, y)) +
  geom_wind_trail(aes(color = hours)) +
  coord_quickmap()
#> Error in wind_trail_layer(mapping, data, stat, GeomWindTrail, position,     seeds, res, fixed_length, length, hours, match.arg(direction),     steps, match.arg(wrap), arrow, na.rm, show.legend, inherit.aes,     ...): Problem while computing aesthetics.
#> ℹ Error occurred in the 1st layer.
#> Caused by error:
#> ! object 'u' not found
```
