# Build a wind connectivity graph

Constructs the directed network that the least-cost functions
([`least_cost()`](https://matthewkling.github.io/windscape/reference/least_cost.md),
[`least_cost_paths()`](https://matthewkling.github.io/windscape/reference/least_cost_paths.md),
and
[`pairwise_least_cost()`](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md))
find paths through. Those functions build the graph from a wind rose
automatically, so you don't normally need to call this function.
Building the graph yourself can save time when making many least-cost
calls on a large grid: pass the graph in place of the rose.

## Usage

``` r
wind_graph(x, direction = "downwind", wrap = NULL)
```

## Arguments

- x:

  A `wind_rose`.

- direction:

  Either `"downwind"` (the default) or `"upwind"`: whether links follow
  the wind, or run against it.

- wrap:

  Logical: join the east and west edges of the grid, linking cells
  across them? The default, `NULL`, does so if `x` is a global grid
  spanning all 360 degrees of longitude, where -180 and 180 are the same
  meridian, and not otherwise. `TRUE` or `FALSE` overrides this; `TRUE`
  on a longitude/latitude grid that isn't global gives a warning.

## Value

A `wind_graph`: a gdistance
[Transition-class](https://AgrDataSci.github.io/gdistance/reference/Transition-classes.html)
object, with its direction recorded.

## Details

The graph links each grid cell to its eight neighbors. Each link's
conductance is the mean, over its two cells, of the rose's conductance
in the link's direction. In a downwind graph, links point the way the
wind carries material; an upwind graph is the same network with every
link reversed, for measuring travel toward a site.

## Examples

``` r
rose <- windscape_example("wind_rose")
graph <- wind_graph(rose)
sites <- cbind(c(-110, -100), c(40, 40))
pairwise_least_cost(graph, sites)
#>          [,1]     [,2]
#> [1,]    0.000 299.9684
#> [2,] 1043.173   0.0000
```
