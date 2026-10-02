# Build a wind connectivity graph

This function constructs a landscape connectivity graph from a windrose
raster stack.

## Usage

``` r
wind_graph(x, direction = "downwind", wrap = FALSE)
```

## Arguments

- x:

  An object of class \`wind_rose\`.

- direction:

  Either "downwind" (the default) or "upwind", indicating whether
  outbound or inbound wind conductance should be computed, respectively.

- wrap:

  Should the left and right edges of the raster be connected? Default is
  FALSE. Set to TRUE if, for example, \`x\` is a global raster where
  -180 and 180 are equivalent longitudes.

## Value

a \`wind_graph\` object consisting of a gdistance
[Transition-class](https://AgrDataSci.github.io/gdistance/reference/Transition-classes.html)
object and additional metadata
