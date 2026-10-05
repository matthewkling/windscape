# Selecting and combining layers of windscape objects

`wind_series`, `wind_field`, and `wind_rose` objects are `SpatRaster`s
whose layers have a fixed structure (u and v layers, or conductance
toward eight neighbors). Selecting layers with `[[` or
[`terra::subset()`](https://rspatial.github.io/terra/reference/subset.html),
or combining objects with [`c()`](https://rdrr.io/r/base/c.html),
generally breaks that structure, so these operations return a plain
`SpatRaster`. To select time steps while keeping a series intact, use
[`subset_series()`](https://matthewkling.github.io/windscape/reference/subset_series.md);
to take one time step as a field, use
[`wind_field()`](https://matthewkling.github.io/windscape/reference/wind_field.md);
and to combine series, use
[`wind_series()`](https://matthewkling.github.io/windscape/reference/wind_series.md)
with a list of series or files.

## Arguments

- x:

  A `wind_series`, `wind_field`, or `wind_rose`.

- i, subset:

  Layers to select, as for a `SpatRaster`.

- j:

  Not used.

- ...:

  Other arguments passed to the `SpatRaster` method, or for
  [`c()`](https://rdrr.io/r/base/c.html), other objects to combine.

## Value

A `SpatRaster`.

## Examples

``` r
series <- windscape_example("wind_series")
class(series[[1:2]]) # a plain SpatRaster
#> [1] "SpatRaster"
#> attr(,"package")
#> [1] "terra"
class(subset_series(series, steps = 1:2)) # a wind_series
#> [1] "wind_series"
#> attr(,"package")
#> [1] "windscape"
```
