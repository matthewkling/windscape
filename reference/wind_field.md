# Create a wind_field

Creates a `wind_field`, the wind across a grid at a single moment, from
a two-layer raster or file, or from one time step of a `wind_series`. To
summarize a whole series as a single field, use
[mean()](https://matthewkling.github.io/windscape/reference/mean-wind_series-method.md)
or
[`net_flow()`](https://matthewkling.github.io/windscape/reference/net_flow.md).

## Usage

``` r
wind_field(x, step = NULL)
```

## Arguments

- x:

  A `SpatRaster` with two layers holding the u and v wind components,
  the path to a raster file with those two layers, or a `wind_series`.
  Data must be on a longitude/latitude grid, with components in m/s,
  oriented to true east and north.

- step:

  For a `wind_series` with more than one time step, the time step to
  use, as an integer index (see
  [`wind_times()`](https://matthewkling.github.io/windscape/reference/wind_times.md)
  for the time of each step).

## Value

A `wind_field` object, which is a particular form of `SpatRaster`.

## See also

[mean()](https://matthewkling.github.io/windscape/reference/mean-wind_series-method.md)
for the time-mean wind of a series.

## Examples

``` r
series <- windscape_example("wind_series")
field <- wind_field(series, step = 1)
```
