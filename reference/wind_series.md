# Create a wind_series

Creates a `wind_series`, a time series of wind fields, from a raster, a
list of rasters, or one or more raster files, such as the monthly files
saved by
[`download_wind_data()`](https://matthewkling.github.io/windscape/reference/download_wind_data.md).
Multiple inputs are combined into one series, in the order given,
without reading the data into memory until needed.

## Usage

``` r
wind_series(x, order = c("uuvv", "uvuv"))
```

## Arguments

- x:

  A multi-layer `SpatRaster` with layers containing u and v wind
  components; a list of them; or a character vector of paths to raster
  files. Each input must contain an equal number of u and v layers,
  arranged as given by `order`, and all inputs must share the same grid.
  Data must be on a longitude/latitude grid with square cells (equal x
  and y resolution in degrees), with u and v components in m/s, oriented
  to true east and north; see
  [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md).

- order:

  Either `"uuvv"`, the default, indicating that each input has all u
  components followed by all v components, or `"uvuv"`, indicating
  alternating u and v components.

## Value

A `wind_series` object, a `SpatRaster` with all the u layers of all
inputs, followed by all the v layers, in time order.

## See also

[`subset_series()`](https://matthewkling.github.io/windscape/reference/subset_series.md)
and
[`wind_times()`](https://matthewkling.github.io/windscape/reference/wind_times.md)
for working with time steps.

## Examples

``` r
series <- windscape_example("wind_series")
series
#> class       : SpatRaster
#> size        : 64, 96, 192  (nrow, ncol, nlyr)
#> resolution  : 0.3157895, 0.3174603  (x, y)
#> extent      : -120.1579, -89.84211, 29.84127, 50.15873  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> sources     : wind_usa.tif
#> names       : u 2000-01-01, u 200~00:00, u 200~00:00, u 200~00:00, u 2000-01-15, u 200~00:00, ...
#> min values  :         -0.4,        -0.4,        -0.6,        -0.6,         -0.4,        -0.5, ...
#> max values  :          0.7,           1,         0.9,           1,          0.7,         1.1, ...

if (FALSE) { # \dontrun{
# combine monthly files from download_wind_data()
series <- wind_series(files)
} # }
```
