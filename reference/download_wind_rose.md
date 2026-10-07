# Download a pre-built wind rose

Downloads a pre-built global wind rose for any combination of years and
calendar months, optionally cropped to a region, so that common analyses
need no wind data downloads or rose building. Roses are built from
hourly 10 m CFSR wind for 1979-2010 with `trans = 1`; see
[`wind_rose_catalog()`](https://matthewkling.github.io/windscape/reference/wind_rose_catalog.md)
for what is available. For other data sets, heights, periods, or `trans`
functions, download wind data with
[`download_wind_data()`](https://matthewkling.github.io/windscape/reference/download_wind_data.md)
and build a rose with
[`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md).

## Usage

``` r
download_wind_rose(
  year = NULL,
  month = NULL,
  ext = NULL,
  source = "cfsr",
  level = "10m",
  cache = is.null(ext),
  quiet = FALSE
)
```

## Arguments

- year:

  Integer vector of years to include. `NULL` (the default) includes all
  years.

- month:

  Integer vector of calendar months (1-12) to include in each year.
  `NULL` (the default) includes all months. For example,
  `year = 1990:1999, month = 6:8` gives a rose for the summers of the
  1990s.

- ext:

  Optional region to crop to: a `SpatExtent`, an object
  [`terra::ext()`](https://rspatial.github.io/terra/reference/ext.html)
  accepts (e.g. a `SpatRaster` or `SpatVector`), or a numeric vector
  `c(xmin, xmax, ymin, ymax)`, in degrees, with longitudes from -180
  to 180. The result includes every cell overlapping `ext`. Regions
  crossing the antimeridian are not supported; for those, use the global
  rose.

- source, level:

  Wind data set and height. Currently only `"cfsr"` and `"10m"`.

- cache:

  Logical. If `TRUE`, whole files are downloaded and kept in a cache on
  disk (see
  [`wind_rose_cache()`](https://matthewkling.github.io/windscape/reference/wind_rose_cache.md)),
  so later requests for them need no download; each file is about 15 MB.
  If `FALSE`, nothing is kept: when `ext` is given, only the cells
  within it are read over the internet, which is much faster than
  downloading whole files for small regions; otherwise files are
  downloaded to a temporary directory. Defaults to `TRUE` when `ext` is
  `NULL` and `FALSE` otherwise.

- quiet:

  Logical. If `TRUE`, progress messages are suppressed.

## Value

A `wind_rose`, with `n_steps` and `trans` recorded so it can be combined
with other roses. Values are wind conductance in 1 / hours.

## Details

Pre-built roses are stored for single months, calendar years, calendar
months across all years, and the full period. A request is met with the
fewest of these, combined with
[`combine_roses()`](https://matthewkling.github.io/windscape/reference/combine_roses.md),
which weights them by their numbers of hours; the result equals the rose
built from all of the requested hours at once.

## See also

[`wind_rose_catalog()`](https://matthewkling.github.io/windscape/reference/wind_rose_catalog.md)
for the available files;
[`wind_rose_cache()`](https://matthewkling.github.io/windscape/reference/wind_rose_cache.md)
to manage the cache.

## Examples

``` r
# \donttest{
# long-term rose for the western US, read remotely
rose <- download_wind_rose(ext = c(-125, -100, 30, 50))
#> reading 1 wind rose file remotely

# summers of the 1990s
summer <- download_wind_rose(year = 1990:1999, month = 6:8, ext = c(-125, -100, 30, 50))
#> reading 30 wind rose files remotely
# }
```
