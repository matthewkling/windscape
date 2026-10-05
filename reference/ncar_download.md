# Download hourly wind data from NCAR

Downloads gridded hourly wind data from the NCAR Geoscience Data
Exchange (GDEX) for a bounding box and a set of months, saving one
GeoTIFF file per month. No account is needed. Data are clipped to the
bounding box on the server, so only the requested region is transferred.

## Usage

``` r
ncar_download(
  source = c("era5", "cfsr", "cfsv2"),
  level = "10m",
  xlim,
  ylim,
  years,
  months = 1:12,
  time_stride = 1,
  dir = tempdir(),
  overwrite = FALSE,
  quiet = FALSE
)
```

## Arguments

- source:

  Data set to download:

  - `"era5"` (the default): ECMWF ERA5 reanalysis, 0.25 degree grid,
    1940 to present (NCAR data set d633000).

  - `"cfsr"`: NCEP Climate Forecast System Reanalysis, ~0.31 degree
    Gaussian grid, 1979 to 2010 (d093001).

  - `"cfsv2"`: NCEP Climate Forecast System version 2, the operational
    continuation of CFSR, ~0.2 degree Gaussian grid, 2011 to present
    (d094001).

- level:

  Height of the wind data. `"10m"` (the default) is wind 10 m above the
  ground, available for all sources. `"100m"` is available for ERA5
  only. Pressure levels `"1000hPa"`, `"850hPa"`, `"700hPa"`, `"500hPa"`,
  and `"200hPa"` are available for CFSR and CFSv2 only.

- xlim:

  Numeric vector of length 2 giving the western and eastern limits of
  the bounding box, in degrees. Longitudes can be given either in the
  -180 to 180 range or the 0 to 360 range, and output files use the same
  convention. Boxes crossing the prime meridian (in -180 to 180
  coordinates, e.g. `c(-10, 30)`) or the antimeridian (in 0 to 360
  coordinates, e.g. `c(170, 200)`) are supported.

- ylim:

  Numeric vector of length 2 giving the southern and northern limits of
  the bounding box, in degrees, in the range -90 to 90.

- years:

  Integer vector of years to download.

- months:

  Integer vector of months (1-12) to download within each year. Defaults
  to all months.

- time_stride:

  Integer. Download every `time_stride`-th hourly time step, starting
  with each month's first time step. The default, 1, downloads all
  hours; for example, 3 downloads every third hour and reduces download
  time and file size about threefold.

- dir:

  Directory where monthly files are saved. Defaults to a temporary
  directory, which is deleted at the end of the R session; set it to a
  permanent location to keep the data and make use of caching across
  sessions.

- overwrite:

  Logical. If `FALSE` (the default), months already downloaded to `dir`
  are skipped.

- quiet:

  Logical. If `TRUE`, progress messages are suppressed.

## Value

A character vector of file paths, one per month, in chronological order
(returned invisibly if all files were already cached).

## Details

Each output file holds one month of data in `wind_series` layout: all u
layers followed by all v layers, named like `"u 2005-08-28 23:00:00"`
(UTC), in m/s, oriented to true east and north, on a longitude/latitude
grid. Combine the files into one series with
[`wind_series()`](https://matthewkling.github.io/windscape/reference/wind_series.md),
which keeps the data on disk, and summarize it with
[`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md).

Downloads are cached: a month whose file already exists in `dir` is not
downloaded again unless `overwrite = TRUE`, so an interrupted download
can be resumed by rerunning the same call. Large requests (many years,
or a large region) can take a long time and use substantial disk space;
thinning hourly data with `time_stride` reduces both.

Requires the ncdf4 package.

## See also

[`wind_series()`](https://matthewkling.github.io/windscape/reference/wind_series.md)
to load the files;
[`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
to summarize them;
[`ncar_land()`](https://matthewkling.github.io/windscape/reference/ncar_land.md)
to download a matching land-water layer.

## Examples

``` r
if (FALSE) { # \dontrun{
# every third hour of 10 m ERA5 wind for summer 2020, Pacific Northwest
files <- ncar_download("era5", xlim = c(-125, -115), ylim = c(42, 49),
                       years = 2020, months = 6:8, time_stride = 3,
                       dir = "~/wind_data")
ws <- wind_series(files)
rose <- wind_rose(files)
} # }
```
