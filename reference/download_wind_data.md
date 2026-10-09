# Download hourly wind data from NCAR

Downloads gridded hourly wind data from the NCAR Geoscience Data
Exchange (GDEX) for a bounding box and a set of months, optionally
limited to particular days and hours, saving one GeoTIFF file per month.
No account is needed. Data are clipped to the bounding box and the
requested times on the server, so only the requested data are
transferred.

## Usage

``` r
download_wind_data(
  source = c("era5", "cfsr", "cfsv2"),
  level = "10m",
  xlim,
  ylim,
  years,
  months = 1:12,
  days = NULL,
  hours = NULL,
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

- days:

  Integer vector of days of the month (1-31) to download, by the UTC
  date of each time step. `NULL` (the default) downloads all days.
  Months with none of the requested days (e.g. day 31 in a 30-day month)
  are skipped.

- hours:

  Integer vector of hours of the day (0-23, UTC) to download. `NULL`
  (the default) downloads all hours; for example,
  `hours = seq(0, 21, 3)` downloads every third hour, which reduces
  download time and file size about threefold. Days and hours refer to
  the UTC date and time of each time step, within each requested month.
  (CFSR and CFSv2 monthly files hold hourly forecasts valid from 01:00
  on the first through 00:00 on the first of the following month. With
  `days` or `hours`, time steps are assigned to their calendar month, so
  00:00 on the first comes from the previous month's file; without them,
  each month's file is downloaded as is, as in the pre-built roses of
  [`download_wind_rose()`](https://matthewkling.github.io/windscape/reference/download_wind_rose.md).)

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
(returned invisibly if all files were already cached). Files for a
subset of days or hours have names ending in a tag recording the
selection (e.g. `_d28_h23`), so they are cached separately from whole
months.

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
selecting hours with `hours` (e.g. every third hour) reduces both. The
server limits the size of each request (about 100 MB), so large requests
are split into several, by time, and reassembled.

Requires the ncdf4 package.

## See also

[`wind_series()`](https://matthewkling.github.io/windscape/reference/wind_series.md)
to load the files;
[`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
to summarize them;
[`download_land_mask()`](https://matthewkling.github.io/windscape/reference/download_land_mask.md)
to download a matching land-water layer;
[`download_wind_rose()`](https://matthewkling.github.io/windscape/reference/download_wind_rose.md)
for pre-built CFSR wind roses, which need no wind data downloads.

## Examples

``` r
if (FALSE) { # \dontrun{
# every third hour of 10 m ERA5 wind for summer 2020, Pacific Northwest
files <- download_wind_data("era5", xlim = c(-125, -115), ylim = c(42, 49),
                            years = 2020, months = 6:8, hours = seq(0, 21, 3),
                            dir = "~/wind_data")
rose <- wind_rose(wind_series(files))

# a single hour: Hurricane Katrina, 2005-08-28 23:00 UTC
f <- download_wind_data("era5", xlim = c(-99, -78), ylim = c(17, 35),
                        years = 2005, months = 8, days = 28, hours = 23)
katrina <- wind_field(wind_series(f))
} # }
```
