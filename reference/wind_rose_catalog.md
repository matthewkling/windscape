# Catalog of pre-built wind roses

Lists the pre-built global wind roses available for download with
[`download_wind_rose()`](https://matthewkling.github.io/windscape/reference/download_wind_rose.md).
Currently these are built from 10 m CFSR wind for 1979-2010, at CFSR's
native resolution (about 0.31 degrees), with `trans = 1`.

## Usage

``` r
wind_rose_catalog()
```

## Value

A data frame with one row per file, and columns:

- `source`, `level`: wind data set and height (as in
  [`download_wind_data()`](https://matthewkling.github.io/windscape/reference/download_wind_data.md)).

- `tier`: the period the file covers: `"month"` (a single year-month),
  `"year"` (a calendar year), `"month_of_year"` (a calendar month across
  all years), or `"all"` (the full period).

- `year`, `month`: the year and month the file covers, or `NA` where it
  spans several.

- `first`, `last`: the first and last year-months covered.

- `n_steps`: the number of hourly time steps summarized.

- `trans`: the wind speed transformation the rose was built with (see
  [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)).

- `file`, `bytes`, `md5`: file name, size in bytes, and MD5 checksum.

- `release`, `url`: the release the file belongs to, and its download
  URL.

## See also

[`download_wind_rose()`](https://matthewkling.github.io/windscape/reference/download_wind_rose.md)

## Examples

``` r
catalog <- wind_rose_catalog()
table(catalog$tier)
#> 
#>           all         month month_of_year          year 
#>             1           384            12            32 
```
