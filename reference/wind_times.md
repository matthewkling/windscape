# Get the time of each step in a wind_series

Returns the time of each step of a `wind_series`, parsed from its layer
names (like `"u 2000-01-01 06:00:00"`, as written by
[`download_wind_data()`](https://matthewkling.github.io/windscape/reference/download_wind_data.md))
or, if the names hold no times, from the times set with
[`terra::time()`](https://rspatial.github.io/terra/reference/time.html).

## Usage

``` r
wind_times(x)
```

## Arguments

- x:

  A `wind_series`.

## Value

A POSIXct vector (UTC), with one element per time step.

## See also

[`subset_series()`](https://matthewkling.github.io/windscape/reference/subset_series.md)
to select time steps.

## Examples

``` r
series <- windscape_example("wind_series")
head(wind_times(series))
#> [1] "2000-01-01 00:00:00 UTC" "2000-01-01 06:00:00 UTC"
#> [3] "2000-01-01 12:00:00 UTC" "2000-01-01 18:00:00 UTC"
#> [5] "2000-01-15 00:00:00 UTC" "2000-01-15 06:00:00 UTC"
```
