# Example wind data sets

Load one of the example data sets shipped with windscape, as a
ready-to-use object.

## Usage

``` r
windscape_example(name = c("wind_series", "wind_rose", "wind_field"))
```

## Source

Climate Forecast System Reanalysis (Saha et al. 2010), via the NSF NCAR
Research Data Archive.

## Arguments

- name:

  Which data set to load: "wind_series", "wind_rose", or "wind_field".

## Value

An object of the class named by `name`.

## Details

All three are derived from 10 m hourly wind from the Climate Forecast
System Reanalysis (CFSR), on CFSR's native grid of approximately 0.32
degrees:

- "wind_series":

  A `wind_series` for the western and central United States (120-90 W,
  30-50 N), with 96 time steps: every sixth hour on the 1st and 15th of
  each month in 2000. Values are rounded to 0.1 m/s. This is for
  demonstrating the `wind_series` -\> `wind_rose` workflow; a real
  analysis would use a longer and denser time series.

- "wind_rose":

  A `wind_rose` (with `trans = 1`) for the same region, built from all
  576 hours on those days. Because it uses six times as many
  observations as `"wind_series"`, it is a better estimate of the
  region's wind regime, and is the better choice for demonstrating
  connectivity analyses.

- "wind_field":

  A `wind_field` for Hurricane Katrina in the Gulf of Mexico, at
  2005-08-28 23:00 UTC (99-78 W, 17-35 N), rounded to 0.1 m/s.

## Examples

``` r
ws <- windscape_example("wind_series")
rose <- windscape_example("wind_rose")
katrina <- windscape_example("wind_field")
```
