# Select time steps from a wind_series

Selects time steps from a `wind_series` by month, hour of day, date
range, or index, keeping the u and v layers of each step together. Use
it to build wind roses for the times that matter for dispersal, such as
a flowering season or the hours when propagules are released.

## Usage

``` r
subset_series(
  x,
  months = NULL,
  hours = NULL,
  start = NULL,
  end = NULL,
  steps = NULL
)
```

## Arguments

- x:

  A `wind_series`, with times in its layer names (see
  [`wind_times()`](https://matthewkling.github.io/windscape/reference/wind_times.md))
  if selecting by `months`, `hours`, `start`, or `end`.

- months:

  Integer vector of months (1-12) to keep.

- hours:

  Integer vector of hours of the day (0-23, UTC) to keep.

- start, end:

  Keep time steps from `start` to `end`, inclusive. Each can be a
  `Date`, a POSIXct time, or a character string like `"2000-06-01"` or
  `"2000-06-01 12:00:00"` (interpreted as UTC). A date without a time
  includes that whole day, so `end = "2000-06-30"` keeps all of June 30.

- steps:

  Time steps to keep, as integer indices or a logical vector with one
  element per time step, for selections the other arguments don't cover.
  Can be combined with them, e.g. `steps = speed > 10`, given a vector
  `speed`.

## Value

A `wind_series` containing the selected time steps, in their original
order.

## Details

Criteria are combined: a time step is kept only if it meets all of them.
Times are in UTC. For hours of the day in local solar time, note that
local time is about UTC plus longitude / 15 hours (e.g. UTC - 7 hours at
105 degrees west), so `hours` for a large region selects different local
times in different places.

## See also

[`wind_times()`](https://matthewkling.github.io/windscape/reference/wind_times.md)

## Examples

``` r
series <- windscape_example("wind_series")

# summer only
summer <- subset_series(series, months = 6:8)
range(wind_times(summer))
#> [1] "2000-06-01 00:00:00 UTC" "2000-08-15 18:00:00 UTC"

# afternoons in the first half of the year
subset_series(series, hours = 18:23, end = "2000-06-30")
#> class       : SpatRaster
#> size        : 64, 96, 24  (nrow, ncol, nlyr)
#> resolution  : 0.3157895, 0.3174603  (x, y)
#> extent      : -120.1579, -89.84211, 29.84127, 50.15873  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> sources     : wind_usa.tif
#> names       : u 200~00:00, u 200~00:00, u 200~00:00, u 200~00:00, u 200~00:00, u 200~00:00, ...
#> min values  :        -0.6,        -0.4,        -0.6,        -0.5,        -0.5,        -0.7, ...
#> max values  :           1,         0.8,         1.1,         1.4,         0.7,         0.8, ...
```
