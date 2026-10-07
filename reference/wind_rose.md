# Build or load a wind rose

A wind rose summarizes a time series of wind fields as the average wind
conductance from each grid cell toward each of its eight neighbors.
Given a `wind_series`, `wind_rose()` builds one; given a raster or file
holding a saved wind rose, it loads it.

## Usage

``` r
wind_rose(
  x,
  trans = 1,
  n_steps = NA_integer_,
  filename = NULL,
  overwrite = FALSE
)
```

## Arguments

- x:

  A `wind_series`, to build a wind rose from; or a saved wind rose to
  load, as an 8-layer `SpatRaster` or the path to a raster file (e.g.
  one written with
  [`terra::writeRaster()`](https://rspatial.github.io/terra/reference/writeRaster.html)).
  A saved rose's layers must hold conductance toward the southwest,
  west, northwest, north, northeast, east, southeast, and south
  neighbors, in that order.

- trans:

  Either a non-negative number indicating the power to raise wind speeds
  to, or an elementwise function of wind speed (it may be applied to
  many cells and time steps at once, so its result for each speed must
  not depend on the others); see details. When loading a saved rose, the
  `trans` it was built with, which is recorded with the rose (e.g. for
  [`combine_roses()`](https://matthewkling.github.io/windscape/reference/combine_roses.md));
  if not given, it is read from the file's metadata where recorded there
  (as in files from
  [`download_wind_rose()`](https://matthewkling.github.io/windscape/reference/download_wind_rose.md)),
  and otherwise defaults to 1.

- n_steps:

  When loading a saved rose, the number of time steps it summarizes,
  needed to combine it with other roses using
  [`combine_roses()`](https://matthewkling.github.io/windscape/reference/combine_roses.md);
  if not given, it is read from the file's metadata where recorded
  there. Ignored when building a rose, which records its number of time
  steps automatically.

- filename:

  When building a rose, an optional file path to write the result to, as
  a raster file (e.g. a GeoTIFF).

- overwrite:

  Logical. Whether to overwrite an existing `filename`.

## Value

A `wind_rose` object. This is an 8-layer raster stack, where each layer
is wind conductance from the focal cell to one of its neighbors
(clockwise starting in the SW). If input windspeeds are in m/s and
`trans = 1`, values are in (1 / hours). When building a rose, cells
missing wind data at any time step are `NA` in all eight layers.

## Details

For each time step, the wind in each cell is divided between the two
neighbors whose directions bracket the wind direction, in proportion to
how closely the wind points toward each, and its transformed speed (see
`trans`) is converted to conductance toward each neighbor; conductance
is then averaged over all time steps. Long series are processed in
chunks of time steps, so a series spanning many files (see
[`wind_series()`](https://matthewkling.github.io/windscape/reference/wind_series.md))
can be summarized without loading it all into memory; the result is
identical to processing it at once. To build a rose from selected time
steps, such as a season or time of day, select them first with
[`subset_series()`](https://matthewkling.github.io/windscape/reference/subset_series.md).
To build many roses (e.g. one per month), run separate `wind_rose()`
calls in parallel processes, e.g. with
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html).

The `trans` parameter defines the transformation function used to
convert wind speed into conductance. If a numeric value is supplied, the
function speed^trans is used. A value of trans = 0 will ignore speed,
assigning weights based on direction only; trans = 1 assumes conductance
is proportional to windspeed, trans = 2 assumes it's proportional to
aerodynamic drag, and trans = 3 assumes it's proportional to force. Any
intermediate value can also be used. Any elementwise function,
transforming each speed independently of the others, can also be
supplied; for example, to model seed dispersal for a species that only
releases seeds when winds exceed 10 m/s, we could specify a threshold
function `trans = function(x){x[x < 10] <- 0; return(x)}`.

Grid geometry: windscape works on longitude/latitude grids with square
cells. Distances and bearings to each cell's neighbors are computed on
the ellipsoid at that cell's latitude, so conductance accounts for the
narrowing of cells toward the poles. Projected grids are not currently
supported. If source data are projected, reproject them to
longitude/latitude before building a `wind_series`, rotating u and v to
true east and north if they are defined relative to the projected grid
(as in some reanalysis products).

## See also

[`combine_roses()`](https://matthewkling.github.io/windscape/reference/combine_roses.md)
to combine roses built from different time periods.

## Examples

``` r
series <- windscape_example("wind_series")
rose <- wind_rose(series)

# save and reload
f <- tempfile(fileext = ".tif")
terra::writeRaster(rose, f)
rose2 <- wind_rose(f, trans = 1, n_steps = rose@n_steps)
```
