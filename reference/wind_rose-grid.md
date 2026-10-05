# Changing the grid of a wind rose

A wind rose's conductances are rates of flow between neighboring cells,
so they depend on cell size. Changing the grid with
[`terra::aggregate()`](https://rspatial.github.io/terra/reference/aggregate.html),
[`terra::disagg()`](https://rspatial.github.io/terra/reference/disaggregate.html),
[`terra::resample()`](https://rspatial.github.io/terra/reference/resample.html),
or
[`terra::project()`](https://rspatial.github.io/terra/reference/project.html)
would average or interpolate conductances without accounting for the new
cell size, giving wrong results, so these operations are not allowed on
a `wind_rose`. To refine a rose's grid, use
[`downscale()`](https://matthewkling.github.io/windscape/reference/downscale.md).
To use a coarser or different grid, change the grid of the `wind_series`
(e.g. with
[`terra::aggregate()`](https://rspatial.github.io/terra/reference/aggregate.html))
and build a new rose from it.

## Arguments

- x:

  A `wind_rose`.

- y, ...:

  Not used.

## Value

These methods signal an error.
