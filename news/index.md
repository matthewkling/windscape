# Changelog

## windscape 2.0.0

This release reorganizes windscape around a simpler workflow: download
wind data or a pre-built wind rose, summarize the wind as a `wind_rose`,
model connectivity with least-cost or random-walk functions that take
the rose directly, and plot each step with ggplot2 layers.

### Breaking changes

#### Wind data

- `cfsr_dl()` and `cfsr_dl_land()` are replaced by
  [`download_wind_data()`](https://matthewkling.github.io/windscape/reference/download_wind_data.md)
  and
  [`download_land_mask()`](https://matthewkling.github.io/windscape/reference/download_land_mask.md),
  which use NCAR’s Geoscience Data Exchange (see New features).
- `cfsr_usa` and `katrina` datasets replaced by
  [`windscape_example()`](https://matthewkling.github.io/windscape/reference/windscape_example.md).
- [`wind_series()`](https://matthewkling.github.io/windscape/reference/wind_series.md)
  and
  [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
  now require lon/lat grids with (near-)square cells.
- Selecting layers of a `wind_series`, `wind_field`, or `wind_rose` with
  `[[` or
  [`subset()`](https://rspatial.github.io/terra/reference/subset.html),
  or combining them with [`c()`](https://rdrr.io/r/base/c.html), now
  returns a plain `SpatRaster`, since the result no longer has the
  object’s structure. Use
  [`subset_series()`](https://matthewkling.github.io/windscape/reference/subset_series.md)
  to select time steps of a series.

#### Wind roses

- `as_wind_rose()` and `read_wind_rose()` are merged into
  [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md),
  which builds a rose from a `wind_series`, or loads a saved rose from
  an 8-layer raster or file path. `rose()` is no longer exported.
- [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
  no longer takes `...` (or `cores`); it gains `filename` and
  `overwrite`. `trans` must be a non-negative number or an elementwise
  function of wind speed.
- [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md):
  cells missing wind data at any time step are now `NA` in all eight
  layers, and calm time steps (zero wind speed) contribute no
  conductance whatever `trans` is (previously, with `trans = 0`, they
  counted toward one fixed direction).
- [`terra::aggregate()`](https://rspatial.github.io/terra/reference/aggregate.html),
  [`disagg()`](https://rspatial.github.io/terra/reference/disaggregate.html),
  [`resample()`](https://rspatial.github.io/terra/reference/resample.html),
  and
  [`project()`](https://rspatial.github.io/terra/reference/project.html)
  now signal an error on a `wind_rose`, since they would average
  conductances without accounting for cell size. Use
  [`downscale()`](https://matthewkling.github.io/windscape/reference/downscale.md),
  or change the grid of the `wind_series` before building the rose.

#### Least-cost model

- `least_cost_surface()` is renamed
  [`least_cost()`](https://matthewkling.github.io/windscape/reference/least_cost.md),
  and `least_cost_distance()` is renamed
  [`pairwise_least_cost()`](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md).
  Both take a `wind_rose` directly (a `wind_graph` built in advance is
  still accepted), and
  [`least_cost()`](https://matthewkling.github.io/windscape/reference/least_cost.md)
  gains a `direction` argument.
- `vrcd()` is removed.
  [`pairwise_least_cost()`](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md)
  now places sites at their actual locations within grid cells
  (`snap = FALSE`, the default), replacing both `vrcd()` and the
  `adjust` argument of `least_cost_distance()`.
- [`wind_graph()`](https://matthewkling.github.io/windscape/reference/wind_graph.md):
  `wrap` now defaults to `NULL`, which wraps global grids (spanning 360
  degrees of longitude) automatically.

#### Random walk model

- [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md):
  `"ratchet"` mode removed; results are now a named list of rasters;
  arguments reordered (general arguments before mode-specific ones). The
  new latitude correction (on by default) changes results on lon/lat
  grids.
- [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md)
  results are now per km^2 by default (`density = TRUE`).
- `particle_flow()` is replaced by
  [`wind_trails()`](https://matthewkling.github.io/windscape/reference/wind_trails.md),
  which traces trails for a given transport time (`hours`) or length
  (`distance`).
- `tesselate()` is removed; global grids now wrap automatically (see
  `wrap` arguments).

#### Requirements

- windscape now requires R \>= 4.0.0 and terra \>= 1.8-42.

### New features

#### Wind data

- [`download_wind_data()`](https://matthewkling.github.io/windscape/reference/download_wind_data.md)
  downloads hourly wind from NCAR’s Geoscience Data Exchange, with no
  account needed: ERA5 (1940 to present; 10 m and 100 m), CFSR
  (1979-2010), and CFSv2 (2011 to present; 10 m and pressure levels).
  Data are clipped to a bounding box on the server, can be limited to
  particular `days` and `hours` (down to a single hour), and are saved
  as one file per month, with downloads cached so interrupted runs
  resume. Large requests are split to fit the server’s size limit.
- [`download_land_mask()`](https://matthewkling.github.io/windscape/reference/download_land_mask.md)
  downloads a land-water layer on the same grid, e.g. for
  [`weight_conductance()`](https://matthewkling.github.io/windscape/reference/weight_conductance.md).
- [`wind_series()`](https://matthewkling.github.io/windscape/reference/wind_series.md)
  combines file paths, rasters, or a list of series into one series that
  stays on disk.
- [`subset_series()`](https://matthewkling.github.io/windscape/reference/subset_series.md)
  selects time steps by month, hour, date range, or index, and
  [`wind_times()`](https://matthewkling.github.io/windscape/reference/wind_times.md)
  returns the time of each step.
- [`wind_field()`](https://matthewkling.github.io/windscape/reference/wind_field.md)
  gains `step`, to take one time step of a series as a field, and
  [`mean()`](https://rspatial.github.io/terra/reference/summarize-generics.html)
  of a `wind_series` gives its mean wind field.
- [`windscape_example()`](https://matthewkling.github.io/windscape/reference/windscape_example.md)
  loads small example data sets.

#### Wind roses

- [`download_wind_rose()`](https://matthewkling.github.io/windscape/reference/download_wind_rose.md)
  provides pre-built global wind roses (10 m CFSR, 1979-2010) for any
  combination of years and calendar months, optionally cropped to a
  region, which is read over the internet without downloading whole
  files.
  [`wind_rose_catalog()`](https://matthewkling.github.io/windscape/reference/wind_rose_catalog.md)
  lists the available files and
  [`wind_rose_cache()`](https://matthewkling.github.io/windscape/reference/wind_rose_cache.md)
  manages downloaded ones.
- [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
  is about five times faster, processes long series in chunks of time
  steps without loading them into memory, and records `n_steps` and
  `trans` in files written with `filename`, so they reload with
  `wind_rose(file)` alone.
- [`combine_roses()`](https://matthewkling.github.io/windscape/reference/combine_roses.md)
  combines roses built from different time periods, weighting them by
  their numbers of time steps.
- [`net_flow()`](https://matthewkling.github.io/windscape/reference/net_flow.md)
  returns a rose’s net flow as a `wind_field`.

#### Least-cost model

- [`least_cost_paths()`](https://matthewkling.github.io/windscape/reference/least_cost_paths.md)
  traces fastest routes between sites and other points.

#### Random walk model

- [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md)
  gains `mode = "stream"` (steady state of continuous release, returning
  residence time and deposition), `half_life` (deposition),
  `latitude_correction`, `direction = "upwind"` (catchments), `flux`
  (net flux of material), `source`, and `wrap`.
- [`pairwise_random_walk()`](https://matthewkling.github.io/windscape/reference/pairwise_random_walk.md)
  computes random-walk connectivity among sets of sites.
- [`random_walk_paths()`](https://matthewkling.github.io/windscape/reference/random_walk_paths.md)
  traces flux streamlines between a source and where its material goes,
  the random-walk analog of
  [`least_cost_paths()`](https://matthewkling.github.io/windscape/reference/least_cost_paths.md).
- [`rw_exit_prob()`](https://matthewkling.github.io/windscape/reference/rw_exit_prob.md)
  gives each cell’s probability of mass leaving the domain.
- [`rw_self_retention()`](https://matthewkling.github.io/windscape/reference/rw_self_retention.md)
  for leave-one-out corrections of stream-mode results.

#### Plotting

- ggplot2 layers for each object:
  [`geom_wind_arrow()`](https://matthewkling.github.io/windscape/reference/geom_wind_arrow.md)
  and
  [`geom_wind_trail()`](https://matthewkling.github.io/windscape/reference/geom_wind_trail.md)
  for wind fields,
  [`geom_wind_rose()`](https://matthewkling.github.io/windscape/reference/geom_wind_rose.md)
  for wind roses, and
  [`geom_wind_path()`](https://matthewkling.github.io/windscape/reference/geom_wind_path.md)
  for paths and trails, with
  [`scale_fill_bearing()`](https://matthewkling.github.io/windscape/reference/scale_fill_bearing.md)
  and
  [`scale_colour_bearing()`](https://matthewkling.github.io/windscape/reference/scale_fill_bearing.md)
  for direction-based colors, and
  [`fortify()`](https://ggplot2.tidyverse.org/reference/fortify.html)
  methods for windscape objects.

### Bug fixes

- [`pairwise_least_cost()`](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md)
  (formerly `least_cost_distance(adjust = TRUE)`) no longer returns
  `NaN` on the diagonal.
- windscape no longer masks
  [`gdistance::geoCorrection()`](https://AgrDataSci.github.io/gdistance/reference/geoCorrection.html).
- Removed hidden dependencies on dplyr and magrittr being attached.

## windscape 1.1.0

- Initial release.
