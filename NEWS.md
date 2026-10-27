# windscape (development version)

# windscape 2.0.0

This release reorganizes windscape around a simpler workflow: download wind data or a pre-built
wind rose, summarize the wind as a `wind_rose`, model connectivity with least-cost or random-walk
functions that take the rose directly, and plot each step with ggplot2 layers.

## Breaking changes

### Wind data

* `cfsr_dl()` and `cfsr_dl_land()` are replaced by `download_wind_data()` and
  `download_land_mask()`, which use NCAR's Geoscience Data Exchange (see New features).
* `cfsr_usa` and `katrina` datasets replaced by `windscape_example()`.
* `wind_series()` and `wind_rose()` now require lon/lat grids with (near-)square cells.
* Selecting layers of a `wind_series`, `wind_field`, or `wind_rose` with `[[` or `subset()`, or
  combining them with `c()`, now returns a plain `SpatRaster`, since the result no longer has
  the object's structure. Use `subset_series()` to select time steps of a series.

### Wind roses

* `as_wind_rose()` and `read_wind_rose()` are merged into `wind_rose()`, which builds a rose
  from a `wind_series`, or loads a saved rose from an 8-layer raster or file path. `rose()` is
  no longer exported.
* `wind_rose()` no longer takes `...` (or `cores`); it gains `filename` and `overwrite`. `trans`
  must be a non-negative number or an elementwise function of wind speed.
* `wind_rose()`: cells missing wind data at any time step are now `NA` in all eight layers, and
  calm time steps (zero wind speed) contribute no conductance whatever `trans` is (previously,
  with `trans = 0`, they counted toward one fixed direction).
* `terra::aggregate()`, `disagg()`, `resample()`, and `project()` now signal an error on a
  `wind_rose`, since they would average conductances without accounting for cell size. Use
  `downscale()`, or change the grid of the `wind_series` before building the rose.

### Least-cost model

* `least_cost_surface()` is renamed `least_cost()`, and `least_cost_distance()` is renamed
  `pairwise_least_cost()`. Both take a `wind_rose` directly (a `wind_graph` built in advance
  is still accepted), and `least_cost()` gains a `direction` argument.
* `vrcd()` is removed. `pairwise_least_cost()` now places sites at their actual locations
  within grid cells (`snap = FALSE`, the default), replacing both `vrcd()` and the `adjust`
  argument of `least_cost_distance()`.
* `wind_graph()`: `wrap` now defaults to `NULL`, which wraps global grids (spanning 360 degrees
  of longitude) automatically.

### Random walk model

* `random_walk()`: `"ratchet"` mode removed; results are now a named list of rasters;
  arguments reordered (general arguments before mode-specific ones). The new latitude
  correction (on by default) changes results on lon/lat grids.
* `random_walk()` results are now per km^2 by default (`density = TRUE`).
* `particle_flow()` is replaced by `wind_trails()`, which traces trails for a given transport
  time (`hours`) or length (`distance`).
* `tesselate()` is removed; global grids now wrap automatically (see `wrap` arguments).

### Requirements

* windscape now requires R >= 4.0.0 and terra >= 1.8-42.

## New features

### Wind data

* `download_wind_data()` downloads hourly wind from NCAR's Geoscience Data Exchange, with no
  account needed: ERA5 (1940 to present; 10 m and 100 m), CFSR (1979-2010), and CFSv2 (2011 to
  present; 10 m and pressure levels). Data are clipped to a bounding box on the server, can be
  limited to particular `days` and `hours` (down to a single hour), and are saved as one file
  per month, with downloads cached so interrupted runs resume. Large requests are split to fit
  the server's size limit.
* `download_land_mask()` downloads a land-water layer on the same grid, e.g. for
  `weight_conductance()`.
* `wind_series()` combines file paths, rasters, or a list of series into one series that
  stays on disk.
* `subset_series()` selects time steps by month, hour, date range, or index, and
  `wind_times()` returns the time of each step.
* `wind_field()` gains `step`, to take one time step of a series as a field, and `mean()` of a
  `wind_series` gives its mean wind field.
* `windscape_example()` loads small example data sets.

### Wind roses

* `download_wind_rose()` provides pre-built global wind roses (10 m CFSR, 1979-2010) for any
  combination of years and calendar months, optionally cropped to a region, which is read over
  the internet without downloading whole files. `wind_rose_catalog()` lists the available files
  and `wind_rose_cache()` manages downloaded ones.
* `wind_rose()` is about five times faster, processes long series in chunks of time steps
  without loading them into memory, and records `n_steps` and `trans` in files written with
  `filename`, so they reload with `wind_rose(file)` alone.
* `combine_roses()` combines roses built from different time periods, weighting them by their
  numbers of time steps.
* `net_flow()` returns a rose's net flow as a `wind_field`.

### Least-cost model

* `least_cost_paths()` traces fastest routes between sites and other points.

### Random walk model

* `random_walk()` gains `mode = "stream"` (steady state of continuous release, returning
  residence time and deposition), `half_life` (deposition), `latitude_correction`,
  `direction = "upwind"` (catchments), `flux` (net flux of material), `source`, and `wrap`.
* `pairwise_random_walk()` computes random-walk connectivity among sets of sites.
* `random_walk_paths()` traces flux streamlines between a source and where its material goes,
  the random-walk analog of `least_cost_paths()`.
* `rw_exit_prob()` gives each cell's probability of mass leaving the domain.
* `rw_self_retention()` for leave-one-out corrections of stream-mode results.

### Plotting

* ggplot2 layers for each object: `geom_wind_arrow()` and `geom_wind_trail()` for wind fields,
  `geom_wind_rose()` for wind roses, and `geom_wind_path()` for paths and trails, with
  `scale_fill_bearing()` and `scale_colour_bearing()` for direction-based colors, and
  `fortify()` methods for windscape objects.

## Bug fixes

* `pairwise_least_cost()` (formerly `least_cost_distance(adjust = TRUE)`) no longer returns
  `NaN` on the diagonal.
* windscape no longer masks `gdistance::geoCorrection()`.
* Removed hidden dependencies on dplyr and magrittr being attached.

# windscape 1.1.0

* Initial release.
