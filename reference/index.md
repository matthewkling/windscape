# Package index

## Wind data

Download hourly wind, and represent it as a time series (`wind_series`)
or a single snapshot (`wind_field`).

- [`download_wind_data()`](https://matthewkling.github.io/windscape/reference/download_wind_data.md)
  : Download hourly wind data from NCAR
- [`download_land_mask()`](https://matthewkling.github.io/windscape/reference/download_land_mask.md)
  : Download a land-water layer from NCAR
- [`wind_series()`](https://matthewkling.github.io/windscape/reference/wind_series.md)
  : Create a wind_series
- [`wind_times()`](https://matthewkling.github.io/windscape/reference/wind_times.md)
  : Get the time of each step in a wind_series
- [`subset_series()`](https://matthewkling.github.io/windscape/reference/subset_series.md)
  : Select time steps from a wind_series
- [`wind_field()`](https://matthewkling.github.io/windscape/reference/wind_field.md)
  : Create a wind_field
- [`mean(`*`<wind_series>`*`)`](https://matthewkling.github.io/windscape/reference/mean-wind_series-method.md)
  : Mean wind of a wind_series

## Wind roses

Summarize a long time series of wind conditions as the conductance from
each cell to its eight neighbors—the basis of both connectivity models.

- [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
  : Build or load a wind rose
- [`download_wind_rose()`](https://matthewkling.github.io/windscape/reference/download_wind_rose.md)
  : Download a pre-built wind rose
- [`wind_rose_catalog()`](https://matthewkling.github.io/windscape/reference/wind_rose_catalog.md)
  : Catalog of pre-built wind roses
- [`wind_rose_cache()`](https://matthewkling.github.io/windscape/reference/wind_rose_cache.md)
  : Manage the wind rose download cache
- [`combine_roses()`](https://matthewkling.github.io/windscape/reference/combine_roses.md)
  : Combine wind roses built from different time periods
- [`weight_conductance()`](https://matthewkling.github.io/windscape/reference/weight_conductance.md)
  : Weight a wind rose's conductances
- [`net_flow()`](https://matthewkling.github.io/windscape/reference/net_flow.md)
  : Net flow of a wind rose, as a wind field
- [`wind_rose-grid`](https://matthewkling.github.io/windscape/reference/wind_rose-grid.md)
  [`aggregate,wind_rose-method`](https://matthewkling.github.io/windscape/reference/wind_rose-grid.md)
  [`disagg,wind_rose-method`](https://matthewkling.github.io/windscape/reference/wind_rose-grid.md)
  [`resample,wind_rose-method`](https://matthewkling.github.io/windscape/reference/wind_rose-grid.md)
  [`project,wind_rose-method`](https://matthewkling.github.io/windscape/reference/wind_rose-grid.md)
  : Changing the grid of a wind rose

## Least-cost model

Travel time along the fastest route through the wind graph.

- [`least_cost()`](https://matthewkling.github.io/windscape/reference/least_cost.md)
  : Least-cost travel time surface
- [`least_cost_paths()`](https://matthewkling.github.io/windscape/reference/least_cost_paths.md)
  : Least-cost paths
- [`pairwise_least_cost()`](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md)
  : Pairwise least-cost travel times among sites
- [`wind_graph()`](https://matthewkling.github.io/windscape/reference/wind_graph.md)
  : Build a wind connectivity graph

## Random walk model

Dispersal of material spreading across the grid in proportion to the
wind.

- [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md)
  : Simulate wind dispersal by random walk
- [`random_walk_paths()`](https://matthewkling.github.io/windscape/reference/random_walk_paths.md)
  : Paths of material through a random walk
- [`pairwise_random_walk()`](https://matthewkling.github.io/windscape/reference/pairwise_random_walk.md)
  : Pairwise random walk connectivity among sites
- [`rw_exit_prob()`](https://matthewkling.github.io/windscape/reference/rw_exit_prob.md)
  : Probability of exiting a random walk's domain
- [`rw_self_retention()`](https://matthewkling.github.io/windscape/reference/rw_self_retention.md)
  : Self-contribution of each source cell in a stream-mode random walk
- [`rw_max_step()`](https://matthewkling.github.io/windscape/reference/rw_max_step.md)
  : Maximum iteration duration for a random walk
- [`iter_length()`](https://matthewkling.github.io/windscape/reference/iter_length.md)
  : Iteration step length of a random walk
- [`check_cell_distance()`](https://matthewkling.github.io/windscape/reference/check_cell_distance.md)
  : Check how grid cells distort distances among sites
- [`downscale()`](https://matthewkling.github.io/windscape/reference/downscale.md)
  : Downscale a wind rose to higher spatial resolution
- [`cell_distance()`](https://matthewkling.github.io/windscape/reference/cell_distance.md)
  : Pairwise distances between cell centroids
- [`ws_summarize()`](https://matthewkling.github.io/windscape/reference/ws_summarize.md)
  : Summary statistics for a windshed

## Hypothesis testing

Test hypotheses about wind connectivity using pairwise matrices of wind,
gene flow, or other relationships among sites.

- [`mantel_test()`](https://matthewkling.github.io/windscape/reference/mantel_test.md)
  : Mantel test
- [`pairwise_ratios()`](https://matthewkling.github.io/windscape/reference/pairwise_ratios.md)
  : Convert asymmetric pairwise matrix to reciprocally symmetrical
  matrix
- [`pairwise_means()`](https://matthewkling.github.io/windscape/reference/pairwise_means.md)
  : Convert asymmetric pairwise matrix to symmetrical pairwise means
- [`point_distance()`](https://matthewkling.github.io/windscape/reference/point_distance.md)
  : Pairwise distances between points

## Plotting

ggplot2 layers for wind data and connectivity results.

### Geoms

- [`geom_wind_rose()`](https://matthewkling.github.io/windscape/reference/geom_wind_rose.md)
  [`stat_wind_rose()`](https://matthewkling.github.io/windscape/reference/geom_wind_rose.md)
  : Wind rose glyphs: a map of local wind roses
- [`geom_wind_arrow()`](https://matthewkling.github.io/windscape/reference/geom_wind_arrow.md)
  [`stat_wind_arrow()`](https://matthewkling.github.io/windscape/reference/geom_wind_arrow.md)
  : Wind field arrows
- [`geom_wind_trail()`](https://matthewkling.github.io/windscape/reference/geom_wind_trail.md)
  [`stat_wind_trail()`](https://matthewkling.github.io/windscape/reference/geom_wind_trail.md)
  : Wind trails and streamlines
- [`geom_wind_path()`](https://matthewkling.github.io/windscape/reference/geom_wind_path.md)
  : Draw existing wind trails or paths
- [`scale_fill_bearing()`](https://matthewkling.github.io/windscape/reference/scale_fill_bearing.md)
  [`scale_colour_bearing()`](https://matthewkling.github.io/windscape/reference/scale_fill_bearing.md)
  [`scale_color_bearing()`](https://matthewkling.github.io/windscape/reference/scale_fill_bearing.md)
  : Color scales for wind direction
- [`fortify(`*`<wind_rose>`*`)`](https://matthewkling.github.io/windscape/reference/fortify.windscape.md)
  [`fortify(`*`<wind_field>`*`)`](https://matthewkling.github.io/windscape/reference/fortify.windscape.md)
  [`fortify(`*`<wind_series>`*`)`](https://matthewkling.github.io/windscape/reference/fortify.windscape.md)
  [`fortify(`*`<random_walk>`*`)`](https://matthewkling.github.io/windscape/reference/fortify.windscape.md)
  : Convert windscape objects to data frames for ggplot2

### Streamlines and particle trails

- [`wind_trails()`](https://matthewkling.github.io/windscape/reference/wind_trails.md)
  : Trace particle trails through a wind field
- [`generate_particles()`](https://matthewkling.github.io/windscape/reference/generate_particles.md)
  : Generate initial particle locations

## Datasets

- [`windscape_example()`](https://matthewkling.github.io/windscape/reference/windscape_example.md)
  : Example wind data sets
- [`birch`](https://matthewkling.github.io/windscape/reference/birch.md)
  : Silver birch landscape genetic data from Tsuda et al. (2017)

## Classes and methods

- [`wind_series-class`](https://matthewkling.github.io/windscape/reference/wind_series-class.md)
  : An S4 object class representing a wind field time series
- [`wind_field-class`](https://matthewkling.github.io/windscape/reference/wind_field-class.md)
  : An S4 object class representing a wind field
- [`wind_rose-class`](https://matthewkling.github.io/windscape/reference/wind_rose-class.md)
  : An S4 object class representing a wind rose
- [`wind_graph-class`](https://matthewkling.github.io/windscape/reference/wind_graph-class.md)
  : An S4 object class representing a wind graph
- [`windscape-layers`](https://matthewkling.github.io/windscape/reference/windscape-layers.md)
  [`[[,wind_series,ANY,ANY-method`](https://matthewkling.github.io/windscape/reference/windscape-layers.md)
  [`[[,wind_field,ANY,ANY-method`](https://matthewkling.github.io/windscape/reference/windscape-layers.md)
  [`[[,wind_rose,ANY,ANY-method`](https://matthewkling.github.io/windscape/reference/windscape-layers.md)
  [`subset,wind_series-method`](https://matthewkling.github.io/windscape/reference/windscape-layers.md)
  [`subset,wind_field-method`](https://matthewkling.github.io/windscape/reference/windscape-layers.md)
  [`subset,wind_rose-method`](https://matthewkling.github.io/windscape/reference/windscape-layers.md)
  [`c,wind_series-method`](https://matthewkling.github.io/windscape/reference/windscape-layers.md)
  [`c,wind_field-method`](https://matthewkling.github.io/windscape/reference/windscape-layers.md)
  [`c,wind_rose-method`](https://matthewkling.github.io/windscape/reference/windscape-layers.md)
  : Selecting and combining layers of windscape objects
- [`geoCorrection(`*`<wind_graph>`*`,`*`<ANY>`*`)`](https://matthewkling.github.io/windscape/reference/geoCorrection-wind_graph-ANY-method.md)
  : geoCorrection is not applicable to wind graphs
