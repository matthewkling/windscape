# Package index

## All functions

- [`as_wind_rose()`](https://matthewkling.github.io/windscape/reference/as_wind_rose.md)
  : Create a wind_rose object from a set of raster layers
- [`birch`](https://matthewkling.github.io/windscape/reference/birch.md)
  : Silver birch landscape genetic data from Tsuda et al. (2017)
- [`cell_distance()`](https://matthewkling.github.io/windscape/reference/cell_distance.md)
  : Pairwise distances between cell centroids
- [`check_cell_distance()`](https://matthewkling.github.io/windscape/reference/check_cell_distance.md)
  : Check how grid cells distort distances among sites
- [`combine_roses()`](https://matthewkling.github.io/windscape/reference/combine_roses.md)
  : Combine wind roses built from different time periods
- [`downscale()`](https://matthewkling.github.io/windscape/reference/downscale.md)
  : Downscale a wind rose to higher spatial resolution
- [`fortify(`*`<wind_rose>`*`)`](https://matthewkling.github.io/windscape/reference/fortify.windscape.md)
  [`fortify(`*`<wind_field>`*`)`](https://matthewkling.github.io/windscape/reference/fortify.windscape.md)
  [`fortify(`*`<wind_series>`*`)`](https://matthewkling.github.io/windscape/reference/fortify.windscape.md)
  [`fortify(`*`<random_walk>`*`)`](https://matthewkling.github.io/windscape/reference/fortify.windscape.md)
  : Convert windscape objects to data frames for ggplot2
- [`generate_particles()`](https://matthewkling.github.io/windscape/reference/generate_particles.md)
  : Generate initial particle locations
- [`geoCorrection(`*`<wind_graph>`*`,`*`<ANY>`*`)`](https://matthewkling.github.io/windscape/reference/geoCorrection-wind_graph-ANY-method.md)
  : geoCorrection is not applicable to wind graphs
- [`geom_wind_arrow()`](https://matthewkling.github.io/windscape/reference/geom_wind_arrow.md)
  [`stat_wind_arrow()`](https://matthewkling.github.io/windscape/reference/geom_wind_arrow.md)
  : Wind field arrows
- [`geom_wind_path()`](https://matthewkling.github.io/windscape/reference/geom_wind_path.md)
  : Draw existing wind trails or paths
- [`geom_wind_rose()`](https://matthewkling.github.io/windscape/reference/geom_wind_rose.md)
  [`stat_wind_rose()`](https://matthewkling.github.io/windscape/reference/geom_wind_rose.md)
  : Wind rose glyphs: a map of local wind roses
- [`geom_wind_trail()`](https://matthewkling.github.io/windscape/reference/geom_wind_trail.md)
  [`stat_wind_trail()`](https://matthewkling.github.io/windscape/reference/geom_wind_trail.md)
  : Wind trails and streamlines
- [`iter_length()`](https://matthewkling.github.io/windscape/reference/iter_length.md)
  : Iteration step length of a random walk
- [`least_cost_paths()`](https://matthewkling.github.io/windscape/reference/least_cost_paths.md)
  : Least-cost paths through a wind graph
- [`least_cost_surface()`](https://matthewkling.github.io/windscape/reference/least_cost_surface.md)
  : Accumulated wind cost surface
- [`mantel_test()`](https://matthewkling.github.io/windscape/reference/mantel_test.md)
  : Mantel test
- [`ncar_download()`](https://matthewkling.github.io/windscape/reference/ncar_download.md)
  : Download hourly wind data from NCAR
- [`ncar_land()`](https://matthewkling.github.io/windscape/reference/ncar_land.md)
  : Download a land-water layer from NCAR
- [`net_flow()`](https://matthewkling.github.io/windscape/reference/net_flow.md)
  : Net flow of a wind rose, as a wind field
- [`pairwise_least_cost()`](https://matthewkling.github.io/windscape/reference/pairwise_least_cost.md)
  : Pairwise least-cost travel times among sites
- [`pairwise_means()`](https://matthewkling.github.io/windscape/reference/pairwise_means.md)
  : Convert asymmetric pairwise matrix to symmetrical matrix of pairwise
  means
- [`pairwise_random_walk()`](https://matthewkling.github.io/windscape/reference/pairwise_random_walk.md)
  : Pairwise random walk connectivity among sites
- [`pairwise_ratios()`](https://matthewkling.github.io/windscape/reference/pairwise_ratios.md)
  : Convert data to reciprocally symmetrical pairwise matrix
- [`point_distance()`](https://matthewkling.github.io/windscape/reference/point_distance.md)
  : Pairwise distances between points
- [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md)
  : Simulate wind dispersal by random walk
- [`read_wind_rose()`](https://matthewkling.github.io/windscape/reference/read_wind_rose.md)
  : Load a wind_rose from a raster file on disk
- [`read_wind_series()`](https://matthewkling.github.io/windscape/reference/read_wind_series.md)
  : Load a wind_series from one or more raster files on disk
- [`rose()`](https://matthewkling.github.io/windscape/reference/rose.md)
  : Calculate 8-neighbor edge loadings from a time series of u and v
  windspeeds
- [`rw_max_step()`](https://matthewkling.github.io/windscape/reference/rw_max_step.md)
  : Maximum iteration duration for a random walk
- [`rw_self_retention()`](https://matthewkling.github.io/windscape/reference/rw_self_retention.md)
  : Self-contribution of each source cell in a stream-mode random walk
- [`scale_fill_bearing()`](https://matthewkling.github.io/windscape/reference/scale_fill_bearing.md)
  [`scale_colour_bearing()`](https://matthewkling.github.io/windscape/reference/scale_fill_bearing.md)
  [`scale_color_bearing()`](https://matthewkling.github.io/windscape/reference/scale_fill_bearing.md)
  : Color scales for wind direction
- [`tesselate()`](https://matthewkling.github.io/windscape/reference/tesselate.md)
  : Extend a global raster with copies of itself
- [`weight_conductance()`](https://matthewkling.github.io/windscape/reference/weight_conductance.md)
  : Weight a windrose conductance raster
- [`wind_field-class`](https://matthewkling.github.io/windscape/reference/wind_field-class.md)
  : An S4 object class representing a wind field
- [`wind_field()`](https://matthewkling.github.io/windscape/reference/wind_field.md)
  : Generate a wind field data set from a pair of rasters.
- [`wind_graph-class`](https://matthewkling.github.io/windscape/reference/wind_graph-class.md)
  : An S4 object class representing a wind graph
- [`wind_graph()`](https://matthewkling.github.io/windscape/reference/wind_graph.md)
  : Build a wind connectivity graph
- [`wind_rose-class`](https://matthewkling.github.io/windscape/reference/wind_rose-class.md)
  : An S4 object class representing a wind rose
- [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
  : Summarize a time series of wind fields into a wind rose
- [`wind_series-class`](https://matthewkling.github.io/windscape/reference/wind_series-class.md)
  : An S4 object class representing a wind field time series
- [`wind_series()`](https://matthewkling.github.io/windscape/reference/wind_series.md)
  : Generate a wind field time series data set from a set of rasters.
- [`wind_trails()`](https://matthewkling.github.io/windscape/reference/wind_trails.md)
  : Trace particle trails through a wind field
- [`windscape_example()`](https://matthewkling.github.io/windscape/reference/windscape_example.md)
  : Example wind data sets
- [`ws_summarize()`](https://matthewkling.github.io/windscape/reference/ws_summarize.md)
  : Summary statistics for a windshed
