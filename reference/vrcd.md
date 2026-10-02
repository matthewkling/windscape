# Variable-resolution costDistance

Increasing the spatial resolution of a wind data set (e.g. by
interpolation using
[downscale](https://matthewkling.github.io/windscape/reference/downscale.md))
can greatly improve wind connectivity estimates among nearby sites, but
can make it computationally impossible to model connectivity over large
geographic regions, a trade-off that presents problems for studies that
include both nearby and distant site pairs. This function gets around
this issue with a hybrid approach, using a broad-scale coarse-resolution
wind grid to model wind cost-distance among distant sites, and a
separate local-scale high-resolution interpolated grid to model
connectivity between each pair of nearby sites.

## Usage

``` r
vrcd(
  rose,
  ll,
  threshold_km = 30,
  pad = 1,
  max_nodes = 1e+06,
  direction = "downwind",
  method = "bilinear"
)
```

## Arguments

- rose:

  wind_rose.

- ll:

  longitude and latitude of site locations, as a two-column matrix.

- threshold_km:

  a positive number representing the distance threshold, in kilometers.
  Site pairs closer together than this distance will get a separate
  local high-resolution connectivity model, while pairs farther apart
  will be modeled at the resolution of `rose`.

- pad:

  a positive number indicating how far beyond a pair of sites the
  modeling domain should extend for local models. This is an expansion
  factor giving the padding as a fraction of the maximum latitudinal or
  longitudinal distance between the two sites. The default of 1 is
  reasonable in most cases; smaller values will increase computational
  speed but may fail to account for wind routes beyond the bounding box
  encompassing the site pair.

- max_nodes:

  an integer representing how finely to interpolate wind grids for local
  site pairs. This represents the number of cells in the high-resolution
  interpolated grid for each site pair; larger values remove artifacts
  of the discrete coarse grid more effectively, but have increased
  computational cost.

- direction:

  either "downwind" or "upwind", indicating whether outbound or inbound
  connectivity should be computed. In this context, changing this
  parameter is equivalent to transposing the resulting wind distance
  matrix.

- method:

  disaggregation method; either "near" or "bilinear'; see
  [disagg](https://rspatial.github.io/terra/reference/disaggregate.html)
  for details.

## Value

A list of square matrices:

- wind_dist:

  Wind cost-distances between site pairs, in hours if `rose` has a p
  = 1. Computed using
  [costDistance](https://AgrDataSci.github.io/gdistance/reference/costDistance-methods.html).

- wind_dist_coarse:

  Wind cost-distances using the coarse input wind grid (`rose`); this
  will differ from `wind_dist` only for site pairs closer than
  `threshold_km`. This is provided for reference, to judge how
  higher-resolution models impact results relative to the coarse grid.

- point_dist:

  Distances between sites, in km.

- cell_dist:

  Distances between the centroids of the grid cells used to model wind
  connectivity, in km. For site pairs where this distance deviates from
  `point_dist` by a substantial percentage, `wind_dist` estimates may
  contain nontrivial rounding noise resulting from the discrete grid.

- cell_dist_coarse:

  Distances between the cell centroids in the coarse-resolution input
  wind dataset (`rose`), in km. This is provided for reference, to judge
  improvements relative to `cell_dist`.
