# Summary statistics for a windshed

A "windshed" is a raster layer representing upwind or downwind
connectivity for a given point. This function computes a range of
statistics summarizing a windshed, such as the size, centroid location,
and isotropy of the wind accessibility surface. It is useful primarily
for comparing the windshed properties of multiple different sites.

## Usage

``` r
ws_summarize(x, origin, radius = NULL)
```

## Arguments

- x:

  SpatRaster representing wind accessibility, with higher values
  indicating greater accessibility. (For example, this could be the
  output of
  [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md),
  or `least_cost_surface(..., rate = TRUE)`.)

- origin:

  Coordinates of center point (matrix with 2 columns and 1 row).

- radius:

  Optional windshed radius, in km; cells farther from the origin than
  this will be excluded from the summary statistics.

## Value

A named vector of summary statistics:

- "centroid_x": Longitude of windshed "centroid". The centroid is the
  mean of all lon-lat coordinates, weighed by the wind accessibility
  values in `x`.

- "centroid_y": Latitude of windshed centroid.

- "centroid_distance": Distance from origin to windshed centroid, in km.

- "centroid_bearing": Bearing from origin to windshed centroid, in
  degrees.

- "windshed_distance: Mean distance from origin, in km. This is the
  average distance from the origin to every grid cell, weighted by
  accessibility.

- "windshed_bearing": Mean bearing from origin, in degrees. This is the
  circular-mean bearing from the origin to every other grid cell,
  weighted by accessibility.

- "windshed_isotropy": Circular standard deviation of the bearing from
  the origin to every other grid cell, weighted by accessibility. Low
  values indicate that accessible sites are concentrated in a narrow
  range of compass directions, while high values indicate a wide
  distribution of directions. The range of possible values is 0-1.

- "windshed_size": Mean accessibility. This is the average of all
  accessibility values, weighted to correct for grid distortions. Higher
  values indicate landscapes with greater overall wind accessibility.

- "windshed_landarea": A realative measure of the land area covered by
  grid cells with nonzero wind accessibility.
