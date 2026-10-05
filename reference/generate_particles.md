# Generate initial particle locations

Initialize a set of particle locations for use in
[`wind_trails()`](https://matthewkling.github.io/windscape/reference/wind_trails.md).

## Usage

``` r
generate_particles(
  x,
  n = 1000,
  sample = "random",
  equalarea = TRUE,
  epsg = 8857,
  expand = 0.25
)
```

## Arguments

- x:

  A wind_field object.

- n:

  Positive integer giving the total number of particles in the grid. If
  `sample = "random"` the result will have this exact number of
  particles, while if `sample = "grid"` it will be approximate.

- sample:

  Sampling scheme, either "random" (the default) to generate points at
  random locations, or "grid" to generate points on a regular grid.

- equalarea:

  Logical indicating whether to generate points on an equal area basis
  (TRUE, the default) or in lon-lat space (FALSE, which over-samples
  higher-latitude regions).

- epsg:

  EPSG code for the projection in which to generate points. Only used if
  `equalarea = TRUE` and `sample = "grid"`. The default is `8857`, the
  Equal Area projection.

- expand:

  Proportion by which to expand the bounding box around `x` after
  projection, to capture the full domain in projected space. Only used
  if `equalarea = TRUE` and `sample = "grid"`.

## Value

A two-column matrix of x (longitude) and y (latitude) values, with a row
for each particle.
