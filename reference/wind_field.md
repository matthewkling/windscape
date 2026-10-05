# Generate a wind field data set from a pair of rasters.

Generate a wind field data set from a pair of rasters.

## Usage

``` r
wind_field(x)
```

## Arguments

- x:

  A `SpatRaster` with two layers representing u and v wind components;
  note that these must be in lat-long coordinates.

## Value

A `wind_field` object, which is a particular form of `SpatRaster`.
