# Generate a wind field time series data set from a set of rasters.

Generate a wind field time series data set from a set of rasters.

## Usage

``` r
wind_series(x, order = c("uuvv", "uvuv"))
```

## Arguments

- x:

  Multi-layer `SpatRaster` with layers containing u and v wind
  components, or an object like a file path that can be converted to a
  `SpatRaster`. Data must be on a longitude/latitude grid with square
  cells (equal x and y resolution in degrees), with u and v components
  oriented to true east and north; see
  [wind_rose](https://matthewkling.github.io/windscape/reference/wind_rose.md).

- order:

  Either `"uuvv"`, the default, indicating \`x\` has all u components
  followed by all v components, or `"uvuv"`, indicating the u and v
  components of \`x\` are alternating.

## Value

A \`wind_series\` object, which is a particular form of `SpatRaster`.
