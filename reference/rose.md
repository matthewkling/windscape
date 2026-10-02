# Calculate 8-neighbor edge loadings from a time series of u and v windspeeds

Calculate 8-neighbor edge loadings from a time series of u and v
windspeeds

## Usage

``` r
rose(x, trans = identity)
```

## Arguments

- x:

  A vector of wind data containing: latitude, resolution, u windspeeds,
  v windspeeds

- trans:

  A function transforming wind speed into conductance (see details)

## Value

A vector of 8 conductance values to neighboring cells, clockwise
starting with the southwest neighbor. If input windspeeds are in m/s,
values are in (1 / hours ^ p)

## Details

`trans` is applied to each wind speed observation. For example,
`function(s) s^0` ignores speed, assigning weights based on direction
only; `identity` (the default) assumes conductance is proportional to
windspeed, `function(s) s^2` assumes it's proportional to aerodynamic
drag, and `function(s) s^3` assumes it's proportional to force. See
[wind_rose](https://matthewkling.github.io/windscape/reference/wind_rose.md).
