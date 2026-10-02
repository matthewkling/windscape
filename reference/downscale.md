# Downscale a wind rose raster data to higher spatial resolution

This function disaggregates a wind_rose raster, interpolating it to
higher spatial resolution and then adjusting values to maintain correct
connectivity units. Downscaling is useful when connectivity is to be
modeled between sites that are separated by few or no grid cells.

## Usage

``` r
downscale(x, fact, method = "bilinear")
```

## Arguments

- x:

  wind_rose

- fact:

  a single integer giving the factor by which to disaggregate the data;
  this is passed to
  [disagg](https://rspatial.github.io/terra/reference/disaggregate.html).

- method:

  disaggregation method; either "near" or "bilinear'; see
  [disagg](https://rspatial.github.io/terra/reference/disaggregate.html)
  for details.

## Value

A higher resolution version of x

## Details

Larger values of `fact` create higher-resolution wind surfaces, which
improves connectivity estimates between nearby sites. The trade-off is
computational cost.
