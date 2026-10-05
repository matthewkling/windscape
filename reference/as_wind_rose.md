# Create a wind_rose object from a set of raster layers

Create a wind_rose object from a set of raster layers

## Usage

``` r
as_wind_rose(x, trans, n_steps = NA_integer_)
```

## Arguments

- x:

  SpatRaster with 8 layers representing flows in each semi-cardinal
  direction, clockwise beginning in southwest.

- trans:

  Transformation to convert speed into conductance. Either a single
  number, or a function. See documentation for
  [wind_rose](https://matthewkling.github.io/windscape/reference/wind_rose.md).

- n_steps:

  An integer giving the number of time steps represented in `x`.

## Value

A `wind_rose` object.
