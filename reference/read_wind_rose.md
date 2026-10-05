# Load a wind_rose from a raster file on disk

Load a wind_rose from a raster file on disk

## Usage

``` r
read_wind_rose(x, trans = 1, n_steps = NA_integer_)
```

## Arguments

- x:

  Path to a raster file containing wind rose data

- trans:

  Transformation to convert speed into conductance. Either a single
  number, or a function. See documentation for
  [wind_rose](https://matthewkling.github.io/windscape/reference/wind_rose.md).

- n_steps:

  An integer giving the number of time steps represented in `x`.

## Value

A `wind_rose` object.
