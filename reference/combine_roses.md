# Combine wind roses built from different time periods

Combines wind roses for the same grid, built from different sets of time
steps (e.g. separate months or years), into a single rose representing
all of those time steps. Because a wind rose is an average over time
steps, the result is the mean of the input roses weighted by their
numbers of time steps, and equals the rose that would be built from all
the time steps at once.

## Usage

``` r
combine_roses(...)
```

## Arguments

- ...:

  Two or more `wind_rose` objects, or a single list of them. They must
  share the same grid and the same `trans` function, and each must have
  a known number of time steps (`n_steps`), as roses built with
  [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
  do.

## Value

A `wind_rose` whose `n_steps` is the total across the inputs.

## See also

[`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md),
which uses this function to build roses from multiple files.
