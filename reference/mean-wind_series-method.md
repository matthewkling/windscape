# Mean wind of a wind_series

Averages the u and v components of a `wind_series` over its time steps,
giving the mean wind vector in each cell as a `wind_field`. Long series
stored in files are processed in blocks, without loading the whole
series into memory.

## Usage

``` r
# S4 method for class 'wind_series'
mean(x, ..., na.rm = FALSE)
```

## Arguments

- x:

  A `wind_series`.

- ...:

  Not used.

- na.rm:

  Logical: ignore missing values when averaging?

## Value

A `wind_field` of the mean u and v components.

## Details

The mean wind vector describes the net drift of the air over the series:
where winds blow from many directions, it is short even if winds are
strong. It is similar to, but not the same as, the
[`net_flow()`](https://matthewkling.github.io/windscape/reference/net_flow.md)
of a wind rose built from the series. Net flow describes the
connectivity model rather than the wind itself: it is shaped by `trans`,
and with `trans = 1` it is typically a few percent weaker than the mean
wind, because the rose divides each wind between two neighbor
directions. Use
[`mean()`](https://rspatial.github.io/terra/reference/summarize-generics.html)
to describe the wind, and
[`net_flow()`](https://matthewkling.github.io/windscape/reference/net_flow.md)
to describe the wind rose used for connectivity modeling.

## See also

[`subset_series()`](https://matthewkling.github.io/windscape/reference/subset_series.md)
to average over selected time steps, such as one season.

## Examples

``` r
series <- windscape_example("wind_series")
prevailing <- mean(series)
```
