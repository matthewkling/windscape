# geoCorrection is not applicable to wind graphs

Wind graph conductances already account for the distances between cell
centers (see
[wind_rose](https://matthewkling.github.io/windscape/reference/wind_rose.md)),
so
[geoCorrection](https://AgrDataSci.github.io/gdistance/reference/geoCorrection.html)
must not be applied to them.

## Usage

``` r
# S4 method for class 'wind_graph,ANY'
geoCorrection(x, type, ...)
```

## Arguments

- x:

  A `wind_graph`.

- type, ...:

  Ignored.

## Value

Always an error.
