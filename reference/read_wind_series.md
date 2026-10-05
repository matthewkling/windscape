# Load a wind_series from one or more raster files on disk

Reads files in `wind_series` layout (all u layers followed by all v
layers), such as the monthly files saved by
[`ncar_download()`](https://matthewkling.github.io/windscape/reference/ncar_download.md),
and combines them into a single `wind_series`. Data are not loaded into
memory until needed.

## Usage

``` r
read_wind_series(x)
```

## Arguments

- x:

  Character vector of file paths. Each file must hold an equal number of
  u and v layers, with u layers first, and all files must share the same
  grid. Time steps are combined in the order the files are given.

## Value

A `wind_series` object whose time steps are those of all the files
combined: all files' u layers, followed by all files' v layers.

## See also

[`ncar_download()`](https://matthewkling.github.io/windscape/reference/ncar_download.md)
