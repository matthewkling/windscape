# Download a land-water layer from NCAR

Downloads a land layer matching the grid of
[`ncar_download()`](https://matthewkling.github.io/windscape/reference/ncar_download.md)
data for the same `source` and bounding box, e.g. for weighting a wind
rose to reduce conductance over water.

## Usage

``` r
ncar_land(source = c("era5", "cfsr", "cfsv2"), xlim, ylim)
```

## Arguments

- source:

  Data set: `"era5"` or `"cfsr"`. (Not yet available for `"cfsv2"`.)

- xlim, ylim:

  Bounding box; see
  [`ncar_download()`](https://matthewkling.github.io/windscape/reference/ncar_download.md).

## Value

A single-layer `SpatRaster` named `"land"`. For ERA5 this is the land
fraction of each cell, from 0 (all water) to 1 (all land); use e.g.
`land >= 0.5` for a binary layer. For CFSR it is binary, 1 for land and
0 for water. (CFSR publishes no land mask, so it is derived from which
cells have soil temperature data.)

## Examples

``` r
if (FALSE) { # \dontrun{
land <- ncar_land("era5", xlim = c(-125, -115), ylim = c(42, 49))
} # }
```
