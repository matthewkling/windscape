# Weight a windrose conductance raster

This function adjusts the conductance values in a windrose raster,
multiplying them by the values of a secondary raster data set `w`, which
can be used to integrate non-wind factors into a wind connectivity
analysis.

## Usage

``` r
weight_conductance(rose, w)
```

## Arguments

- rose:

  A raster stack created using
  [`windrose_rasters()`](https://matthewkling.github.io/windscape/reference/wind_rose.md).

- w:

  A raster layer with values to be multiplied by `rose`. This layer must
  have the same spatial properties as the windrose raster.

## Value

A version of `rose`, with conductance values weighted by `w`.

## Details

As a hypothetical example, to incorporate a decreased but nonzero
likelihood of dispersal over inhospitable areas, `w` could be a raster
layer with 0.1 indicating water or mountains and 1.0 elsewhere. This
would have the effect of down-weighting conductance over water or
mountains by 90%.
