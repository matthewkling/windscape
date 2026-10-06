# Weight a wind rose's conductances

This function adjusts the conductance values in a wind rose, multiplying
them by the values of a secondary raster data set `w`, which can be used
to integrate non-wind factors into a wind connectivity analysis.

## Usage

``` r
weight_conductance(rose, w)
```

## Arguments

- rose:

  A `wind_rose`.

- w:

  A single-layer `SpatRaster` with values to be multiplied by `rose`, on
  the same grid as `rose`.

## Value

A `wind_rose`, with conductance values weighted by `w`.

## Details

As a hypothetical example, to incorporate a decreased but nonzero
likelihood of dispersal over inhospitable areas, `w` could be a raster
layer with 0.1 indicating water or mountains and 1.0 elsewhere. This
would have the effect of down-weighting conductance over water or
mountains by 90%.
