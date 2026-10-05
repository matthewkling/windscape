# Color scales for wind direction

Cyclic color scales for compass bearings (degrees clockwise from north),
using the same hues as the windscape direction color wheel. Bearings
outside 0-360 are wrapped, so 0 and 360 (or -90 and 270) get the same
color. By default the legend shows the eight compass directions; pass
`breaks` and `labels` to change this.

## Usage

``` r
scale_fill_bearing(..., chroma = 90, luminance = 65, aesthetics = "fill")

scale_colour_bearing(..., chroma = 90, luminance = 65, aesthetics = "colour")

scale_color_bearing(..., chroma = 90, luminance = 65, aesthetics = "colour")
```

## Arguments

- ...:

  Other arguments passed to
  [`ggplot2::continuous_scale()`](https://ggplot2.tidyverse.org/reference/continuous_scale.html),
  such as `name` or `guide`.

- chroma, luminance:

  HCL chroma and luminance of the colors. Defaults match the windscape
  rose plots.

- aesthetics:

  The aesthetics to which the scale applies.

## Value

A ggplot2 scale.

## Examples

``` r
library(ggplot2)
d <- data.frame(x = 1:8, bearing = seq(0, 315, 45))
ggplot(d, aes(x, 1, fill = bearing)) + geom_tile() + scale_fill_bearing()
```
