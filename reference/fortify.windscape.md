# Convert windscape objects to data frames for ggplot2

[`fortify()`](https://ggplot2.tidyverse.org/reference/fortify.html)
methods that convert windscape objects to tidy data frames, so they can
be passed directly to
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
or a layer's `data` argument. They can also be called directly, to
inspect or modify the data before plotting.

## Usage

``` r
# S3 method for class 'wind_rose'
fortify(model, data, na.rm = TRUE, ...)

# S3 method for class 'wind_field'
fortify(model, data, na.rm = TRUE, ...)

# S3 method for class 'wind_series'
fortify(model, data, na.rm = TRUE, ...)

# S3 method for class 'random_walk'
fortify(model, data, na.rm = TRUE, ...)
```

## Arguments

- model:

  A windscape object: a `wind_rose`, `wind_field`, `wind_series`, or the
  result of
  [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md).

- data:

  Not used.

- na.rm:

  Logical: drop grid cells with missing values? Default `TRUE`.

- ...:

  Not used.

## Value

A data frame with grid cell center coordinates `x` and `y`, plus:

- `wind_rose`: one column of conductance per direction, `SW`, `W`, `NW`,
  `N`, `NE`, `E`, `SE`, `S`, and `total`, their sum. For
  longitude/latitude roses, also per-cell flow summaries matching those
  computed by
  [`geom_wind_rose()`](https://matthewkling.github.io/windscape/reference/geom_wind_rose.md):
  `speed`, the mean wind speed in km/h (the sum of the flows toward the
  eight neighbors, if `trans = 1` and wind speeds are in m/s); `net`,
  the speed of net flow in km/h; `bearing`, its direction in degrees
  clockwise from north; and `consistency`, `net / speed`, the steadiness
  of wind direction over time (near 1 where wind blows predominantly one
  way, near 0 where it has no net direction).

- `wind_field`: wind components `u` and `v`; `speed`, in the units of
  `u` and `v`; and `bearing`, the direction the wind blows toward, in
  degrees clockwise from north.

- `wind_series`: one row per grid cell and time step, with `step` (the
  time step's index), `time` (parsed from layer names where possible,
  otherwise `NA`), and `u`, `v`, `speed`, and `bearing` as for a
  `wind_field`.

- [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md)
  result, pulse mode: one row per grid cell and recorded iteration, with
  `iteration`, `hours` (elapsed time), `airborne`, and `deposition`.

- [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md)
  result, stream mode: `residence`, `deposition`, and for upwind walks,
  `origin` (the `flux` element is a `wind_field`; fortify it
  separately).

## Examples

``` r
rose <- windscape_example("wind_rose")
head(ggplot2::fortify(rose))
#>           x  y          SW          W         NW          N         NE
#> 1 -120.0000 50 0.001837510 0.01303786 0.04689814 0.03452373 0.02867673
#> 2 -119.6842 50 0.003325526 0.01249226 0.04783026 0.04005871 0.02012151
#> 3 -119.3684 50 0.003556589 0.01323038 0.04368901 0.04622269 0.01887328
#> 4 -119.0526 50 0.003929742 0.01285706 0.03766154 0.05081339 0.02046587
#> 5 -118.7368 50 0.002442006 0.01549426 0.03567446 0.05171587 0.02091837
#> 6 -118.4211 50 0.002181489 0.02034521 0.03644594 0.05170206 0.01930627
#>            E         SE           S     total    speed      net  bearing
#> 1 0.13035671 0.02295375 0.006273328 0.2845578 8.895296 3.901069 44.44692
#> 2 0.11763944 0.02839319 0.004975851 0.2748367 8.714891 3.427445 42.85106
#> 3 0.10748533 0.03674658 0.003969136 0.2737730 8.817199 3.257523 45.75774
#> 4 0.09756345 0.04342743 0.003780435 0.2704989 8.848722 3.172106 50.06727
#> 5 0.08117093 0.04646795 0.005902223 0.2597861 8.643859 2.865233 48.87757
#> 6 0.06930302 0.04561561 0.007976120 0.2528757 8.475258 2.498716 43.23467
#>   consistency
#> 1   0.4385542
#> 2   0.3932861
#> 3   0.3694510
#> 4   0.3584819
#> 5   0.3314762
#> 6   0.2948248
```
