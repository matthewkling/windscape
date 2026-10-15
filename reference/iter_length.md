# Iteration step length of a random walk

Iteration step length of a random walk

## Usage

``` r
iter_length(x)
```

## Arguments

- x:

  The list returned by
  [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md),
  or one of its elements.

## Value

Length of one iteration (time step), in hours if `trans = 1` and wind
speeds are in m/s.

## Examples

``` r
rose <- windscape_example("wind_rose")
rw <- random_walk(rose, cbind(-100, 40), iter = 10)
#>  iteration timestep: ~0.935
#>  simulation duration: ~9.35
#> (these values are in hours, IF trans == 1 and wind_field units are m/s)
iter_length(rw)
#> [1] 0.9349788
```
