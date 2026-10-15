# Pairwise distances between points

A wrapper around [distm](https://rdrr.io/pkg/geosphere/man/distm.html),
returning distances in km.

## Usage

``` r
point_distance(ll)
```

## Arguments

- ll:

  two-column matrix of point coordinates

## Value

Pairwise distances between the points, in km

## Examples

``` r
sites <- cbind(c(-110, -100, -95), c(35, 40, 45))
point_distance(sites)
#>          [,1]      [,2]      [,3]
#> [1,]    0.000 1042.9516 1690.1820
#> [2,] 1042.952    0.0000  690.6346
#> [3,] 1690.182  690.6346    0.0000
```
