# Pairwise distances between cell centroids

Calculate the pairwise distances between the centroids of the raster
grid cells that a set of points fall into. This is useful as an
alternative to the distances between the points themselves, if distances
are going to be compared to other metrics that are based on cell centers
(e.g. wind connectivity).

## Usage

``` r
cell_distance(x, ll)
```

## Arguments

- x:

  SpatRaster

- ll:

  two-column matrix of site coordinates

## Value

Pairwise distances between the cell centroids, in km

## Examples

``` r
rose <- windscape_example("wind_rose")
sites <- cbind(c(-110, -100, -95), c(35, 40, 45))
cell_distance(rose, sites)
#>          [,1]      [,2]      [,3]
#> [1,]    0.000 1013.6924 1668.6407
#> [2,] 1013.692    0.0000  700.7436
#> [3,] 1668.641  700.7436    0.0000
```
