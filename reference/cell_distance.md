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
