# Convert asymmetric pairwise matrix to symmetrical pairwise means

This function produces a symmetrical matrix with values representing the
pairwise means of directional flows in an asymmetric input matrix. The
output can be used in tests about the **strength** of direction-agnostic
connectivity among nodes in a data set (e.g. the "isolation" hypothesis;
Kling and Ackerly (2021)).

## Usage

``` r
pairwise_means(x)
```

## Arguments

- x:

  A square matrix with entries representing directional edge weights
  between pairs of nodes.

## Value

A symmetrical square matrix the same size as `x`, with values
representing the pairwise means of directional flows in `x`.

## References

Kling and Ackerly (2021). Global wind patterns shape genetic
differentiation, asymmetric gene flow, and genetic diversity in trees.
Proceedings of the National Academy of Sciences 118: e2017317118.
[doi:10.1073/pnas.2017317118](https://doi.org/10.1073/pnas.2017317118)

## Examples

``` r
# directional flows among three sites, from row to column
x <- matrix(c(0, 2, 6,
              4, 0, 1,
              2, 3, 0), 3, byrow = TRUE)
pairwise_means(x)
#>      [,1] [,2] [,3]
#> [1,]    0    3    4
#> [2,]    3    0    2
#> [3,]    4    2    0
```
