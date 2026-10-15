# Convert asymmetric pairwise matrix to reciprocally symmetrical matrix

This function produces a "pairwise ratio matrix" that can be used in
tests about asymmetric flows among nodes in a data set (e.g. the
"asymmetry" and "diversity" hypotheses; Kling and Ackerly (2021)). The
input data set can be an asymmetric matrix with entries representing
directed edge weights (e.g. directional flows of wind or genes
connecting each pair of populations or "nodes"), or it can be a vector
of node attributes (e.g. properties of a site or population, such as
genetic diversity). The result is a matrix that represents each pair's
relative asymmetry, removing information about their absolute
magnitudes.

## Usage

``` r
pairwise_ratios(x, log = TRUE)
```

## Arguments

- x:

  Either a square matrix with entries representing directional edge
  weights between pairs of nodes, or a numeric vector with a value for
  each node.

- log:

  Should ratios be log-transformed? Default is `TRUE`, since this
  generates a symmetrical frequency distribution desirable for most
  analyses.

## Value

A square matrix the same size as `x`, with ones in the diagonal (zeros
if `log = TRUE`) and reciprocal values in the upper and lower triangles.
The diagonal is set regardless of the diagonal of `x`, which for a
matrix of flows is typically zero. A pair with zero flow in one
direction gives a ratio of 0 or `Inf` (`-Inf` or `Inf` if logged), and a
pair with no flow in either direction gives `NaN`. If `x` is a matrix,
the output represents the ratios of above-diagonal to below-diagonal
values in `x`. If `x` is a vector, the output represents the ratios of
every pairwise combination of `x` values; e.g. a matrix in which the
value of cell `[2, 4]` equals `x[2] / x[4]`.

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
pairwise_ratios(x)
#>            [,1]       [,2]      [,3]
#> [1,]  0.0000000 -0.6931472  1.098612
#> [2,]  0.6931472  0.0000000 -1.098612
#> [3,] -1.0986123  1.0986123  0.000000

# ratios of a site attribute, such as genetic diversity
pairwise_ratios(c(0.2, 0.4, 0.8))
#>           [,1]       [,2]       [,3]
#> [1,] 0.0000000 -0.6931472 -1.3862944
#> [2,] 0.6931472  0.0000000 -0.6931472
#> [3,] 1.3862944  0.6931472  0.0000000
```
