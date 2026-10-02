# Convert asymmetric pairwise matrix to symmetrical matrix of pairwise means

This function produces a symmetrical matrix with values representing the
pairwise means of directional flows in an asymmetric in put matrix. The
output can be used in tests about the strength of direction-agnostic
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
Proceedings of the National Academy of Sciences.
https://doi.org/10.1073/pnas.2017317118
