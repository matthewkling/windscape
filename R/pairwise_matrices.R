#' Convert data to reciprocally symmetrical pairwise matrix
#'
#' This function produces a "pairwise ratio matrix" that can be used in tests about asymmetric flows among nodes in a data set
#' (e.g. the "asymmetry" and "diversity" hypotheses; Kling and Ackerly (2021)). The input data set can be an asymmetric matrix with entries representing directed edge weights (e.g. directional flows of wind or genes connecting each pair of populations or "nodes"),
#' or it can be a vector of node attributes (e.g. properties of a site or population, such as genetic diversity).
#'
#' @param x Either a square matrix with entries representing directional edge weights between pairs of nodes,
#' or a numeric vector with a value for each node.
#' @param log Should ratios be log-transformed? Default is \code{TRUE}, since this generates a symmetrical frequency distribution desirable for most analyses.
#' @return A square matrix the same size as \code{x}, with ones in the diagonal and reciprocal values in the upper and lower triangles.
#' If \code{x} is a matrix, the output represents the ratios of above-diagonal to below-diagonal values in \code{x}.
#' If \code{x} is a vector, the output represents the ratios of every pairwise combination of \code{x} values; e.g. a matrix in which the value of cell [2,4] equals x[2] / x[4].
#' @references Kling and Ackerly (2021). Global wind patterns shape genetic differentiation, asymmetric gene flow, and genetic diversity in trees. Proceedings of the National Academy of Sciences. https://doi.org/10.1073/pnas.2017317118
#' @export
pairwise_ratios <- function(x, log = TRUE){
      if(!inherits(x, "matrix")) x <- matrix(rep(x, length(x)), length(x))
      x <- x / t(x)
      if(log) x <- log(x)
      return(x)
}

#' Convert asymmetric pairwise matrix to symmetrical matrix of pairwise means
#'
#' This function produces a symmetrical matrix with values representing the pairwise means of directional flows in an asymmetric in put matrix.
#' The output can be used in tests about the strength of direction-agnostic connectivity among nodes in a data set (e.g. the "isolation" hypothesis; Kling and Ackerly (2021)).
#'
#' @param x A square matrix with entries representing directional edge weights between pairs of nodes.
#' @return A symmetrical square matrix the same size as \code{x}, with values representing the pairwise means of directional flows in \code{x}.
#' @references Kling and Ackerly (2021). Global wind patterns shape genetic differentiation, asymmetric gene flow, and genetic diversity in trees. Proceedings of the National Academy of Sciences. https://doi.org/10.1073/pnas.2017317118
#' @export
pairwise_means <- function(x){
      x <- (x + t(x)) / 2
      return(x)
}
