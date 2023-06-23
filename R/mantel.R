#' Mantel test
#'
#' This version of the Mantel test supports asymmetric matrices and multiple partial-Mantel controls,
#' a combination of features that is not possible using other widely available Mantel functions such
#' as \link[vegan]{mantel.partial}.
#'
#' @param x,y Square matrices (can be asymmetric).
#' @param z Optional list of one or more square control matrices if partial Mantel is desired.
#' @param nperm Number of random permutations to perform (positive integer)
#' @param method Correlation method (see \link[stats]{cor} for details)
#' @return A list with the following components:
#' \itemize{
#'  \item{`stat`: }{The correlation coefficient between \code{x} and \code{y}, measured on the data provided. This will be a partial correlation if \code{z} is specified.}
#'  \item{`quantile`: }{The quantile of the observed stat in the null distribution.}
#' }
#' @examples
#' # A convenience function to generate symmetric matrices:
#' sym_matrix <- function(x){
#'    n <- which(cumsum(1:100) == length(x)) + 1
#'    m <- matrix(0, n, n)
#'    m[lower.tri(m)] <- x
#'    m[upper.tri(m)] <- t(m)[upper.tri(t(m))]
#'    m
#' }
#'
#' # Simulate some example data:
#' g <- runif(45)
#' w <- runif(45)
#' d <- runif(45)
#' e <- sample(0:1, 45, T)
#' G <- sym_matrix(g) # e.g. a genetic matrix
#' W <- sym_matrix(w) # e.g. a wind matrix
#' D <- sym_matrix(d) # e.g. a distance matrix
#' E <- sym_matrix(e) # e.g. an environment matrix
#'
#' # Basic function usage:
#' mantel_test(G, W)                      # simple first-order Mantel
#' mantel_test(G, W, list(D, E))          # partial Mantel
#' mantel_test(G, W, method = "kendall")  # use of alternative test statistic
#'
#' # Demonstrate that the function matches vegan::mantel.partial output:
#' # (using only symmetric matrices, and just one control variable, since
#' # vegan does not work with asymmetric matrices or multiple controls)
#' a <- mantel_test(G, W, list(D), nperm = 99999)
#' b <- vegan::mantel.partial(G, W, D, permutations = 99999)
#' c(a$stat, b$statistic)
#' c(1 - a$quantile, b$signif) # will differ slightly due to randomization
#' @export
mantel_test <- function(x, y, # square matrices (can be asymmetric)
                        z = NULL, # optional: list of square control matrices if partial Mantel is desired
                        nperm = 999,
                        method = "pearson"){
      diag(x) <- NA
      diag(y) <- NA

      tri <- upper.tri(x) | lower.tri(x)
      R1 <- x
      R2 <- y

      if(!is.null(z)){
            z <- lapply(z, function(m){diag(m) <- NA; return(m)})
            R1[tri] <- residuals(lm(as.vector(x) ~ sapply(z, as.vector)))
            R2[tri] <- residuals(lm(as.vector(y) ~ sapply(z, as.vector)))
      }

      stat <- cor(R1[tri], R2[tri], method = method)

      perm <- rep(NA, nperm)
      for(i in 1:nperm){
            p <- sample(ncol(R1), ncol(R1))
            perm[i] <- cor(R1[p, p][tri], R2[tri], method = method)
      }

      list(stat = stat,
           quantile = mean(stat > perm))
}
