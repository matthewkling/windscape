% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mantel.R
\name{mantel_test}
\alias{mantel_test}
\title{Mantel test}
\usage{
mantel_test(
  x,
  y,
  z = NULL,
  nperm = 999,
  method = "pearson",
  alternative = c("two.sided", "less", "greater")
)
}
\arguments{
\item{x, y}{Square matrices (can be asymmetric).}

\item{z}{Optional list of one or more square control matrices if partial Mantel is desired.}

\item{nperm}{Number of random permutations to perform (positive integer).}

\item{method}{Correlation method (see \link[stats]{cor} for details).}

\item{alternative}{A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".}
}
\value{
A list with the following components:
\itemize{
 \item{`stat`: }{The correlation coefficient between \code{x} and \code{y}, measured on the data provided. This will be a partial correlation if \code{z} is specified.}
 \item{`quantile`: }{The quantile of the observed stat in the null distribution.}
 \item{`p.value`: }{Significance, derived from the combination of \code{quantile} and \code{alternative}.}
 \item{`perm`: }{A vector of null correlations based on permuted data.}
}
}
\description{
This version of the Mantel test supports asymmetric matrices and multiple partial-Mantel controls,
a combination of features that is not possible using other widely available Mantel functions such
as \link[vegan]{mantel.partial}.
}
\examples{
# A convenience function to generate symmetric matrices:
sym_matrix <- function(x){
   n <- which(cumsum(1:100) == length(x)) + 1
   m <- matrix(0, n, n)
   m[lower.tri(m)] <- x
   m[upper.tri(m)] <- t(m)[upper.tri(t(m))]
   m
}

# Simulate some example data:
g <- runif(45)
w <- runif(45)
d <- runif(45)
e <- sample(0:1, 45, T)
G <- sym_matrix(g) # e.g. a genetic matrix
W <- sym_matrix(w) # e.g. a wind matrix
D <- sym_matrix(d) # e.g. a distance matrix
E <- sym_matrix(e) # e.g. an environment matrix

# Basic function usage:
mantel_test(G, W)                      # simple first-order Mantel
mantel_test(G, W, list(D, E))          # partial Mantel
mantel_test(G, W, method = "kendall")  # use of alternative test statistic

# Demonstrate that the function matches vegan::mantel.partial output:
# (using only symmetric matrices, and just one control variable, since
# vegan does not work with asymmetric matrices or multiple controls)
a <- mantel_test(G, W, list(D), nperm = 99999)
b <- vegan::mantel.partial(G, W, D, permutations = 99999)
c(a$stat, b$statistic)
c(1 - a$quantile, b$signif) # will differ slightly due to randomization
}
