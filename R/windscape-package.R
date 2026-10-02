#' @keywords internal
"_PACKAGE"

#' @useDynLib windscape, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import terra
#' @import methods
#' @importFrom magrittr %>%
#' @importFrom stats complete.cases cor lm residuals runif setNames
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importClassesFrom gdistance TransitionLayer
#' @importMethodsFrom gdistance geoCorrection
NULL
