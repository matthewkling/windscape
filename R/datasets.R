#' Silver birch landscape genetic data from Tsuda et al. (2017)
#'
#' An example landscape genetic dataset for the species silver birch (Betula pendula) across 30 sampling
#' sites in Asia, originally published by Tsuda et al. (2017).
#'
#' @format `birch`
#' A list with 4 entries:
#' \describe{
#'   \item{sites}{A matrix with two columns giving the longitude and latitude of each site.}
#'   \item{div}{A numeric vector giving the allelic richness sampled at each site.}
#'   \item{mig}{A square, asymmetric matrix with estimated gene flow rates for each pair of sites.}
#'   \item{fst}{A square, symmetric matrix with Fst values for each pair of sites.}
#' }
#' @source Y. Tsuda, V. Semerikov, F. Sebastiani, G. G. Vendramin, M. Lascoux, Multispecies genetic
#'    structure and hybridization in the Betula genus across Eurasia. Molecular Ecology 26, 589-605 (2017).
#'    <https://doi.org/10.1111/mec.13885>
"birch"



#' Example wind data sets
#'
#' Load one of the example data sets shipped with windscape, as a ready-to-use object.
#'
#' All three are derived from 10 m hourly wind from the Climate Forecast System Reanalysis (CFSR),
#' on CFSR's native grid of approximately 0.32 degrees:
#' \describe{
#'   \item{"wind_series"}{A \code{wind_series} for the western and central United States
#'     (120-90 W, 30-50 N), with 96 time steps: every sixth hour on the 1st and 15th of each month
#'     in 2000. Values are rounded to 0.1 m/s. This is for demonstrating the
#'     \code{wind_series} -> \code{wind_rose} workflow; a real analysis would use a longer
#'     and denser time series.}
#'   \item{"wind_rose"}{A \code{wind_rose} (with \code{trans = 1}) for the same region, built from
#'     all 576 hours on those days. Because it uses six times as many observations as
#'     \code{"wind_series"}, it is a better estimate of the region's wind regime, and is the
#'     better choice for demonstrating connectivity analyses.}
#'   \item{"wind_field"}{A \code{wind_field} for Hurricane Katrina in the Gulf of Mexico, at
#'     2005-08-28 23:00 UTC (99-78 W, 17-35 N), rounded to 0.1 m/s.}
#' }
#'
#' @param name Which data set to load: "wind_series", "wind_rose", or "wind_field".
#' @return An object of the class named by \code{name}.
#' @source Climate Forecast System Reanalysis (Saha et al. 2010), via the NSF NCAR Research
#'    Data Archive.
#' @examples
#' ws <- windscape_example("wind_series")
#' rose <- windscape_example("wind_rose")
#' katrina <- windscape_example("wind_field")
#' @export
windscape_example <- function(name = c("wind_series", "wind_rose", "wind_field")){
      name <- match.arg(name)
      file <- c(wind_series = "wind_usa.tif", wind_rose = "rose_usa.tif",
                wind_field = "katrina.tif")[[name]]
      x <- terra::rast(system.file("extdata", file, package = "windscape", mustWork = TRUE))
      switch(name,
             wind_series = wind_series(x),
             wind_rose = as_wind_rose(x, trans = 1, n_steps = 576),
             wind_field = wind_field(x))
}
