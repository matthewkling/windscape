
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



#' Hurricane Katrina wind field data
#'
#' An example of the u and v components of a wind field raster, from CFSR for 2005-08-28 23:00 UTC,
#' when Hurricane Katrina was in the Gulf of Mexico.
#'
#' @format `katrina`
#' A SpatRaster with layers for u and v wind speeds, in m/s. This is a "wrapped" SpatRaster, and needs to
#' be unpacked with \code{unwrap(katrina)} to be used.
#' @source Climate Forecast System Reanalysis
"katrina"



#' Wind time series for the contiguous US
#'
#' An example of a wind dataset that can be used to create a \code{wind_series}. This is CFSR data for part
#' of the central United States for a smattering of days in year 2000. A real analysis would likely require
#' a denser or longer time series, so this is useful only as an example.
#'
#' @format `cfsr_usa`
#' A SpatRaster with layers for u and v wind speeds for multiple time steps, in m/s. This is a "wrapped"
#' SpatRaster, and needs to be unpacked with \code{unwrap(cfsr_usa)} to be used.
#' @source Climate Forecast System Reanalysis
"cfsr_usa"
