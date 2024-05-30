#' costDistance with local interpolation
#'
#' Increasing the spatial resolution of a wind data set (e.g. by interpolation using \link{downscale})
#' can greatly improve wind connectivity estimates among nearby sites, but can make it computationally
#' impossible to model connectivity over large geographic regions, a trade-off that presents problems
#' for studies that include both nearby and distant site pairs.
#'
#' This function gets around this issue with a hybrid approach, using a broad-scale coarse-resolution
#' wind grid to model connectivity among distant sites, and a separate local-scale high-resolution
#' interpolated grid to model connectivity between each pair of nearby sites.
#'
#' @param rose wind_rose.
#' @param ll longitude and latitude of site locations, as a two-column matrix.
#' @param threshold_km a positive number representing the distance threshold, in kilometers. Site pairs
#' closer together than this distance will get a separate local high-resolution connectivity model, while pairs
#' farther apart will be modeled at the resolution of \code{rose}.
#' @param pad a positive number indicating how far beyond a pair of sites the modeling domain should extend
#' for local models. This is an expansion factor giving the padding as a fraction of the maximum latitudinal or
#' longitudinal distance between the two sites. The default of 1 is reasonable in most cases; smaller values will
#' increase computational speed but may fail to account for wind routes beyond the bounding box encompassing the
#' site pair.
#' @param max_nodes an integer representing how finely to interpolate wind grids for local site pairs. This represents
#' the number of cells in the high-resolution interpolated grid for each site pair; larger values remove artifacts of
#' the discrete coarse grid more effectively, but have increased computational cost.
#' @param direction either "downwind" or "upwind", indicating whether outbound or inbound connectivity should be computed.
#' In this context, changing this parameter is equivalent to transposing the resulting wind distance matrix.
#' @param method disaggregation method; either "near" or "bilinear'; see \link[terra]{disagg} for details.
#' @return A list of square matrices:
#' \describe{
#'  \item{wind_dist}{Wind cost-distances between site pairs, in hours if \code{rose} has a p = 1. Computed using \link[gdistance]{costDistance}.}
#'  \item{wind_dist_coarse}{Wind cost-distances using the coarse input wind grid (\code{rose}); this will differ from \code{wind_dist} only for site pairs closer than \code{threshold_km}. This is provided for reference, to judge how higher-resolution models impact results relative to the coarse grid.}
#'  \item{point_dist}{Distances between sites, in km.}
#'  \item{cell_dist}{Distances between the centroids of the grid cells used to model wind connectivity, in km. For site pairs where this distance deviates from \code{point_dist} by a substantial percentage, \code{wind_dist} estimates may contain nontrivial rounding noise resulting from the discrete grid.}
#'  \item{cell_dist_coarse}{Distances between the cell centroids in the coarse-resolution input wind dataset (\code{rose}), in km. This is provided for reference, to judge improvements relative to \code{cell_dist}.}
#' }
#' @export
cdli <- function(rose, ll, threshold_km = 30, pad = 1, max_nodes = 1e6, direction = "downwind", method = "bilinear"){

      # point distances
      pdist <- geosphere::distm(ll) / 1000

      # global connectivity and cell distances
      message("... computing global connectivity ...")
      wdist <- wdist0 <- rose %>%
            wind_graph(direction) %>% # convert to connectivity graph
            costDistance(ll) # calculate travel times between sites
      cdist <- cdist0 <- cell_distance(rose, ll)

      # bookkeeping
      npw <- sum(pdist[upper.tri(pdist)] <= threshold_km)
      if(npw > 0){
            message("... computing pairwise connectivity for ", npw, " site pairs closer than ", threshold_km, " km ...")
            pb <- txtProgressBar(min = 0, max = npw, initial = 0, style = 3)
            pwi <- 0

            # high-res connectivity for nearby site pairs
            for(i in 1:nrow(ll)){
                  for(j in 1:nrow(ll)){
                        if(j >= i) next()
                        if(pdist[i, j] > threshold_km) next()

                        # bookkeeping
                        setTxtProgressBar(pb, pwi)
                        pwi <- pwi + 1
                        wi <- c(i, j)
                        lli <- ll[wi,]

                        # crop wind rose to local neighborhood of at least a few cells
                        expn <- lli %>% apply(2, range) %>% apply(2, diff) %>% max() %>% "*"(pad)
                        exp <- max(expn, res(rose)[2]*2)
                        xt <- extent(t(apply(lli, 2, range))) + c(-exp, exp, -exp, exp)
                        ri <- rose %>% crop(xt)

                        # interpolate
                        fact <- floor(sqrt(max_nodes / prod(dim(ri[[1]]))))
                        ri <- ri %>% downscale(fact, method = method)

                        # crop to immediate region around focal sites
                        exp <- max(expn, res(ri)[2]*2)
                        xt <- extent(t(apply(lli, 2, range))) + c(-exp, exp, -exp, exp)
                        ri <- ri %>% crop(xt)

                        # wind connectivity, cell distance, and pairwise flag
                        wdist[wi, wi] <- ri %>% wind_graph(direction) %>% costDistance(lli)
                        cdist[wi, wi] <- cell_distance(ri, lli)
                  }
            }
            setTxtProgressBar(pb, npw)
            close(pb)
      }

      return(list(wind_dist = wdist,
                  wind_dist_coarse = wdist0,
                  point_dist = pdist,
                  cell_dist = cdist,
                  cell_dist_coarse = cdist0))
}
