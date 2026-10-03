#' Pairwise random walk connectivity among sites
#'
#' For every pair of sites, computes how strongly particles released at one site reach the other
#' under a stream-mode random walk (see [random_walk()]): by default, the probability that a
#' particle released at the first site is deposited at the second. This is the random walk
#' counterpart to [pairwise_least_cost()].
#'
#' @param rose A `wind_rose`.
#' @param sites A two-column matrix (or data frame) of site coordinates, or a `SpatVector` of
#'    points.
#' @param half_life Half-life of airborne mass, in hours of transport time; see [random_walk()].
#'    Required for `value = "deposition"`, which is zero without deposition.
#' @param value `"deposition"` (the default) or `"residence"`; see Value.
#' @param density Logical: divide by the area of the destination site's grid cell, giving values
#'    per km^2? Default `TRUE`. Without this, values depend on cell size: a larger cell receives
#'    more of the particles passing nearby. On a longitude/latitude grid, cell area shrinks
#'    toward the poles (by about 25 percent from 30 to 50 degrees latitude), so raw values are
#'    biased toward lower-latitude destinations even within a single analysis. Dividing by area
#'    removes that bias, and also makes values comparable across grid resolutions (though see
#'    Details).
#' @param timescale,latitude_correction See [random_walk()]. Results do not depend on
#'    `timescale`.
#' @param chunk Number of sites to solve for at once. Larger values are faster but use more
#'    memory.
#'
#' @return A square matrix with one row and column per site, giving connectivity from the row's
#'    site (the origin) to the column's site (the destination). Large values mean strong
#'    connectivity, the opposite orientation from the travel times of [pairwise_least_cost()].
#'    With `value = "deposition"`, entries are the probability that a particle released at the
#'    origin is deposited in the destination's grid cell (per km^2 if `density = TRUE`). With
#'    `value = "residence"`, they are the time, in hours, that a unit release at the origin
#'    spends airborne over the destination's cell (per km^2 if `density = TRUE`). Because wind
#'    is directional, the matrix is generally asymmetric. The diagonal is each site's
#'    self-connectivity: release that is deposited (or stays airborne) in its own cell, which
#'    can be large; see [rw_self_retention()]. Sites in the same grid cell get identical rows and
#'    columns.
#'
#' @details
#' Each origin's row is the destination values of a stream-mode random walk released from that
#' origin, so it equals the result of `random_walk(rose, origin, mode = "stream", ...)` read at
#' each destination. All origins share one matrix factorization, so the cost grows slowly with
#' the number of sites.
#'
#' `density = TRUE` removes the direct effect of cell size on how much a destination receives,
#' but the random walk's spread itself depends somewhat on grid resolution (see
#' [random_walk()]), so values are not fully resolution-independent.
#'
#' @examples
#' rose <- windscape_example("wind_rose")
#' sites <- cbind(c(-110, -105, -100, -95), c(40, 42, 38, 44))
#' pairwise_random_walk(rose, sites, half_life = 24)
#' @export
pairwise_random_walk <- function(rose, sites, half_life = NULL,
                                 value = c("deposition", "residence"), density = TRUE,
                                 timescale = 1, latitude_correction = TRUE, chunk = 200){
      if(!inherits(rose, "wind_rose")) stop("`rose` must be a wind_rose")
      value <- match.arg(value)
      if(inherits(sites, "SpatVector")) sites <- terra::crds(sites)
      sites <- as.matrix(sites)
      if(ncol(sites) != 2 || !is.numeric(sites)) stop("`sites` must be a two-column matrix of coordinates")
      if(is.null(half_life)){
            if(value == "deposition") stop("`half_life` is required for `value = \"deposition\"`")
            half_life <- Inf
      }
      if(!(timescale > 0 && timescale <= 1)) stop("'timescale' must be greater than 0 and less than or equal to 1.")

      if(latitude_correction) rose <- rw_latitude_correction(rose)
      t <- rw_max_step(rose) * timescale
      lambda <- rw_decay(half_life, t)
      if(value == "deposition" && lambda == 0) stop("`value = \"deposition\"` requires a finite `half_life`")
      cells <- terra::cellFromXY(rose, sites)
      if(anyNA(cells)) stop("some `sites` fall outside the extent of `rose`")
      P <- rw_matrix(rw_prob(rose, t))
      valid <- attr(P, "valid")
      if(any(!valid[cells])) stop("some `sites` fall in grid cells that are NA in `rose`")
      if(lambda == 0) rw_check_drainage(P)

      # one stream solve per distinct origin cell, sharing a factorization: column s of
      # (I - (1 - lambda) P')^-1 is the per-step steady state from a unit release at s
      N <- nrow(P)
      lu <- Matrix::lu(Matrix::Diagonal(N) - (1 - lambda) * Matrix::t(P))
      origins <- unique(cells)
      M <- matrix(0, length(origins), length(origins))
      for(s in split(seq_along(origins), ceiling(seq_along(origins) / chunk))){
            B <- Matrix::sparseMatrix(i = origins[s], j = seq_along(s), x = 1, dims = c(N, length(s)))
            X <- as.matrix(Matrix::solve(lu, B))
            M[s, ] <- t(X[origins, , drop = FALSE])
      }
      M <- M * if(value == "deposition") lambda else (1 - lambda) * t

      if(density){
            area <- terra::values(rw_cell_area(rose))[cells, 1]
      }
      idx <- match(cells, origins)
      out <- M[idx, idx, drop = FALSE]
      if(density) out <- sweep(out, 2, area, "/")
      if(!is.null(rownames(sites))) dimnames(out) <- list(rownames(sites), rownames(sites))
      out
}
