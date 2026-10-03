#' Trace particle trails through a wind field
#'
#' Traces the paths of particles carried by a wind field, upwind and/or downwind from a set of
#' starting points. Each step, particles move in the direction of the local wind (interpolated
#' bilinearly), either by the distance the wind carries them in that time (`hours`), or by a
#' fixed distance (`distance`), giving trails that show direction only. Trails through a single
#' wind field are streamlines.
#'
#' @param x A `wind_field`, created with [wind_field()].
#' @param seeds Starting points: either a two-column matrix of longitude and latitude, or a
#'    single integer giving the number of points to generate with [generate_particles()], in
#'    which case further arguments to that function can be passed in `...`.
#' @param hours Transport time in hours spanned by each trail, assuming wind speeds in `x` are
#'    in m/s. Particles move at the local wind speed, so trail length is proportional to speed.
#'    Supply exactly one of `hours` or `distance`.
#' @param distance Length in km of each trail. Particles move at the same speed regardless of
#'    wind speed, so trails show only direction. Supply exactly one of `hours` or `distance`.
#' @param steps Number of integration steps per trail. More steps give smoother and more
#'    accurate trails. Default 100.
#' @param direction Trace trails `"both"` ways from each seed (the default; half of `hours` or
#'    `distance` upwind and half downwind), only `"downwind"`, or only `"upwind"`.
#' @param wrap Wrap particles that leave the field across its edges back in on the opposite
#'    side: `"neither"` (the default), `"horizontal"` (e.g. for global fields), `"vertical"`, or
#'    `"both"`. Otherwise, trails end where they leave the field.
#' @param sf Logical: return trails as an `sf` object of linestrings instead of a data frame?
#'    Requires the sf package. Default `FALSE`.
#' @param ... Further arguments to [generate_particles()], used only if `seeds` is an integer.
#'
#' @return A data frame with one row per particle position, ordered along each trail from upwind
#'    to downwind:
#' * `trail`: trail id. Each particle's path is one trail, unless it wraps across the field's
#'   edge, which starts a new trail so that lines don't cross the map.
#' * `particle`: particle id, the row of `seeds` it started from.
#' * `step`: integration step, negative upwind of the seed and positive downwind.
#' * `hours` (if `hours` was given): signed transport time from the seed. Or `km` (if
#'   `distance` was given): signed distance from the seed.
#' * `x`, `y`: longitude and latitude.
#' * `speed`: wind speed at the particle's position, in the units of `x`.
#'
#' With `sf = TRUE`, an `sf` object with one linestring per trail, with `trail` and `particle`
#' columns and the elapsed `hours` or `km` as the linestrings' M coordinate.
#'
#' @examples
#' katrina <- windscape_example("wind_field")
#' seeds <- cbind(c(-92, -86), c(24, 30))
#'
#' # where does air at these points come from and go to, over 12 hours?
#' tr <- wind_trails(katrina, seeds, hours = 12)
#' head(tr)
#'
#' library(ggplot2)
#' ggplot(tr, aes(x, y)) +
#'   geom_wind_trail(aes(colour = speed)) +
#'   coord_quickmap()
#' @export
wind_trails <- function(x, seeds, hours = NULL, distance = NULL, steps = 100,
                        direction = c("both", "downwind", "upwind"),
                        wrap = c("neither", "horizontal", "vertical", "both"), sf = FALSE, ...){

      if(!inherits(x, "wind_field")) stop("`x` must be a `wind_field` object.")
      direction <- match.arg(direction)
      wrap <- match.arg(wrap)
      if(is.null(hours) == is.null(distance)) stop("supply exactly one of `hours` or `distance`")
      total <- if(is.null(hours)) distance else hours
      if(length(total) != 1 || !is.finite(total) || total <= 0)
            stop("`hours` or `distance` must be a single positive number")
      if(length(steps) != 1 || !is.finite(steps) || steps < 1) stop("`steps` must be a number >= 1")
      fixed <- !is.null(distance)

      if(length(seeds) == 1 && is.numeric(seeds)) seeds <- generate_particles(x, n = seeds, ...)
      seeds <- as.matrix(seeds)
      if(ncol(seeds) != 2) stop("`seeds` must be a two-column matrix of coordinates, or an integer")

      dirs <- if(direction == "both") c("downwind", "upwind") else direction
      n_step <- max(1, round(steps / length(dirs))) # steps per direction
      dt <- total / length(dirs) / n_step           # hours, or km, per step

      xlim <- as.vector(terra::ext(x))[1:2]
      ylim <- as.vector(terra::ext(x))[3:4]
      out <- vector("list", length(dirs))

      for(di in seq_along(dirs)){
            s <- if(dirs[di] == "downwind") 1 else -1
            p <- seeds
            wraps <- rep(0, nrow(p))
            rec <- vector("list", n_step + 1)
            for(i in 1:(n_step + 1)){
                  uv <- as.matrix(terra::extract(x, p, method = "bilinear"))
                  spd <- sqrt(rowSums(uv^2))
                  rec[[i]] <- data.frame(particle = seq_len(nrow(p)), wraps = wraps,
                                         step = (i - 1) * s, x = p[, 1], y = p[, 2], speed = spd)
                  if(i > n_step) break
                  # displacement in km: wind speed (m/s) x time, or a fixed distance
                  km <- if(fixed) uv / ifelse(spd > 0, spd, Inf) * dt else uv * 3.6 * dt
                  km <- km * s
                  p <- p + cbind(km[, 1] / KM_LAT * aspect(p[, 2]), km[, 2] / KM_LAT)

                  if(wrap %in% c("horizontal", "both")){
                        lo <- which(p[, 1] < xlim[1])
                        hi <- which(p[, 1] > xlim[2])
                        p[lo, 1] <- p[lo, 1] + diff(xlim)
                        p[hi, 1] <- p[hi, 1] - diff(xlim)
                        wraps[c(lo, hi)] <- wraps[c(lo, hi)] + s
                  }
                  if(wrap %in% c("vertical", "both")){
                        lo <- which(p[, 2] < ylim[1])
                        hi <- which(p[, 2] > ylim[2])
                        p[lo, 2] <- p[lo, 2] + diff(ylim)
                        p[hi, 2] <- p[hi, 2] - diff(ylim)
                        wraps[c(lo, hi)] <- wraps[c(lo, hi)] + s
                  }
            }
            out[[di]] <- do.call(rbind, rec)
      }

      b <- do.call(rbind, out)
      b <- b[!(duplicated(b[, c("particle", "step")])), ] # the seed is recorded in both directions
      b <- b[stats::complete.cases(b), ]
      b <- b[order(b$particle, b$step), ]
      b$trail <- as.integer(factor(paste(b$particle, b$wraps), levels = unique(paste(b$particle, b$wraps))))
      b[[if(fixed) "km" else "hours"]] <- b$step * dt
      b <- b[, c("trail", "particle", "step", if(fixed) "km" else "hours", "x", "y", "speed")]
      rownames(b) <- NULL

      if(sf) b <- trails_to_sf(b)
      b
}




#' Convert particle trails data frame to sf
#'
#' @param x Data frame generated by [wind_trails()].
#' @noRd
trails_to_sf <- function(x){
      if(!requireNamespace("sf", quietly = TRUE)) stop("the sf package is required for `sf = TRUE`")
      m <- intersect(c("hours", "km"), names(x))
      trails <- split(x, x$trail)
      lines <- lapply(trails, function(d) sf::st_linestring(as.matrix(d[, c("x", "y", m)]), dim = "XYM"))
      sf::st_sf(trail = as.integer(names(trails)),
                particle = vapply(trails, function(d) d$particle[1], integer(1)),
                geometry = sf::st_sfc(lines, crs = 4326))
}


#' Generate spatial grid of points
#'
#' @param x A raster or extent
#' @param n Number of sample points (approximate)
#' @noRd
geo_grid <- function(x, n){
      dx <- xmax(x) - xmin(x)
      dy <- ymax(x) - ymin(x)
      ar <- dy / dx
      nx <- round(sqrt(n / ar))
      ny <- round(ar * sqrt(n / ar))
      px <- rep(seq(xmin(x), xmax(x), length.out = nx + 1)[1:nx] + dx / nx / 2, ny)
      py <- rep(seq(ymin(x), ymax(x), length.out = ny + 1)[1:ny] + dy / ny / 2, each = nx)
      cbind(x = px, y = py)
}

#' Generate initial particle locations
#'
#' Initialize a set of particle locations for use in [wind_trails()].
#'
#' @param x A wind_field object.
#' @param n Positive integer giving the total number of particles in the grid. If \code{sample = "random"} the result
#'    will have this exact number of particles, while if \code{sample = "grid"} it will be approximate.
#' @param sample Sampling scheme, either "random" (the default) to generate points at random locations, or "grid" to
#'    generate points on a regular grid.
#' @param equalarea Logical indicating whether to generate points on an equal area basis (TRUE, the default) or in
#'    lon-lat space (FALSE, which over-samples higher-latitude regions).
#' @param epsg EPSG code for the projection in which to generate points. Only used if \code{equalarea = TRUE} and
#' \code{sample = "grid"}. The default is \code{8857}, the Equal Area projection.
#' @param expand Proportion by which to expand the bounding box around \code{x} after projection, to capture the full
#'    domain in projected space. Only used if \code{equalarea = TRUE} and \code{sample = "grid"}.
#' @return A two-column matrix of x (longitude) and y (latitude) values, with a row for each particle.
#' @export
generate_particles <- function(x, n = 1000, sample = "random", equalarea = TRUE,
                               epsg = 8857, expand = .25){

      if(!equalarea & sample == "grid"){
            p <- geo_grid(x, n)
      }

      if(equalarea & sample == "grid"){
            if(!requireNamespace("sf", quietly = TRUE))
                  stop("the sf package is required for `sample = \"grid\", equalarea = TRUE`")
            g <- as.data.frame(geo_grid(x, n))
            g <- sf::st_as_sf(g, coords = 1:2, crs = sf::st_crs(4326))
            g <- sf::st_transform(g, crs = sf::st_crs(epsg))
            h <- sf::st_convex_hull(sf::st_union(g))
            p <- as.data.frame(geo_grid(terra::ext(g) * (1 + expand), n))
            p <- sf::st_as_sf(p, coords = 1:2, crs = sf::st_crs(epsg))
            p <- p[sf::st_contains(h, p)[[1]],]
            p <- sf::st_coordinates(sf::st_transform(p, crs = sf::st_crs(4326)))
            colnames(p) <- c("x", "y")
      }


      if(equalarea & sample == "random"){
            # calculate number of samples needed
            y <- seq(ymin(x), ymax(x), length.out = 1000)
            ws <- 1 / aspect(y)
            ws <- ws / max(ws)
            ns <- round(n / mean(ws)) + 1000 # safety margin of 1000

            # weighted sample
            px <- runif(ns, xmin(x), xmax(x))
            py <- runif(ns, ymin(x), ymax(x))
            w <- 1 / aspect(py)
            s <- sample(ns, ns, prob = w / max(ws))
            p <- cbind(x = px[s][1:n], y = py[s][1:n])
      }

      if(!equalarea & sample == "random"){
            p <- cbind(runif(n, xmin(x), xmax(x)),
                       runif(n, ymin(x), ymax(x)))
      }

      p
}


# # experimental
# particle_walk <- function(x, p, n_iter = 100, scale = .1,
#                           wrap = "neither",
#                           direction = "both", #ignore_speed = FALSE,
#                           sf = FALSE){
#
#       if(direction == "both") direction <- c("downwind", "upwind")
#       xlim <- ext(x)[1:2]
#       xrng <- diff(xlim)
#       ylim <- ext(x)[3:4]
#       yrng <- diff(ylim)
#
#       p0 <- p
#       b <- data.frame()
#
#       # step <- matrix(c(0, 0,  -1, -1,  -1, 0,  -1, 1,  0, 1,
#       #                  1, 1,  1, 0,  1, -1,  0, -1),
#       #                byrow = T, ncol = 2) * res(x)
#       step <- matrix(c(-1, -1,  -1, 0,  -1, 1,  0, 1,
#                        1, 1,  1, 0,  1, -1,  0, -1),
#                      byrow = T, ncol = 2) * res(x)
#
#       for(drn in direction){
#             s <- ifelse(drn == "downwind", 1, -1)
#             p <- p0
#             t <- array(NA, dim = c(nrow(p), 2, n_iter + 1))
#             for(i in 1:(n_iter + 1)){
#
#                   d <- as.matrix(terra::extract(x, p, method = "bilinear"))
#                   # d <- t(apply(d, 1, function(x) step[sample.int(9, 1, prob = x, useHash = FALSE),]))
#                   d <- t(apply(d, 1, function(x){
#                         ss <- sample.int(8, 1, useHash = TRUE)
#                         step[ss,] * x[ss]
#                   }))
#
#                   t[, 1:2, i] <- p
#                   p <- p + d * s * scale
#
#                   if(wrap %in% c("horizontal", "both")){
#                         z <- which(p[, 1] < xlim[1])
#                         p[z, 1] <- p[z, 1] + xrng
#                         z <- which(p[,1] > xlim[2])
#                         p[z, 1] <- p[z, 1] - xrng
#                   }
#                   if(wrap %in% c("vertical", "both")){
#                         z <- which(p[, 2] < ylim[1])
#                         p[z, 2] <- p[z, 2] + yrng
#                         z <- which(p[,2] > ylim[2])
#                         p[z, 2] <- p[z, 2] - yrng
#                   }
#             }
#
#             d <- lapply(1:(n_iter + 1), function(i) data.frame(p = 1:nrow(p),
#                                                                x = t[,1,i],
#                                                                y = t[,2,i],
#                                                                t = (i - 1) * s))
#             d <- do.call("rbind", d)
#             b <- rbind(b, d)
#       }
#
#       b <- distinct(b)
#       b <- b[order(b$p, b$t),]
#       b <- na.omit(b)
#
#       if(sf){
#             bb <- b %>%
#                   st_as_sf(coords = c(2, 3, 5), dim = "XYM",
#                            crs = sf::st_crs(4326)) %>%
#                   group_by(p) %>%
#                   summarize(do_union = FALSE) %>%
#                   st_cast("LINESTRING")
#       }
#
#       b
# }
