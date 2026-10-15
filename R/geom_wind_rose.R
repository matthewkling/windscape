#' Wind rose glyphs: a map of local wind roses
#'
#' Draws a star glyph for each block of grid cells in a wind rose, showing the strength of flow
#' toward each of the eight neighbors. `stat_wind_rose()` does the computation: converting
#' conductance to flow, averaging over blocks, and building glyph polygons. `geom_wind_rose()`
#' draws them, coloring each glyph by its net flow direction and graying it out where flow
#' has little net direction.
#'
#' @param mapping Aesthetic mappings created by [ggplot2::aes()]. `x` and `y` (cell center
#'    coordinates, in degrees longitude and latitude) must be mapped, usually in the plot's
#'    main call: `ggplot(rose, aes(x, y))`. The eight conductance columns (`SW`, `W`, `NW`, `N`,
#'    `NE`, `E`, `SE`, `S`) are mapped automatically to columns of the same names, as produced
#'    by `fortify()` on a `wind_rose` (see [fortify.windscape]). If your columns are named
#'    differently, map them in this layer's `mapping` rather than in [ggplot2::ggplot()]: the
#'    automatic mappings are part of the layer, so they take precedence over plot-level ones.
#' @param data A data frame, or a `wind_rose`, which is converted with `fortify()`. Default
#'    is to inherit the plot's data.
#' @param stat,geom Use to override the default pairing of `stat_wind_rose()` and
#'    `geom_wind_rose()`.
#' @param position Position adjustment; see [ggplot2::layer()].
#' @param ... Other arguments passed to [ggplot2::layer()], such as fixed aesthetics like
#'    `color = "black"`.
#' @param res Approximate number of glyphs along the longer side of the data's extent. Blocks
#'    of grid cells are sized to be approximately square in km. The same blocks are used in
#'    every panel, so faceted glyphs are directly comparable. Default 15.
#' @param scale Multiplier for the size of all glyphs. At the default of 1, the longest ray of a
#'    typical large glyph (the 90th percentile across all glyphs in the layer, including all
#'    panels) is half the spacing between glyphs.
#' @param saturation Consistency (see Computed variables) at and above which glyph fill is shown
#'    at full strength; glyphs with lower consistency are blended toward gray. Set to 0 to
#'    disable blending. Default 0.5.
#' @param normalize Logical: scale each ray by the range of wind directions its neighbor
#'    receives? If `FALSE` (the default), rays show flow toward each neighbor, so the glyph
#'    shows how the grid's eight connections divide the wind. Away from the equator, neighbor
#'    bearings are unevenly spaced (the north and south neighbors lie closer to the diagonals),
#'    so even an isotropic wind draws as an east-west peanut. If `TRUE`, each ray is divided by
#'    the width of the range of directions allocated to that neighbor, in units of 45 degrees,
#'    giving flow per 45 degrees of wind direction: an isotropic wind draws as a circle at any
#'    latitude, which suits maps of wind regimes. Either way, directions are resolved only to
#'    the eight neighbor bearings, and the computed variables (`speed`, `net`, `bearing`,
#'    `consistency`) describe the unnormalized flows.
#' @param smooth Logical: join the ray tips with a smooth curve (the default) rather than
#'    straight lines? See Details.
#' @param center Logical: draw a point at each glyph's center? Default `TRUE`.
#' @param bearing_scale Logical: add [scale_fill_bearing()] to the plot? Default `TRUE` if
#'    `fill` is neither mapped nor set, otherwise `FALSE`.
#' @param na.rm,show.legend,inherit.aes See [ggplot2::layer()].
#'
#' @return A ggplot2 layer (with `bearing_scale = TRUE`, a list of a layer and a scale).
#'
#' @section Computed variables:
#' * `bearing`: direction of net flow, in degrees clockwise from north.
#' * `consistency`: net flow divided by total flow (mean resultant length); near 1 where wind
#'   blows predominantly one way, near 0 where it has no net direction.
#' * `speed`: total flow, the mean wind speed in km/h (if `trans = 1` in [wind_rose()] and
#'   wind speeds are in m/s).
#' * `net`: net flow speed, in km/h.
#' * `n`: number of grid cells in the glyph's block.
#' * `x0`, `y0`: glyph center, the mean coordinates of the cells in its block.
#'
#' By default, `fill` is mapped to `after_stat(bearing)` and the `consistency` aesthetic to
#' `after_stat(consistency)`.
#'
#' @details
#' Each glyph's eight rays point toward the eight neighbors, with length proportional to flow
#' in that direction: conductance times the distance to the neighbor, averaged over the block,
#' in km/h. Unlike conductance, which is a rate per grid cell, flow does not depend on cell
#' size, so it can be averaged across cells and compared across latitudes. The flows sum to
#' the mean wind speed, and their vector sum is the net flow (with `normalize = TRUE`, these
#' hold only approximately for the drawn rays).
#'
#' The outline passes through the tips of the eight rays. With `smooth = TRUE`, they are joined
#' by a smooth curve: radius is
#' interpolated as a function of bearing with a periodic, shape-preserving cubic (PCHIP). Between
#' two adjacent rays the radius stays within the range of their lengths, so the outline adds no
#' lobes that aren't in the rose. The curve between rays is a visual connector: the rose holds
#' flow toward the eight neighbors only. With `smooth = FALSE`, the tips are joined by straight
#' lines, giving a star-shaped polygon.
#'
#' Glyphs are built in longitude/latitude coordinates, scaled so that their shapes are correct
#' in km at each glyph's latitude. Use [ggplot2::coord_quickmap()] or [ggplot2::coord_sf()]
#' with a longitude/latitude CRS so that the map's aspect ratio matches; in projected
#' coordinates, glyph angles are distorted.
#'
#' @examples
#' library(ggplot2)
#' rose <- windscape_example("wind_rose")
#' ggplot(rose, aes(x, y)) +
#'   geom_wind_rose() +
#'   coord_quickmap()
#'
#' # a background of mean wind speed: the glyphs use the fill scale, so show the background
#' # with alpha (or use a second fill scale, e.g. from the ggnewscale package)
#' ggplot(rose, aes(x, y)) +
#'   geom_raster(aes(alpha = speed), fill = "gray20") +
#'   geom_wind_rose() +
#'   coord_quickmap()
#'
#' # faceting: glyphs use the same blocks and size scale in every panel
#' d <- rbind(cbind(fortify(rose), version = "original"),
#'            cbind(fortify(rose), version = "copy"))
#' ggplot(d, aes(x, y)) +
#'   geom_wind_rose(res = 10) +
#'   facet_wrap(~version) +
#'   coord_quickmap()
#' @name geom_wind_rose
NULL

ROSE_DIRS <- c("SW", "W", "NW", "N", "NE", "E", "SE", "S")

# default mapping of aesthetics to columns (by default, same-named), unless overridden in
# `mapping`
default_mapping <- function(mapping, aesthetics, columns = aesthetics){
      default <- do.call(ggplot2::aes, lapply(stats::setNames(columns, aesthetics), as.name))
      if(is.null(mapping)) return(default)
      out <- c(default[setdiff(names(default), names(mapping))], mapping)
      class(out) <- class(mapping)
      out
}

# Blend colors toward gray in proportion to (1 - consistency / saturation), capped at 0. NA
# consistency, or saturation <= 0, leaves colors unchanged.
blend_gray <- function(cols, consistency, saturation = 0.5){
      if(saturation <= 0 || all(is.na(consistency))) return(cols)
      w <- pmin(1, consistency / saturation)
      w[is.na(w)] <- 1
      ok <- !is.na(cols)
      if(!any(ok)) return(cols)
      rgb <- grDevices::col2rgb(cols[ok]) * rep(w[ok], each = 3) +
            grDevices::col2rgb("#A0A0A0")[, rep(1, sum(ok)), drop = FALSE] * rep(1 - w[ok], each = 3)
      cols[ok] <- grDevices::rgb(t(rgb), maxColorValue = 255)
      cols
}

# columns that are constant within a group, to carry through aggregation
constant_cols <- function(df, exclude){
      keep <- vapply(df, function(z) length(unique(z)) == 1, logical(1))
      df[1, setdiff(names(df)[keep], exclude), drop = FALSE]
}

#' @rdname geom_wind_rose
#' @format NULL
#' @usage NULL
#' @export
StatWindRose <- ggplot2::ggproto("StatWindRose", ggplot2::Stat,
                                 required_aes = c("x", "y", ROSE_DIRS),
                                 extra_params = c("na.rm", "res", "scale", "normalize", "smooth"),
                                 default_aes = ggplot2::aes(fill = ggplot2::after_stat(bearing),
                                                            consistency = ggplot2::after_stat(consistency)),

                                 compute_layer = function(self, data, params, layout){
                                       data <- data[stats::complete.cases(data[, c("x", "y", ROSE_DIRS)]), , drop = FALSE]
                                       if(nrow(data) == 0) return(data.frame())
                                       spec <- block_spec(data, params$res)
                                       cell <- mean(c(spec$dx, spec$dy)) # as in wind_rose()

                                       # flow (km/h) toward each neighbor: conductance x distance to that neighbor
                                       data[, ROSE_DIRS] <- rose_flows(data[, ROSE_DIRS], data$y, cell)

                                       # ray lengths: flows, or flows per 45 degrees of direction,
                                       # normalized cell by cell before block averaging
                                       normalize <- isTRUE(params$normalize)
                                       RAYS <- paste0(".ray_", ROSE_DIRS)
                                       data[, RAYS] <- if(normalize){
                                             data[, ROSE_DIRS] * 45 / rose_hat_widths(data$y, cell)
                                       }else data[, ROSE_DIRS]

                                       # block means within each panel and group
                                       pieces <- split(data, list(data$PANEL, data$group), drop = TRUE)
                                       glyphs <- lapply(pieces, function(d){
                                             b <- block_means(d, c(ROSE_DIRS, RAYS), spec = spec)
                                             if(nrow(b) == 0) return(NULL)
                                             cc <- constant_cols(d, c("x", "y", ROSE_DIRS, RAYS))
                                             cc <- cc[rep(1, nrow(b)), , drop = FALSE]
                                             rownames(cc) <- NULL
                                             cbind(b, cc)
                                       })
                                       g <- do.call(rbind, glyphs)
                                       if(is.null(g) || nrow(g) == 0) return(data.frame())
                                       fl <- as.matrix(g[, ROSE_DIRS])
                                       rays <- as.matrix(g[, RAYS])
                                       g <- g[, setdiff(names(g), RAYS), drop = FALSE]

                                       g <- cbind(g, flow_stats(fl, g$y, cell))

                                       # one size reference across the whole layer, so panels are comparable
                                       ref <- stats::quantile(apply(rays, 1, max), 0.9, names = FALSE)
                                       if(!is.finite(ref) || ref <= 0) ref <- 1
                                       k <- params$scale * 0.5 * spec$spacing_km / ref # km of glyph per km/h of flow

                                       order <- match(c("N", "NE", "E", "SE", "S", "SW", "W", "NW"), ROSE_DIRS)
                                       other <- setdiff(names(g), c("x", "y", ROSE_DIRS, "group"))
                                       n <- if(isFALSE(params$smooth)) 1 else 10 # outline points per gap
                                       out <- do.call(rbind, lapply(seq_len(nrow(g)), function(i){
                                             d <- rw_neighbor_displacements(g$y[i], cell)[order, ]
                                             b <- (atan2(d[, 1], d[, 2]) * 180 / pi) %% 360 # N first, clockwise
                                             o <- rose_outline(b, rays[i, order] * k, n)     # km
                                             v <- data.frame(x = g$x[i] + o$x / km_lon(g$y[i]),
                                                             y = g$y[i] + o$y / KM_LAT,
                                                             x0 = g$x[i], y0 = g$y[i], group = i)
                                             o <- g[rep(i, nrow(v)), other, drop = FALSE]
                                             rownames(o) <- NULL
                                             cbind(v, o)
                                       }))
                                       rownames(out) <- NULL
                                       out
                                 },

                                 compute_panel = function(data, scales, ...) data
)

#' @rdname geom_wind_rose
#' @format NULL
#' @usage NULL
#' @export
GeomWindRose <- ggplot2::ggproto("GeomWindRose", ggplot2::GeomPolygon,
                                 default_aes = ggplot2::aes(colour = "gray15", fill = "gray60", linewidth = 0.25,
                                                            linetype = 1, alpha = NA, subgroup = NULL, consistency = NA),

                                 draw_panel = function(self, data, panel_params, coord, saturation = 0.5, center = TRUE,
                                                       lineend = "butt", linejoin = "round", linemitre = 10){
                                       data$fill <- blend_gray(data$fill, data$consistency, saturation)
                                       poly <- ggplot2::GeomPolygon$draw_panel(data, panel_params, coord, lineend = lineend,
                                                                               linejoin = linejoin, linemitre = linemitre)
                                       if(!center) return(poly)
                                       ctr <- unique(data[, c("x0", "y0", "group")])
                                       ctr <- coord$transform(data.frame(x = ctr$x0, y = ctr$y0), panel_params)
                                       pts <- grid::pointsGrob(ctr$x, ctr$y, pch = 19, size = grid::unit(0.6, "mm"),
                                                               gp = grid::gpar(col = "black"))
                                       grid::grobTree(poly, pts)
                                 }
)

#' @rdname geom_wind_rose
#' @export
geom_wind_rose <- function(mapping = NULL, data = NULL, stat = StatWindRose,
                           position = "identity", ..., res = 15, scale = 1, normalize = FALSE,
                           smooth = TRUE, saturation = 0.5, center = TRUE, bearing_scale = NULL,
                           na.rm = FALSE, show.legend = NA, inherit.aes = TRUE){
      wind_rose_layer(mapping, data, stat, GeomWindRose, position, res, scale, normalize,
                      smooth, saturation, center, bearing_scale, na.rm, show.legend, inherit.aes,
                      ...)
}

#' @rdname geom_wind_rose
#' @export
stat_wind_rose <- function(mapping = NULL, data = NULL, geom = GeomWindRose,
                           position = "identity", ..., res = 15, scale = 1, normalize = FALSE,
                           smooth = TRUE, saturation = 0.5, center = TRUE, bearing_scale = NULL,
                           na.rm = FALSE, show.legend = NA, inherit.aes = TRUE){
      wind_rose_layer(mapping, data, StatWindRose, geom, position, res, scale, normalize,
                      smooth, saturation, center, bearing_scale, na.rm, show.legend, inherit.aes,
                      ...)
}

wind_rose_layer <- function(mapping, data, stat, geom, position, res, scale, normalize,
                            smooth, saturation, center, bearing_scale, na.rm, show.legend, inherit.aes,
                            ...){
      if(length(scale) != 1 || !is.finite(scale) || scale <= 0) stop("`scale` must be a positive number")
      if(length(res) != 1 || !is.finite(res) || res < 1) stop("`res` must be a number >= 1")
      if(!isTRUE(normalize) && !isFALSE(normalize)) stop("`normalize` must be TRUE or FALSE")
      if(!isTRUE(smooth) && !isFALSE(smooth)) stop("`smooth` must be TRUE or FALSE")
      mapping <- default_mapping(mapping, ROSE_DIRS)
      dots <- list(...)
      if(is.null(bearing_scale)) bearing_scale <- !"fill" %in% c(names(mapping), names(dots))
      l <- ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = geom,
                          position = position, show.legend = show.legend,
                          inherit.aes = inherit.aes,
                          params = c(list(res = res, scale = scale, normalize = normalize,
                                          smooth = smooth, saturation = saturation,
                                          center = center, na.rm = na.rm), dots))
      if(bearing_scale) list(l, scale_fill_bearing(name = "net flow\ndirection")) else l
}


# Width (degrees) of the range of wind directions that wind_rose() allocates to each neighbor
# (see rose_accumulate() in src/rose.cpp): the integral of the neighbor's tent-shaped share,
# which is half the sum of the gaps to the two adjacent neighbor bearings. Returns a matrix
# with one row per value of `lat` and columns in ROSE_DIRS order. ROSE_DIRS runs clockwise, so
# adjacent columns are adjacent neighbors. Widths sum to 360; all are 45 on square-km cells.
rose_hat_widths <- function(lat, cell){
      lats <- unique(lat)
      w <- t(vapply(lats, function(l){
            d <- rw_neighbor_displacements(l, cell)
            b <- atan2(d[, 1], d[, 2]) * 180 / pi
            gap <- (b[c(2:8, 1)] - b) %% 360 # to the next neighbor clockwise
            (gap + gap[c(8, 1:7)]) / 2
      }, numeric(8)))
      colnames(w) <- ROSE_DIRS
      w[match(lat, lats), , drop = FALSE]
}

# Smooth glyph outline through ray tips: radii `r` at increasing compass bearings `b`
# (degrees, one full turn), interpolated in bearing with a periodic, shape-preserving cubic
# Hermite spline (PCHIP; Fritsch-Butland slopes, as in scipy's PchipInterpolator). Slopes are
# zero at local extrema and between same-signed secants are a weighted harmonic mean, so
# between adjacent rays the radius stays within the range of their lengths. Returns x (east)
# and y (north) coordinates of `n` points per gap, starting at each ray tip, so vertices
# 1, n + 1, 2n + 1, ... are the tips themselves; n = 1 gives just the tips (a star polygon).
rose_outline <- function(b, r, n = 10){
      k <- length(b)
      h <- diff(c(b, b[1] + 360))               # gap after each ray
      s <- (r[c(2:k, 1)] - r) / h               # secant slope over each gap
      hp <- h[c(k, 1:(k - 1))]                  # gap before each ray
      sp <- s[c(k, 1:(k - 1))]
      m <- numeric(k)
      ok <- sp * s > 0
      w1 <- 2 * h + hp
      w2 <- h + 2 * hp
      m[ok] <- (w1[ok] + w2[ok]) / (w1[ok] / sp[ok] + w2[ok] / s[ok])

      t <- (seq_len(n) - 1) / n                 # position within each gap
      i <- rep(seq_len(k), each = n)
      j <- c(2:k, 1)[i]
      tt <- rep(t, k)
      hh <- h[i]
      rad <- (2 * tt^3 - 3 * tt^2 + 1) * r[i] + (tt^3 - 2 * tt^2 + tt) * hh * m[i] +
            (-2 * tt^3 + 3 * tt^2) * r[j] + (tt^3 - tt^2) * hh * m[j]
      rad <- pmax(rad, 0)
      th <- (b[i] + tt * hh) * pi / 180
      data.frame(x = rad * sin(th), y = rad * cos(th))
}
