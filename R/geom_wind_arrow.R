#' Wind field arrows
#'
#' Draws an arrow for each block of grid cells in a wind field, showing the direction and speed
#' of the block's mean wind vector. `stat_wind_arrow()` averages wind components over blocks
#' and computes arrow positions; `geom_wind_arrow()` draws them.
#'
#' @param mapping Aesthetic mappings created by [ggplot2::aes()]. `x` and `y` (cell center
#'    coordinates, in degrees longitude and latitude) must be mapped, usually in the plot's
#'    main call: `ggplot(field, aes(x, y))`. The wind components `u` (eastward) and `v`
#'    (northward) are mapped automatically to columns of the same names, as produced by
#'    `fortify()` on a `wind_field` (see [fortify.windscape]). If your columns are named
#'    differently, map them in this layer's `mapping` rather than in [ggplot2::ggplot()].
#' @param data A data frame, or a `wind_field`, which is converted with `fortify()`. Default is
#'    to inherit the plot's data.
#' @param stat,geom Use to override the default pairing of `stat_wind_arrow()` and
#'    `geom_wind_arrow()`.
#' @param position Position adjustment; see [ggplot2::layer()].
#' @param ... Other arguments passed to [ggplot2::layer()], such as fixed aesthetics like
#'    `color = "white"` or `linewidth = 1`.
#' @param res Approximate number of arrows along the longer side of the data's extent. Blocks of
#'    grid cells are sized to be approximately square in km, and the same blocks are used in
#'    every panel. Default 20.
#' @param scale Multiplier for the length of all arrows. At the default of 1, an arrow at the
#'    90th percentile of speed across the layer (including all panels) is 0.9 times the spacing
#'    between arrows, or with `fixed_length = TRUE`, every arrow is.
#' @param fixed_length Logical: give all arrows the same length, so that length shows only
#'    direction? Default `FALSE` (length proportional to speed). Use with a speed mapping such as
#'    `aes(color = after_stat(speed))` or `aes(linewidth = after_stat(speed))` to encode speed
#'    another way. Blocks with zero mean wind get no arrow.
#' @param pivot Position of the arrow relative to its block center, as the fraction of its
#'    length that lies behind (upwind of) the center: 0.5 (the default) centers the arrow on the
#'    block, 0 starts it at the center.
#' @param arrow Arrowhead specification created by [grid::arrow()], or `NULL` for none.
#' @param saturation If a `consistency` aesthetic is mapped (e.g.
#'    `aes(consistency = after_stat(consistency))`), arrow colors are blended toward gray where
#'    consistency is below this value. Default 0.5; set to 0 to disable blending.
#' @param center Logical: draw a point at each block center? Default `FALSE`. With `pivot = 0`
#'    and `arrow = NULL`, this gives a point-and-spoke plot.
#' @param na.rm,show.legend,inherit.aes See [ggplot2::layer()].
#'
#' @return A ggplot2 layer.
#'
#' @section Computed variables:
#' * `bearing`: direction of the block's mean wind vector, in degrees clockwise from north.
#' * `speed`: speed of the mean wind vector, in the units of `u` and `v`. This sets arrow
#'   length.
#' * `mean_speed`: mean of the wind speeds of the block's cells.
#' * `consistency`: `speed / mean_speed`, the spatial coherence of wind direction within the
#'   block; 1 where all cells' winds point the same way, lower where they diverge, converge, or
#'   rotate (as in a cyclone's eye), meaning the arrow understates the block's winds.
#' * `n`: number of grid cells in the block.
#' * `x0`, `y0`: block center, the mean coordinates of its cells.
#'
#' @details
#' Wind components are averaged as vectors, so the arrow shows the block's net wind. Arrows are
#' built in longitude/latitude coordinates, scaled so that their directions and relative lengths
#' are correct in km at each block's latitude. Use [ggplot2::coord_quickmap()] or
#' [ggplot2::coord_sf()] with a longitude/latitude CRS.
#'
#' @examples
#' library(ggplot2)
#' katrina <- windscape_example("wind_field")
#' # map background variables in the raster layer, not in ggplot(), so that the arrow layer
#' # doesn't inherit them
#' ggplot(katrina, aes(x, y)) +
#'   geom_raster(aes(fill = speed)) +
#'   scale_fill_gradient(low = "gray95", high = "gray40", name = "speed") +
#'   geom_wind_arrow() +
#'   coord_quickmap()
#'
#' # fixed-length arrows, with speed shown by color
#' ggplot(katrina, aes(x, y)) +
#'   geom_wind_arrow(aes(color = after_stat(speed)), fixed_length = TRUE) +
#'   scale_color_viridis_c() +
#'   coord_quickmap()
#'
#' # arrows colored by direction
#' ggplot(katrina, aes(x, y)) +
#'   geom_wind_arrow(aes(color = after_stat(bearing)), linewidth = 0.8) +
#'   scale_color_bearing() +
#'   coord_quickmap()
#' @name geom_wind_arrow
NULL

#' @rdname geom_wind_arrow
#' @format NULL
#' @usage NULL
#' @export
StatWindArrow <- ggplot2::ggproto("StatWindArrow", ggplot2::Stat,
                                  required_aes = c("x", "y", "u", "v"),
                                  extra_params = c("na.rm", "res", "scale", "pivot", "fixed_length"),

                                  compute_layer = function(self, data, params, layout){
                                        data <- data[stats::complete.cases(data[, c("x", "y", "u", "v")]), , drop = FALSE]
                                        if(nrow(data) == 0) return(data.frame())
                                        spec <- block_spec(data, params$res)
                                        data$cell_speed <- sqrt(data$u^2 + data$v^2)

                                        pieces <- split(data, list(data$PANEL, data$group), drop = TRUE)
                                        blocks <- lapply(pieces, function(d){
                                              b <- block_means(d, c("u", "v", "cell_speed"), spec = spec)
                                              if(nrow(b) == 0) return(NULL)
                                              cc <- constant_cols(d, c("x", "y", "u", "v", "cell_speed"))
                                              cc <- cc[rep(1, nrow(b)), , drop = FALSE]
                                              rownames(cc) <- NULL
                                              cbind(b, cc)
                                        })
                                        b <- do.call(rbind, blocks)
                                        if(is.null(b) || nrow(b) == 0) return(data.frame())

                                        b$speed <- sqrt(b$u^2 + b$v^2)
                                        b$mean_speed <- b$cell_speed
                                        b$cell_speed <- NULL
                                        b$bearing <- (atan2(b$u, b$v) * 180 / pi) %% 360
                                        b$consistency <- ifelse(b$mean_speed > 0, b$speed / b$mean_speed, 0)

                                        # arrow length in km: proportional to speed, with one reference across the whole
                                        # layer so panels are comparable, or the same for every arrow
                                        if(isTRUE(params$fixed_length)){
                                              len <- ifelse(b$speed > 0, params$scale * 0.9 * spec$spacing_km, 0)
                                              ux <- ifelse(b$speed > 0, b$u / b$speed, 0)
                                              uy <- ifelse(b$speed > 0, b$v / b$speed, 0)
                                        }else{
                                              ref <- stats::quantile(b$speed, 0.9, names = FALSE)
                                              if(!is.finite(ref) || ref <= 0) ref <- 1
                                              len <- b$speed * params$scale * 0.9 * spec$spacing_km / ref
                                              ux <- ifelse(b$speed > 0, b$u / b$speed, 0)
                                              uy <- ifelse(b$speed > 0, b$v / b$speed, 0)
                                        }
                                        dx <- ux * len / km_lon(b$y) # arrow vector, in degrees
                                        dy <- uy * len / KM_LAT
                                        b$x0 <- b$x
                                        b$y0 <- b$y
                                        b$x <- b$x0 - params$pivot * dx
                                        b$y <- b$y0 - params$pivot * dy
                                        b$xend <- b$x + dx
                                        b$yend <- b$y + dy
                                        b$u <- b$v <- NULL
                                        b$group <- seq_len(nrow(b))
                                        b
                                  },

                                  compute_panel = function(data, scales, ...) data
)

#' @rdname geom_wind_arrow
#' @format NULL
#' @usage NULL
#' @export
GeomWindArrow <- ggplot2::ggproto("GeomWindArrow", ggplot2::GeomSegment,
                                  default_aes = ggplot2::aes(colour = "gray10", linewidth = 0.4, linetype = 1, alpha = NA,
                                                             consistency = NA),

                                  draw_panel = function(self, data, panel_params, coord, arrow = NULL, saturation = 0.5,
                                                        center = FALSE, lineend = "butt", linejoin = "mitre"){
                                        data$colour <- blend_gray(data$colour, data$consistency, saturation)
                                        seg <- ggplot2::GeomSegment$draw_panel(data, panel_params, coord, arrow = arrow,
                                                                               arrow.fill = data$colour, lineend = lineend,
                                                                               linejoin = linejoin)
                                        if(!center) return(seg)
                                        ctr <- coord$transform(data.frame(x = data$x0, y = data$y0), panel_params)
                                        pts <- grid::pointsGrob(ctr$x, ctr$y, pch = 19, size = grid::unit(0.8, "mm"),
                                                                gp = grid::gpar(col = data$colour))
                                        grid::grobTree(seg, pts)
                                  }
)

#' @rdname geom_wind_arrow
#' @export
geom_wind_arrow <- function(mapping = NULL, data = NULL, stat = StatWindArrow,
                            position = "identity", ..., res = 20, scale = 1, fixed_length = FALSE, pivot = 0.5,
                            arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed"),
                            saturation = 0.5, center = FALSE, na.rm = FALSE, show.legend = NA,
                            inherit.aes = TRUE){
      wind_arrow_layer(mapping, data, stat, GeomWindArrow, position, res, scale, fixed_length, pivot, arrow,
                       saturation, center, na.rm, show.legend, inherit.aes, ...)
}

#' @rdname geom_wind_arrow
#' @export
stat_wind_arrow <- function(mapping = NULL, data = NULL, geom = GeomWindArrow,
                            position = "identity", ..., res = 20, scale = 1, fixed_length = FALSE, pivot = 0.5,
                            arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed"),
                            saturation = 0.5, center = FALSE, na.rm = FALSE, show.legend = NA,
                            inherit.aes = TRUE){
      wind_arrow_layer(mapping, data, StatWindArrow, geom, position, res, scale, fixed_length, pivot, arrow,
                       saturation, center, na.rm, show.legend, inherit.aes, ...)
}

wind_arrow_layer <- function(mapping, data, stat, geom, position, res, scale, fixed_length, pivot, arrow,
                             saturation, center, na.rm, show.legend, inherit.aes, ...){
      if(length(scale) != 1 || !is.finite(scale) || scale <= 0) stop("`scale` must be a positive number")
      if(length(res) != 1 || !is.finite(res) || res < 1) stop("`res` must be a number >= 1")
      if(length(pivot) != 1 || !is.finite(pivot)) stop("`pivot` must be a number")
      mapping <- default_mapping(mapping, c("u", "v"))
      ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = geom,
                     position = position, show.legend = show.legend, inherit.aes = inherit.aes,
                     params = c(list(res = res, scale = scale, fixed_length = fixed_length,
                                     pivot = pivot, arrow = arrow,
                                     saturation = saturation, center = center, na.rm = na.rm),
                                list(...)))
}
