#' Wind trails and streamlines
#'
#' `geom_wind_trail()` draws particle trails: paths through a wind field, with an arrowhead at
#' the downwind end. `stat_wind_trail()` computes trails from a wind field with [wind_trails()],
#' seeding them on a regular grid (or at points you supply) and tracing them upwind and
#' downwind; trails traced through a single wind field are streamlines. To draw trails from
#' [wind_trails()] output, use `geom_wind_trail()` with its default identity stat.
#'
#' @param mapping Aesthetic mappings created by [ggplot2::aes()]. `x` and `y` (coordinates, in
#'    degrees longitude and latitude) must be mapped, usually in the plot's main call.
#'    \cr For `geom_wind_trail()` (drawing existing trails), `group` and `t` are mapped
#'    automatically to the columns `trail` and `step` of [wind_trails()] output: `group`
#'    identifies each trail, and `t` orders the points along it, from upwind to downwind. Map
#'    them in this layer if your columns are named differently.
#'    \cr For `stat_wind_trail()` (computing trails), the wind components `u` and `v` are mapped
#'    automatically to columns of the same names, as produced by `fortify()` on a `wind_field`
#'    (see [fortify.windscape]).
#' @param data A data frame, or for `stat_wind_trail()`, a `wind_field`, which is converted with
#'    `fortify()`. Default is to inherit the plot's data.
#' @param stat,geom Use to override the default pairings: `geom_wind_trail()` uses the identity
#'    stat, and `stat_wind_trail()` uses `geom_wind_trail()`.
#' @param position Position adjustment; see [ggplot2::layer()].
#' @param ... Other arguments passed to [ggplot2::layer()], such as fixed aesthetics like
#'    `colour = "white"` or `linewidth = 0.8`.
#' @param arrow Arrowhead drawn at the downwind end of each trail, created by [grid::arrow()],
#'    or `NULL` for none.
#' @param seeds `stat_wind_trail()`: starting points for trails. `NULL` (the default) seeds
#'    one trail at the center of each block of grid cells (see `res`). Alternatively, a
#'    two-column matrix of longitude and latitude, e.g. to trace trails from study sites; the
#'    same seeds are used in every panel.
#' @param res `stat_wind_trail()`: approximate number of seeds along the longer side of the
#'    data's extent, when `seeds = NULL`. Seeds are placed at the centers of blocks of grid cells
#'    that are approximately square in km. Also sets the spacing used by `length`. Default 20.
#' @param fixed_length `stat_wind_trail()`: logical: give all trails the same length, so they
#'    show only direction (streamlines)? Default `FALSE`: trails span the same transport time,
#'    so their length is proportional to wind speed, as for [geom_wind_arrow()]. In fields with
#'    very uneven speeds, `TRUE` shows flow structure more evenly, since trails in calm areas
#'    don't shrink.
#' @param length `stat_wind_trail()`: length of each trail, as a multiple of the spacing between
#'    seeds set by `res`. With `fixed_length = FALSE`, this is the length of a trail at the 90th
#'    percentile of wind speed, unless `hours` is given. Default 1.5.
#' @param hours `stat_wind_trail()`, with `fixed_length = FALSE`: transport time in hours
#'    spanned by each trail, in place of `length`, so that trails show actual transport
#'    distance (assuming wind speeds in m/s). Default `NULL`.
#' @param direction `stat_wind_trail()`: trace trails `"both"` ways from each seed (the
#'    default), only `"downwind"`, or only `"upwind"`.
#' @param steps `stat_wind_trail()`: number of integration steps per trail; more gives smoother
#'    trails. Default 20.
#' @param wrap `stat_wind_trail()`: wrap trails across the field's edges; see [wind_trails()].
#'    Default `"neither"`.
#' @param na.rm,show.legend,inherit.aes See [ggplot2::layer()].
#'
#' @return A ggplot2 layer.
#'
#' @section Computed variables (`stat_wind_trail()`):
#' * `speed`: wind speed at each point along the trail, in the units of `u` and `v`.
#' * `t`: integration step, negative upwind of the seed and positive downwind.
#' * `hours` (with `fixed_length = FALSE`) or `km` (with `fixed_length = TRUE`): signed
#'   transport time or distance from the seed.
#' * `progress`: position along the trail from its upwind end (0) to its downwind end (1).
#'
#' @details
#' Trails are traced with [wind_trails()], using bilinear interpolation of the wind field, and
#' stop where they leave the field (unless wrapped), so trails seeded near the edges can be
#' shorter; extend the field beyond the area of interest to avoid this. They are computed in
#' longitude/latitude coordinates with correct directions and lengths in km; use
#' [ggplot2::coord_quickmap()] or [ggplot2::coord_sf()] with a longitude/latitude CRS. Seeds and
#' trail lengths are set once for the whole layer, so faceted panels are comparable.
#'
#' @examples
#' library(ggplot2)
#' katrina <- windscape_example("wind_field")
#'
#' # trails over wind speed, with length proportional to speed
#' ggplot(katrina, aes(x, y)) +
#'   geom_raster(aes(fill = speed)) +
#'   stat_wind_trail(colour = "white") +
#'   coord_quickmap()
#'
#' # streamlines: equal-length trails showing direction only
#' ggplot(katrina, aes(x, y)) +
#'   stat_wind_trail(fixed_length = TRUE) +
#'   coord_quickmap()
#'
#' # 6 hours of transport from two sites
#' sites <- cbind(c(-92, -86), c(24, 30))
#' ggplot(katrina, aes(x, y)) +
#'   geom_raster(aes(fill = speed)) +
#'   stat_wind_trail(seeds = sites, hours = 6, steps = 50, colour = "white") +
#'   coord_quickmap()
#'
#' # trails from wind_trails()
#' trails <- wind_trails(katrina, sites, hours = 12)
#' ggplot(trails, aes(x, y)) +
#'   geom_wind_trail(aes(colour = speed)) +
#'   coord_quickmap()
#' @name geom_wind_trail
NULL

#' @rdname geom_wind_trail
#' @format NULL
#' @usage NULL
#' @export
GeomWindTrail <- ggplot2::ggproto("GeomWindTrail", ggplot2::GeomPath,
                                  default_aes = ggplot2::aes(colour = "gray10", linewidth = 0.4, linetype = 1, alpha = NA),
                                  optional_aes = "t",

                                  setup_data = function(data, params){
                                        # order points along each trail, so paths run from upwind to downwind
                                        if("t" %in% names(data)) data <- data[order(data$PANEL, data$group, data$t), , drop = FALSE]
                                        data
                                  }
)

#' @rdname geom_wind_trail
#' @format NULL
#' @usage NULL
#' @export
StatWindTrail <- ggplot2::ggproto("StatWindTrail", ggplot2::Stat,
                                  required_aes = c("x", "y", "u", "v"),
                                  extra_params = c("na.rm", "seeds", "res", "fixed_length", "length", "hours", "direction",
                                                   "steps", "wrap"),

                                  compute_layer = function(self, data, params, layout){
                                        data <- data[stats::complete.cases(data[, c("x", "y", "u", "v")]), , drop = FALSE]
                                        if(nrow(data) == 0) return(data.frame())
                                        spec <- block_spec(data, params$res)
                                        pieces <- split(data, list(data$PANEL, data$group), drop = TRUE)

                                        # trail extent, set once for the whole layer so panels are comparable
                                        km <- params$length * spec$spacing_km
                                        if(isTRUE(params$fixed_length)){
                                              extent <- list(distance = km)
                                        }else if(!is.null(params$hours)){
                                              extent <- list(hours = params$hours)
                                        }else{
                                              ref <- stats::quantile(sqrt(data$u^2 + data$v^2), 0.9, names = FALSE)
                                              if(!is.finite(ref) || ref <= 0) ref <- 1
                                              extent <- list(hours = km / (ref * 3.6)) # time for the 90th-percentile wind to travel km
                                        }

                                        out <- lapply(seq_along(pieces), function(i){
                                              d <- pieces[[i]]
                                              seeds <- params$seeds
                                              if(is.null(seeds)) seeds <- as.matrix(block_means(d, c("u", "v"), spec = spec)[, c("x", "y")])
                                              if(nrow(seeds) == 0) return(NULL)
                                              r <- terra::rast(d[, c("x", "y", "u", "v")], type = "xyz", crs = "EPSG:4326")
                                              tr <- do.call(wind_trails, c(list(x = wind_field(r), seeds = seeds,
                                                                                steps = params$steps,
                                                                                direction = params$direction,
                                                                                wrap = params$wrap), extent))
                                              if(nrow(tr) == 0) return(NULL)
                                              rng <- tapply(tr$step, tr$trail, range)
                                              lo <- vapply(rng, `[`, numeric(1), 1)[as.character(tr$trail)]
                                              hi <- vapply(rng, `[`, numeric(1), 2)[as.character(tr$trail)]
                                              cc <- constant_cols(d, c("x", "y", "u", "v", "group"))
                                              cc <- cc[rep(1, nrow(tr)), , drop = FALSE]
                                              rownames(cc) <- NULL
                                              elapsed <- tr[, intersect(c("hours", "km"), names(tr)), drop = FALSE]
                                              cbind(data.frame(x = tr$x, y = tr$y, speed = tr$speed, t = tr$step,
                                                               progress = ifelse(hi > lo, (tr$step - lo) / (hi - lo), 0),
                                                               trail = paste(i, tr$trail)), elapsed, cc)
                                        })
                                        out <- do.call(rbind, out)
                                        if(is.null(out) || nrow(out) == 0) return(data.frame())
                                        out$group <- as.integer(factor(out$trail, levels = unique(out$trail)))
                                        out$trail <- NULL
                                        out
                                  },

                                  compute_panel = function(data, scales, ...) data
)

#' @rdname geom_wind_trail
#' @export
geom_wind_trail <- function(mapping = NULL, data = NULL, stat = "identity",
                            position = "identity", ...,
                            arrow = grid::arrow(length = grid::unit(0.1, "cm"), type = "closed"),
                            na.rm = FALSE, show.legend = NA, inherit.aes = TRUE){
      if(identical(stat, "identity")) mapping <- default_mapping(mapping, c("group", "t"), c("trail", "step"))
      ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = GeomWindTrail,
                     position = position, show.legend = show.legend, inherit.aes = inherit.aes,
                     params = c(list(arrow = arrow, na.rm = na.rm), list(...)))
}

#' @rdname geom_wind_trail
#' @export
stat_wind_trail <- function(mapping = NULL, data = NULL, geom = GeomWindTrail,
                            position = "identity", ..., seeds = NULL, res = 20,
                            fixed_length = FALSE, length = 1.5, hours = NULL,
                            direction = c("both", "downwind", "upwind"), steps = 20,
                            wrap = c("neither", "horizontal", "vertical", "both"),
                            arrow = grid::arrow(length = grid::unit(0.1, "cm"), type = "closed"),
                            na.rm = FALSE, show.legend = NA, inherit.aes = TRUE){
      direction <- match.arg(direction)
      wrap <- match.arg(wrap)
      if(!is.null(seeds)){
            seeds <- as.matrix(seeds)
            if(ncol(seeds) != 2) stop("`seeds` must be a two-column matrix of coordinates")
      }
      if(!is.null(hours) && (base::length(hours) != 1 || !is.finite(hours) || hours <= 0))
            stop("`hours` must be a positive number")
      if(length(res) != 1 || !is.finite(res) || res < 1) stop("`res` must be a number >= 1")
      if(base::length(length) != 1 || !is.finite(length) || length <= 0)
            stop("`length` must be a positive number")
      if(base::length(steps) != 1 || !is.finite(steps) || steps < 1) stop("`steps` must be a number >= 1")
      mapping <- default_mapping(mapping, c("u", "v"))
      ggplot2::layer(data = data, mapping = mapping, stat = StatWindTrail, geom = geom,
                     position = position, show.legend = show.legend, inherit.aes = inherit.aes,
                     params = c(list(seeds = seeds, res = res, fixed_length = fixed_length,
                                     length = length, hours = hours, direction = direction,
                                     steps = steps, wrap = wrap, arrow = arrow, na.rm = na.rm),
                                list(...)))
}
