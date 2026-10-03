# ggplot2 support: fortify methods, block aggregation for stats, and bearing color scales.


#' Convert windscape objects to data frames for ggplot2
#'
#' `fortify()` methods that convert windscape objects to tidy data frames, so they can be
#' passed directly to [ggplot2::ggplot()] or a layer's `data` argument. They can also be called
#' directly, to inspect or modify the data before plotting.
#'
#' @param model A windscape object: a `wind_rose`, `wind_field`, `wind_series`, or the result
#'    of [random_walk()].
#' @param data Not used.
#' @param na.rm Logical: drop grid cells with missing values? Default `TRUE`.
#' @param ... Not used.
#' @return A data frame with grid cell center coordinates `x` and `y`, plus:
#' * `wind_rose`: one column of conductance per direction, `SW`, `W`, `NW`, `N`, `NE`, `E`,
#'   `SE`, `S`, and `total`, their sum. For longitude/latitude roses, also per-cell flow
#'   summaries matching those computed by [geom_wind_rose()]: `speed`, the mean wind speed in
#'   km/h (the sum of the flows toward the eight neighbors, if `trans = 1` and wind speeds are
#'   in m/s); `net`, the speed of net flow in km/h; `bearing`, its direction in degrees
#'   clockwise from north; and `consistency`, `net / speed`, the steadiness of wind direction
#'   over time (near 1 where wind blows predominantly one way, near 0 where it has no net
#'   direction).
#' * `wind_field`: wind components `u` and `v`; `speed`, in the units of `u` and `v`; and
#'   `bearing`, the direction the wind blows toward, in degrees clockwise from north.
#' * `wind_series`: one row per grid cell and time step, with `step` (the time step's index),
#'   `time` (parsed from layer names where possible, otherwise `NA`), and `u`, `v`, `speed`,
#'   and `bearing` as for a `wind_field`.
#' * `random_walk()` result, pulse mode: one row per grid cell and recorded iteration, with
#'   `iteration`, `hours` (elapsed time), `airborne`, and `deposition`.
#' * `random_walk()` result, stream mode: `residence`, `deposition`, and for upwind walks,
#'   `origin` (the `flux` element is a `wind_field`; fortify it separately).
#' @examples
#' rose <- windscape_example("wind_rose")
#' head(ggplot2::fortify(rose))
#' @name fortify.windscape
NULL

#' @rdname fortify.windscape
#' @importFrom ggplot2 fortify
#' @export
fortify.wind_rose <- function(model, data, na.rm = TRUE, ...){
      d <- terra::as.data.frame(as(model, "SpatRaster"), xy = TRUE, na.rm = na.rm)
      cond <- as.matrix(d[, ROSE_DIRS])
      d$total <- rowSums(cond)
      if(isTRUE(terra::is.lonlat(model, perhaps = TRUE, warn = FALSE))){
            cell <- mean(terra::res(model)) # as in wind_rose()
            d <- cbind(d, flow_stats(rose_flows(cond, d$y, cell), d$y, cell))
      }
      d
}

#' @rdname fortify.windscape
#' @export
fortify.wind_field <- function(model, data, na.rm = TRUE, ...){
      d <- terra::as.data.frame(as(model, "SpatRaster"), xy = TRUE, na.rm = na.rm)
      names(d) <- c("x", "y", "u", "v")
      d$speed <- sqrt(d$u^2 + d$v^2)
      d$bearing <- (atan2(d$u, d$v) * 180 / pi) %% 360
      d
}

#' @rdname fortify.windscape
#' @export
fortify.wind_series <- function(model, data, na.rm = TRUE, ...){
      n <- model@n_steps
      d <- terra::as.data.frame(as(model, "SpatRaster"), xy = TRUE, na.rm = na.rm)
      out <- data.frame(x = rep(d$x, n), y = rep(d$y, n),
                        step = rep(seq_len(n), each = nrow(d)),
                        u = unlist(d[2 + seq_len(n)], use.names = FALSE),
                        v = unlist(d[2 + n + seq_len(n)], use.names = FALSE))
      out$time <- rep(layer_times(names(model)[seq_len(n)]), each = nrow(d))
      out$speed <- sqrt(out$u^2 + out$v^2)
      out$bearing <- (atan2(out$u, out$v) * 180 / pi) %% 360
      out[, c("x", "y", "step", "time", "u", "v", "speed", "bearing")]
}

#' @rdname fortify.windscape
#' @export
fortify.random_walk <- function(model, data, na.rm = TRUE, ...){
      first <- model[[1]]
      d <- lapply(model, function(r) terra::as.data.frame(as(r, "SpatRaster"), xy = TRUE, na.rm = FALSE))
      if(first@mode == "stream"){
            out <- data.frame(x = d[[1]]$x, y = d[[1]]$y)
            for(nm in intersect(c("residence", "deposition", "origin"), names(model))) out[[nm]] <- d[[nm]][[nm]]
      }else{
            layers <- names(first)
            iteration <- as.integer(sub("^iter", "", layers))
            nc <- nrow(d[[1]])
            out <- data.frame(x = rep(d[[1]]$x, length(layers)), y = rep(d[[1]]$y, length(layers)),
                              iteration = rep(iteration, each = nc),
                              hours = rep(iteration * first@iter_length, each = nc),
                              airborne = unlist(d$airborne[layers], use.names = FALSE),
                              deposition = unlist(d$deposition[layers], use.names = FALSE))
      }
      if(na.rm) out <- out[stats::complete.cases(out), ]
      rownames(out) <- NULL
      out
}


# Parse times from wind_series layer names like "u 2000-01-01 06:00:00" (midnight is often
# written without a time). Returns POSIXct (UTC), with NA where parsing fails.
layer_times <- function(nm){
      s <- sub("^[uv] ", "", nm)
      s <- ifelse(grepl("^\\d{4}-\\d{2}-\\d{2}$", s), paste(s, "00:00:00"), s)
      as.POSIXct(s, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
}


# Block layout for gridded data in degrees: blocks approximately square in km, with about
# `res` blocks along the longer side of the data's extent (as for the star plot). Cell size is
# inferred from coordinate spacing. Blocks are anchored at the top-left corner of the grid.
block_spec <- function(df, res){
      if(length(res) != 1 || !is.finite(res) || res < 1) stop("`res` must be a number >= 1")
      spacing <- function(z){
            u <- sort(unique(z))
            if(length(u) < 2) return(NA_real_)
            min(diff(u)[diff(u) > 1e-9 * max(abs(u), 1)])
      }
      dx <- spacing(df$x)
      dy <- spacing(df$y)
      if(is.na(dx)) dx <- dy
      if(is.na(dy)) dy <- dx
      if(is.na(dx)) dx <- dy <- 1
      xmin <- min(df$x) - dx / 2
      xmax <- max(df$x) + dx / 2
      ymin <- min(df$y) - dy / 2
      ymax <- max(df$y) + dy / 2
      mid <- (ymin + ymax) / 2
      block_km <- max((xmax - xmin) * km_lon(mid), (ymax - ymin) * KM_LAT) / res
      fx <- max(1, round(block_km / (dx * km_lon(mid))))
      fy <- max(1, round(block_km / (dy * KM_LAT)))
      list(fx = fx, fy = fy, dx = dx, dy = dy, xmin = xmin, ymax = ymax,
           spacing_km = min(fx * dx * km_lon(mid), fy * dy * KM_LAT))
}

KM_LAT <- 110.57                                         # km per degree latitude, approximately
km_lon <- function(lat) 111.32 * cos(lat * pi / 180)     # km per degree longitude, approximately


# Average gridded values over the blocks defined by `spec` (see block_spec()). Returns one row
# per block: block centroid `x`, `y` (mean coordinates of its cells), the mean of each column
# in `cols`, and the number of cells `n`. Rows with any NA in `cols` are dropped first.
block_means <- function(df, cols, res = NULL, spec = block_spec(df, res)){
      df <- df[stats::complete.cases(df[, c("x", "y", cols)]), , drop = FALSE]
      if(nrow(df) == 0) return(cbind(df[0, c("x", "y", cols)], n = integer(0)))
      bx <- floor((df$x - spec$xmin) / (spec$fx * spec$dx) + 1e-9)
      by <- floor((spec$ymax - df$y) / (spec$fy * spec$dy) + 1e-9)
      key <- paste(by, bx)
      out <- stats::aggregate(df[, c("x", "y", cols)], by = list(key = key), FUN = mean)
      out$n <- as.vector(table(key)[out$key])
      out <- out[order(as.numeric(sub(" .*", "", out$key)), as.numeric(sub(".* ", "", out$key))), ]
      out$key <- NULL
      rownames(out) <- NULL
      attr(out, "block") <- c(fx = spec$fx, fy = spec$fy, dx = spec$dx, dy = spec$dy,
                              spacing_km = spec$spacing_km)
      out
}


#' Color scales for wind direction
#'
#' Cyclic color scales for compass bearings (degrees clockwise from north), using the same hues
#' as the windscape direction color wheel. Bearings outside 0-360 are wrapped, so 0 and 360 (or
#' -90 and 270) get the same color. By default the legend shows the eight compass directions;
#' pass `breaks` and `labels` to change this.
#'
#' @param ... Other arguments passed to [ggplot2::continuous_scale()], such as `name` or
#'    `guide`.
#' @param chroma,luminance HCL chroma and luminance of the colors. Defaults match the windscape
#'    rose plots.
#' @param aesthetics The aesthetics to which the scale applies.
#' @return A ggplot2 scale.
#' @examples
#' library(ggplot2)
#' d <- data.frame(x = 1:8, bearing = seq(0, 315, 45))
#' ggplot(d, aes(x, 1, fill = bearing)) + geom_tile() + scale_fill_bearing()
#' @export
scale_fill_bearing <- function(..., chroma = 90, luminance = 65, aesthetics = "fill"){
      bearing_scale(aesthetics, chroma, luminance, ...)
}

#' @rdname scale_fill_bearing
#' @export
scale_colour_bearing <- function(..., chroma = 90, luminance = 65, aesthetics = "colour"){
      bearing_scale(aesthetics, chroma, luminance, ...)
}

#' @rdname scale_fill_bearing
#' @export
scale_color_bearing <- scale_colour_bearing

bearing_scale <- function(aesthetics, chroma, luminance, ...){
      args <- list(aesthetics = aesthetics,
                   palette = function(x) grDevices::hcl(h = x * 360, c = chroma, l = luminance),
                   limits = c(0, 360),
                   oob = function(x, range = c(0, 360)) x %% 360,
                   breaks = seq(0, 315, 45),
                   labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))
      args <- utils::modifyList(args, list(...)) # user arguments (e.g. breaks) take precedence
      # `scale_name` is required before ggplot2 3.5.0 and deprecated after
      if(utils::packageVersion("ggplot2") < "3.5.0") args$scale_name <- "bearing"
      do.call(ggplot2::continuous_scale, args)
}


# Flow (km/h) toward each of the eight neighbors: conductance times the distance to that
# neighbor, at each row's latitude. `cond` is a matrix with columns in ROSE_DIRS order, `lat`
# the latitude of each row, and `cell` the cell size in degrees.
rose_flows <- function(cond, lat, cell){
      lats <- unique(lat)
      dist <- t(vapply(lats, function(l) sqrt(rowSums(rw_neighbor_displacements(l, cell)^2)),
                       numeric(8)))
      out <- as.matrix(cond) * dist[match(lat, lats), , drop = FALSE]
      colnames(out) <- ROSE_DIRS
      out
}

# Summaries of flow matrices (rows in ROSE_DIRS column order) at latitudes `lat`: total
# `speed`, `net` flow speed, its `bearing`, and `consistency` (net / speed).
flow_stats <- function(flows, lat, cell){
      lats <- unique(lat)
      U <- lapply(lats, function(l){ d <- rw_neighbor_displacements(l, cell); d / sqrt(rowSums(d^2)) })
      ux <- t(vapply(U, function(u) u[, 1], numeric(8)))[match(lat, lats), , drop = FALSE]
      uy <- t(vapply(U, function(u) u[, 2], numeric(8)))[match(lat, lats), , drop = FALSE]
      nx <- rowSums(flows * ux)
      ny <- rowSums(flows * uy)
      speed <- rowSums(flows)
      net <- sqrt(nx^2 + ny^2)
      data.frame(speed = speed, net = net, bearing = (atan2(nx, ny) * 180 / pi) %% 360,
                 consistency = ifelse(speed > 0, net / speed, 0))
}
