# Synthetic wind roses for tests. All are built with rose(), so they share the package's
# conductance geometry; no downloaded wind data are needed.

# Rose from a function giving u and v time series for a cell at (x, y).
build_rose <- function(nr, nc, uv_fun, xmin = -100, ymin = 30, res = 1,
                       crs = "EPSG:4326", lat_fixed = NULL, n_steps = 50){
      g <- terra::rast(nrows = nr, ncols = nc, xmin = xmin, xmax = xmin + nc * res,
                       ymin = ymin, ymax = ymin + nr * res, crs = crs)
      xy <- terra::crds(g)
      vals <- t(apply(xy, 1, function(p){
            uv <- uv_fun(p[1], p[2])
            lat <- if(is.null(lat_fixed)) p[2] else lat_fixed
            rose(c(lat, res, uv$u, uv$v))
      }))
      r <- terra::rast(g, nlyrs = 8)
      terra::values(r) <- vals
      as_wind_rose(r, trans = 1, n_steps = n_steps)
}

# Spatially varying, noisy lon/lat rose: westerlies strengthening northward.
noisy_rose <- function(nr = 12, nc = 15, seed = 1){
      set.seed(seed)
      build_rose(nr, nc, function(x, y){
            list(u = 3 + 0.1 * (y - 30) + stats::rnorm(50, 0, 4),
                 v = sin(x / 5) + stats::rnorm(50, 0, 4))
      })
}

# Spatially uniform rose on a planar (square-cell) grid. Directions are spread evenly
# around a mean wind (u, v), so the rose is smooth.
uniform_rose <- function(nr = 21, nc = 21, u = 0, v = 0, spread = 2){
      th <- (1:72 - 0.5) / 72 * 2 * pi
      build_rose(nr, nc, function(x, y) list(u = u + spread * sin(th), v = v + spread * cos(th)),
                 xmin = 0, ymin = 0, crs = "local", lat_fixed = 0)
}

# Run an expression silently (random_walk prints messages and a progress bar).
quietly <- function(expr){
      out <- NULL
      utils::capture.output(out <- suppressMessages(suppressWarnings(expr)))
      out
}

# Single-layer raster with given values at given cells, zero elsewhere.
point_raster <- function(rose, cells, values = 1){
      r <- terra::rast(rose, nlyrs = 1, vals = 0)
      r[cells] <- values
      r
}

# Mass-weighted mean row and column of a single-layer raster.
centroid <- function(r){
      v <- terra::values(r)[, 1]
      rc <- terra::rowColFromCell(r, seq_along(v))
      c(row = sum(rc[, 1] * v) / sum(v), col = sum(rc[, 2] * v) / sum(v))
}

vals <- function(x) terra::values(x)[, 1]

# u and v layers (uuvv order) on a small grid
uv_raster <- function(nr = 4, nc = 5, u = 5, v = 0, n_steps = 3, xmin = -100, ymin = 30,
                      xres = 1, yres = xres, crs = "EPSG:4326"){
      r <- terra::rast(nrows = nr, ncols = nc, xmin = xmin, xmax = xmin + nc * xres,
                       ymin = ymin, ymax = ymin + nr * yres, crs = crs, nlyrs = 2 * n_steps)
      terra::values(r) <- matrix(rep(c(rep(u, n_steps), rep(v, n_steps)), each = nr * nc), nr * nc)
      r
}

# Wind rose from a uv_raster(), via the full wind_series -> wind_rose pipeline.
uv_rose <- function(...) wind_rose(wind_series(uv_raster(...)), trans = 1)

# Conductance-free check of rose() geometry: geodesic distance (m) from a cell center at
# latitude `lat` to each of its 8 neighbors, in rose() layer order (SW, W, NW, N, NE, E, SE, S).
neighbor_distances <- function(lat, res = 1){
      dx <- c(-1, -1, -1, 0, 1, 1, 1, 0) * res
      dy <- c(-1, 0, 1, 1, 1, 0, -1, -1) * res
      geosphere::distGeo(c(0, lat), cbind(dx, lat + dy))
}

# rose() input vector for wind blowing toward compass bearing(s) `toward` at speed(s) `speed`.
rose_input <- function(lat, toward, speed = 5, res = 1){
      th <- toward * pi / 180
      c(lat, res, speed * sin(th), speed * cos(th))
}
