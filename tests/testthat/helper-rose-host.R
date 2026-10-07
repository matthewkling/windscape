# A fake host of pre-built wind roses, for testing download_wind_rose() offline: small global
# roses for two years, in the same tiers, names, and catalog layout as the real release (see
# data-raw/cfsr_roses.R), served from a local directory via file:// URLs. Sets the windscape
# options pointing to the host, its catalog, and a fresh cache, for the calling test.
fake_rose_host <- function(env = parent.frame()){
      host <- withr::local_tempdir(.local_envir = env)
      release <- "test-roses"
      dir.create(file.path(host, release))
      rows <- list()
      put <- function(rose, tier, year = NA, month = NA){
            name <- paste0(paste(c("cfsr_10m", tier, if(!is.na(year)) year,
                                   if(!is.na(month)) sprintf("%02d", month)), collapse = "_"), ".tif")
            f <- file.path(host, release, name)
            terra::writeRaster(rose, f, datatype = "FLT4S", overwrite = TRUE)
            rows[[length(rows) + 1]] <<- data.frame(
                  source = "cfsr", level = "10m", tier = tier, year = year, month = month,
                  first = NA, last = NA, n_steps = rose@n_steps, trans = 1, file = name,
                  bytes = file.size(f), md5 = unname(tools::md5sum(f)), release = release)
            as_wind_rose(terra::rast(f), trans = 1, n_steps = rose@n_steps)
      }

      g <- terra::rast(nrows = 18, ncols = 36, xmin = -180, xmax = 180, ymin = -90, ymax = 90,
                       nlyrs = 8, crs = "EPSG:4326")
      ym <- expand.grid(month = 1:12, year = 2001:2002)
      ym$n_steps <- days_in_month(ym$year, ym$month) * 24
      set.seed(10)
      monthly <- lapply(seq_len(nrow(ym)), function(i){
            r <- g
            terra::values(r) <- stats::runif(terra::ncell(r) * 8)
            put(as_wind_rose(r, trans = 1, n_steps = ym$n_steps[i]), "month", ym$year[i], ym$month[i])
      })
      annual <- lapply(unique(ym$year), function(y)
            put(combine_roses(monthly[ym$year == y]), "year", y))
      for(m in 1:12) put(combine_roses(monthly[ym$month == m]), "month_of_year", NA, m)
      put(combine_roses(annual), "all")

      catalog <- file.path(host, "catalog.csv")
      utils::write.csv(do.call(rbind, rows), catalog, row.names = FALSE)
      withr::local_options(list(windscape.rose_host = paste0("file://", host),
                                windscape.rose_catalog = catalog,
                                windscape.cache_dir = file.path(host, "cache")),
                           .local_envir = env)
      list(host = host, release = release, ym = ym, monthly = monthly)
}
