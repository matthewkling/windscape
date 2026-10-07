# Pre-built global wind roses, hosted as files on a GitHub release (built by
# data-raw/cfsr_roses.R) and listed in a catalog shipped with the package.

# Where the files live: <host>/<release>/<file>. Options, so tests and mirrors can redirect them.
rose_host <- function() getOption("windscape.rose_host",
                                  "https://github.com/matthewkling/windscape/releases/download")
rose_cache_dir <- function() getOption("windscape.cache_dir",
                                        tools::R_user_dir("windscape", "cache"))


#' Catalog of pre-built wind roses
#'
#' Lists the pre-built global wind roses available for download with [download_wind_rose()].
#' Currently these are built from 10 m CFSR wind for 1979-2010, at CFSR's native resolution
#' (about 0.31 degrees), with `trans = 1`.
#'
#' @return A data frame with one row per file, and columns:
#'   * `source`, `level`: wind data set and height (as in [download_wind_data()]).
#'   * `tier`: the period the file covers: `"month"` (a single year-month), `"year"` (a calendar
#'     year), `"month_of_year"` (a calendar month across all years), or `"all"` (the full period).
#'   * `year`, `month`: the year and month the file covers, or `NA` where it spans several.
#'   * `first`, `last`: the first and last year-months covered.
#'   * `n_steps`: the number of hourly time steps summarized.
#'   * `trans`: the wind speed transformation the rose was built with (see [wind_rose()]).
#'   * `file`, `bytes`, `md5`: file name, size in bytes, and MD5 checksum.
#'   * `release`, `url`: the release the file belongs to, and its download URL.
#' @seealso [download_wind_rose()]
#' @examples
#' catalog <- wind_rose_catalog()
#' table(catalog$tier)
#' @export
wind_rose_catalog <- function(){
      f <- getOption("windscape.rose_catalog",
                     system.file("extdata", "wind_rose_catalog.csv", package = "windscape"))
      x <- utils::read.csv(f, stringsAsFactors = FALSE)
      x$url <- paste(rose_host(), x$release, x$file, sep = "/")
      x
}


#' Download a pre-built wind rose
#'
#' Downloads a pre-built global wind rose for any combination of years and calendar months,
#' optionally cropped to a region, so that common analyses need no wind data downloads or
#' rose building. Roses are built from hourly 10 m CFSR wind for 1979-2010 with `trans = 1`; see
#' [wind_rose_catalog()] for what is available. For other data sets, heights, periods, or `trans`
#' functions, download wind data with [download_wind_data()] and build a rose with [wind_rose()].
#'
#' Pre-built roses are stored for single months, calendar years, calendar months across all
#' years, and the full period. A request is met with the fewest of these, combined with
#' [combine_roses()], which weights them by their numbers of hours; the result equals the rose
#' built from all of the requested hours at once.
#'
#' @param year Integer vector of years to include. `NULL` (the default) includes all years.
#' @param month Integer vector of calendar months (1-12) to include in each year. `NULL` (the
#'   default) includes all months. For example, `year = 1990:1999, month = 6:8` gives a rose for
#'   the summers of the 1990s.
#' @param ext Optional region to crop to: a `SpatExtent`, an object `terra::ext()` accepts (e.g. a
#'   `SpatRaster` or `SpatVector`), or a numeric vector `c(xmin, xmax, ymin, ymax)`, in degrees,
#'   with longitudes from -180 to 180. The result includes every cell overlapping `ext`. Regions
#'   crossing the antimeridian are not supported; for those, use the global rose.
#' @param source,level Wind data set and height. Currently only `"cfsr"` and `"10m"`.
#' @param cache Logical. If `TRUE`, whole files are downloaded and kept in a cache on disk (see
#'   [wind_rose_cache()]), so later requests for them need no download; each file is about 15
#'   MB. If `FALSE`, nothing is kept: when `ext` is given, only the cells within it are read
#'   over the internet, which is much faster than downloading whole files for small regions;
#'   otherwise files are downloaded to a temporary directory. Defaults to `TRUE` when `ext` is
#'   `NULL` and `FALSE` otherwise.
#' @param quiet Logical. If `TRUE`, progress messages are suppressed.
#' @return A `wind_rose`, with `n_steps` and `trans` recorded so it can be combined with other
#'   roses. Values are wind conductance in 1 / hours.
#' @seealso [wind_rose_catalog()] for the available files; [wind_rose_cache()] to manage the cache.
#' @examples
#' \donttest{
#' # long-term rose for the western US, read remotely
#' rose <- download_wind_rose(ext = c(-125, -100, 30, 50))
#'
#' # summers of the 1990s
#' summer <- download_wind_rose(year = 1990:1999, month = 6:8, ext = c(-125, -100, 30, 50))
#' }
#' @export
download_wind_rose <- function(year = NULL, month = NULL, ext = NULL,
                               source = "cfsr", level = "10m",
                               cache = is.null(ext), quiet = FALSE){
      catalog <- wind_rose_catalog()
      cat_sl <- catalog[catalog$source == source & catalog$level == level, , drop = FALSE]
      if(nrow(cat_sl) == 0){
            have <- unique(paste0('source = "', catalog$source, '", level = "', catalog$level, '"'))
            stop("no pre-built roses for source = \"", source, "\", level = \"", level,
                 "\"; available: ", paste(have, collapse = "; "), call. = FALSE)
      }
      files <- rose_files(cat_sl, year, month)
      if(!is.null(ext)) ext <- rose_ext(ext)
      remote <- !is.null(ext) && !isTRUE(cache)
      dir <- if(isTRUE(cache)) rose_cache_dir() else file.path(tempdir(), "windscape_roses")

      if(!quiet) message(if(remote) "reading " else "getting ", nrow(files), " wind rose file",
                         if(nrow(files) > 1) "s", if(remote) " remotely")
      roses <- vector("list", nrow(files))
      for(i in seq_len(nrow(files))){
            f <- files[i, ]
            r <- NULL
            if(remote){
                  r <- tryCatch(suppressWarnings(read_rose_remote(f$url, ext)), error = function(e){
                        if(!quiet) message("remote read of ", f$file, " failed (",
                                           conditionMessage(e), "); downloading the whole file")
                        NULL
                  })
            }
            if(is.null(r)){
                  r <- terra::rast(fetch_rose(f, dir, quiet))
                  if(!is.null(ext)) r <- terra::crop(r, ext, snap = "out")
            }
            roses[[i]] <- as_wind_rose(r, trans = f$trans, n_steps = f$n_steps)
      }
      if(length(roses) == 1) roses[[1]] else combine_roses(roses)
}


#' Manage the wind rose download cache
#'
#' Lists, or deletes, the pre-built wind rose files that [download_wind_rose()] has saved to the
#' cache. The cache is in the user cache directory given by `tools::R_user_dir("windscape",
#' "cache")`; set `options(windscape.cache_dir = ...)` to use a different location.
#'
#' @param clear Logical. If `TRUE`, delete all cached wind rose files.
#' @return A data frame of cached files, with columns `file`, `release`, `bytes`, and `path`
#'   (invisibly, listing the deleted files, if `clear = TRUE`).
#' @seealso [download_wind_rose()]
#' @examples
#' wind_rose_cache()
#' @export
wind_rose_cache <- function(clear = FALSE){
      dir <- rose_cache_dir()
      paths <- if(dir.exists(dir)) list.files(dir, pattern = "\\.tif(\\.part)?$", recursive = TRUE,
                                              full.names = TRUE) else character()
      out <- data.frame(file = basename(paths), release = basename(dirname(paths)),
                        bytes = file.size(paths), path = paths, stringsAsFactors = FALSE)
      if(!isTRUE(clear)) return(out)
      unlink(paths)
      for(d in unique(dirname(paths))){
            empty <- length(list.files(d, all.files = TRUE, no.. = TRUE)) == 0
            if(empty) unlink(d, recursive = TRUE)
      }
      message("deleted ", nrow(out), " file(s), ", format_mb(sum(out$bytes)))
      invisible(out)
}


# Helpers ----------------------------------------------------------------------------------------

# Choose the fewest catalog files covering the requested years x months: the full-period file,
# month-of-year files (all years), annual files (all months), or single months.
rose_files <- function(catalog, year, month){
      yrs <- sort(unique(catalog$year[catalog$tier == "month"]))
      if(is.null(year)) year <- yrs
      if(is.null(month)) month <- 1:12
      if(length(year) == 0 || !is.numeric(year) || anyNA(year) || any(year %% 1 != 0))
            stop("`year` must be a vector of whole years", call. = FALSE)
      if(length(month) == 0 || !is.numeric(month) || anyNA(month) || any(!month %in% 1:12))
            stop("`month` must be a vector of integers from 1 to 12", call. = FALSE)
      bad <- setdiff(year, yrs)
      if(length(bad)) stop("pre-built ", catalog$source[1], " roses are available for ", min(yrs),
                           "-", max(yrs), "; requested years outside this range: ",
                           paste(sort(bad), collapse = ", "), call. = FALSE)
      year <- sort(unique(year))
      month <- sort(unique(month))
      all_y <- setequal(year, yrs)
      all_m <- setequal(month, 1:12)
      tier <- catalog$tier
      if(all_y && all_m){
            sel <- tier == "all"
            n <- 1
      }else if(all_y){
            sel <- tier == "month_of_year" & catalog$month %in% month
            n <- length(month)
      }else if(all_m){
            sel <- tier == "year" & catalog$year %in% year
            n <- length(year)
      }else{
            sel <- tier == "month" & catalog$year %in% year & catalog$month %in% month
            n <- length(year) * length(month)
      }
      out <- catalog[sel, , drop = FALSE]
      if(nrow(out) != n)
            stop("the wind rose catalog is missing files for this request", call. = FALSE)
      out
}

# Validate and standardize a cropping extent
rose_ext <- function(ext){
      e <- if(is.numeric(ext)){
            if(length(ext) != 4 || anyNA(ext))
                  stop("a numeric `ext` must be c(xmin, xmax, ymin, ymax)", call. = FALSE)
            terra::ext(ext)
      }else{
            tryCatch(terra::ext(ext), error = function(e)
                  stop("`ext` must be a SpatExtent, an object with an extent, or ",
                       "c(xmin, xmax, ymin, ymax)", call. = FALSE))
      }
      v <- as.vector(e)
      if(v[1] < -180 || v[2] > 180 || v[3] < -90 || v[4] > 90)
            stop("`ext` must be in degrees, with longitudes from -180 to 180 and latitudes from ",
                 "-90 to 90", call. = FALSE)
      e
}

# Read the cells of a hosted file within `ext`, over HTTP via GDAL's /vsicurl/, into memory
read_rose_remote <- function(url, ext){
      cfg <- c(GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR", GDAL_HTTP_MAX_RETRY = "3",
               GDAL_HTTP_RETRY_DELAY = "1")
      old <- vapply(names(cfg), terra::getGDALconfig, character(1))
      on.exit(for(k in names(cfg)) terra::setGDALconfig(k, old[[k]]), add = TRUE)
      for(k in names(cfg)) terra::setGDALconfig(k, cfg[[k]])
      x <- terra::crop(terra::rast(paste0("/vsicurl/", url)), ext, snap = "out")
      out <- terra::rast(x)
      terra::values(out) <- terra::values(x)
      out
}

# Path to a local copy of a catalog file in `dir`, downloading it if needed
fetch_rose <- function(f, dir, quiet){
      d <- file.path(dir, f$release)
      dir.create(d, showWarnings = FALSE, recursive = TRUE)
      path <- file.path(d, f$file)
      if(file.exists(path) && isTRUE(file.size(path) == f$bytes)) return(path)
      if(!quiet) message("downloading ", f$file, " (", format_mb(f$bytes), ")")
      tmp <- paste0(path, ".part")
      on.exit(unlink(tmp), add = TRUE)
      op <- options(timeout = max(600, getOption("timeout")))
      on.exit(options(op), add = TRUE)
      tryCatch(utils::download.file(f$url, tmp, mode = "wb", quiet = TRUE),
               error = function(e) stop("failed to download ", f$url, ": ", conditionMessage(e),
                                        call. = FALSE),
               warning = function(w) stop("failed to download ", f$url, ": ", conditionMessage(w),
                                          call. = FALSE))
      if(!identical(unname(tools::md5sum(tmp)), f$md5))
            stop("downloaded file ", f$file, " does not match its checksum; try again",
                 call. = FALSE)
      file.rename(tmp, path)
      path
}

format_mb <- function(bytes) paste0(format(round(bytes / 1e6, 1), nsmall = 1), " MB")
