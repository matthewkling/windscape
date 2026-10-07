# Manage the wind rose download cache

Lists, or deletes, the pre-built wind rose files that
[`download_wind_rose()`](https://matthewkling.github.io/windscape/reference/download_wind_rose.md)
has saved to the cache. The cache is in the user cache directory given
by `tools::R_user_dir("windscape", "cache")`; set
`options(windscape.cache_dir = ...)` to use a different location.

## Usage

``` r
wind_rose_cache(clear = FALSE)
```

## Arguments

- clear:

  Logical. If `TRUE`, delete all cached wind rose files.

## Value

A data frame of cached files, with columns `file`, `release`, `bytes`,
and `path` (invisibly, listing the deleted files, if `clear = TRUE`).

## See also

[`download_wind_rose()`](https://matthewkling.github.io/windscape/reference/download_wind_rose.md)

## Examples

``` r
wind_rose_cache()
#> [1] file    release bytes   path   
#> <0 rows> (or 0-length row.names)
```
