# windscape (development version)

## Breaking changes

* `random_walk()`: `"ratchet"` mode removed; results are now a named list of rasters;
  arguments reordered (general arguments before mode-specific ones). The new latitude
  correction (on by default) changes results on lon/lat grids.
* `cfsr_usa` and `katrina` datasets replaced by `windscape_example()`.
* `wind_series()` and `wind_rose()` now require lon/lat grids with (near-)square cells.

## New features

* `random_walk()` gains `mode = "stream"` (steady state of continuous release, returning
  residence time and deposition), `half_life` (deposition), and `latitude_correction`.
* `rw_self_retention()` for leave-one-out corrections of stream-mode results.
* `windscape_example()` loads small example data sets.

## Bug fixes

* `particle_flow()`: east-west displacement now accounts for latitude, and
  `ignore_speed = TRUE` no longer applies `scale` twice.
* `as_wind_rose()` and `read_wind_rose()` now work with numeric `trans`.
* `least_cost_distance(adjust = TRUE)` no longer returns `NaN` on the diagonal.
* windscape no longer masks `gdistance::geoCorrection()`.
* Removed hidden dependencies on dplyr and magrittr being attached.

# windscape 1.1.0

* Initial release.
