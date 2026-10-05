# Changelog

## windscape (development version)

### Breaking changes

- [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md):
  `"ratchet"` mode removed; results are now a named list of rasters;
  arguments reordered (general arguments before mode-specific ones). The
  new latitude correction (on by default) changes results on lon/lat
  grids.
- `cfsr_usa` and `katrina` datasets replaced by
  [`windscape_example()`](https://matthewkling.github.io/windscape/reference/windscape_example.md).
- [`wind_series()`](https://matthewkling.github.io/windscape/reference/wind_series.md)
  and
  [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)
  now require lon/lat grids with (near-)square cells.

### New features

- [`random_walk()`](https://matthewkling.github.io/windscape/reference/random_walk.md)
  gains `mode = "stream"` (steady state of continuous release, returning
  residence time and deposition), `half_life` (deposition), and
  `latitude_correction`.
- [`rw_self_retention()`](https://matthewkling.github.io/windscape/reference/rw_self_retention.md)
  for leave-one-out corrections of stream-mode results.
- [`windscape_example()`](https://matthewkling.github.io/windscape/reference/windscape_example.md)
  loads small example data sets.

### Bug fixes

- `particle_flow()`: east-west displacement now accounts for latitude,
  and `ignore_speed = TRUE` no longer applies `scale` twice.
- `as_wind_rose()` and `read_wind_rose()` now work with numeric `trans`.
- `least_cost_distance(adjust = TRUE)` no longer returns `NaN` on the
  diagonal.
- windscape no longer masks
  [`gdistance::geoCorrection()`](https://AgrDataSci.github.io/gdistance/reference/geoCorrection.html).
- Removed hidden dependencies on dplyr and magrittr being attached.

## windscape 1.1.0

- Initial release.
