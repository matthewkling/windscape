# Download land-water raster from the Climate Forecast System Reanalysis (CFSR)

This function downloads a CFSR land-water raster layer from the NCAR
THREDDS server, for an area specified by a spatial bounding box.

## Usage

``` r
cfsr_dl_land(xlim = c(260, 270), ylim = c(40, 50))
```

## Arguments

- xlim:

  Vector of length 2, giving longitudinal limits of the data to
  download, in the range 0 to 360.

- ylim:

  Vector of length 2, giving latitudinal limits of the data to download,
  in the range -90 to 90.

## Value

A SpatRaster layer with 1 indicating land and 0 indicating water.

## Details

(Actually, since CFSR does not publish a land layer per se, this
function downloads a soil temperature layer for an arbitrary date and
classifies areas with data as land and no data as water.)
