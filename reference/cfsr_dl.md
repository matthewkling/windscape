# Download wind data from the Climate Forecast System Reanalysis (CFSR)

This function downloads hourly CFSR wind data from the NCAR THREDDS
server, for one or more historic dates, within a spatial bounding box.
Data for all factorial combinations of 'years', 'months', 'days', and
'hours' are downloaded and combined.

## Usage

``` r
cfsr_dl(
  variable = "wnd10m",
  years = 1979,
  months = 1,
  days = 1,
  hlim = c(0, 23),
  xlim = c(260, 270),
  ylim = c(40, 50)
)
```

## Arguments

- variable:

  Name of the CFSR wind variable to download. Options include "wnd10m"
  (the default, representing wind speed 10 m above the ground), as well
  as "wnd1000", "wnd850", "wnd700", "wnd500", and "wnd200" (representing
  higher altitudes, with the numbers indicating pressure levels in hPa).

- years:

  Integer vector representing the historic year(s) for which data are to
  be accessed, in the range 1979-2010.

- months:

  Integer vector representing the month(s) for which data are to be
  accessed, in the range 1-12.

- days:

  Integer vector representing day(s) for which data is to be requested,
  in the range 0-31.

- hlim:

  Vector of length 2, containing integers between 0 and 23, giving the
  range of hours (in UTC time zone) to download.

- xlim:

  Vector of length 2, giving longitudinal limits of the data to
  download, in the range 0 to 360.

- ylim:

  Vector of length 2, giving latitudinal limits of the data to download,
  in the range -90 to 90.

## Value

A list containing u and v components of the wind field, each as a
SpatRaster with a layer for each hour of each requested date.
