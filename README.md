windscape: time-integrated wind dispersal modeling in R
================

The `windscape` R package contains functions to model landscape
connectivity by wind dispersal. While prior implementations like the
`rWind` package support connectivity modeling based on a snapshot of
wind conditions, the functions implemented here are designed to
integrate over short-term or long-term time series of wind conditions
across a landscape, which can be used to predict wind connectivity
across scales ranging from days to millennia. These models estimate
directional wind travel times between any pair of sites, and can
therefore be used to generate upwind immigration and downwind emigration
accessibility surfaces for a focal site such as in the example below.

![](https://matthewkling.github.io/img/images/windscape_hawaii.png)

### Installation

`devtools::install_github("matthewkling/windscape")`

### Use

The core workflow currently supported by the package begins with a
raster stack representing a time series of windfields, summarizes this
into a set of time-integrated wind conductance values between every grid
cell and its eight neighbors, and then converts this into an asymmetric
transition object (defined in the `gdistance` package) representing
either upwind or downwind connectivity.

The example code below shows a workflow that models and visualizes
upwind and downwind accessibility for a focal site in North America,
based on wind data for a few days in the year 2000. A more detailed set
of vignettes is in development.

``` r
library(windscape)
library(raster)
library(tidyverse)
library(gdistance)

# Download some hourly wind data from the Climate Forecast System Reanalysis.
# We'll grab data for a chunk of North America for 24 days distributed across a single year,
# though a proper study would want to use a denser sampling of data
wind <- cfsr_dl(years = 2000, months = 1:12, days = c(1, 15),
                xlim = c(-120, -90) + 360, # shift longitudes to be in [0, 360] range for CFSR
                ylim = c(30, 50))

# convert to raster::stack and shift longitude to the standard [-180, 180] range
wind <- shift(stack(wind), dx = -360)

# summarize wind time series (n = 1152 hourly layers) into 
# a wind conductance or "windrose" raster (n = 8 directional layers)
conductance <- windrose_rasters(wind, order = "uuvv", p = 1)

# generate downwind and upwind dispersal surfaces for one focal site 
site <- matrix(c(-105, 40), ncol = 2)
downwind <- build_wind_graph(conductance, "downwind") %>%
      accCost(site)
upwind <- build_wind_graph(conductance, "upwind") %>%
      accCost(site)

# restructure data and plot
d <- stack(downwind, upwind) %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      gather(direction, wind_hours, layer.1, layer.2) %>%
      mutate(direction = recode(direction, 
                                layer.1 = "downwind (outbound)", 
                                layer.2 = "upwind (inbound)"))
ggplot(d) +
      facet_wrap(~direction) +
      geom_raster(aes(x, y, fill = wind_hours)) +
      geom_contour(aes(x, y, z = wind_hours), bins = 20, color = "white", linewidth = .25) +
      geom_point(data = as.data.frame(site), aes(V1, V2)) +
      coord_fixed(ratio = 1.2) +
      scale_fill_gradientn(colors=c("yellow", "red", "blue", "black")) +
      theme_void() +
      theme(legend.position="top",
            strip.text=element_text(size=15))
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

### Web tool

This R package powers the [WINDSCAPE WEB
APP](http://matthewkling.net/shiny/windscape/) (still under development
and in beta) that models and visualizes the wind dispersal landscape for
any location on earth. You just click the map to select a site, choose
between “inbound” (upwind) or “outbound” (downwind) dispersal, and the
app generates a global map of relative accessibility by wind. The app
utilizes decades of hourly wind data from the [Climate Forecast System
Renanalysis](https://cfs.ncep.noaa.gov/cfsr/) to estimate long-term
average wind travel times between locations.

### Publications

Publications that have used `windscape` include:

- **Kling, M.**, and D. Ackerly. (2021) Global wind patterns shape
  genetic differentiation, asymmetric gene flow, and genetic diversity
  in trees. Proceedings of the National Academy of Sciences, 118(17)
  \[<https://doi.org/10.1073/pnas.2017317118>\]
- **Kling, M.**, and D. Ackerly. (2020) Global wind patterns and the
  vulnerability of wind-dispersed species to climate change. Nature
  Climate Change, 10: 868-875
  \[<https://doi.org/10.1038/s41558-020-0848-3>\]
