# windscape: time-integrated wind dispersal modeling in R

The `windscape` R package contains functions to model landscape connectivity by wind dispersal. While prior implementations like the `rWind` package support connectivity modeling based on a snapshot of wind conditions, the functions implemented here are designed to integrate over short-term or long-term time series of wind conditions across a landscape, which can be used to predict wind connectivity across scales ranging from days to millennia. These models predict directional wind travel times between any pair of sites, and can therefore be used to generate upwind immigration and downwind emigration accessibility surfaces for a focal site such as in the example below.

![windsheds](https://github.com/matthewkling/matthewkling.github.io/blob/608ef040272ef1f44efc3a8519cbf8e99e39c9ef/images/windscape_demo.png)

### Installation

`devtools::install_github("matthewkling/windscape")`

### Use

The core workflow currently supported by the package begins with a raster stack representing a time series of windfields, summarizes this into a set of time-integrated wind conductance values between every grid cell and its eight neighbors, and then converts this into an asymmetric transition object (defined in the `gdistance` package) representing either upwind or downwind connectivity. The example code below demonstrates a workflow that reproduces the figure above.

```r
library(windscape)
library(rWind)
library(raster)
library(tidyverse)
library(gdistance)

# download some wind time series data using the rWind package
# we'll grab 24 days distributed across a single year
times <- expand.grid(yyyy = 2015, mm = 1:12, dd = c(1, 15), tt = 12)
wd <- pmap(times, wind.dl, lon1 = -120, lon2 = -90, lat1 = 30, lat2 = 50)

# convert to raster stack
u <- wd %>%
      map(select, lon, lat, ugrd10m) %>%
      map(rasterFromXYZ) %>%
      stack()
v <- wd %>%
      map(select, lon, lat, vgrd10m) %>%
      map(rasterFromXYZ) %>%
      stack()
wind <- stack(u, v)

# convert wind time series into a wind conductance raster
conductance <- wind %>%
   windrose_rasters(outfile="windrose.tif", order="uuvv", p=1) %>% 
   add_coords()

# generate downwind and upwind dispersal surfaces for one focal site 
origin <- matrix(c(-105, 40), ncol=2)
downwind <- conductance %>% 
   transition_stack(windflow, directions=8, symm=F, direction="downwind") %>%
   accCost(origin)
upwind <- conductance %>% 
   transition_stack(windflow, directions=8, symm=F, direction="upwind") %>%
   accCost(origin)

# convert seconds to hours, resturcture data, and plot
d <- stack(downwind / 3600, upwind / 3600) %>%
   rasterToPoints() %>%
   as.data.frame() %>%
   rename(downwind=layer.1, upwind=layer.2) %>%
   gather(direction, wind_hours, downwind, upwind)
ggplot() +
   facet_wrap(~direction) +
   geom_raster(data=d, aes(x, y, fill=wind_hours)) +
   geom_point(data=as.data.frame(origin), aes(V1, V2)) +
   scale_fill_gradientn(colors=c("yellow", "red", "blue", "black")) +
   theme_void() +
   theme(legend.position="top",
         strip.text=element_text(size=15))
```
