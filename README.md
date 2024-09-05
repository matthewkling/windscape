windscape: landscape connectivity modeling by wind
================

![](https://matthewkling.github.io/img/images/windscape_hawaii.png)

The `windscape` R package models landscape connectivity by wind
dispersal. It includes functions to help you import wind data, model
wind diffusion times across a region, visualize the results, and test
statistical relationships between wind and your ecological data.

### Installation

`devtools::install_github("matthewkling/windscape")`

### Usage

Here’s a demonstration of a workflow to create a map of estimated wind
travel times from a particular location. See the package vignette for
details and other functionality.

``` r
library(windscape)
library(tidyverse)

site <- vect(matrix(c(-105, 40), 1)) # lat-lon coordinates of a focal site

rast(system.file("extdata/wind.tif", package = "windscape")) %>% # load wind time series rasters
      wind_field(order = "uuvv") %>% # convert to a formal wind field object
      wind_rose() %>% # summarize into wind rose conductance object
      wind_graph("downwind") %>% # format as connectivity graph
      least_cost_surface(site) %>% # calculate least cost path from coordinates
      plot(main = "outbound wind travel time (hours)") # plot map of travel times

points(site, col = "red") # add origin location to map
```

![](man/figures/example-1.png)<!-- -->

<!-- ### In the wild -->
<!-- Publications that have used the `windscape` framework include: -->
<!-- -   **Kling, M.**, and D. Ackerly. (2021) Global wind patterns shape genetic differentiation, asymmetric gene flow, and genetic diversity in trees. Proceedings of the National Academy of Sciences, 118(17) [<https://doi.org/10.1073/pnas.2017317118>] -->
<!-- -   **Kling, M.**, and D. Ackerly. (2020) Global wind patterns and the vulnerability of wind-dispersed species to climate change. Nature Climate Change, 10: 868-875 [<https://doi.org/10.1038/s41558-020-0848-3>] -->
