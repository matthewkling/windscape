#' setGeneric("plot", function(x) standardGeneric("plot"))
#'
#' #' Plot a wind_rose
#' #'
#' #' @param x A wind_rose object.
#' #' @include wind_rose.R
#' #' @return A ggplot.
#' plot.wind_rose <- function(x){
#'       require(tidyverse)
#'       if(nlyr(x) != 8) return(terra::plot(as(x, "SpatRaster")))
#'       as.data.frame(x, xy = T) %>%
#'             as_tibble() %>%
#'             gather(layer, value, -x, -y) %>%
#'             mutate(layer = factor(layer, levels = c("NW", "N", "NE", "W", "", "E", "SW", "S", "SE"))) %>%
#'             ggplot(aes(x, y, fill = value)) +
#'             facet_wrap(~ layer, nrow = 3, drop = F) +
#'             geom_raster() +
#'             scale_fill_viridis_c() +
#'             theme_void() +
#'             labs(fill = "conductance")
#' }
#' setMethod("plot", "wind_rose", plot.wind_rose)
