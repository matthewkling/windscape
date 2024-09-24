## code to prepare `katrina` dataset

library(windscape)
library(tidyverse)

katrina <- cfsr_dl(years = 2005, months = 8, days = 28, hlim = c(23, 23),
             xlim = c(-99, -78) + 360, ylim = c(17, 35)) %>%
      shift(dx = -360) %>%
      wrap()

usethis::use_data(katrina, overwrite = TRUE)
