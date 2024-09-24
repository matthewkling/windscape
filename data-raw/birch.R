
## code to prepare `birch` dataset, sourced from Tsuda et al. 2017

library(tidyverse)
library(diveRsity)

sites <- system.file("data-raw/tsuda2017_populations.csv", package = "windscape") %>%
      read_csv() %>%
      select(x = longitude, y = latitude) %>%
      as.matrix()

genes <- system.file("data-raw/tsuda2017_genepop.txt", package = "windscape")

mig <- divMigrate(genes, stat = "d")$dRelMig
div <- divBasic(genes)[["Allelic_richness"]] %>% tail(1) %>% as.vector()
fst <- diffCalc(genes, fst = TRUE, pairwise = TRUE)$pairwise$Fst
fst[upper.tri(fst)] <- t(fst)[upper.tri(fst)]

birch <- list(sites = sites,
              div = div,
              mig = mig,
              fst = fst)

usethis::use_data(birch, overwrite = TRUE)
