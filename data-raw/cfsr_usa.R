## code to prepare `cfsr_usa` dataset

cfsr_usa <- cfsr_dl(variable = "wnd10m", years = 2000, months = 6:9, days = seq(1, 28, 3),
                    hlim = c(7, 18) + 5, # 0700-1800 US central time, converted to UTC
                    xlim = c(-120, -90) + 360, # shift longitudes to be in [0, 360] range for CFSR
                    ylim = c(30, 50)) %>%
      shift(dx = -360) # shift longitude back to the standard [-180, 180] range

writeRaster(cfsr_usa, "inst/extdata/cfsr_usa.tif", overwrite = TRUE)

cfsr_usa <- wrap(cfsr_usa)

usethis::use_data(cfsr_usa, overwrite = TRUE)
