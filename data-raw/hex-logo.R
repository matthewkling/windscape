## Code to create the package hex logo, man/figures/logo.png.
##
## The logo shows real model output from the example rose: a focal point in the southern Plains
## (100 W, 36 N), with its upwind catchment (teal, from the south, where the Great Plains
## low-level jet comes from) and downwind windshed (amber, to the north), from stream-mode random
## walks with a 36-hour half-life. Plumes are drawn as many nested contours of residence time,
## each enclosing a fixed fraction of the airborne mass, giving a soft glow. Lines are flux
## streamlines from random_walk_paths(), each carrying an equal share of material. The background
## "stars" are arrowheads pointing along the rose's net flow at a random subset of cells, jittered
## within their cells.
##
## Run from the package root. The title uses Helvetica Neue Condensed Bold (macOS), falling back
## to Helvetica Neue Bold or Nimbus Sans Narrow Bold if it isn't installed. Everything else is
## deterministic, so reruns give the same logo.
##
## The image bleeds to the hexagon's edges, so instead of hexSticker (which doesn't clip its
## subplot to the hexagon), the plot is rendered with ragg and masked to the hexagon with magick,
## giving a transparent background outside it.

library(ggplot2)
devtools::load_all()

out_file <- tempfile(fileext = ".png")   # full-resolution render; use_logo() resizes it
dpi <- 600

# design parameters ---------------------------------------------------------------------------

site <- cbind(-100, 36)                  # focal point (lon, lat)
half_life <- 36                          # hours
n_paths <- 160                           # streamlines per direction
R <- 520                                 # hexagon circumradius, km
inch_per_R <- 1.176                      # print scale: inches of logo per R; sets the
# sticker size (sqrt(3) x 2 R) that sizes below assume
star_frac <- 0.4                         # fraction of cells with a star
seed <- 7

col <- list(bg = "#14213D",              # night sky
            border = "#0B1426",
            up = "#35C2CF", line_up = "#C9F4F6",   # upwind: teal
            dn = "#F2A33A", line_dn = "#FFE1A8")   # downwind: amber

# title font: first available of these faces
font_choices <- list(c("Helvetica Neue", "Condensed Bold"),
                     c("Helvetica Neue", "Bold"),
                     c("Nimbus Sans Narrow", "Bold"))

# model output --------------------------------------------------------------------------------

rose <- windscape_example("wind_rose")
lon0 <- site[1]
lat0 <- site[2]

# local km coordinates around the focal point, with north up and true shapes at its latitude
km_lat <- 110.57
km_lon <- 111.32 * cos(lat0 * pi / 180)
to_km <- function(d){
      d$X <- (d$x - lon0) * km_lon
      d$Y <- (d$y - lat0) * km_lat
      d
}

# residence time surface, scaled to a maximum of 1 and smoothed by bilinear disaggregation
plume <- function(direction){
      r <- random_walk(rose, site, mode = "stream", direction = direction,
                       half_life = half_life)$residence
      r <- r / max(terra::values(r), na.rm = TRUE)
      r <- terra::disagg(r, 6, method = "bilinear")
      d <- to_km(as.data.frame(r, xy = TRUE))
      names(d)[3] <- "v"
      d
}
up <- plume("upwind")
dn <- plume("downwind")

paths <- function(direction){
      to_km(random_walk_paths(rose, site, n = n_paths, direction = direction,
                              half_life = half_life))
}
paths_up <- paths("upwind")
paths_dn <- paths("downwind")

# contour levels enclosing fractions p of the total mass (highest-density regions)
mass_levels <- function(d, p){
      v <- sort(d$v, decreasing = TRUE)
      cm <- cumsum(v) / sum(v)
      vapply(p, function(q) v[which(cm >= q)[1]], numeric(1))
}
p_levels <- seq(0.97, 0.03, length.out = 24)

# stars: net flow arrowheads at a jittered random subset of cells
set.seed(seed)
nf <- to_km(fortify(net_flow(rose)))
cell <- terra::res(rose)
keep <- runif(nrow(nf)) < star_frac
nf$X <- nf$X + runif(nrow(nf), -0.35, 0.35) * cell[1] * km_lon
nf$Y <- nf$Y + runif(nrow(nf), -0.35, 0.35) * cell[2] * km_lat
stars <- nf[keep & abs(nf$X) < 1.2 * R & abs(nf$Y) < 1.2 * R, ]
rel <- stars$speed / max(stars$speed)
stars$alpha <- 0.3 + 0.6 * rel
stars$dx <- stars$u / stars$speed
stars$dy <- stars$v / stars$speed
stars$len <- 0.035 + 0.06 * rel           # arrowhead length, cm

# title font ----------------------------------------------------------------------------------

fonts <- systemfonts::system_fonts()
family <- NULL
for(fc in font_choices){
      i <- which(fonts$family == fc[1] & fonts$style == fc[2])
      if(length(i)){
            systemfonts::register_font("windscape-logo",
                                       plain = list(fonts$path[i[1]], fonts$index[i[1]]))
            family <- "windscape-logo"
            message("title font: ", fc[1], " ", fc[2])
            break
      }
}
if(is.null(family)) stop("none of the title fonts in `font_choices` is installed")

# plot ----------------------------------------------------------------------------------------

hex_angles <- (90 + 60 * 0:5) * pi / 180  # pointy-top
hex <- data.frame(X = R * cos(hex_angles), Y = R * sin(hex_angles))

contour_layers <- function(d, fill){
      lapply(mass_levels(d, p_levels), function(lev)
            geom_contour_filled(data = d, aes(X, Y, z = v), breaks = c(lev, 2),
                                fill = fill, alpha = 0.05))
}

star_layers <- lapply(split(stars, cut(stars$len, 3)), function(s)
      geom_segment(data = s, aes(X - dx * 4, Y - dy * 4, xend = X, yend = Y, alpha = alpha),
                   color = "white", linewidth = 0.15,
                   arrow = grid::arrow(length = grid::unit(mean(s$len), "cm"),
                                       type = "closed", angle = 18)))

p <- ggplot() +
      geom_polygon(data = hex, aes(X, Y), fill = col$bg) +
      star_layers +
      contour_layers(up, col$up) +
      contour_layers(dn, col$dn) +
      geom_path(data = paths_up, aes(X, Y, group = trail),
                color = col$line_up, alpha = 0.22, linewidth = 0.22) +
      geom_path(data = paths_dn, aes(X, Y, group = trail),
                color = col$line_dn, alpha = 0.22, linewidth = 0.22) +
      annotate("point", 0, 0, size = 1.6, color = "white") +
      annotate("text", 0, -0.47 * R, label = "windscape", color = "white",
               size = 6.6, family = family) +
      scale_alpha_identity() +
      coord_equal(xlim = c(-1, 1) * sqrt(3) / 2 * R, ylim = c(-1, 1) * R, expand = FALSE) +
      theme_void() +
      theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0),
            plot.background = element_blank(), panel.background = element_blank())

# render, mask to the hexagon, and add a border -----------------------------------------------

width <- sqrt(3) * inch_per_R             # hexagon width sqrt(3) R, height 2 R
height <- 2 * inch_per_R

# The border is a filled ring of uniform width between two hexagons, drawn entirely inside the
# image, so its width is the same on every side and nothing is clipped at the edges. The outer
# hexagon is inset slightly so its points aren't cut off by the image boundary.
border_in <- 0.05                         # border width, inches (perpendicular to each side)
r_out <- inch_per_R * (1 - 0.004)         # outer hexagon circumradius, inches
r_in <- r_out - border_in / (sqrt(3) / 2) # inner: apothem smaller by exactly border_in
hex_npc <- function(r) list(x = 0.5 + r * cos(hex_angles) / width,
                            y = 0.5 + r * sin(hex_angles) / height)
outer <- hex_npc(r_out)
inner <- hex_npc(r_in)

render <- function(draw){
      f <- tempfile(fileext = ".png")
      ragg::agg_png(f, width = width, height = height, units = "in", res = dpi,
                    background = "transparent")
      grid::grid.newpage()
      draw()
      invisible(dev.off())
      magick::image_read(f)
}
art <- render(function() grid::grid.draw(ggplotGrob(p)))
mask <- render(function() grid::grid.polygon(outer$x, outer$y,
                                             gp = grid::gpar(col = NA, fill = "black")))
border <- render(function() grid::grid.path(c(outer$x, inner$x), c(outer$y, inner$y),
                                            id = rep(1:2, each = 6), rule = "evenodd",
                                            gp = grid::gpar(col = NA, fill = col$border)))
logo <- magick::image_composite(art, mask, operator = "DstIn")
logo <- magick::image_composite(logo, border, operator = "Over")
magick::image_write(logo, out_file)

# use_logo() resizes to fit `geometry` (default 240x278); twice that keeps fine detail sharp on
# high-DPI screens, since pkgdown and the README display the logo much smaller than this
usethis::use_logo(out_file, geometry = "480x556")
