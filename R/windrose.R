
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


# direction from u&v components -- in degrees, clockwise from 12:00
direction <- function(x, y) atan2(x, y) * 180 / pi
spin90 <- function(x){
      x <- x - 90
      x[x<(-180)] <- x[x<(-180)] + 360
      x
}

# convert an angle [0-45] to a loading [0-1] in the 45-degree direction
loading <- function(theta){ # theta in [0,45]
      sin((4*(theta-22.5))*pi/180)/2+.5
}

#velocity-weighted frequency of wind in each quadrant
windrose <- function(x, weighting="velocity"){
      #x0 <- s[277,349]

      wfx <- switch(weighting,
                    time = function(x) 1,
                    velocity = function(x) sqrt(sum(x^2)),
                    force = function(x) sqrt(sum(x^2))^3)

      # restructure: row=timestep, col=u&v components
      m <- matrix(x, ncol=2, byrow=F)

      # calcualte force and direction
      frc <- apply(m, 1, wfx) # fun should be force or velocity
      dir <- apply(m, 1, function(x) spin90(direction(x[2], -1*x[1])))

      # octant bounds for each vector
      d1 <- plyr::round_any(dir, 45, floor)
      d1[d1 == -180] <- 180
      d2 <- plyr::round_any(dir, 45, ceiling)
      d2[d2 == -180] <- 180

      # loadings on each octant bound, based on angle
      theta <- dir %% 45
      l1 <- loading(theta)
      l2 <- 1-l1

      # add some dummy data to ensure all quadrants are represented
      frc <- c(frc, rep(0,8))
      l1 <- c(l1, rep(0,8))
      l2 <- c(l2, rep(0,8))
      d1 <- c(d1, c(-135, -90, -45, 0, 45, 90, 135, 180))
      d2 <- c(d2, c(-135, -90, -45, 0, 45, 90, 135, 180))

      # repackage
      d <- c(d1, d2)
      l <- c(l1, l2)
      frc <- c(frc, frc)

      # summarize
      d <- split(frc*l, d)
      d <- unlist(lapply(d, sum))
      return(as.numeric(d))
}
