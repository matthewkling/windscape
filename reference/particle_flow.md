# Simulate the flow of particles across a wind field

This function simulates the movement of particles over a wind field. At
each iteration, particles are displaced in accordance with the direction
and speed of the local wind conditions.

## Usage

``` r
particle_flow(
  x,
  p,
  n_iter = 100,
  scale = 0.01,
  wrap = "neither",
  direction = "both",
  ignore_speed = FALSE,
  sf = FALSE,
  ...
)
```

## Arguments

- x:

  A `wind_field` object created by
  [wind_field](https://matthewkling.github.io/windscape/reference/wind_field.md).

- p:

  Either a two-column lon-lat matrix of initial particle locations (e.g.
  produced by
  [generate_particles](https://matthewkling.github.io/windscape/reference/generate_particles.md)),
  or an integer giving the number of particles to initialize (in which
  case, provide additional arguments to
  [generate_particles](https://matthewkling.github.io/windscape/reference/generate_particles.md)
  in `...`).

- n_iter:

  Positive integer giving the number of iterations to simulate.

- scale:

  A multiplier to convert step sizes from velocity (in whatever then
  native units of `x` are) to degrees lon-lat. Smaller scale values
  yield smoother but shorter particle trails over a given number of
  iterations.

- wrap:

  Either "vertical", "horizontal", "both", or "neither" (the default).

- direction:

  Either "downwind", "upwind", or "both" (the default).

- ignore_speed:

  Logical: if TRUE, only direction is considered and all particles are
  displaced the same distance. Default is FALSE.

- sf:

  Logical indicating whether the output shoudl be break_trailted as an
  \`sf\` object instead of the default data frame.

- ...:

  Further arguments to
  [generate_particles](https://matthewkling.github.io/windscape/reference/generate_particles.md),
  used only if `p` is an integer.

## Value

If \`sf = FALSE\`, a data frame representing particle positions over
time. Variables include "p" (the id number of the particle), "x"
(longitude), "y" (latitude), "v" (velocity, in the same units as the
input wind data), and "t" (the iteration number, which is positive for
downwind and negative for upwind; arrange by this variable when
plotting). If \`sf = TRUE\`, a simple feature collection of linestrings
(in this case, velocity is not included, and the iteration index is
included as the third coordination dimension).
