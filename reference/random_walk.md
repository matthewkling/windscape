# Simulate wind dispersal by random walk

Estimates the dispersal of particles across a windscape with a Markov
random walk on the grid. Each iteration, the "particle mass" in each
grid cell is dispersed within its 9-cell neighborhood in proportion to
wind conductance. Depending on the application, this mass could
represent probability, numbers of seeds or spores, etc. Particles can be
released once (`mode = "pulse"`) or continuously (`mode = "stream"`),
and can be deposited along the way (`half_life`).

## Usage

``` r
random_walk(
  rose,
  init,
  mode = c("pulse", "stream"),
  direction = c("downwind", "upwind"),
  half_life = Inf,
  timescale = 1,
  latitude_correction = TRUE,
  density = TRUE,
  iter = 100,
  record = iter,
  method = c("auto", "solve", "iterate"),
  tol = 1e-08,
  max_iter = 1e+05,
  flux = FALSE,
  source = NULL,
  progress = FALSE
)
```

## Arguments

- rose:

  A `wind_rose`.

- init:

  Where particles are released: either a two-column matrix of
  coordinates, which places a mass of 1 in each cell containing a
  coordinate, or a single-layer SpatRaster of non-negative mass values
  with the same geometry as `rose`. In stream mode, `init` can be read
  either as a release rate (mass per hour) or as a one-time release; see
  Value. With `direction = "upwind"`, `init` instead gives the receptor
  locations (and weights) whose sources are traced.

- mode:

  Release mode. `"pulse"` (the default) releases `init` once and tracks
  the drifting, spreading cloud through time. `"stream"` releases
  particles continuously and returns the steady state, where release is
  balanced by deposition and loss across domain edges; equivalently,
  this is the all-time total for a single release (see Value). It is
  computed directly rather than by simulating forward.

- direction:

  `"downwind"` (the default) simulates where particles released at
  `init` go: their outbound windshed. `"upwind"` simulates where
  particles reaching `init` come from: its inbound windshed, or
  catchment. See Value and Details.

- half_life:

  Half-life of airborne mass, in hours of transport time (assuming
  `trans = 1` and wind speeds in m/s; see
  [`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)).
  Particles are deposited at the continuous rate
  `k = log(2) / half_life`, in the cell where they are airborne. Default
  `Inf`: no deposition.

- timescale:

  A value between 0 and 1 that scales the time step length (see
  [`rw_max_step()`](https://matthewkling.github.io/windscape/reference/rw_max_step.md)).
  The default of 1 uses the longest possible step, advancing the
  simulation farthest per iteration. Smaller values give more accurate
  pulse-mode dynamics or set the step to a desired duration. Stream-mode
  results do not depend on `timescale`.

- latitude_correction:

  Logical: correct the compression of east-west spread on
  longitude/latitude grids, where cells narrow toward the poles? Default
  `TRUE`. This adds conductance to east and west edges without changing
  drift, which shortens the time step (pulse mode needs about 1.3 times
  as many iterations at 45 degrees, 1.9 at 60). Has no effect on
  projected rasters. See Details.

- density:

  Logical: express results per km^2 rather than per grid cell? Default
  `TRUE`, matching
  [`pairwise_random_walk()`](https://matthewkling.github.io/windscape/reference/pairwise_random_walk.md).
  On a longitude/latitude grid, cell area shrinks toward the poles, so
  per-cell values are biased toward lower latitudes; per-km^2 values
  remove that bias and are comparable across grid resolutions. Downwind,
  each cell's values are divided by its own area. Upwind, values are per
  particle released from each origin, so they are instead divided by the
  area of the receptor (each receptor's own, if several); with a single
  receptor this rescales the whole map by a constant. `flux` is divided
  by each cell's area in either direction. Use `density = FALSE` for
  per-cell masses, which sum to totals (e.g. in the pulse-mode mass
  balance); equivalently, multiply per-km^2 results by
  [`terra::cellSize()`](https://rspatial.github.io/terra/reference/cellSize.html).

- iter:

  Pulse mode only: number of iterations (time steps) to simulate.

- record:

  Pulse mode only: integer vector of iterations to return, between 0
  (the initial state) and `iter`. Default is the final iteration only.

- method:

  Stream mode only: how to compute the steady state. `"solve"` is exact,
  using a sparse LU solve. `"iterate"` repeats the stream update rule
  until the estimated error is below `tol`; it uses less memory for very
  large grids. `"auto"` (the default) uses `"solve"` for grids of up to
  5e5 cells and `"iterate"` otherwise.

- tol:

  Stream mode with `method = "iterate"` only: convergence tolerance, as
  estimated L1 error relative to total airborne mass. For finite
  `half_life` the error estimate is a rigorous bound; for
  `half_life = Inf` it is extrapolated from the observed convergence
  rate.

- max_iter:

  Stream mode with `method = "iterate"` only: maximum number of
  iterations, after which a warning is given.

- flux:

  Stream mode only: logical: also return the net flux of material, the
  direction and rate at which it moves through each cell? Default
  `FALSE`. For upwind walks, requires a finite `half_life`. See Value.

- source:

  Upwind stream mode only: where material is released, used for `origin`
  and `flux`: either a two-column matrix of coordinates (a unit release
  in each cell containing a coordinate) or a single-layer SpatRaster of
  non-negative release per grid cell, like `init`. Default `NULL`
  releases uniformly per km^2 (each cell in proportion to its area). To
  use a map of release per km^2, multiply it by
  [`terra::cellSize()`](https://rspatial.github.io/terra/reference/cellSize.html).

- progress:

  Pulse mode only: logical: show a progress bar while iterating? Default
  `FALSE`.

## Value

A named list of two `wind_walk` rasters: `airborne` and `deposition` in
pulse mode, or `residence` and `deposition` in stream mode. Masses are
in the units of `init`, per km^2 by default or per grid cell with
`density = FALSE` (see `density`), and hours assume `trans = 1` and wind
speeds in m/s. Cells that are NA in `rose` are NA in the output.

**Pulse mode.** Each element has one layer per iteration in `record`
(named `iter0`, `iter1`, ...):

- `airborne`: the mass airborne over each cell at that iteration.

- `deposition`: the cumulative mass deposited in each cell up to that
  iteration. Zero everywhere if `half_life = Inf`.

At every iteration, airborne mass plus deposition plus mass lost across
domain edges equals the released mass (summing per-cell values, i.e.
with `density = FALSE`). Use
[`iter_length()`](https://matthewkling.github.io/windscape/reference/iter_length.md)
to convert iterations to hours.

**Stream mode.** Each element is a single layer, and neither depends on
`timescale`. The two differ in units by a factor of time:

- `residence`: airborne mass integrated over time, in mass-hours (units
  of `init` times hours).

- `deposition`: deposited mass, in units of `init`. It equals
  `k * residence`, where `k = log(2) / half_life` is the deposition rate
  per hour, so it is zero everywhere if `half_life = Inf`.

Stream results have two equivalent readings. They are the steady state
when `init` is released continuously as a rate (mass per hour), in which
case `residence` is the steady-state airborne mass and `deposition` the
steady-state deposition rate (mass per hour). They are also the all-time
totals for a single release of `init`, in which case `residence` is the
total mass-hours spent airborne over each cell and `deposition` the
total mass eventually deposited; for a unit release from one cell,
`deposition` is the probability distribution of where a particle lands
(a probability density per km^2 by default, or a probability per cell
with `density = FALSE`).

**Net flux.** With `flux = TRUE`, the list also contains `flux`, a
[`wind_field()`](https://matthewkling.github.io/windscape/reference/wind_field.md)
whose `u` and `v` layers are the eastward and northward components of
the net transport of material through each cell: the flow from the cell
to each neighbor minus the flow back, times the distance to that
neighbor, halved (each flow is shared between two cells). Units are mass
times km per hour (if `init` is a release rate, `trans = 1`, and wind
speeds are in m/s); direction and relative magnitude are usually what
matter. Draw it with
[`geom_wind_arrow()`](https://matthewkling.github.io/windscape/reference/geom_wind_arrow.md),
or
[`stat_wind_trail()`](https://matthewkling.github.io/windscape/reference/geom_wind_trail.md)
with `fixed_length = TRUE`. Unlike the gradient of `residence` or
`deposition`, which describes the shape of the windshed, net flux shows
the transport that produces it: across a plume's flanks, for example,
material moves mostly downwind, not sideways down the gradient. Each
cell's net outflow equals its release minus its deposition.

For upwind walks, `flux` is the transport of the material that ends up
at the receptor: material released per `source` (by default, uniformly,
one unit per km^2 per hour) and eventually deposited at the receptor
(weighted by `init`). It shows the routes by which the receptor's
material arrives. Each cell's net outflow equals its release of such
material (its release times its probability of eventual deposition at
the receptor) minus the amount deposited in the cell, which is nonzero
only at the receptor.

**Relationship between modes.** For the same `rose`, `init`,
`half_life`, and `latitude_correction`, stream results are the totals
that a pulse walk accumulates over all time, at any `timescale`. Stream
`deposition` equals pulse `deposition` once the pulse has run until
negligible mass remains airborne, and stream `residence` equals pulse
`airborne` mass integrated over time (exactly, `(1 - lambda) * t` times
its sum over iterations 0, 1, 2, ...; see Details). This holds because
the walk is linear (particles move independently, so their contributions
add) and time-invariant (the same `rose` applies at every time step),
which is also why continuous and single releases give the same numbers.

**Upwind walks.** With `direction = "upwind"`, each cell's values
describe particles released from that cell and their contribution to the
receptor at `init` (weighted by its values, if a raster):

- Stream `deposition`: the probability that a particle released from the
  cell is deposited at the receptor. This maps where the receptor's
  deposited material comes from.

- Stream `residence`: the time, in hours, that a unit release from the
  cell spends airborne over the receptor.

- Pulse `airborne` at iteration `k`: the probability that a particle
  released from the cell is airborne over the receptor `k` time steps
  later; and pulse `deposition`: the probability it has been deposited
  there by then.

- Stream `origin` (with a finite `half_life`): the probability that a
  particle deposited at the receptor was released from the cell. Where
  `deposition` asks "if a particle were released here, would it reach
  the receptor?", `origin` asks "of the particles reaching the receptor,
  what share came from here?" By Bayes' rule, `origin` is `deposition`
  times release (see `source`), normalized to sum to 1 over cells, or
  with `density = TRUE`, to integrate to 1 per km^2. With the default
  uniform release per km^2, `origin` with `density = TRUE` is
  proportional to `deposition`. Shares are among sources inside the
  domain: material from beyond its edges is not represented. For a
  receptor spanning several cells, `origin` covers particles landing
  anywhere in it.

## Details

### Transition probabilities

The wind rose is converted to a simplex of nine probabilities for each
grid cell, giving the rates at which particles remain in the cell or
move to each of its eight neighbors. Probabilities of moving to a
neighbor are proportional to conductance, and retention probabilities
are scaled so that the cell with the highest total conductance has zero
retention. This keeps conductance proportional across cells while
maximizing the dispersal occurring at each iteration.

This is uniformization of the continuous-time Markov chain defined by
the conductances, which are rates (in units of 1 / hours, if `trans = 1`
and wind speeds are in m/s). With time step `t` and rate matrix `Q`, the
transition matrix is `P = I + tQ`. Pulse iterations are therefore a
first-order approximation of the continuous-time process at times `t`,
`2t`, ...; reducing `timescale` improves their accuracy.

### Deposition and the two release modes

Each time step, airborne mass disperses and a fraction
`lambda = kt / (1 + kt)` of it is deposited, where
`k = log(2) / half_life`. This is the uniformization of transport and
deposition together, so stream-mode results are exact solutions of the
continuous-time process and independent of `timescale`. In pulse mode,
mass halves after approximately `half_life` hours (exactly in the limit
of small `timescale`).

The pulse update rule is `n <- (1 - lambda) * disperse(n)`, with
`lambda * n` deposited in place each step before dispersal. The stream
update rule is `n <- n0 + (1 - lambda) * disperse(n)`, with fixed point
`n* = (I - (1 - lambda) P')^-1 n0`; residence is
`(1 - lambda) * t * n*`. The two modes are linked exactly: stream-mode
residence equals `(1 - lambda) * t` times the sum of the pulse-mode
surfaces over all time steps, starting from step 0.

### Upwind walks

Write the downwind stream solution as `n* = G n0`, where entry `G[r, s]`
is the airborne mass at cell `r` per unit released at cell `s`. A
downwind walk from a source `s` gives column `s` of `G`: everywhere that
source's particles go. An upwind walk to a receptor `r` gives row `r`:
every source's contribution to that receptor. It is computed with the
adjoint process, using the transpose of the downwind transition matrix
(`P` in place of `P'`), not by reversing the wind, so for any source `s`
and receptor `r`, the upwind value at `s` equals the downwind value at
`r`. The same holds for each pulse iteration.

### Domain edges

Domain edges are absorbing: mass that disperses off the grid, or into NA
cells, is lost and never returns. Values near edges are therefore biased
low, because they receive no inflow from beyond the edge. In stream mode
with finite `half_life`, the fraction of released mass lost across edges
is reported as a message; the domain should be buffered by several decay
lengths beyond the area of interest, until this fraction is small or the
bias in the area of interest is acceptable. With `half_life = Inf`,
edges are the only place stream-mode mass can leave, so the result
measures connectivity within the chosen domain and depends on its
extent; every valid cell must have a path to an edge or NA cell, or an
error is raised.

### Latitude correction

Conductances account for latitude (see
[`wind_rose()`](https://matthewkling.github.io/windscape/reference/wind_rose.md)),
so the speed at which mass drifts is correct at all latitudes for winds
aligned with a neighbor direction. However, the spread of mass in a
nearest-neighbor walk is numerical diffusion whose magnitude scales with
hop length. Because longitude/latitude cells narrow east-west toward the
poles, east-west hops are shorter than north-south hops, and without
correction east-west spread is compressed: for isotropic wind, to
roughly 0.86, 0.70, and 0.49 times north-south spread at 30, 45, and 60
degrees latitude.

The correction scales each cell's east-west spread rate up to what the
same rose would produce on square cells of the same north-south size, by
adding equal conductance to its east and west edges. Because these two
neighbors are exactly opposite, drift is unchanged. In tests against the
same winds on square cells, the correction brings east-west spread to
within about 8 percent of the square-cell value across a range of wind
regimes at 45 and 60 degrees. It does not correct the orientation
(east-west/north-south covariance) of spread for oblique winds, which
remains underestimated at high latitude, nor the small drift bias for
winds blowing between neighbor directions (up to about 8 percent at the
equator, more at high latitude). The added conductance shortens the time
step by an amount set by the highest-latitude cells in the domain.
Spread remains dependent on grid resolution with or without the
correction.
