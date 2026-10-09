# Maximum iteration duration for a random walk

Calculate the maximum possible iteration length for a random walk
simulation, which is the residence time of the grid cell with the
greatest total conductance. Walks can proceed no faster than this.

## Usage

``` r
rw_max_step(rose)
```

## Arguments

- rose:

  A `wind_rose`.

## Value

A number: the residence time of the cell with the greatest total
conductance, in hours for a rose built from wind speeds in m/s with
`trans = 1`.
