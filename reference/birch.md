# Silver birch landscape genetic data from Tsuda et al. (2017)

An example landscape genetic dataset for the species silver birch
(Betula pendula) across 30 sampling sites in Asia, originally published
by Tsuda et al. (2017).

## Usage

``` r
birch
```

## Format

\`birch\` A list with 4 entries:

- sites:

  A matrix with two columns giving the longitude and latitude of each
  site.

- div:

  A numeric vector giving the allelic richness sampled at each site.

- mig:

  A square, asymmetric matrix with estimated gene flow rates for each
  pair of sites.

- fst:

  A square, symmetric matrix with Fst values for each pair of sites.

## Source

Y. Tsuda, V. Semerikov, F. Sebastiani, G. G. Vendramin, M. Lascoux,
Multispecies genetic structure and hybridization in the Betula genus
across Eurasia. Molecular Ecology 26, 589-605 (2017).
\<https://doi.org/10.1111/mec.13885\>
