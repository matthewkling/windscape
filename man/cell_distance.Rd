% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cell_distance.R
\name{cell_distance}
\alias{cell_distance}
\title{Pairwise distances between cell centroids}
\usage{
cell_distance(x, ll)
}
\arguments{
\item{x}{SpatRaster}

\item{ll}{two-column matrix of site coordinates}
}
\value{
Pairwise distances between the cell centroids, in km
}
\description{
Calculate the pairwise distances between the centroids of the raster grid cells that a set of points fall into.
This is useful as an alternative to the distances between the points themselves, if distances are going to be
compared to other metrics that are based on cell centers (e.g. wind connectivity).
}
