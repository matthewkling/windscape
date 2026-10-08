#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Wind roses are built from a cells x (2 * k) matrix `m` holding k time steps of wind data for
// every cell: k columns of u (eastward) components followed by k columns of v (northward).

// Wind speed at each cell and time step, as a cells x k matrix.
// [[Rcpp::export]]
NumericMatrix wind_speeds(NumericMatrix m, int k) {
      const int nc = m.nrow();
      NumericMatrix s(nc, k);
      for(int t = 0; t < k; ++t) {
            for(int c = 0; c < nc; ++c) {
                  const double u = m(c, t), v = m(c, k + t);
                  s(c, t) = std::sqrt(u * u + v * v);
            }
      }
      return s;
}

// Add each time step's wind to running per-cell sums of loadings toward the 8 neighbors, in
// place. The wind at each step is split between the two neighbors whose bearings bracket its
// direction, in proportion to how closely it points toward each, weighted by its transformed
// speed: `w` (a cells x k matrix), or speed ^ `p` when `w` has no rows. Calm steps (u and v
// both zero) have no direction and add nothing, whatever the weight (e.g. 0 ^ 0 = 1).
//
// nb: one row per grid row, holding the rhumb-line bearings (degrees) to the N, NE, E, SE, S, SW,
//     W, and NW neighbors, then 360.
// row: each cell's (0-based) row in `nb`.
// acc: cells x 8 running sums, in the same neighbor order as `nb`. A cell with a missing value
//      at any time step becomes NA in all 8.
// [[Rcpp::export]]
void rose_accumulate(NumericMatrix m, int k, NumericMatrix w, double p,
                     NumericMatrix nb, IntegerVector row, NumericMatrix acc) {
      const int nc = m.nrow();
      const bool use_w = w.nrow() > 0;
      for(int t = 0; t < k; ++t) {
            for(int c = 0; c < nc; ++c) {
                  const double u = m(c, t), v = m(c, k + t);
                  const double wt = use_w ? w(c, t) : std::pow(std::sqrt(u * u + v * v), p);
                  if(std::isnan(u) || std::isnan(v) || std::isnan(wt)) {
                        for(int j = 0; j < 8; ++j) acc(c, j) = NA_REAL;
                        continue;
                  }
                  if(u == 0.0 && v == 0.0) continue;   // calm

                  // compass bearing the wind blows toward, in (0, 360]
                  double b = std::atan2(v, -u) * 180.0 / M_PI - 90.0;
                  if(b < -180.0) b += 360.0;
                  if(b < 0.0) b += 360.0;
                  if(b == 0.0) b = 360.0;

                  // the bracketing pair of neighbors, and the share toward each
                  const int r = row[c];
                  int ni = 0;
                  for(int j = 0; j < 9; ++j) {
                        if(b > nb(r, j)) ni = j;
                  }
                  const double prop = (b - nb(r, ni)) / (nb(r, ni + 1) - nb(r, ni));
                  acc(c, ni) += (1.0 - prop) * wt;
                  acc(c, (ni + 1) % 8) += prop * wt;
            }
      }
}
