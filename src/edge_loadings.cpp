
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector edge_loadings(NumericVector b, NumericVector w, NumericVector nb) {
      int n = b.size();
      NumericMatrix loadings(n, 9);
      int ni = 0;
      for(int i = 0; i < n; ++i) {
            for(int j = 0; j < 9; ++j) {
                  if(b[i] > nb[j]) {
                        ni = j;
                  }
            }
            double prop = (b[i] - nb[ni]) / (nb[ni+1] - nb[ni]);
            loadings(i, ni) = (1 - prop) * w[i];
            loadings(i, ni+1) = prop * w[i];
            loadings(i, 0) = loadings(i, 0) + loadings(i, 8);
      }
      loadings = loadings(_, Range(0, 7));
      NumericVector out(8);
      for(int j = 0; j < 8; ++j) {
            out[j] = sum(loadings(_, j));
      }
      return out;
}
