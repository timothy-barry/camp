#include <Rcpp.h>
using namespace Rcpp;
#include <algorithm>

// side: -1 is left, 0 is both, 1 is right
// [[Rcpp::export]]
double compute_empirical_p_value(const std::vector<double>& null_statistics, double z_orig, int side) {
  double p;
  if (side == -1 || side == 1) { // left or right
    double counter = 0;
    double B = (double) null_statistics.size();
    if (side == -1) { // left
      for (int i = 0; i < null_statistics.size(); i ++) if (z_orig >= null_statistics[i]) counter ++;
    } else { // right
      for (int i = 0; i < null_statistics.size(); i ++) if (z_orig <= null_statistics[i]) counter ++;
    }
    p = (1.0 + counter)/(1.0 + B);
  } else { // two-sided
    p = std::min(1.0, 2 * std::min(compute_empirical_p_value(null_statistics, z_orig, -1),
                                   compute_empirical_p_value(null_statistics, z_orig, 1)));
  }
  return p;
}
