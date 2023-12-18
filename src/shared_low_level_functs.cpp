// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
using namespace Rcpp;
#include <boost/math/distributions/skew_normal.hpp>
using boost::math::skew_normal;
#include <cmath>
#include <math.h>
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
    p = 2 * std::min(compute_empirical_p_value(null_statistics, z_orig, -1),
                     compute_empirical_p_value(null_statistics, z_orig, 1));
  }
  return p;
}


// [[Rcpp::export]]
std::vector<double> fit_skew_normal_funct(const std::vector<double>& y) {
  // initialize variables
  int n = y.size();
  double MAX_GAMMA_1 = 0.995;
  double n_doub = (double) n, s_1 = 0, s_2 = 0, s_3 = 0;

  // compute mean and standard deviation
  for (int i = 0; i < n; i ++) {
    s_1 += y[i];
    s_2 += y[i] * y[i];
  }
  double m_y = s_1/n_doub;
  double sd_y = sqrt(s_2/n_doub - m_y * m_y);

  // compute gamma1
  for (int i = 0; i < n; i ++) {
    s_3 += pow(y[i] - m_y, 3);
  }
  double gamma1 = s_3/(n_doub * pow(sd_y, 3));
  if (gamma1 > MAX_GAMMA_1) gamma1 = 0.9 * MAX_GAMMA_1;

  // reparameterize solution
  double b = sqrt(2.0/M_PI);
  double r = std::copysign(1.0, gamma1) * pow(2 * std::abs(gamma1)/(4 - M_PI), 1.0/3.0);
  double delta = r/(b * sqrt(1 + r * r));
  double alpha = delta/sqrt(1 - delta * delta);
  double mu_z = b * delta;
  double sd_z = sqrt(1 - mu_z * mu_z);
  double omega = sd_y/sd_z;
  double xi = m_y - omega * mu_z;

  return std::vector<double> {xi, omega, alpha, m_y, sd_y};
}


// [[Rcpp::export]]
std::vector<double> fit_and_evaluate_skew_normal(double z_orig, std::vector<double>& null_statistics, int side_code) {
  // 0. define variables
  double p = -1.0;
  bool finite_params = true;

  // 1. fit the skew normal
  std::vector<double> fitted_params = fit_skew_normal_funct(null_statistics);

  // 2. verify that the fitted parameters are finite
  for (int i = 0; i < fitted_params.size(); i ++) {
    if (!std::isfinite(fitted_params[i])) finite_params = false;
  }

  if (finite_params) {
    skew_normal dist(fitted_params[0], fitted_params[1], fitted_params[2]);
    if (side_code == 0) { // two-tailed
      p = 2.0 * std::min(cdf(complement(dist, z_orig)), cdf(dist, z_orig));
    } else if (side_code == 1) { // right-tailed
      p = cdf(complement(dist, z_orig));
    } else { // left-tailed
      p = cdf(dist, z_orig);
    }
    if (p <= 1.0e-250) p = 1.0e-250;
  }
  // 9. return p, xi, omega, alpha
  return std::vector<double> {fitted_params[0], fitted_params[1], fitted_params[2], p};
}
