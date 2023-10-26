#include <boost/random.hpp>
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

void draw_wor_sample(boost::random::mt19937& generator, boost::random::uniform_real_distribution<double>& distribution, const std::vector<double>& i_doub_array, std::vector<int>& x, int n, int s) {
  for (int i = 0; i < n; i ++) x[i] = i;
  double n_tot_doub = (double) n, u, temp;
  int pos;
  // perform the swap
  for (int i = 0; i < s; i ++) {
    u = distribution(generator);
    pos = floor((n_tot_doub - i_doub_array[i]) * u);
    temp = x[pos];
    x[pos] = x[n - i - 1];
    x[n - i - 1] = temp;
  }
  return;
}


//' @title Permute a Bernoulli treatment vector
//' @description This function permutes a Bernoulli treatment vector of length n with s entries equal to 1.
//' @param n the sample size
//' @param s the number of entries equal to 1
//' @param B the number of WOR samples to generate
//' @noRd
// [[Rcpp::export]]
SEXP fisher_yates_samlper(int n, int s, int B) {
   // initialize output vector, x vector shared across samples
   std::vector<std::vector<int>>* synth_idx_list = new std::vector<std::vector<int>>(B);
   std::vector<int> x(n);

   // initialize the random number generator
   boost::random::mt19937 generator(4);
   boost::random::uniform_real_distribution<double> distribution(0, 1);

   // initialize array of i doubles
   std::vector<double> i_doub_array(s);
   for (int i = 0; i < s; i ++) i_doub_array[i] = (double) i;

   // loop from 0 to B, generating WOR samples
   for (int j = 0; j < B; j ++) {
     std::vector<int> v(s);
     // perform the swap
     draw_wor_sample(generator, distribution, i_doub_array, x, n, s);
     // load the s WOR samples into v
     for (int i = n - s; i < n; i ++) v[i - n + s] = x[i];
     // store v within synth_idx_list
     (*synth_idx_list)[j] = v;
   }
   Rcpp::XPtr<std::vector<std::vector<int>>> ptr(synth_idx_list);
   return ptr;
}
