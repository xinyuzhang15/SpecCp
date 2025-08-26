#include <Rcpp.h>
using namespace Rcpp;

// A Rcpp function of truncate a vector by keeping its leading s elements (in absolute value).

// [[Rcpp::export]]
NumericVector mytruncate(NumericVector x, int s) {
  int n = x.size();
  NumericVector xtruncate(n, 0.0);
  std::vector<std::pair<double, int>> abs_x(n);

  for (int i = 0; i < n; ++i) {
    abs_x[i] = std::make_pair(std::abs(x[i]), i);
  }

  std::partial_sort(abs_x.begin(), abs_x.begin() + s, abs_x.end(), std::greater<std::pair<double, int>>());

  for (int i = 0; i < s; ++i) {
    xtruncate[abs_x[i].second] = x[abs_x[i].second];
  }

  return xtruncate;
}
