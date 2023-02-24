//This file is copied and modified from the version of code.cpp published in the package "qad" (Version 1.0.4) on CRAN
//by Thimo Kasper, Florian Griessenberger, Robert R. Junker, Valentin Petzel and Wolfgang Trutschnig
//The original file can be found under https://cran.r-project.org/web/packages/qad/index.html and https://github.com/griefl/qad
//The original file is distributed under the license GPL-2
//The changes are (20.01.2023):
//-An adapted version of build_checkerboard_weights which considers different resolutions in asymmetric_checkerboard_mass
//-An adapted version of build_checkerboard_weights which calculates a single entry A_kl instead of the entire matrix A in asymmetric_checkerboard_index

#include <Rcpp.h>
#include <cmath>
#include <cstdio>
using namespace Rcpp;

namespace std
{
template<>
struct hash<pair<double, double> >
{
  size_t operator()(pair<double, double> const& p) const
  {
    auto hash1 = hash<double>{}(p.first);
    auto hash2 = hash<double>{}(p.second);
    return hash1*1000000000 + hash2;
  }
};
}


std::unordered_map<std::pair<double, double>, R_xlen_t>* pair_range (const NumericVector& x, const NumericVector& y, const NumericVector& x_range, const NumericVector& y_range) {
  std::unordered_map<std::pair<double, double>, R_xlen_t>* result = new std::unordered_map<std::pair<double, double>, R_xlen_t>();
  for (R_xlen_t i = 0; i < x.length(); i++) {
    if(i%100000 == 0) Rcpp::checkUserInterrupt();
    if (x_range[(R_xlen_t)x[i]-1] > 1 && y_range[(R_xlen_t)y[i]-1] > 1) {
      std::pair<double, double> key(x[i], y[i]);
      if ( result->find(key) == result->end() ) {
        (*result)[key] = 0;
      }
      (*result)[key] = (*result)[key] + 1;
    }
  }


  return result;
}


NumericVector range(const NumericVector& x) {
  NumericVector result(x.length());

  for (R_xlen_t i = 0; i < x.length(); i++) {
    if(i%100000 == 0) Rcpp::checkUserInterrupt();
    result[(R_xlen_t) x[i] - 1] = result[(R_xlen_t) x[i] - 1] + 1;
  }

  return result;
}



// [[Rcpp::export]]
NumericMatrix asymmetric_checkerboard_mass (const NumericVector& X, const NumericVector& Y, R_xlen_t resolution1, R_xlen_t resolution2) {
  R_xlen_t sample_size = std::min(X.length(), Y.length());
  R_xlen_t x_upper;
  R_xlen_t x_lower;
  R_xlen_t y_upper;
  R_xlen_t y_lower;
  NumericMatrix result(resolution1, resolution2);
  NumericVector x_range = range(X);
  NumericVector y_range = range(Y);
  R_xlen_t rx;
  R_xlen_t ry;
  R_xlen_t rp;
  double lambda_x;
  double lambda_y;
  std::unordered_map<std::pair<double, double>, R_xlen_t>* p_range = pair_range(X, Y, x_range, y_range);

  for (R_xlen_t i = 0; i < sample_size; i++) {
    if(i%100000 == 0) Rcpp::checkUserInterrupt();
    rx = x_range[(R_xlen_t)X[i] - 1];
    ry = y_range[(R_xlen_t)Y[i] - 1];
    if (rx > 1 && ry > 1) {
      std::pair<double, double> key(X[i], Y[i]);
      rp = (*p_range)[key];
      (*p_range)[key] = 0;
    } else {
      rp = 1;
    }

    if (rp != 0) {

      x_upper = std::ceil((X[i]) / (double)sample_size * resolution1);
      x_lower = std::max(std::ceil((X[i] - (double)rx) / (double)sample_size * resolution1), 1.);
      y_upper = std::ceil((Y[i]) / (double)sample_size * resolution2);
      y_lower = std::max(std::ceil((Y[i] - (double)ry) / (double)sample_size * resolution2), 1.);

      for (R_xlen_t x = x_lower; x <= x_upper; x++) {
        lambda_x = std::min(X[i], (double)x / (double)resolution1 * sample_size) - std::max(X[i] - (double)rx, (double)(x - 1) / (double)resolution1 * sample_size);
        for (R_xlen_t y = y_lower; y <= y_upper; y++) {
          lambda_y = std::min(Y[i], (double)y / (double)resolution2 * sample_size) - std::max(Y[i] - (double)ry, (double)(y - 1) / (double)resolution2 * sample_size);
          result(x-1, y-1) = result(x-1, y-1) + (lambda_x * lambda_y * rp)/(double)(sample_size * rx * ry);
        }
      }

    }
  }

  delete p_range;

  return result;
}


// [[Rcpp::export]]
double asymmetric_checkerboard_index (const NumericVector& X, const NumericVector& Y, R_xlen_t k, R_xlen_t l, R_xlen_t resolution1, R_xlen_t resolution2) {
   R_xlen_t sample_size = std::min(X.length(), Y.length());
   R_xlen_t x_upper;
   R_xlen_t x_lower;
   R_xlen_t y_upper;
   R_xlen_t y_lower;
   double result = 0;
   NumericVector x_range = range(X);
   NumericVector y_range = range(Y);
   R_xlen_t rx;
   R_xlen_t ry;
   R_xlen_t rp;
   double lambda_x;
   double lambda_y;
   std::unordered_map<std::pair<double, double>, R_xlen_t>* p_range = pair_range(X, Y, x_range, y_range);

   for (R_xlen_t i = 0; i < sample_size; i++) {
     if(i%100000 == 0) Rcpp::checkUserInterrupt();
     rx = x_range[(R_xlen_t)X[i] - 1];
     ry = y_range[(R_xlen_t)Y[i] - 1];
     if (rx > 1 && ry > 1) {
       std::pair<double, double> key(X[i], Y[i]);
       rp = (*p_range)[key];
       (*p_range)[key] = 0;
     } else {
       rp = 1;
     }

     if (rp != 0) {

       x_upper = std::ceil((X[i]) / (double)sample_size * resolution1);
       x_lower = std::max(std::ceil((X[i] - (double)rx) / (double)sample_size * resolution1), 1.);
       y_upper = std::ceil((Y[i]) / (double)sample_size * resolution2);
       y_lower = std::max(std::ceil((Y[i] - (double)ry) / (double)sample_size * resolution2), 1.);

       if(x_lower > k || x_upper < k || y_lower > l || y_upper < l) {
         continue;
       }

       lambda_x = std::min(X[i], (double)k / (double)resolution1 * sample_size) - std::max(X[i] - (double)rx, (double)(k - 1) / (double)resolution1 * sample_size);
       lambda_y = std::min(Y[i], (double)l / (double)resolution2 * sample_size) - std::max(Y[i] - (double)ry, (double)(l - 1) / (double)resolution2 * sample_size);
       result += (lambda_x * lambda_y * rp)/(double)(sample_size * rx * ry);
     }
   }
   delete p_range;
   return result;
}
