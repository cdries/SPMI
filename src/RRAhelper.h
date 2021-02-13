#ifndef RRAHELPER_H
#define RRAHELPER_H

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


List RRAfnV(arma::cube x, int N, int M);

arma::cube RRArearrange(arma::cube x, arma::umat pimat, int N, int D, int M);

arma::umat RRArandomperm(int N, int D);

#endif
