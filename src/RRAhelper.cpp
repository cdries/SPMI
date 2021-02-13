#include "RcppArmadillo.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


List RRAfnV(arma::cube x, int N, int M) {
  // x        : cube of dimension 3 for which the objective value has to be computed
  //
  // output:
  // V        : objective value
  // rS       : matrix with rowsums
  //
  // author: Dries Cornilly
  
  double n(N);
  double m(M);
  arma::mat rS = arma::sum(x, 1);                 // rowSums per matrix
  arma::mat S = arma::repmat(arma::sum(rS, 0) / n, N, 1); // matrix sums divided by number of rows
  double V = arma::sum(arma::sum(arma::square(rS - S))) / (n * m); // average variance of rowsums
  
  List out;
  out["V"] = V;
  out["rS"] = rS;
  
  return out;
}


// [[Rcpp::export]]
arma::cube RRArearrange(arma::cube x, arma::umat pimat, int N, int D, int M) {
  // x        : cube of dimension 3 to rearrange
  // pimat    : matrix with the permuatation
  // N        : number of rows
  // D        : number of columns
  // M        : number of matrices
  //
  // output:
  // y        : rearranged matrix
  //
  // author: Dries Cornilly
  
  arma::cube y(x);
  for (int ii = 0; ii < D; ii++) {
    arma::mat xt = x.subcube(0, ii, 0, N - 1, ii, M - 1);
    for (int kk = 0; kk < M; kk++) {
      arma::vec xtt = xt.col(kk);
      y.subcube(0, ii, kk, N - 1, ii, kk) = xtt(pimat.col(ii));
    }
  }
  
  return y;
}


inline int randWrapper(const int n) { return floor(unif_rand() * n); }


arma::umat RRArandomperm(int N, int D) {
  // N        : number of rows
  // D        : number of columns
  //
  // output:
  // pimat    : permutation matrix
  //
  // author: Dries Cornilly

  arma::umat pimat(N, D);
  
  for (int jj = 0; jj < D; jj++) {
    
    arma::uvec b = arma::regspace<arma::uvec>(0, N - 1);
    std::random_shuffle(b.begin(), b.end(), randWrapper);
    
    pimat.col(jj) = b;
  }
  
  return pimat;
}

