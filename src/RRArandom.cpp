#include "RcppArmadillo.h"
#include "RRAhelper.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List RRArandom(Rcpp::NumericVector xx, int nsteps, double eps, int N, int D, int M) {
  // Rearrangement of random matrices with a randomized algorithm
  // 
  // arguments:
  // xx       : R array of dimension 3 with the matrices to rearrange 
  //            over the 3rd dimension
  // nsteps   : maximum number of attempted swaps before termination
  // eps      : terminate if V < eps
  // N        : number of rows
  // D        : number of columns
  // M        : number of matrices
  //
  // output:
  // x        : rearranged array
  // V        : vector with objective values
  // pimat    : matrix with optimal permutation
  // rS       : matrix with rowsums per matrix
  //
  // author: Dries Cornilly
  
  // set data into cube format, more efficient that reading in using arma::cube x
  arma::cube x(xx.begin(), N, D, M, false);       // y as cube of xx input
  
  // set initial variables
  arma::vec V = -arma::ones(nsteps + 1);          // stores the objective value over the iterations
  List tmpFn = RRAfnV(x, N, M);                   // get rowsums and objective
  V(0) = tmpFn["V"];
  arma::umat pimat = arma::repelem(arma::regspace<arma::uvec>(0, N - 1), 1, D); // initial permutation
  
  // iterate over the possible swaps, stop when termination criteria are met
  int iter = 0;
  while(iter < nsteps) {
    
    arma::umat pimat_try = RRArandomperm(N, D);
    arma::cube x_new = RRArearrange(x, pimat_try, N, D, M);
    double Vtemp = RRAfnV(x_new, N, M)["V"];
    
    if (Vtemp < V(iter)) {
      V(iter + 1) = Vtemp;
      pimat = pimat_try;
    } else {
      V(iter + 1) = V(iter);
    }
    
    iter++;
  }
  arma::mat rS = RRAfnV(RRArearrange(x, pimat, N, D, M), N, M)["rS"];
  
  // output
  List out;
  out["x"] = x;
  out["V"] = V;
  out["pimat"] = pimat + 1;
  out["rS"] = rS;
  
  return out;
}
