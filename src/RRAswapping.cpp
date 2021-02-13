#include "RcppArmadillo.h"
#include "RRAhelper.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List RRAswapping(Rcpp::NumericVector xx, int nsteps, int maxnoimprove, double eps,
                 int N, int D, int M) {
  // Rearrangement of random matrices with the swapping algorithm
  // 
  // arguments:
  // xx       : R array of dimension 3 with the matrices to rearrange 
  //            over the 3rd dimension
  // nsteps   : maximum number of attempted swaps before termination
  // maxnoimprove : maximum number of tried swaps without improvement
  //                before termination
  // eps      : terminate if V < eps
  // N        : number of rows
  // D        : number of columns
  // M        : number of matrices
  //
  // output:
  // x        : rearranged array
  // V        : vector with objective values
  // pimat    : matrix with optimal permutation
  // ksteps   : number of steps
  // status   : 0 (V < eps); 1 (no improvement in maxnoimprove steps); 2 (nsteps reached)
  //
  // author: Dries Cornilly
  
  // set data into cube format, more efficient that reading in using arma::cube x
  double n(N);
  double m(M);
  arma::cube x(xx.begin(), N, D, M, false);       // y as cube of xx input

  // set initial variables
  arma::vec V = -arma::ones(nsteps + 1);          // stores the objective value over the iterations
  List tmpFn = RRAfnV(x, N, M);                   // get rowsums and objective
  V(0) = tmpFn["V"];
  arma::mat rS = tmpFn["rS"];
  arma::mat pimat = arma::repelem(arma::regspace(0, N - 1), 1, D); // initial permutation
  
  // iterate over the possible swaps, stop when termination criteria are met
  int iter = 0;
  int iternoimprove = 0;
  bool converged = false;
  int status = 2;
  while (iter < nsteps && !converged) {
    
    // random column and row indices
    int j = arma::randi(1, arma::distr_param(0, D - 1))(0);
    int i1 = arma::randi(1, arma::distr_param(0, N - 1))(0);
    int i2 = arma::randi(1, arma::distr_param(0, N - 2))(0);
    if (i2 >= i1) i2++;                           // i2 always different than i1
    int pi_i1 = pimat(i1, j);
    int pi_i2 = pimat(i2, j);
    
    // improvement criterion
    arma::rowvec dk = x.tube(pi_i2, j) - x.tube(pi_i1, j);
    double ic = arma::sum(arma::square(rS.row(i1) + dk) + arma::square(rS.row(i2) - dk) -
                          arma::square(rS.row(i1)) - arma::square(rS.row(i2))) / (n * m);
    
    // update matrix and rowsums if there is improvement; update objective value
    if (ic < 0) {
      pimat(i1, j) = pi_i2;
      pimat(i2, j) = pi_i1;
      rS.row(i1) += dk;
      rS.row(i2) -= dk;
      
      V(iter + 1) = V(iter) + ic;
      iternoimprove = 0;
      
    } else {
      V(iter + 1) = V(iter);
      iternoimprove++;
    }
    
    // update convergence criterion
    if (V(iter + 1) < eps) {
      converged = true;
      status = 0;
    }
    if (iternoimprove >= maxnoimprove) {
      converged = true;
      status = 1;
    }
    
    iter++;
  }
  
  // output
  List out;
  out["x"] = x;
  out["V"] = V;
  out["pimat"] = pimat + 1;
  out["ksteps"] = iter;
  out["status"] = status;
  out["rS"] = rS;
  
  return out;
}
