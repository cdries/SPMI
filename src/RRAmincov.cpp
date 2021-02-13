#include "RcppArmadillo.h"
#include "RRAhelper.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List RRAmincov(Rcpp::NumericVector xx, arma::mat beta, int nsteps, 
               int maxnoimprove, double eps, int N, int D, int M) {
  // Rearrangement of random matrices with the improved (minimum covariance) algorithm
  // 
  // arguments:
  // xx       : R array of dimension 3 with the matrices to rearrange 
  //            over the 3rd dimension
  // beta     : matrix (dimension M x D) containing beta values, 1 row per matrix
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
  arma::mat S = arma::repmat(arma::sum(rS, 0) / n, N, 1);
  arma::umat pimat = arma::repelem(arma::regspace<arma::uvec>(0, N - 1), 1, D); // initial permutation
  
  // iterate over the possible swaps, stop when termination criteria are met
  int iter = 0;
  int iternoimprove = 0;
  bool converged = false;
  int status = 2;
  while (iter < nsteps && !converged) {
    
    // random column index
    int j = arma::randi(1, arma::distr_param(0, D - 1))(0);
    arma::uvec pi_j = pimat.col(j);
    
    // construct variable to which to be anti-comonotonic
    arma::mat jcol = x.tube(0, j, N - 1, j);
    jcol = jcol.rows(pi_j);
    arma::mat L = rS - jcol;
    arma::vec bL = L * beta.col(j);
    
    // increasing order of bL and decreasing order of column j of first matrix
    arma::uvec orderbL = arma::sort_index(bL);
    arma::uvec orderjcol = arma::sort_index(jcol.col(0), "descend");
    
    // // reorder as to make column j anti-comonotonic with bL
    arma::uvec pj = pimat.col(j);
    pj(orderbL) = pi_j(orderjcol);
    pimat.col(j) = pj;
    
    arma::mat jz = x.tube(0, j, N - 1, j);
    rS = L + jz.rows(pimat.col(j));
    
    V(iter + 1) = arma::sum(arma::sum(arma::square(rS - S))) / (n * m);
    if (V(iter + 1) < V(iter)) {
      iternoimprove = 0;
    } else {
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
