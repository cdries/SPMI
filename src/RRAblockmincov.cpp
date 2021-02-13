#include "RcppArmadillo.h"
#include "RRAhelper.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List RRAblockmincov(Rcpp::NumericVector xx, int nsteps, int maxnoimprove, 
                    double eps, int N, int D, int M, arma::vec psample) {
  // Rearrangement of random matrices with the improved block (minimum covariance) algorithm
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
  // psample  : vector with expected block size as percentage (in (0, 1)) of D
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

    // random column and row indices
    arma::vec randomcolumns = arma::randu(D);     // sample in [0, 1]
    arma::uvec jvec = randomcolumns < psample(iter); // which columns to use
    if (arma::sum(jvec) == 0) {                   // at least 1 column selected
      int j = arma::randi(1, arma::distr_param(0, D - 1))(0);
      jvec(j) = 1;
    } else if (arma::sum(jvec) == D) {            // at least 1 columns should not be selected
      int j = arma::randi(1, arma::distr_param(0, D - 1))(0);
      jvec(j) = 0;
    }
    int i1 = arma::randi(1, arma::distr_param(0, N - 1))(0);
    int i2 = arma::randi(1, arma::distr_param(0, N - 2))(0);
    if (i2 >= i1) i2++;                           // i2 always different than i1

    // construct variable to which to be anti-comonotonic
    arma::mat jcol = arma::zeros(N, M);
    for (int jj = 0; jj < D; jj++) {
      if (jvec(jj) == 1) {
        arma::mat tmp = x.tube(0, jj, N - 1, jj);
        jcol += tmp.rows(pimat.col(jj));
      }
    }
    arma::mat L = rS - jcol;

    // construct beta and finalise variable to which to be anti-comonotonic
    arma::vec beta = arma::ones(M);
    double tmpVar = arma::sum(arma::square(jcol.col(0) - arma::sum(jcol.col(0)) / n)) / n;
    for (int jj = 1; jj < M; jj++) {
      beta(jj) = (arma::sum((jcol.col(jj) - arma::sum(jcol.col(jj)) / n) %
        (jcol.col(0) - arma::sum(jcol.col(0)) / n)) / n) / tmpVar;
    }
    arma::vec bL = L * beta;

    // increasing order of bL and decreasing order of column j of first matrix
    arma::uvec orderbL = arma::sort_index(bL);
    arma::uvec orderjcol = arma::sort_index(jcol.col(0), "descend");

    // reorder as to make column j anti-comonotonic with bL
    rS = L;
    for (int jj = 0; jj < D; jj++) {
      if (jvec(jj) == 1) {
        arma::uvec pj = pimat.col(jj);
        pj(orderbL) = pj(orderjcol);
        pimat.col(jj) = pj;

        arma::mat jz = x.tube(0, jj, N - 1, jj);
        rS += jz.rows(pimat.col(jj));
      }
    }

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
