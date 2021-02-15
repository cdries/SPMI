#' Rearrange multiple matrices
#'
#' wrapper function for rearrangement of multiple matrices
#'
#' 
#' There are currently four algorithms implemented: 1. swap 2. blockswap 3. custom and
#' 4. random. Each algorithm tries to minimize the average variance of the row-wise sums
#' across the matrices.
#' Control parameters include 'maxiter' for maximum number of iterations (default 1e5) and
#' 'eps' the tolerance to stop when V < eps (default 1e-6).
#'
#' @name rra
#' @encoding UTF-8
#' @concept rra
#' @param x Initial matrices to rearrange, of dimension (n, d, m)
#' @param algo Algorithm, one of ('swap', 'blockswap', 'custom', 'random')
#' @param maxiter Maximum number of iterations to run the algorithm
#' @param maxnoimprove If no immprovement after this many iterations, stop early
#' @param eps Stop if the function value is lower than this
#' @param p_block Percentage of rows in a block, either a numeric or a vector
#' @param n_preprocess Number of pre-processing steps for the custom algorithm
#' @param p_vec Vector with block-sizes (as percentage of total number of rows) for
#' the custom algorithm
#' @author Dries Cornilly
#' @references
#' Cornilly, D., Puccetti, G., Rüschendorf, L., & Vanduffel, S. (2021). 
#' On a synchronization problem with multiple instances.
#'
#' @import Rcpp
#' @useDynLib SPMI
#' @export rra
rra <- function(x, algo='swap', maxiter=1e5, maxnoimprove=1e3, eps=1e-6, 
                p_block=0.5, n_preprocess=3, p_vec=NULL) {
  
  # initialize properties
  n_rows <- dim(x)[1]
  d_cols <- dim(x)[2]
  m_matrices <- dim(x)[3]
  
  # call the requested algorithm
  if (algo == 'swap') {
    out <- swapping_wrapper(x, n_rows, d_cols, m_matrices, maxiter, maxnoimprove, eps)
  } else if (algo == 'blockswap') {
    if (length(p_block) < maxiter) {
      p_block <- rep(p_block, maxiter)
    }
    out <- blockswapping_wrapper(
      x, n_rows, d_cols, m_matrices, maxiter, maxnoimprove, eps, p_vec = p_block
    )
  } else if (algo == 'custom') {
    out <- custom_wrapper(
      x, n_rows, d_cols, m_matrices, maxiter, maxnoimprove, eps, 
      n_preprocess=n_preprocess, p_vec=p_vec
    )
  } else if (algo == 'random') {
    out <- random_wrapper(x, n_rows, d_cols, m_matrices, maxiter, maxnoimprove, eps)
  } else {
    warning('Chosen algorithm not implemented.')
  }
  
  return (out) 
}
