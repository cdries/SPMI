#' Rearrange multiple matrices
#'
#' wrapper function for rearrangement of multiple matrices
#'
#' TODO
#'
#' @name rra
#' @encoding UTF-8
#' @concept rra
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
    out <- blockswapping_wrapper(
      x, n_rows, d_cols, m_matrices, maxiter, maxnoimprove, eps, p_block = p_block
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
