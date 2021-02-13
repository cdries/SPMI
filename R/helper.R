
swapping_wrapper <- function(
  x, n_rows, d_cols, m_matrices, maxiter, maxnoimprove, eps
  ) {
  
  # call swapping
  out <- SPMI:::RRAswapping(x, maxiter, maxnoimprove, eps, n_rows, d_cols, m_matrices)
  
  return (out)
}


blockswapping_wrapper <- function(
  x, n_rows, d_cols, m_matrices, maxiter, maxnoimprove, eps, p_block = 0.5
  ) {
  
  # call blockswapping
  out <- SPMI:::RRAblockswapping(
    x, maxiter, maxnoimprove, eps, n_rows, d_cols, m_matrices, p_block
    )
  
  return (out)
}


custom_wrapper <- function(
  x, n_rows, d_cols, m_matrices, maxiter, maxnoimprove, eps, n_preprocess = 3, p_vec = NULL
  ) {
  
  # Default percentage vec
  if (is.null(p_vec)) {
    n1 <- round(n_rows * log(5 + d_cols))
    alpha <- c(rep(0.5, n1), rep(0.25, n1), rep(0.125, n1), rep(1 / (2 * d_cols), maxiter))
  } else {
    alpha <- p_vec
  }
  
  # call custom algorithm
  out_custom1 <- SPMI:::RRAblockmincov(
    x, n_preprocess, n_preprocess, eps, n_rows, d_cols, m_matrices, rep(0.5, 3)
    )
  out <- SPMI:::RRAblockswapping(
    SPMI:::RRArearrange(
      x, out_custom1$pimat - 1, n_rows, d_cols, m_matrices), 
    niter - n_preprocess, maxnoimprove, eps, n_rows, d_cols, m_matrices, alpha
    )
  
  # Adapt output
  # TODO
  
  return (out)
}


random_wrapper <- function(x, n_rows, d_cols, m_matrices, maxiter, maxnoimprove, eps) {
  
  # call random algorithm
  out <- SPMI:::RRArandom(x, maxiter, maxnoimprove, eps, n_rows, d_cols, m_matrices)
  # TODO: adapt C++ code for maxnoimprove usage
  
  return (out)
}
