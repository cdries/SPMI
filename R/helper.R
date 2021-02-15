
swapping_wrapper <- function(
  x, n_rows, d_cols, m_matrices, maxiter, maxnoimprove, eps
) {
  
  # call swapping
  out <- RRAswapping(x, maxiter, maxnoimprove, eps, n_rows, d_cols, m_matrices)
  
  # Add output
  out$Vopt <- out$V[1 + out$ksteps]
  out$V <- out$V[1:(1 + out$ksteps)]
  
  return (out)
}


blockswapping_wrapper <- function(
  x, n_rows, d_cols, m_matrices, maxiter, maxnoimprove, eps, p_vec
) {
  
  # call blockswapping
  out <- RRAblockswapping(
    x, maxiter, maxnoimprove, eps, n_rows, d_cols, m_matrices, p_vec
  )
  
  # Add output
  out$Vopt <- out$V[1 + out$ksteps]
  out$V <- out$V[1:(1 + out$ksteps)]
  
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
  out_preprocess <- RRAblockmincov(
    x, n_preprocess, n_preprocess, eps, n_rows, d_cols, m_matrices, rep(0.5, 3)
  )
  out <- RRAblockswapping(
    RRArearrange(
      x, out_preprocess$pimat - 1, n_rows, d_cols, m_matrices
    ), 
    maxiter - n_preprocess, maxnoimprove, eps, n_rows, d_cols, m_matrices, alpha
  )
  
  # Adapt output
  out$Vopt <- out$V[1 + out$ksteps]
  out$Vpreprop <- out_preprocess$V[1 + n_preprocess]
  out$V <- c(out_preprocess$V, out$V[2:(1 + out$ksteps)])
  out$ksteps <- out$ksteps + n_preprocess

  return (out)
}


random_wrapper <- function(x, n_rows, d_cols, m_matrices, maxiter, maxnoimprove, eps) {
  
  # call random algorithm
  out <- RRArandom(x, maxiter, maxnoimprove, eps, n_rows, d_cols, m_matrices)
  
  # Add output
  out$Vopt <- out$V[1 + out$ksteps]
  out$V <- out$V[1:(1 + out$ksteps)]
  
  return (out)
}
