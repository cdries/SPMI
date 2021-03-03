#' Simulation scenario used in the paper
#'
#' simulate multiple matrices
#'
#' @name simulate
#' @encoding UTF-8
#' @concept simulate
#' @param n number of rows
#' @param d number of columns
#' @param m number of matrices
#' @param rho correlation parameter 
#' @param normalise whether to ensure a V=0 solution exists
#' @param tdist whether to use marginals that are t-distributed or Gaussian
#' @author Dries Cornilly
#' @references
#' Cornilly, D., Puccetti, G., Rüschendorf, L., & Vanduffel, S. (2021). 
#' On a synchronization problem with multiple instances.
#'
#' @importFrom stats rnorm
#' @importFrom stats pnorm
#' @importFrom stats qt
#' @export simulate
simulate <- function(n, d = 2, m = 2, rho = 0.5, normalise=FALSE, tdist=FALSE) {
  
  # cross-matrix covariance matrix
  Sigma <- matrix(rho, nrow = m, ncol = m)
  diag(Sigma) <- 1
  
  B <- chol(Sigma)
  
  x <- array(NA, dim = c(n, d, m))
  
  for (ii in 1:d) {
    y <- matrix(stats::rnorm(n * m), ncol = m) %*% B
    x[, ii,] <- y
  }
  
  if (normalise) {
    # make sure V=0 solution exists
    for (ii in 1:m) {
      x[,, ii] <- x[,, ii] - rowMeans(x[,, ii])
    }
    
    # re-randomise
    for (ii in 1:d) {
      x[, ii,] <- x[sample(n), ii,]
    }
  }
  
  if (tdist) {
    nu <- 5
    x <- stats::qt(stats::pnorm(x), nu)
  }
  
  return (x)
}
