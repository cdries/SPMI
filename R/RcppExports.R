# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

RRAblockmincov <- function(xx, nsteps, maxnoimprove, eps, N, D, M, psample) {
    .Call('_SPMI_RRAblockmincov', PACKAGE = 'SPMI', xx, nsteps, maxnoimprove, eps, N, D, M, psample)
}

RRAblockswapping <- function(xx, nsteps, maxnoimprove, eps, N, D, M, psample) {
    .Call('_SPMI_RRAblockswapping', PACKAGE = 'SPMI', xx, nsteps, maxnoimprove, eps, N, D, M, psample)
}

RRArearrange <- function(x, pimat, N, D, M) {
    .Call('_SPMI_RRArearrange', PACKAGE = 'SPMI', x, pimat, N, D, M)
}

RRAmincov <- function(xx, beta, nsteps, maxnoimprove, eps, N, D, M) {
    .Call('_SPMI_RRAmincov', PACKAGE = 'SPMI', xx, beta, nsteps, maxnoimprove, eps, N, D, M)
}

RRArandom <- function(xx, nsteps, maxnoimprove, eps, N, D, M) {
    .Call('_SPMI_RRArandom', PACKAGE = 'SPMI', xx, nsteps, maxnoimprove, eps, N, D, M)
}

RRAswapping <- function(xx, nsteps, maxnoimprove, eps, N, D, M) {
    .Call('_SPMI_RRAswapping', PACKAGE = 'SPMI', xx, nsteps, maxnoimprove, eps, N, D, M)
}

