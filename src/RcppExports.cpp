// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// RRAblockmincov
List RRAblockmincov(Rcpp::NumericVector xx, int nsteps, int maxnoimprove, double eps, int N, int D, int M, arma::vec psample);
RcppExport SEXP _SPMI_RRAblockmincov(SEXP xxSEXP, SEXP nstepsSEXP, SEXP maxnoimproveSEXP, SEXP epsSEXP, SEXP NSEXP, SEXP DSEXP, SEXP MSEXP, SEXP psampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< int >::type maxnoimprove(maxnoimproveSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type psample(psampleSEXP);
    rcpp_result_gen = Rcpp::wrap(RRAblockmincov(xx, nsteps, maxnoimprove, eps, N, D, M, psample));
    return rcpp_result_gen;
END_RCPP
}
// RRAblockswapping
List RRAblockswapping(Rcpp::NumericVector xx, int nsteps, int maxnoimprove, double eps, int N, int D, int M, arma::vec psample);
RcppExport SEXP _SPMI_RRAblockswapping(SEXP xxSEXP, SEXP nstepsSEXP, SEXP maxnoimproveSEXP, SEXP epsSEXP, SEXP NSEXP, SEXP DSEXP, SEXP MSEXP, SEXP psampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< int >::type maxnoimprove(maxnoimproveSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type psample(psampleSEXP);
    rcpp_result_gen = Rcpp::wrap(RRAblockswapping(xx, nsteps, maxnoimprove, eps, N, D, M, psample));
    return rcpp_result_gen;
END_RCPP
}
// RRArearrange
arma::cube RRArearrange(arma::cube x, arma::umat pimat, int N, int D, int M);
RcppExport SEXP _SPMI_RRArearrange(SEXP xSEXP, SEXP pimatSEXP, SEXP NSEXP, SEXP DSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type pimat(pimatSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(RRArearrange(x, pimat, N, D, M));
    return rcpp_result_gen;
END_RCPP
}
// RRAmincov
List RRAmincov(Rcpp::NumericVector xx, arma::mat beta, int nsteps, int maxnoimprove, double eps, int N, int D, int M);
RcppExport SEXP _SPMI_RRAmincov(SEXP xxSEXP, SEXP betaSEXP, SEXP nstepsSEXP, SEXP maxnoimproveSEXP, SEXP epsSEXP, SEXP NSEXP, SEXP DSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< int >::type maxnoimprove(maxnoimproveSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(RRAmincov(xx, beta, nsteps, maxnoimprove, eps, N, D, M));
    return rcpp_result_gen;
END_RCPP
}
// RRArandom
List RRArandom(Rcpp::NumericVector xx, int nsteps, int maxnoimprove, double eps, int N, int D, int M);
RcppExport SEXP _SPMI_RRArandom(SEXP xxSEXP, SEXP nstepsSEXP, SEXP maxnoimproveSEXP, SEXP epsSEXP, SEXP NSEXP, SEXP DSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< int >::type maxnoimprove(maxnoimproveSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(RRArandom(xx, nsteps, maxnoimprove, eps, N, D, M));
    return rcpp_result_gen;
END_RCPP
}
// RRAswapping
List RRAswapping(Rcpp::NumericVector xx, int nsteps, int maxnoimprove, double eps, int N, int D, int M);
RcppExport SEXP _SPMI_RRAswapping(SEXP xxSEXP, SEXP nstepsSEXP, SEXP maxnoimproveSEXP, SEXP epsSEXP, SEXP NSEXP, SEXP DSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< int >::type maxnoimprove(maxnoimproveSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(RRAswapping(xx, nsteps, maxnoimprove, eps, N, D, M));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SPMI_RRAblockmincov", (DL_FUNC) &_SPMI_RRAblockmincov, 8},
    {"_SPMI_RRAblockswapping", (DL_FUNC) &_SPMI_RRAblockswapping, 8},
    {"_SPMI_RRArearrange", (DL_FUNC) &_SPMI_RRArearrange, 5},
    {"_SPMI_RRAmincov", (DL_FUNC) &_SPMI_RRAmincov, 8},
    {"_SPMI_RRArandom", (DL_FUNC) &_SPMI_RRArandom, 7},
    {"_SPMI_RRAswapping", (DL_FUNC) &_SPMI_RRAswapping, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_SPMI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
