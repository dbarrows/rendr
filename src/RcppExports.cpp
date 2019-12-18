// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// nsm_cpp
Rcpp::List nsm_cpp(SEXP rnet_xptr, SEXP volume_xptr, arma::vec D, arma::vec tspan);
RcppExport SEXP _reactor_nsm_cpp(SEXP rnet_xptrSEXP, SEXP volume_xptrSEXP, SEXP DSEXP, SEXP tspanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type volume_xptr(volume_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tspan(tspanSEXP);
    rcpp_result_gen = Rcpp::wrap(nsm_cpp(rnet_xptr, volume_xptr, D, tspan));
    return rcpp_result_gen;
END_RCPP
}
// ssa_cpp
Rcpp::DataFrame ssa_cpp(SEXP network_xptr, arma::vec y, arma::vec tspan);
RcppExport SEXP _reactor_ssa_cpp(SEXP network_xptrSEXP, SEXP ySEXP, SEXP tspanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type network_xptr(network_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tspan(tspanSEXP);
    rcpp_result_gen = Rcpp::wrap(ssa_cpp(network_xptr, y, tspan));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_volume_cpp();

static const R_CallMethodDef CallEntries[] = {
    {"_reactor_nsm_cpp", (DL_FUNC) &_reactor_nsm_cpp, 4},
    {"_reactor_ssa_cpp", (DL_FUNC) &_reactor_ssa_cpp, 3},
    {"_rcpp_module_boot_volume_cpp", (DL_FUNC) &_rcpp_module_boot_volume_cpp, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_reactor(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
