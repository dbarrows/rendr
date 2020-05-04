// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// issa_cpp
Rcpp::List issa_cpp(SEXP rnet_xptr, arma::vec D, SEXP volume_xptr, double T, bool verbose);
RcppExport SEXP _rendr_issa_cpp(SEXP rnet_xptrSEXP, SEXP DSEXP, SEXP volume_xptrSEXP, SEXP TSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type D(DSEXP);
    Rcpp::traits::input_parameter< SEXP >::type volume_xptr(volume_xptrSEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(issa_cpp(rnet_xptr, D, volume_xptr, T, verbose));
    return rcpp_result_gen;
END_RCPP
}
// nsm_cpp
Rcpp::List nsm_cpp(SEXP rnet_xptr, arma::vec D, SEXP volume_xptr, double T, bool verbose);
RcppExport SEXP _rendr_nsm_cpp(SEXP rnet_xptrSEXP, SEXP DSEXP, SEXP volume_xptrSEXP, SEXP TSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type D(DSEXP);
    Rcpp::traits::input_parameter< SEXP >::type volume_xptr(volume_xptrSEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(nsm_cpp(rnet_xptr, D, volume_xptr, T, verbose));
    return rcpp_result_gen;
END_RCPP
}
// ssa_cpp
Rcpp::DataFrame ssa_cpp(SEXP rnet_xptr, arma::vec y, double T, Rcpp::Nullable<arma::vec> k_vec, bool record_all);
RcppExport SEXP _rendr_ssa_cpp(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP k_vecSEXP, SEXP record_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type k_vec(k_vecSEXP);
    Rcpp::traits::input_parameter< bool >::type record_all(record_allSEXP);
    rcpp_result_gen = Rcpp::wrap(ssa_cpp(rnet_xptr, y, T, k_vec, record_all));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rendr_issa_cpp", (DL_FUNC) &_rendr_issa_cpp, 5},
    {"_rendr_nsm_cpp", (DL_FUNC) &_rendr_nsm_cpp, 5},
    {"_rendr_ssa_cpp", (DL_FUNC) &_rendr_ssa_cpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_rendr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
