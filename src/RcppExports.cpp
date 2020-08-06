// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// issa_cpp
Rcpp::List issa_cpp(SEXP rnet_xptr, arma::vec D, SEXP volume_xptr, double T, int length_out, bool all_out, bool verbose);
RcppExport SEXP _rendr_issa_cpp(SEXP rnet_xptrSEXP, SEXP DSEXP, SEXP volume_xptrSEXP, SEXP TSEXP, SEXP length_outSEXP, SEXP all_outSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type D(DSEXP);
    Rcpp::traits::input_parameter< SEXP >::type volume_xptr(volume_xptrSEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type length_out(length_outSEXP);
    Rcpp::traits::input_parameter< bool >::type all_out(all_outSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(issa_cpp(rnet_xptr, D, volume_xptr, T, length_out, all_out, verbose));
    return rcpp_result_gen;
END_RCPP
}
// nsm_cpp
Rcpp::List nsm_cpp(SEXP rnet_xptr, arma::vec D, SEXP volume_xptr, double T, int length_out, bool all_out, bool verbose);
RcppExport SEXP _rendr_nsm_cpp(SEXP rnet_xptrSEXP, SEXP DSEXP, SEXP volume_xptrSEXP, SEXP TSEXP, SEXP length_outSEXP, SEXP all_outSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type D(DSEXP);
    Rcpp::traits::input_parameter< SEXP >::type volume_xptr(volume_xptrSEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type length_out(length_outSEXP);
    Rcpp::traits::input_parameter< bool >::type all_out(all_outSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(nsm_cpp(rnet_xptr, D, volume_xptr, T, length_out, all_out, verbose));
    return rcpp_result_gen;
END_RCPP
}
// ssa_cpp
Rcpp::DataFrame ssa_cpp(SEXP rnet_xptr, arma::vec y, double T, int length_out, bool all_out, Rcpp::Nullable<arma::vec> k_vec);
RcppExport SEXP _rendr_ssa_cpp(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP length_outSEXP, SEXP all_outSEXP, SEXP k_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type length_out(length_outSEXP);
    Rcpp::traits::input_parameter< bool >::type all_out(all_outSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type k_vec(k_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(ssa_cpp(rnet_xptr, y, T, length_out, all_out, k_vec));
    return rcpp_result_gen;
END_RCPP
}
// ssa_cpp_pest
Rcpp::NumericVector ssa_cpp_pest(SEXP rnet_xptr, arma::vec y, double T, int trajectories, Rcpp::Nullable<arma::vec> k_vec);
RcppExport SEXP _rendr_ssa_cpp_pest(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP trajectoriesSEXP, SEXP k_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type trajectories(trajectoriesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type k_vec(k_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(ssa_cpp_pest(rnet_xptr, y, T, trajectories, k_vec));
    return rcpp_result_gen;
END_RCPP
}
// tauleap_cpp
Rcpp::DataFrame tauleap_cpp(SEXP rnet_xptr, arma::vec y, double T, arma::vec hors, arma::vec hots, arma::vec reverse, int length_out, bool all_out, Rcpp::Nullable<arma::vec> k_vec, bool verbose);
RcppExport SEXP _rendr_tauleap_cpp(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP horsSEXP, SEXP hotsSEXP, SEXP reverseSEXP, SEXP length_outSEXP, SEXP all_outSEXP, SEXP k_vecSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hors(horsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hots(hotsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type reverse(reverseSEXP);
    Rcpp::traits::input_parameter< int >::type length_out(length_outSEXP);
    Rcpp::traits::input_parameter< bool >::type all_out(all_outSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type k_vec(k_vecSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(tauleap_cpp(rnet_xptr, y, T, hors, hots, reverse, length_out, all_out, k_vec, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rendr_issa_cpp", (DL_FUNC) &_rendr_issa_cpp, 7},
    {"_rendr_nsm_cpp", (DL_FUNC) &_rendr_nsm_cpp, 7},
    {"_rendr_ssa_cpp", (DL_FUNC) &_rendr_ssa_cpp, 6},
    {"_rendr_ssa_cpp_pest", (DL_FUNC) &_rendr_ssa_cpp_pest, 5},
    {"_rendr_tauleap_cpp", (DL_FUNC) &_rendr_tauleap_cpp, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_rendr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
