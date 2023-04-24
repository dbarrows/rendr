// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/rendr.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// issa_cpp
Rcpp::List issa_cpp(SEXP rnet_xptr, arma::vec D, SEXP volume_xptr, double T, int length_out, bool all_out, bool verbose, Rcpp::Nullable<arma::vec> k_vec);
RcppExport SEXP _rendr_issa_cpp(SEXP rnet_xptrSEXP, SEXP DSEXP, SEXP volume_xptrSEXP, SEXP TSEXP, SEXP length_outSEXP, SEXP all_outSEXP, SEXP verboseSEXP, SEXP k_vecSEXP) {
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
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type k_vec(k_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(issa_cpp(rnet_xptr, D, volume_xptr, T, length_out, all_out, verbose, k_vec));
    return rcpp_result_gen;
END_RCPP
}
// nsm_cpp
Rcpp::List nsm_cpp(SEXP rnet_xptr, arma::vec D, SEXP volume_xptr, double T, int length_out, bool all_out, bool verbose, Rcpp::Nullable<arma::vec> k_vec);
RcppExport SEXP _rendr_nsm_cpp(SEXP rnet_xptrSEXP, SEXP DSEXP, SEXP volume_xptrSEXP, SEXP TSEXP, SEXP length_outSEXP, SEXP all_outSEXP, SEXP verboseSEXP, SEXP k_vecSEXP) {
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
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type k_vec(k_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(nsm_cpp(rnet_xptr, D, volume_xptr, T, length_out, all_out, verbose, k_vec));
    return rcpp_result_gen;
END_RCPP
}
// nsm_cpp_pest
Rcpp::List nsm_cpp_pest(SEXP rnet_xptr, arma::vec D, SEXP volume_xptr, double T, int trajectories, Rcpp::Nullable<arma::vec> k_vec);
RcppExport SEXP _rendr_nsm_cpp_pest(SEXP rnet_xptrSEXP, SEXP DSEXP, SEXP volume_xptrSEXP, SEXP TSEXP, SEXP trajectoriesSEXP, SEXP k_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type D(DSEXP);
    Rcpp::traits::input_parameter< SEXP >::type volume_xptr(volume_xptrSEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type trajectories(trajectoriesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type k_vec(k_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(nsm_cpp_pest(rnet_xptr, D, volume_xptr, T, trajectories, k_vec));
    return rcpp_result_gen;
END_RCPP
}
// ssa_cpp
Rcpp::DataFrame ssa_cpp(SEXP rnet_xptr, arma::vec y, double T, int length_out, bool all_out, Rcpp::Nullable<arma::vec> k_vec, Rcpp::Nullable<int> seed_val);
static SEXP _rendr_ssa_cpp_try(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP length_outSEXP, SEXP all_outSEXP, SEXP k_vecSEXP, SEXP seed_valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type length_out(length_outSEXP);
    Rcpp::traits::input_parameter< bool >::type all_out(all_outSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type k_vec(k_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type seed_val(seed_valSEXP);
    rcpp_result_gen = Rcpp::wrap(ssa_cpp(rnet_xptr, y, T, length_out, all_out, k_vec, seed_val));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rendr_ssa_cpp(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP length_outSEXP, SEXP all_outSEXP, SEXP k_vecSEXP, SEXP seed_valSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rendr_ssa_cpp_try(rnet_xptrSEXP, ySEXP, TSEXP, length_outSEXP, all_outSEXP, k_vecSEXP, seed_valSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// ssa_cpp_pest
Rcpp::List ssa_cpp_pest(SEXP rnet_xptr, arma::vec y, double T, int trajectories, Rcpp::Nullable<arma::vec> k_vec, Rcpp::Nullable<int> seed_val);
static SEXP _rendr_ssa_cpp_pest_try(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP trajectoriesSEXP, SEXP k_vecSEXP, SEXP seed_valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type trajectories(trajectoriesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type k_vec(k_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type seed_val(seed_valSEXP);
    rcpp_result_gen = Rcpp::wrap(ssa_cpp_pest(rnet_xptr, y, T, trajectories, k_vec, seed_val));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rendr_ssa_cpp_pest(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP trajectoriesSEXP, SEXP k_vecSEXP, SEXP seed_valSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rendr_ssa_cpp_pest_try(rnet_xptrSEXP, ySEXP, TSEXP, trajectoriesSEXP, k_vecSEXP, seed_valSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// ssa_cpp_trajest
Rcpp::List ssa_cpp_trajest(SEXP rnet_xptr, arma::vec y, double T, int trajectories, int length_out, Rcpp::Nullable<arma::vec> k_vec);
static SEXP _rendr_ssa_cpp_trajest_try(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP trajectoriesSEXP, SEXP length_outSEXP, SEXP k_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type trajectories(trajectoriesSEXP);
    Rcpp::traits::input_parameter< int >::type length_out(length_outSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type k_vec(k_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(ssa_cpp_trajest(rnet_xptr, y, T, trajectories, length_out, k_vec));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rendr_ssa_cpp_trajest(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP trajectoriesSEXP, SEXP length_outSEXP, SEXP k_vecSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rendr_ssa_cpp_trajest_try(rnet_xptrSEXP, ySEXP, TSEXP, trajectoriesSEXP, length_outSEXP, k_vecSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// ssa_cpp_count
Rcpp::List ssa_cpp_count(SEXP rnet_xptr, arma::vec y, double T, Rcpp::Nullable<arma::vec> k_vec, Rcpp::Nullable<int> seed_val);
static SEXP _rendr_ssa_cpp_count_try(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP k_vecSEXP, SEXP seed_valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type k_vec(k_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type seed_val(seed_valSEXP);
    rcpp_result_gen = Rcpp::wrap(ssa_cpp_count(rnet_xptr, y, T, k_vec, seed_val));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rendr_ssa_cpp_count(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP k_vecSEXP, SEXP seed_valSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rendr_ssa_cpp_count_try(rnet_xptrSEXP, ySEXP, TSEXP, k_vecSEXP, seed_valSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// prop_px
double prop_px(SEXP rnet_xptr, arma::vec x, int pi, int xi);
static SEXP _rendr_prop_px_try(SEXP rnet_xptrSEXP, SEXP xSEXP, SEXP piSEXP, SEXP xiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type pi(piSEXP);
    Rcpp::traits::input_parameter< int >::type xi(xiSEXP);
    rcpp_result_gen = Rcpp::wrap(prop_px(rnet_xptr, x, pi, xi));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rendr_prop_px(SEXP rnet_xptrSEXP, SEXP xSEXP, SEXP piSEXP, SEXP xiSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rendr_prop_px_try(rnet_xptrSEXP, xSEXP, piSEXP, xiSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// tauleap_cpp
Rcpp::DataFrame tauleap_cpp(SEXP rnet_xptr, arma::vec y, double T, arma::mat hots, bool use_implicit, int length_out, bool all_out, Rcpp::Nullable<arma::vec> k_vec, bool verbose);
RcppExport SEXP _rendr_tauleap_cpp(SEXP rnet_xptrSEXP, SEXP ySEXP, SEXP TSEXP, SEXP hotsSEXP, SEXP use_implicitSEXP, SEXP length_outSEXP, SEXP all_outSEXP, SEXP k_vecSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rnet_xptr(rnet_xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type hots(hotsSEXP);
    Rcpp::traits::input_parameter< bool >::type use_implicit(use_implicitSEXP);
    Rcpp::traits::input_parameter< int >::type length_out(length_outSEXP);
    Rcpp::traits::input_parameter< bool >::type all_out(all_outSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type k_vec(k_vecSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(tauleap_cpp(rnet_xptr, y, T, hots, use_implicit, length_out, all_out, k_vec, verbose));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _rendr_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("Rcpp::DataFrame(*ssa_cpp)(SEXP,arma::vec,double,int,bool,Rcpp::Nullable<arma::vec>,Rcpp::Nullable<int>)");
        signatures.insert("Rcpp::List(*ssa_cpp_pest)(SEXP,arma::vec,double,int,Rcpp::Nullable<arma::vec>,Rcpp::Nullable<int>)");
        signatures.insert("Rcpp::List(*ssa_cpp_trajest)(SEXP,arma::vec,double,int,int,Rcpp::Nullable<arma::vec>)");
        signatures.insert("Rcpp::List(*ssa_cpp_count)(SEXP,arma::vec,double,Rcpp::Nullable<arma::vec>,Rcpp::Nullable<int>)");
        signatures.insert("double(*prop_px)(SEXP,arma::vec,int,int)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _rendr_RcppExport_registerCCallable() { 
    R_RegisterCCallable("rendr", "_rendr_ssa_cpp", (DL_FUNC)_rendr_ssa_cpp_try);
    R_RegisterCCallable("rendr", "_rendr_ssa_cpp_pest", (DL_FUNC)_rendr_ssa_cpp_pest_try);
    R_RegisterCCallable("rendr", "_rendr_ssa_cpp_trajest", (DL_FUNC)_rendr_ssa_cpp_trajest_try);
    R_RegisterCCallable("rendr", "_rendr_ssa_cpp_count", (DL_FUNC)_rendr_ssa_cpp_count_try);
    R_RegisterCCallable("rendr", "_rendr_prop_px", (DL_FUNC)_rendr_prop_px_try);
    R_RegisterCCallable("rendr", "_rendr_RcppExport_validate", (DL_FUNC)_rendr_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_rendr_issa_cpp", (DL_FUNC) &_rendr_issa_cpp, 8},
    {"_rendr_nsm_cpp", (DL_FUNC) &_rendr_nsm_cpp, 8},
    {"_rendr_nsm_cpp_pest", (DL_FUNC) &_rendr_nsm_cpp_pest, 6},
    {"_rendr_ssa_cpp", (DL_FUNC) &_rendr_ssa_cpp, 7},
    {"_rendr_ssa_cpp_pest", (DL_FUNC) &_rendr_ssa_cpp_pest, 6},
    {"_rendr_ssa_cpp_trajest", (DL_FUNC) &_rendr_ssa_cpp_trajest, 6},
    {"_rendr_ssa_cpp_count", (DL_FUNC) &_rendr_ssa_cpp_count, 5},
    {"_rendr_prop_px", (DL_FUNC) &_rendr_prop_px, 4},
    {"_rendr_tauleap_cpp", (DL_FUNC) &_rendr_tauleap_cpp, 9},
    {"_rendr_RcppExport_registerCCallable", (DL_FUNC) &_rendr_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_rendr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
