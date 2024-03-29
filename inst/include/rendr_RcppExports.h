// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_rendr_RCPPEXPORTS_H_GEN_
#define RCPP_rendr_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace rendr {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("rendr", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("rendr", "_rendr_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in rendr");
            }
        }
    }

    inline Rcpp::DataFrame ssa_cpp(SEXP rnet_xptr, arma::vec y, double T, int length_out = 100, bool all_out = false, Rcpp::Nullable<arma::vec> k_vec = R_NilValue, Rcpp::Nullable<int> seed_val = R_NilValue) {
        typedef SEXP(*Ptr_ssa_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_ssa_cpp p_ssa_cpp = NULL;
        if (p_ssa_cpp == NULL) {
            validateSignature("Rcpp::DataFrame(*ssa_cpp)(SEXP,arma::vec,double,int,bool,Rcpp::Nullable<arma::vec>,Rcpp::Nullable<int>)");
            p_ssa_cpp = (Ptr_ssa_cpp)R_GetCCallable("rendr", "_rendr_ssa_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_ssa_cpp(Shield<SEXP>(Rcpp::wrap(rnet_xptr)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(T)), Shield<SEXP>(Rcpp::wrap(length_out)), Shield<SEXP>(Rcpp::wrap(all_out)), Shield<SEXP>(Rcpp::wrap(k_vec)), Shield<SEXP>(Rcpp::wrap(seed_val)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::DataFrame >(rcpp_result_gen);
    }

    inline Rcpp::List ssa_cpp_pest(SEXP rnet_xptr, arma::vec y, double T, int trajectories = 1, Rcpp::Nullable<arma::vec> k_vec = R_NilValue, Rcpp::Nullable<int> seed_val = R_NilValue) {
        typedef SEXP(*Ptr_ssa_cpp_pest)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_ssa_cpp_pest p_ssa_cpp_pest = NULL;
        if (p_ssa_cpp_pest == NULL) {
            validateSignature("Rcpp::List(*ssa_cpp_pest)(SEXP,arma::vec,double,int,Rcpp::Nullable<arma::vec>,Rcpp::Nullable<int>)");
            p_ssa_cpp_pest = (Ptr_ssa_cpp_pest)R_GetCCallable("rendr", "_rendr_ssa_cpp_pest");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_ssa_cpp_pest(Shield<SEXP>(Rcpp::wrap(rnet_xptr)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(T)), Shield<SEXP>(Rcpp::wrap(trajectories)), Shield<SEXP>(Rcpp::wrap(k_vec)), Shield<SEXP>(Rcpp::wrap(seed_val)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline Rcpp::List ssa_cpp_trajest(SEXP rnet_xptr, arma::vec y, double T, int trajectories = 1, int length_out = 100, Rcpp::Nullable<arma::vec> k_vec = R_NilValue) {
        typedef SEXP(*Ptr_ssa_cpp_trajest)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_ssa_cpp_trajest p_ssa_cpp_trajest = NULL;
        if (p_ssa_cpp_trajest == NULL) {
            validateSignature("Rcpp::List(*ssa_cpp_trajest)(SEXP,arma::vec,double,int,int,Rcpp::Nullable<arma::vec>)");
            p_ssa_cpp_trajest = (Ptr_ssa_cpp_trajest)R_GetCCallable("rendr", "_rendr_ssa_cpp_trajest");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_ssa_cpp_trajest(Shield<SEXP>(Rcpp::wrap(rnet_xptr)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(T)), Shield<SEXP>(Rcpp::wrap(trajectories)), Shield<SEXP>(Rcpp::wrap(length_out)), Shield<SEXP>(Rcpp::wrap(k_vec)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline Rcpp::List ssa_cpp_count(SEXP rnet_xptr, arma::vec y, double T, Rcpp::Nullable<arma::vec> k_vec = R_NilValue, Rcpp::Nullable<int> seed_val = R_NilValue) {
        typedef SEXP(*Ptr_ssa_cpp_count)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_ssa_cpp_count p_ssa_cpp_count = NULL;
        if (p_ssa_cpp_count == NULL) {
            validateSignature("Rcpp::List(*ssa_cpp_count)(SEXP,arma::vec,double,Rcpp::Nullable<arma::vec>,Rcpp::Nullable<int>)");
            p_ssa_cpp_count = (Ptr_ssa_cpp_count)R_GetCCallable("rendr", "_rendr_ssa_cpp_count");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_ssa_cpp_count(Shield<SEXP>(Rcpp::wrap(rnet_xptr)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(T)), Shield<SEXP>(Rcpp::wrap(k_vec)), Shield<SEXP>(Rcpp::wrap(seed_val)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline double prop_px(SEXP rnet_xptr, arma::vec x, int pi, int xi) {
        typedef SEXP(*Ptr_prop_px)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_prop_px p_prop_px = NULL;
        if (p_prop_px == NULL) {
            validateSignature("double(*prop_px)(SEXP,arma::vec,int,int)");
            p_prop_px = (Ptr_prop_px)R_GetCCallable("rendr", "_rendr_prop_px");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_prop_px(Shield<SEXP>(Rcpp::wrap(rnet_xptr)), Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(pi)), Shield<SEXP>(Rcpp::wrap(xi)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

}

#endif // RCPP_rendr_RCPPEXPORTS_H_GEN_
