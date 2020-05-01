#include "core/volume.h"

RCPP_EXPOSED_CLASS(volume)
RCPP_MODULE(volume_cpp) {
    Rcpp::class_<core::volume>("volume_cpp")
        .constructor<arma::uvec, double>()
        .constructor<arma::uvec, double, arma::vec>()
        .field_readonly("h", &core::volume::h, "Voxel side length")
        .method("set", &core::volume::set, "Set the state at a given index")
        .method("get", &core::volume::get, "Get the state at a given index")
        .property("dims", &core::volume::dims, "Volume dimensions")
        .property("xptr", &core::volume::xptr, "External pointer");
}
