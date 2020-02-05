#include "core/volume.h"

RCPP_EXPOSED_CLASS(volume)
RCPP_MODULE(volume_cpp) {
    Rcpp::class_<volume>("volume_cpp")
        .constructor<arma::uvec, double>()
        .constructor<arma::uvec, double, arma::vec>()
        .field_readonly("h", &volume::h, "Voxel side length")
        .method("set", &volume::set, "Set the state at a given index")
        .method("get", &volume::get, "Get the state at a given index")
        .property("dims", &volume::dims, "Volume dimensions")
        .property("xptr", &volume::xptr, "External pointer");
}
