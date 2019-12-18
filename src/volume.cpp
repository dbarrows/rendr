#include "core/volume.h"

RCPP_EXPOSED_CLASS(volume)
RCPP_MODULE(volume_cpp) {
    Rcpp::class_<volume>("volume_cpp")
        .constructor<uvec>()
        .constructor<uvec, vec>()
        .method("set", &volume::set, "Set the state at a given index")
        .method("get", &volume::get, "Get the state at a given index")
        .property("dims", &volume::dims, "Volume dimensions")
        .property("xptr", &volume::xptr, "External pointer");
}
