// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
using namespace Rcpp;



MergeMat::MergeMat(int bk, int bl, double bv, arma::mat mat){
    k=bk;
    l=bl;
    merge_mat=mat;
    v=bv;
};
