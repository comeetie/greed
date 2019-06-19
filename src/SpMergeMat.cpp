// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "SpMergeMat.h"
using namespace Rcpp;



SpMergeMat::SpMergeMat(int bk, int bl, double bv, arma::sp_mat mat){
    k=bk;
    l=bl;
    merge_mat=mat;
    v=bv;
};
