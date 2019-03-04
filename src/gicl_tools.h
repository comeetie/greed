#ifndef GICL_TOOLS
#define GICL_TOOLS

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat submatcross(int oldcl,int newcl,int K);

arma::vec count(arma::vec cl,int K);

arma::mat gsum_mat(arma::vec cl,const arma::sp_mat& x, int K);

arma::mat gsum_mm(arma::vec cl,const arma::sp_mat& x, int K);

arma::mat gsum_col(arma::vec cl,const arma::sp_mat& x, int i, int K);

arma::mat update_count(arma::mat counts,int oldcl, int newcl);

arma::vec sum_cols(const arma::sp_mat& x);

#endif

