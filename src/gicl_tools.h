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

arma::mat update_count(arma::vec counts,int oldcl, int newcl);

arma::vec sum_cols(const arma::sp_mat& x);
arma::vec sum_cols(const arma::mat x);
arma::vec sum_rows(const arma::mat x);

double sum_lfact(const arma::sp_mat& x);

List lm_post(const arma::mat X,const arma::colvec& y,double regu, double a0, double b0);

List lm_post_add(List current, const arma::mat X,const arma::colvec& y,double regu, double a0, double b0);

List lm_post_del(List current, const arma::mat X,const arma::colvec& y,double regu, double a0, double b0);

List lm_post_add1(List current, const arma::rowvec X,double y,double regu, double a0, double b0);

List lm_post_del1(List current, const arma::rowvec X,double y,double regu, double a0, double b0);


List lm_post_merge(List current_k,List current_l,double regu, double a0, double b0);

#endif

