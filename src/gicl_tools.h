#ifndef GICL_TOOLS
#define GICL_TOOLS

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat submatcross(int oldcl,int newcl,int K);

arma::vec count(arma::vec cl,int K);

arma::mat gsum_mat(arma::vec cl,const arma::sp_mat& x, int K);

arma::sp_mat gsum_mm(arma::vec cl,const arma::sp_mat& x, int K);

arma::sp_mat gsum_col(arma::vec cl,const arma::sp_mat& x, int i, int K);

arma::mat update_count(arma::vec counts,int oldcl, int newcl);

arma::sp_mat delcol(const arma::sp_mat & a, int ci);
arma::sp_mat delrowcol(const arma::sp_mat & a, int ci);
arma::sp_mat add_sppat(const arma::sp_mat & a, const arma::sp_mat & b);



List lm_post(const arma::mat X,const arma::colvec& y,double regu, double a0, double b0);

List lm_post_add(List current, const arma::mat X,const arma::colvec& y,double regu, double a0, double b0);

List lm_post_del(List current, const arma::mat X,const arma::colvec& y,double regu, double a0, double b0);

List lm_post_add1(List current, const arma::rowvec X,double y,double regu, double a0, double b0);

List lm_post_del1(List current, const arma::rowvec X,double y,double regu, double a0, double b0);

List lm_post_merge(List current_k,List current_l,double regu, double a0, double b0);

List mvlm_post(const arma::mat X,const arma::mat Y,double alpha, double N0);

List mvlm_post_add1(List current, const arma::rowvec X,const arma::rowvec Y,double alpha, double N0);

List mvlm_post_del1(List current, const arma::rowvec X,const arma::rowvec Y,double alpha, double N0);

List mvlm_post_merge(List current1,List current2,double alpha, double N0);

arma::uvec possible_moves(int k,arma::sp_mat & move_mat);
#endif

