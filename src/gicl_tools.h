#ifndef GICL_TOOLS
#define GICL_TOOLS

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


arma::vec count(const arma::uvec & cl,int K);
arma::vec update_count(const arma::vec  & counts,int oldcl, int newcl);
arma::uvec to_zero_based(const arma::uvec & cl);
arma::umat table_count(arma::uvec cl, arma::uvec x, int K, int nbmod);
arma::mat gsum_mat(arma::uvec cl,const arma::sp_mat& x, int K);
arma::umat submatcross(int oldcl,int newcl,int K);
arma::sp_mat gsum_col(arma::uvec cl,const arma::sp_mat& x, int i, int K);
arma::cube gsum_cube(arma::uvec cl,const arma::cube& x, int K);
arma::sp_mat gsum_mm(arma::uvec cl,const arma::sp_mat& x, int K);
arma::mat gsum_bimat(arma::uvec clr,arma::uvec clc, const arma::sp_mat& x,int K);


arma::sp_mat sp_cross(arma::sp_mat colvec,arma::sp_mat rowvec,int self, int oldcl, int newcl, int K);






arma::sp_mat gsum_mat_sp(arma::vec cl,const arma::sp_mat& x, int K);




arma::sp_mat delcol(const arma::sp_mat & a, int ci);
void delrowcol(arma::sp_mat & a, int ci);
arma::sp_mat delrowcol_copy(const arma::sp_mat & a, int ci);
arma::sp_mat add_sppat(const arma::sp_mat & a, const arma::sp_mat & b);
arma::sp_mat add_spmatpat(const arma::sp_mat & a, const arma::sp_mat & b);
arma::sp_mat which_spmatpat(const arma::sp_mat & a, const arma::sp_mat & b);

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
List mvlm_post_comp(const arma::mat X,const arma::mat Y,const arma::mat K,const arma::mat M, const arma::mat S0, double N0);
List mvlm_post_add1_comp(List current, const arma::rowvec X,const arma::rowvec Y,const arma::mat K,const arma::mat M, const arma::mat S0, double N0);
List mvlm_post_del1_comp(List current, const arma::rowvec X,const arma::rowvec Y,const arma::mat K,const arma::mat M, const arma::mat S0, double N0);
List mvlm_post_merge_comp(List current1,List current2,const arma::mat K,const arma::mat M, const arma::mat S0, double N0);

List gmm_marginal(const arma::mat X,double tau,int N0, const arma::mat epsilon, const arma::rowvec mu);
List gmm_marginal_add1(List current, const arma::rowvec X,double tau,int N0, const arma::mat epsilon, const arma::rowvec mu);
List gmm_marginal_del1(List current, const arma::rowvec X,double tau,int N0, const arma::mat epsilon, const arma::rowvec mu);
List gmm_marginal_merge(List current1, List current2,double tau,int N0, const arma::mat epsilon, const arma::rowvec mu);

List gmm_marginal_spherical(const arma::mat X,double kappa,double tau,double beta, const arma::rowvec mu);
List gmm_marginal_spherical_add1(List current, const arma::rowvec X,double kappa,double tau,double beta, const arma::rowvec mu);
List gmm_marginal_spherical_del1(List current, const arma::rowvec X,double kappa,double tau,double beta, const arma::rowvec mu);
List gmm_marginal_spherical_merge(List current1, List current2,double kappa,double tau,double beta, const arma::rowvec mu);

arma::uvec possible_moves(int k,arma::sp_mat & move_mat);


#endif

