#ifndef MM
#define MM

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "IclModel.h"

using namespace Rcpp;

class Mm : public IclModel
{
public:
  Mm(arma::sp_mat& xp,S4 model,arma::vec& cl,bool verb=false);
  void set_cl(arma::vec cl);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  arma::mat delta_swap(int i,arma::uvec iclust);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  void merge_update(int k, int l);
  List get_obs_stats();
private:
  arma::sp_mat xt;
  // matrix of observed counts for each clusters
  arma::sp_mat x_counts;
  // col_sums of x_counts
  arma::rowvec col_sums;  
  double beta;
  double norm_fact;
};

#endif

