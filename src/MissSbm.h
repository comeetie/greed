#ifndef MISSSBM
#define MISSSBM

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"
using namespace Rcpp;

class MissSbm : public IclModel
{
public:
  MissSbm(arma::sp_mat& xp,arma::sp_mat& xpobs,double alpha,double a0,double b0,double a0obs,double b0obs,arma::vec& cl,bool verb=false);
  void set_cl(arma::vec clt);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  arma::mat delta_swap(int i,arma::uvec iclust);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
  MergeMat delta_merge(arma::mat delta, int obk, int obl, const List & old_stats);
  void merge_update(int k, int l);
  List get_obs_stats();
protected:
  arma::sp_mat x;
  arma::sp_mat xt;
  arma::sp_mat xobs;
  arma::sp_mat xtobs;
  // matrix of observed counts for each clusters
  arma::mat x_counts;
  arma::mat x_counts_obs;
  double a0;
  double b0;
  double a0obs;
  double b0obs;
};

#endif

