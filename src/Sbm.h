#ifndef SBM
#define SBM

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"
using namespace Rcpp;

class Sbm : public IclModel
{
public:
  Sbm(arma::sp_mat& xp,int K,double alpha,double a0,double b0, bool verb=false);
  Sbm(arma::sp_mat& xp,int K,double alpha,double a0,double b0,arma::vec& cl,bool verb=false);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  arma::mat delta_swap(int i,arma::uvec iclust);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  void merge_update(int k, int l);
  List get_obs_stats();
private:
  arma::sp_mat x;
  arma::sp_mat xt;
  // matrix of observed counts for each clusters
  arma::mat x_counts;
  double a0;
  double b0;
};

#endif

