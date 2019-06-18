#ifndef MM
#define MM

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"

using namespace Rcpp;

class Mm : public IclModel
{
public:
  Mm(arma::sp_mat& xp,int K,double alpha,double beta,arma::vec& cl,bool verb=false);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  arma::mat delta_swap(int i);
  void swap_update(int i, int newcl);
  MergeMat delta_merge();
  MergeMat delta_merge(arma::mat delta, int obk, int obl);
  void merge_update(int k, int l);
  List get_obs_stats();
private:
  arma::sp_mat x;
  arma::sp_mat xt;
  // matrix of observed counts for each clusters
  arma::mat x_counts;
  double beta;
  double norm_fact;
};

#endif

