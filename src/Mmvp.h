#ifndef MMVP
#define MMVP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"

using namespace Rcpp;

class Mmvp : public IclModel
{
public:
  Mmvp(arma::sp_mat& xp,int K,double alpha);
  Mmvp(arma::sp_mat& xp,int K,double alpha,arma::vec& cl);
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
  double beta;
  arma::vec p;
};

#endif

