#ifndef MULTSBM
#define MULTSBM

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"
using namespace Rcpp;

class MultSbm : public IclModel
{
public:
  MultSbm(const arma::cube& xp,S4 model,arma::vec& cl,bool verb=false);
  void set_cl(arma::vec clt);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  arma::mat delta_swap(int i,arma::uvec iclust);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
  void merge_update(int k, int l);
  List get_obs_stats();
protected:
  int M;
  double beta;
  // matrix of observed counts for each clusters
  arma::cube x;
  arma::cube x_counts;
  double cst; 
};

#endif
