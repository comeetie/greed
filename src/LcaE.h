#ifndef LCAE
#define LCAE

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModelEmission.h"
using namespace Rcpp;

class LcaE : public IclModelEmission
{
public:
  LcaE(arma::mat& x,S4 model,bool verb=false);
  void set_cl(arma::vec clt);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  arma::mat delta_swap(int i, int K,const arma::vec cl,arma::uvec iclust);
  void swap_update(int i,const arma::vec cl,bool dead_cluster, int newcl);
  double delta_merge(int k, int l);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){return 0;};
  void merge_update(int k, int l);
  List get_obs_stats();
protected:
  arma::mat x;
  arma::vec nbmod;
  arma::vec counts;
  std::vector<arma::mat> x_counts;
  double beta;
  int K;
};

#endif

