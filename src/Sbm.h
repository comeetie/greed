#ifndef SBM
#define SBM

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include "IclModelEmission.h"
using namespace Rcpp;



class Sbm : public IclModelEmission
{
public:
  Sbm(const arma::sp_mat & xp,S4 model,bool verb=false);
  void set_cl(arma::uvec clt);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl,bool dead_cluster);
  arma::vec delta_swap(const int i,arma::uvec & cl, bool almost_dead_cluster, arma::uvec iclust, int K);
  void swap_update(const int i,arma::uvec &  cl,bool dead_cluster,const int newcl);
  double delta_merge(int k, int l);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
  void merge_update(const int k,const int l);
  List get_obs_stats();
protected:
  arma::sp_mat  x;
  arma::sp_mat  xt;
  // matrix of observed counts for each clusters
  arma::mat x_counts;
  double cst; 
  arma::vec counts;
  int K;
  int N;
  double a0;
  double b0;
};

#endif

