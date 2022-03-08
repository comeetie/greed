#ifndef DCLBM
#define DCLBM

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include "IclModelEmission.h"
using namespace Rcpp;



class DcLbm : public IclModelEmission
{
public:
  DcLbm(const arma::sp_mat& xp,int Nr, int Nc,S4 model,bool verb=false);
  void set_cl(arma::uvec clt);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl,bool dead_cluster);
  arma::vec delta_swap(const int i,arma::uvec & cl, bool almost_dead_cluster, arma::uvec iclust, int K);
  void swap_update(const int i,arma::uvec &  cl,bool dead_cluster,const int newcl);
  double delta_merge(int k, int l);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
  void merge_update(const int k,const int l);
  List get_obs_stats();
  List get_obs_stats_cst();
protected:
  arma::sp_mat x;
  arma::sp_mat xt;
  arma::mat x_counts;
  arma::vec counts;
  double p;
  double alpha;
  double cst; 
  arma::vec dr;
  arma::vec dc;
  arma::uvec row_clusts;
  arma::uvec col_clusts;
  arma::vec clusttypes;
  int N;
  int K;
  int Kr;
  int Kc;
  int Nr;
  int Nc;
};

#endif

