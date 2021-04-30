#ifndef CODCSBM
#define CODCSBM

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"
using namespace Rcpp;

class CoDcSbm : public IclModel
{
public:
  CoDcSbm(const arma::sp_mat& xp,int Nri, int Nci,S4 model,arma::vec& cl,bool verb=false);
  void set_cl(arma::vec clt);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  double icl(const List & obs_stats);
  double icl(const List & obs_stats,int oldcl,int newcl);
  arma::mat delta_swap(int i,arma::uvec iclust);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
  void merge_update(int k, int l);
  List get_obs_stats();
  List get_obs_stats_cst();
private:
  arma::sp_mat  x;
  // matrix of observed counts for each clusters
  arma::mat x_counts;
  double p;
  double cst; 
  arma::vec dr;
  arma::vec dc;
  arma::vec row_clusts;
  arma::vec col_clusts;
  arma::vec clusttypes;
  int Kr;
  int Kc;
  int Nr;
  int Nc;
};

#endif

