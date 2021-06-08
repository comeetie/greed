#ifndef COMBINEDICLMODEL
#define COMBINEDICLMODEL

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"
#include "IclModelEmission.h"

using namespace Rcpp;

class CombinedIclModel : public IclModel
{
public:
  CombinedIclModel(std::vector<IclModelEmission*> IclModelsi, S4 model,arma::vec cli,bool verb=false);
  void set_cl(arma::vec cli);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  double icl_prop(arma::vec counts);
  double icl_prop(arma::vec counts,int oldcl,int newcl);
  arma::mat delta_prop_swap(int i,arma::uvec iclust);
  arma::mat delta_swap(int i,arma::uvec iclust);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
  void merge_update(int k, int l);
  List get_obs_stats();
private:
  std::vector<IclModelEmission*> IclModels;
};

#endif

