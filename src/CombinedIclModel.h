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
  CombinedIclModel(std::vector<IclModelEmission*> IclModelsi, S4 model,arma::uvec cli,bool verb=false);
  void set_cl(arma::uvec cli);
  double icl(const List & obs_stats);
  double icl(const List & obs_stats,int oldcl,int newcl);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  double icl_prop(arma::vec counts);
  double icl_prop(arma::vec counts,int oldcl,int newcl);
  arma::vec delta_prop_swap(int i,arma::uvec iclust);
  arma::vec delta_swap(int i,arma::uvec iclust);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
  double delta_merge_correction_prop(int k,int l,int obk,int obl,const List & old_stats);
  void merge_update(int k, int l);
  List get_obs_stats();
  S4 get_model();
  ~CombinedIclModel(){
    for(int m=0;m<IclModels.size();m++){
      delete IclModels[m];
    }
  };
private:
  std::vector<IclModelEmission*> IclModels;
  CharacterVector components_names;
};

#endif

