#ifndef SIMPLEICLMODEL
#define SIMPLEICLMODEL

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"
#include "IclModelEmission.h"

using namespace Rcpp;

class SimpleIclModel : public IclModel
{
public:
  SimpleIclModel(IclModelEmission * emission_modeli, S4 model,arma::uvec cli,bool verb=false);
  void set_cl(arma::uvec cli);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  double icl_prop(arma::vec counts);
  double icl_prop(arma::vec counts,int oldcl,int newcl);
  arma::vec delta_prop_swap(int i,arma::uvec iclust);
  arma::vec delta_swap(int i,arma::uvec iclust);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
  void merge_update(int k, int l);
  S4 get_model(){return emission_model->get_model();};
  List get_obs_stats();
  List get_obs_stats_cst(){return emission_model->get_obs_stats_cst();};
  ~SimpleIclModel(){
    delete emission_model;
  };
private:
  IclModelEmission * emission_model;
};

#endif

