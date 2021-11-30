#ifndef SIMPLEICLCOMODEL
#define SIMPLEICLCOMODEL

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"
#include "IclModelEmission.h"

using namespace Rcpp;

class SimpleIclCoModel : public IclModel
{
public:
  SimpleIclCoModel(IclModelEmission * emission_modeli, S4 model,arma::uvec cli,int Nri, int Nci,bool verb=false);
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
  S4 get_model(){return emission_model->get_model();};
  List get_obs_stats();
  List get_obs_stats_cst(){return emission_model->get_obs_stats_cst();};
  ~SimpleIclCoModel(){
    delete emission_model;
  };
private:
  IclModelEmission * emission_model;
  arma::uvec row_clusts;
  arma::uvec col_clusts;
  arma::vec clusttypes;
  int Kr;
  int Kc;
  int Nr;
  int Nc;
};

#endif

