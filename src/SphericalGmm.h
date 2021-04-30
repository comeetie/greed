#ifndef SPHERICALGMM
#define SPHERICALGMM

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"

using namespace Rcpp;

class SphericalGmm : public IclModel
{
public:
  SphericalGmm(const arma::mat & X,S4 model,arma::vec& cl,bool verb=false);
  void set_cl(arma::vec cli);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  arma::mat delta_swap(int i,arma::uvec iclust);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  void merge_update(int k, int l);
  List get_obs_stats();
private:
  arma::mat X;
  List regs;
  double tau;
  double kappa;
  double beta;
  arma::rowvec mu;
};

#endif

