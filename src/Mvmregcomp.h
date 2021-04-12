#ifndef MREGCOMP
#define MREGCOMP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"

using namespace Rcpp;

class Mvmregcomp : public IclModel
{
public:
  Mvmregcomp(const arma::mat & X,const arma::mat & Y,double alpha,double beta, double N0,arma::vec& cl,bool verb=false);
  void set_cl(arma::vec clt);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  arma::mat delta_swap(int i,arma::uvec iclust);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  void merge_update(int k, int l);
  List get_obs_stats();
private:
  arma::mat X;
  arma::mat Y;
  List regs;
  arma::mat M;
  arma::mat Kp;
  arma::mat S0;
  double N0;
};

#endif

