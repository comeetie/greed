#ifndef MREG
#define MERG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"

using namespace Rcpp;

class Mvmreg : public IclModel
{
public:
  Mvmreg(const arma::mat & X,const arma::mat & Y,int K,double alpha,double beta, double N0,arma::vec& cl,bool verb=false);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  arma::mat delta_swap(int i);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  void merge_update(int k, int l);
  List get_obs_stats();
private:
  arma::mat X;
  arma::mat Y;
  List regs;
  double beta;
  double N0;
};

#endif

