#ifndef MREG
#define MERG

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"

using namespace Rcpp;

class Mreg : public IclModel
{
public:
  Mreg(const arma::mat & X,const arma::colvec & y, int K,double alpha,double reg, double a0, double b0,bool verb=false);
  Mreg(const arma::mat & X,const arma::colvec & y,int K,double alpha,double reg, double a0, double b0,arma::vec& cl,bool verb=false);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  arma::mat delta_swap(int i,arma::uvec iclust);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  void merge_update(int k, int l);
  List get_obs_stats();
private:
  arma::mat X;
  arma::colvec y;
  List regs;
  double reg;
  double a0;
  double b0;
};

#endif

