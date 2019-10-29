#ifndef GMM
#define GMM

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"

using namespace Rcpp;

class Gmm : public IclModel
{
public:
  Gmm(const arma::mat & X,double alpha,double tau,int N0,arma::mat epsilon, arma::rowvec mu,arma::vec& cl,bool verb=false);
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
  arma::mat S;
  double normfact;
  List regs;
  double tau;
  int N0;
  arma::mat epsilon;
  arma::rowvec mu;
};

#endif

