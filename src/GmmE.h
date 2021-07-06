#ifndef GMME
#define GMME

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "IclModelEmission.h"
using namespace Rcpp;



class GmmE : public IclModelEmission
{
public:
  GmmE(const arma::mat & X,S4 model,bool verb=false);
  void set_cl(arma::vec clt);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  arma::mat delta_swap(const int i,int K, Partition clp, arma::uvec iclust);
  void swap_update(const int i,Partition clp,bool dead_cluster,const int newcl);
  double delta_merge(int k, int l);
  void merge_update(const int k,const int l);
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
  arma::vec counts;
  int N;
  int K;
};

#endif

