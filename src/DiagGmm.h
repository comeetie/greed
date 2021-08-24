#ifndef DIAGGMM
#define DIAGGMM

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "IclModelEmission.h"
using namespace Rcpp;



class DiagGmm : public IclModelEmission
{
public:
  DiagGmm(const arma::mat & X,S4 model,bool verb=false);
  void set_cl(arma::uvec clt);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl,bool dead_cluster);
  arma::vec delta_swap(const int i,arma::uvec & cl, bool almost_dead_cluster, arma::uvec iclust, int K);
  void swap_update(const int i,arma::uvec &  cl,bool dead_cluster,const int newcl);
  double delta_merge(int k, int l);
  void merge_update(const int k,const int l);
  List get_obs_stats();
private:
  arma::mat X;
  double tau;
  double kappa;
  double beta;
  arma::rowvec mu;
  int K;
  List regs;
};

#endif

