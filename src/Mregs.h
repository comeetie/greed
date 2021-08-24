#ifndef MREGS
#define MREGS

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include "IclModelEmission.h"
using namespace Rcpp;



class Mregs : public IclModelEmission
{
public:
  Mregs(const arma::mat & X,const arma::mat & Y,S4 model,bool verb=false);
  void set_cl(arma::uvec clt);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl,bool dead_cluster);
  arma::vec delta_swap(const int i,arma::uvec & cl, bool almost_dead_cluster, arma::uvec iclust, int K);
  void swap_update(const int i,arma::uvec &  cl,bool dead_cluster,const int newcl);
  double delta_merge(int k, int l);
  void merge_update(const int k,const int l);
  List get_obs_stats();
protected:
  arma::mat X;
  arma::mat Y;
  List regs;
  arma::mat M;
  arma::mat Kp;
  arma::mat S0;
  double N0;
  int N;
  int K;
};

#endif

