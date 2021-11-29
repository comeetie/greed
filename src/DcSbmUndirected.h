#ifndef DCSBMUNDIRECTED
#define DCSBMUNDIRECTED

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "DcSbm.h"
using namespace Rcpp;

class DcSbmUndirected : public DcSbm
{
  using DcSbm::DcSbm;
public:
  void set_cl(arma::uvec clt);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl,bool dead_cluster);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
};

#endif
