#ifndef SBMUNDIRECTED
#define SBMUNDIRECTED

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"
#include "Sbm.h"
using namespace Rcpp;

class SbmUndirected : public Sbm
{
  using Sbm::Sbm;
public:
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl,bool dead_cluster);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
};

#endif
