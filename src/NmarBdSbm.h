#ifndef NMARBDSBM
#define NMARBDSBM

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"
#include "MissSbm.h"
using namespace Rcpp;

class NmarBdSbm : public MissSbm
{
using MissSbm::MissSbm;
public:
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
};

#endif

