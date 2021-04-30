#ifndef MARSBM
#define MARSBM

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"
#include "MissSbm.h"
using namespace Rcpp;

class MarSbm : public MissSbm
{
public:
  MarSbm(arma::sp_mat& xp,arma::sp_mat& xpobs,S4 model,arma::vec& clt,bool verb);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
private:
  double cst;
};

#endif

