#ifndef NMARBCSBMUNDIRECTED
#define NMARBCSBMUNDIRECTED

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "NmarBcSbm.h"
using namespace Rcpp;

class NmarBcSbmUndirected : public NmarBcSbm 
{
public:
  NmarBcSbmUndirected(arma::sp_mat& xp,arma::sp_mat& xpobs,S4 model,arma::vec& cl,bool verb=false);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
};

#endif

