#ifndef MARSBMUNDIRECTED
#define MARSBMUNDIRECTED

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "IclModel.h"
#include "MissSbm.h"
using namespace Rcpp;

class MarSbmUndirected : public MissSbm
{
public:
  MarSbmUndirected(arma::sp_mat& xp,arma::sp_mat& xpobs,double alphai,double a0i,double b0i,double a0obsi,double b0obsi,arma::vec& clt,bool verb);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
private:
  double cst;
};

#endif

