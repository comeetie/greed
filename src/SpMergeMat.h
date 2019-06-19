#ifndef SPMERGEMAT
#define SPMERGEMAT

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;



class SpMergeMat
{
public :
  SpMergeMat(int bk, int bl, double bv, arma::sp_mat mat);
  int getK(){return k;};
  int getL(){return l;};
  arma::sp_mat getMergeMat(){return merge_mat;};
  double getValue(){return v;};
private :
  int k;
  int l;
  double v;
  arma::sp_mat merge_mat;
};

#endif

