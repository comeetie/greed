#ifndef MERGEMAT
#define MERGEMAT

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;



class MergeMat
{
public :
  MergeMat(int bk, int bl, double bv, arma::mat mat);
  int getK(){return k;};
  int getL(){return l;};
  arma::mat getMergeMat(){return merge_mat;};
  double getValue(){return v;};
private :
  int k;
  int l;
  double v;
  arma::mat merge_mat;
};

#endif

