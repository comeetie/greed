#ifndef GAUSSMERGE
#define GAUSSMERGE

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


double log_prob_gauss(arma::vec x,arma::vec mu,arma::mat iS,arma::mat S);
List gauss_merge(List g1, List g2,arma::mat iSigma_prior);


#endif

