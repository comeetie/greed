#ifndef glm_H
#define glm_H

#include "family.h"
#include "linkfunctions.h"
using namespace Rcpp;
// taken from https://github.com/dirkschumacher/rcppglm for testing purpose
List glm_fit(const arma::mat& X, const arma::colvec& y, const double lambda, int maxit, double tol);
List glm_fit_hot(const arma::mat& X, const arma::colvec& y, arma::colvec s,const double lambda, int maxit, double tol);
arma::colvec glm_log_lik(const arma::mat& X, const arma::colvec& y, List fit);
double delta_merge_post(List fit1,List fit2,double lambda);
double log_mvn_pdf(arma::colvec x,arma::colvec mu, arma::mat S);
#endif