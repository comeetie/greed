// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double log_prob_gauss(arma::vec x,arma::vec mu,arma::mat iS){
  double H = arma::as_scalar((x-mu).t()*iS*(x-mu));
  int d = iS.n_cols;
  return d*0.5*log(det(iS)) -d*0.5*log(2*M_PI) -0.5*H;
}

// [[Rcpp::export]]
double log_evi_gauss(arma::vec mu,arma::mat iS){
  int d = iS.n_cols;
  return -d*0.5*log(det(iS)) +d*0.5*log(2*M_PI);
}


// [[Rcpp::export]]
List GaussMerge(List current_k,List current_l,arma::mat iSigma_prior) {
  int n = as<int>(current_k["n"])+as<int>(current_l["n"]);

  arma::mat iS1 = as<arma::mat>(current_k["iS"]);
  arma::mat iS2 = as<arma::mat>(current_l["iS"]);
  int d = iS1.n_cols;

  arma::mat iS = iS1 + iS2 - iSigma_prior;
  arma::mat S = inv_sympd(iS);

  arma::vec mu1 = as<arma::vec>(current_k["mu"]);
  arma::vec mu2 = as<arma::vec>(current_l["mu"]);

  mu2 = inv_sympd(iS2-iSigma_prior)*iS2*mu2;

  arma::mat mu = S*(iS1*mu1+(iS2-iSigma_prior)*mu2);
  
  double le1 = d*0.5*log(det(iS1)) -d*0.5*log(2*M_PI);
  double le2 = d*0.5*log(det(iS2-iSigma_prior)) -d*0.5*log(2*M_PI);
  double lemerge =d*0.5*log(det(iS)) -d*0.5*log(2*M_PI);
  double log_evidence = lemerge-le1-le2;
  return List::create(Named("S")  = S,
                      Named("mu") = mu,
                      Named("iS") = iS,
                      Named("n") = n,
                      Named("log_evidence")=log_evidence);
  
}

