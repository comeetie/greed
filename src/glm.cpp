#include "glm.h"
#include "family.h"
using namespace Rcpp;
// Based on the algorithm on page 137:
// Arnold, T., Kane, M., & Lewis, B. W. (2019). A Computational Approach to Statistical Learning. CRC Press.
// taken and adpated from https://github.com/dirkschumacher/rcppglm for testing purpose
// [[Rcpp::export]]
List glm_fit_hot(const arma::mat& X, const arma::colvec& y, arma::colvec s, const double lambda= 0.01, int maxit=100, double tol=1e-6) {
  
  std::unique_ptr<Link::LinkFunction> ptr(new Link::Log());
  Family::Poisson * family_pt = new Family::Poisson(ptr);
  const Family::ExponentialFamily& family = *family_pt;
  
  const int n_cols = X.n_cols;
  const int n_rows = X.n_rows;

  
  arma::colvec s_old;
  arma::colvec eta = arma::ones<arma::colvec>(n_rows);
  arma::mat Lambdas = lambda*arma::diagmat(arma::vec(n_cols, arma::fill::ones)); 
  arma::mat H;
  arma::mat iH;
  int i;
  for (i = 0; i < maxit; i++) {
    s_old = s;
    arma::colvec eta = X*s;
    const arma::colvec mu = family.link_inverse(eta);
    const arma::colvec mu_p = family.link_mu_eta(eta);
    const arma::colvec z = eta + (y - mu) / mu_p;
    const arma::colvec W = arma::square(mu_p) / family.variance(mu);
    H  = X.t() * (X.each_col() % W)  + Lambdas;
    iH = arma::inv_sympd(H);
    s = iH * X.t() * ( W % z);
    const bool is_converged = std::sqrt(arma::accu(arma::square(s - s_old))) < tol;
    if (is_converged) break;
  }
  // compute LL
  // compute Laplace approx
  return List::create(Named("beta", s),Named("lambda", lambda),Named("H", H),Named("iH", iH),Named("N",n_rows),Named("nb_iter",i));
}
// [[Rcpp::export]]
List glm_fit(const arma::mat& X, const arma::colvec& y, const double lambda= 0.01, int maxit=100, double tol=1e-6){
  arma::colvec s = arma::zeros<arma::colvec>(X.n_cols);
  return glm_fit_hot(X,y,s,lambda,maxit,tol);
}

arma::colvec glm_log_lik(const arma::mat& X, const arma::colvec& y, List fit){
  const arma::colvec beta  = as<arma::colvec>(fit["beta"]); 

  std::unique_ptr<Link::LinkFunction> ptr(new Link::Log());
  Family::Poisson * family_pt = new Family::Poisson(ptr);
  const Family::ExponentialFamily& family = *family_pt;
  
  arma::colvec eta = X*beta;
  const arma::colvec mu = family.link_inverse(eta);
  
  return exp(mu%y);
  
  
}
// [[Rcpp::export]]
double delta_merge_post(List fit1,List fit2,double lambda){
  int D = as<int>(fit1["beta"]);
  arma::mat Lambdas = lambda*arma::diagmat(arma::vec(D, arma::fill::ones)); 
  arma::mat iLambdas = (1/lambda)*arma::diagmat(arma::vec(D, arma::fill::ones)); 
  const arma::mat H1 = as<arma::mat>(fit1["H"]);
  const arma::mat H2 = as<arma::mat>(fit2["H"]);
  const arma::mat iH1 = as<arma::mat>(fit1["iH"]);
  const arma::mat iH2 = as<arma::mat>(fit2["iH"]);
  const arma::colvec beta1 = as<arma::colvec>(fit1["beta"]);
  const arma::colvec beta2 = as<arma::colvec>(fit2["beta"]);
  const arma::mat H12 = H1+H2;
  const arma::mat iH12 = arma::inv_sympd(H12);
  const arma::colvec beta12 = iH12*(H1*beta1+H2*beta2);
  const arma::mat S = arma::inv_sympd(iH1+iH2);
  double qdets = -arma::log_det_sympd(iLambdas-iH12)+std::pow(1/lambda,D);
  double delta = qdets + log_mvn_pdf(beta1,beta2,S)-log_mvn_pdf(arma::colvec(D, arma::fill::zeros),beta12,iLambdas-iH12);
  return delta;
}

// [[Rcpp::export]]
double log_mvn_pdf(arma::colvec x,arma::colvec mu, arma::mat S){
  int D = x.n_rows; 
  double ll = -0.5*D*log(2*M_PI)-0.5*arma::log_det_sympd(S)-0.5*((x-mu).t()*arma::inv_sympd(S)*(x-mu)).eval()(0,0);
  return ll;
}

