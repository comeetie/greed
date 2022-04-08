
#include <RcppArmadillo.h>
#include "family.h"
// taken from https://github.com/dirkschumacher/rcppglm for testing purpose
namespace Family {

arma::colvec Family::Binomial::log_lik(const arma::colvec& eta,const arma::colvec& y) const {
  return y;
}
arma::colvec Family::Binomial::log_lik_norm(const arma::colvec& eta,const arma::colvec& y) const {
  return y;
}

arma::colvec Family::Binomial::variance(const arma::colvec& mu) const {
  return mu % (1.0 - mu);
}

arma::colvec Family::Poisson::log_lik(const arma::colvec& eta,const arma::colvec& y) const {
  return eta % y - exp(eta);
}
arma::colvec Family::Poisson::log_lik_norm(const arma::colvec& eta,const arma::colvec& y) const {
  return eta % y - exp(eta)-lgamma(y+1);
}

arma::colvec Family::Poisson::variance(const arma::colvec& mu) const {
  return mu;
}

arma::colvec Family::Gamma::log_lik(const arma::colvec& eta,const arma::colvec& y) const {
  return eta % y - exp(eta);
}
arma::colvec Family::Gamma::log_lik_norm(const arma::colvec& eta,const arma::colvec& y) const {
  return eta % y - exp(eta)-lgamma(y+1);
}

arma::colvec Family::Gamma::variance(const arma::colvec& mu) const {
  return arma::square(mu);
}

}