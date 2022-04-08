#include <RcppArmadillo.h>
#include "linkfunctions.h"
// taken from https://github.com/dirkschumacher/rcppglm for testing purpose
arma::colvec Link::Identity::link(const arma::colvec& mu) const {
  return mu;
}

arma::colvec Link::Identity::inverse(const arma::colvec& eta) const {
  return eta;
}

arma::colvec Link::Identity::mu_eta(const arma::colvec& eta) const {
  return arma::ones<arma::colvec>(eta.n_elem);
}


arma::colvec Link::Logit::link(const arma::colvec& mu) const {
  return arma::log(mu / (1.0 - mu));
}

arma::colvec Link::Logit::inverse(const arma::colvec& eta) const {
  return (arma::exp(eta) / (1.0 + arma::exp(eta)));
}

arma::colvec Link::Logit::mu_eta(const arma::colvec& eta) const {
  return arma::exp(eta) / arma::square(arma::exp(eta) + 1.0);
}


arma::colvec Link::Probit::link(const arma::colvec& mu) const {
  arma::colvec mu_local = mu;
  // TODO: Use Boost???
  mu_local.for_each([](arma::colvec::elem_type& val) {
    val = R::qnorm5(val, 0.0, 1.0, 0, 0);
  });
  return mu_local;
}

arma::colvec Link::Probit::inverse(const arma::colvec& eta) const {
  return arma::normcdf(eta);
}

arma::colvec Link::Probit::mu_eta(const arma::colvec& eta) const {
  return arma::normpdf(eta);
}


arma::colvec Link::Log::link(const arma::colvec& mu) const {
  return arma::log(mu);
}

arma::colvec Link::Log::inverse(const arma::colvec& eta) const {
  return arma::exp(eta);
}

arma::colvec Link::Log::mu_eta(const arma::colvec& eta) const {
  return arma::exp(eta);
}


arma::colvec Link::Inverse::link(const arma::colvec& mu) const {
  return 1.0 / mu;
}

arma::colvec Link::Inverse::inverse(const arma::colvec& eta) const {
  return 1.0 / eta;
}

arma::colvec Link::Inverse::mu_eta(const arma::colvec& eta) const {
  return -1.0 / arma::square(eta);
}


arma::colvec Link::Sqrt::link(const arma::colvec& mu) const {
  return arma::sqrt(mu);
}

arma::colvec Link::Sqrt::inverse(const arma::colvec& eta) const {
  return arma::square(eta);
}

arma::colvec Link::Sqrt::mu_eta(const arma::colvec& eta) const {
  return 2.0 * eta;
}