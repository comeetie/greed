#ifndef family_H
#define family_H

#include <RcppArmadillo.h>
#include "linkfunctions.h"
// taken from https://github.com/dirkschumacher/rcppglm for testing purpose
namespace Family {

class ExponentialFamily {
  // TODO: make variance a parameter as well
public:
  arma::colvec link_fun(const arma::colvec& eta) const {
    return link_function->link(eta);
  }
  arma::colvec link_inverse(const arma::colvec& eta) const {
    return link_function->inverse(eta);
  }
  arma::colvec link_mu_eta(const arma::colvec& eta) const {
    return link_function->mu_eta(eta);
  }
  virtual arma::colvec variance(const arma::colvec& mu) const = 0;
  ExponentialFamily(std::unique_ptr<Link::LinkFunction>& link) : link_function(std::move(link)) {
  }
  virtual arma::colvec log_lik(const arma::colvec& eta,const arma::colvec& y) const = 0;
  virtual arma::colvec log_lik_norm(const arma::colvec& eta,const arma::colvec& y) const = 0;
  virtual ~ExponentialFamily() {}
private:
  std::unique_ptr<Link::LinkFunction> link_function;
};


class Binomial : public ExponentialFamily {
public:
  arma::colvec variance(const arma::colvec& mu) const;
  arma::colvec log_lik(const arma::colvec& eta,const arma::colvec& y) const;
  arma::colvec log_lik_norm(const arma::colvec& eta,const arma::colvec& y) const;
  Binomial(std::unique_ptr<Link::LinkFunction>& link) : ExponentialFamily(link) {}
};

class Poisson : public ExponentialFamily {
public:
  arma::colvec variance(const arma::colvec& mu) const;
  arma::colvec log_lik(const arma::colvec& eta,const arma::colvec& y) const;
  arma::colvec log_lik_norm(const arma::colvec& eta,const arma::colvec& y) const;
  Poisson(std::unique_ptr<Link::LinkFunction>& link) : ExponentialFamily(link) {}
};

class Gamma : public ExponentialFamily {
public:
  arma::colvec variance(const arma::colvec& mu) const;
  arma::colvec log_lik(const arma::colvec& eta,const arma::colvec& y) const;
  arma::colvec log_lik_norm(const arma::colvec& eta,const arma::colvec& y) const;
  Gamma(std::unique_ptr<Link::LinkFunction>& link) : ExponentialFamily(link) {}
};

}
#endif