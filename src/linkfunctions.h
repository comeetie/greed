#ifndef linkfunctions_H
#define linkfunctions_H
// taken from https://github.com/dirkschumacher/rcppglm for testing purpose
#include <RcppArmadillo.h>

namespace Link {

class LinkFunction {
public:
  virtual arma::colvec link(const arma::colvec& eta) const = 0;
  virtual arma::colvec inverse(const arma::colvec& eta) const = 0;
  virtual arma::colvec mu_eta(const arma::colvec& eta) const = 0;
  virtual ~LinkFunction() {}
};

class Identity : public LinkFunction {
public:
  arma::colvec link(const arma::colvec& eta) const;
  arma::colvec inverse(const arma::colvec& eta) const;
  arma::colvec mu_eta(const arma::colvec& eta) const;
};

class Logit : public LinkFunction {
public:
  arma::colvec link(const arma::colvec& eta) const;
  arma::colvec inverse(const arma::colvec& eta) const;
  arma::colvec mu_eta(const arma::colvec& eta) const;
};

class Probit : public LinkFunction {
public:
  arma::colvec link(const arma::colvec& eta) const;
  arma::colvec inverse(const arma::colvec& eta) const;
  arma::colvec mu_eta(const arma::colvec& eta) const;
};

class Log : public LinkFunction {
public:
  arma::colvec link(const arma::colvec& eta) const;
  arma::colvec inverse(const arma::colvec& eta) const;
  arma::colvec mu_eta(const arma::colvec& eta) const;
};

class Inverse : public LinkFunction {
public:
  arma::colvec link(const arma::colvec& eta) const;
  arma::colvec inverse(const arma::colvec& eta) const;
  arma::colvec mu_eta(const arma::colvec& eta) const;
};

class Sqrt : public LinkFunction {
public:
  arma::colvec link(const arma::colvec& eta) const;
  arma::colvec inverse(const arma::colvec& eta) const;
  arma::colvec mu_eta(const arma::colvec& eta) const;
};

}

#endif
