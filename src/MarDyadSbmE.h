#ifndef MARDYADSBME
#define MARDYADSBME

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "IclModelEmission.h"
using namespace Rcpp;

class MarDyadSbmE : public IclModelEmission
{
public:
  MarDyadSbmE(arma::sp_mat& xpobs,S4 model,bool verb){
    List sampling_priors = as<List>(model.slot("sampling_priors"));
    double a0obs = sampling_priors["a0obs"];
    double b0obs = sampling_priors["b0obs"];
    int nbobs = arma::accu(xpobs);
    int nbdyads = xpobs.n_cols*xpobs.n_rows;
    cst = lgamma(a0obs+nbobs)+lgamma(nbdyads-nbobs+b0obs)+lgamma(a0obs+b0obs)- lgamma(a0obs) - lgamma(b0obs) - lgamma(nbdyads+a0obs+b0obs);
  };
  void set_cl(arma::vec cl){};
  double icl_emiss(const List & obs_stats){return cst;};
  double icl_emiss(const List & obs_stats,int oldcl,int newcl){return 0;};
  arma::mat delta_swap(int i,int K, const arma::vec cl,arma::uvec iclust){
    arma::vec delta(K);
    delta.fill(0);
    return delta;
  };
  void swap_update(int i,const arma::vec cl,int newcl){};
  double delta_merge(int k, int l){return 0;};
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){return 0;};
  void merge_update(int k, int l){};
  List get_obs_stats(){return List::create(Named("cst", cst));};
private:
  double cst;
};

#endif

