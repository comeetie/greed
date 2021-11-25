#ifndef ICLMODELEMISSION
#define ICLMODELEMISSION

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;



// Main class with only virtual method for models to implement / expected to work with CombinedIclModel
class IclModelEmission
{
public:
  virtual void set_cl(arma::uvec clt){};
  // compute icl
  virtual double icl_emiss(const List & obs_stats){return 0;};
  // virtual methods to be implemented by models to compute log(p(X|Z)) optimized for deltas
  virtual double icl_emiss(const List & obs_stats,int oldcl,int newcl,bool dead_cluster){return 0;};
  // compute the delta for each possible swap of a node
  virtual arma::vec delta_swap(const int i,arma::uvec & cl,bool dead_cluster, arma::uvec iclust,int K){return NULL;};
  // update the stats when a node is swapped
  virtual void swap_update(const int i,arma::uvec & cl,bool dead_cluster,const int newcl){};
  // virtual methods to be implemented by models to compute merge deltas
  virtual double delta_merge(int k, int l){return 0;};
  // update the stats when two clusters are merged
  virtual void merge_update(const int k,const int l){};
  // compute correction if needed to merge matrix
  virtual double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){return 0;};
  // methods  to compute merge matrix deltas 
  virtual List get_obs_stats(){return List::create();};
  virtual List get_obs_stats_cst(){return List::create();};
  virtual ~IclModelEmission(){};
  S4 get_model(){return model;};
protected:
  bool verbose;
  // priors params
  S4 model;
};

#endif

