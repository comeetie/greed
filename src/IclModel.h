#ifndef ICLMODEL
#define ICLMODEL

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "MergeMat.h"
using namespace Rcpp;



// Main class with algorithmic logic and virtual method for models to implement
class IclModel
{
public:
  // compute icl
  double icl(const List & obs_stats);
  // compute icl optimized for deltas
  double icl(const List & obs_stats,int oldcl,int newcl);
  // virtual methods to be implemented by models to compute log(p(X|Z))
  virtual double icl_emiss(const List & obs_stats){};
  // virtual methods to be implemented by models to compute log(p(X|Z)) optimized for deltas
  virtual double icl_emiss(const List & obs_stats,int oldcl,int newcl){};
  // compute the delta for each possible swap of a node
  virtual arma::mat delta_swap(const int i){};
  // update the stats when a node is swaped
  virtual void swap_update(const int i,const int newcl){};
  // main method for greedy swaping
  void greedy_swap(int nbpassmax);
  // virtual methods to be implemented by models to compute merge matrix deltas
  virtual MergeMat delta_merge(){};
  // virtual methods to be implemented by models to compute merge matrix deltas update version
  virtual MergeMat delta_merge(arma::mat delta, int obk, int obl){};
  // update the stats when two clusters are merged
  virtual void merge_update(const int k,const int l){};
  // main method for greedy swaping
  void greedy_merge();
  // ain method for greedy merge
  List greedy_merge_path();
  // accessors
  virtual List get_obs_stats(){};
  arma::vec get_cl(){return cl;};
  arma::vec get_counts(){return counts;};
  arma::mat get_x_counts(){return x_counts;};
protected:
  // prior 
  double alpha;
  // nb clusters
  int K;
  // problem dim
  int N;
  // vector of cluster labels
  arma::vec cl;
  // stats
  // vector of cluster counts
  arma::vec counts;
  // matrix of observed counts for each clusters
  arma::mat x_counts;
  double icl_value;
};

#endif

