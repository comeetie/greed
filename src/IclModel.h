#ifndef ICLMODEL
#define ICLMODEL

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "SpMergeMat.h"
using namespace Rcpp;



// Main class with algorithmic logic and virtual method for models to implement
class IclModel
{
public:
  virtual void set_cl(arma::vec clt){};
  // compute icl
  double icl(const List & obs_stats);
  // compute icl optimized for deltas
  double icl(const List & obs_stats,int oldcl,int newcl);
  // compute icl optimized for deltas and sparse updates
  double icl(const List & obs_stats,const List & up_stats,int oldcl,int newcl);
  // virtual methods to be implemented by models to compute log(p(X|Z))
  virtual double icl_emiss(const List & obs_stats){return 0;};
  // virtual methods to be implemented by models to compute log(p(X|Z)) optimized for deltas
  virtual double icl_emiss(const List & obs_stats,int oldcl,int newcl){return 0;};
  // virtual methods to be implemented by models to compute log(p(X|Z)) optimized for deltas and sparse updates
  virtual double icl_emiss(const List & obs_stats,const List & up_stats,int oldcl,int newcl){return 0;};
  // compute the delta for each possible swap of a node
  virtual arma::mat delta_swap(const int i,arma::uvec iclust){return NULL;};
  // update the stats when a node is swaped
  virtual void swap_update(const int i,const int newcl){};
  // main method for greedy swaping
  void greedy_swap(int nbpassmax, arma::vec workingset,arma::uvec iclust);
  // main method for greedy swaping with move constraints KxK sparse matrix
  void greedy_swap(int nbpassmax, arma::vec workingset,arma::sp_mat & move_mat);
  // virtual methods to be implemented by models to compute merge deltas
  virtual double delta_merge(int k, int l){return 0;};
  // update the stats when two clusters are merged
  virtual void merge_update(const int k,const int l){};
  // compute correction if needed to merge matrix
  virtual double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){return 0;};
  // methods  to compute merge matrix deltas 
  MergeMat delta_merge();
  // update version
  MergeMat delta_merge(arma::mat delta, int obk, int obl,const List & old_stats);
  // method  to compute merge matrix deltas with constraints on possible merge
  SpMergeMat nasty_delta_merge(const arma::sp_mat & merge_graph);
  SpMergeMat delta_merge(const arma::sp_mat & merge_graph);
  // method  to compute merge matrix deltas with constraints on possible merge (update version)
  SpMergeMat delta_merge(const arma::sp_mat & merge_graph, int obk, int obl,const List & old_stats);
  // main method for greedy swaping
  void greedy_merge();
  arma::sp_mat greedy_merge(const arma::sp_mat & merge_graph);
  arma::sp_mat batch_greedy_merge(const arma::sp_mat & merge_graph,int nb_try,double reduction_factor);
  // main method for greedy merge
  List greedy_merge_path();
  // get posterior probs p(Zi|X,Z-i)
  arma::mat get_probs();
  // accessors
  virtual List get_obs_stats(){return List::create();};
  arma::vec get_cl(){return cl;};
  arma::vec get_counts(){return counts;};
  int get_K(){return K;};
  virtual ~IclModel(){};
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

  double icl_value;
  // verbose ?
  bool verbose;
};

#endif

