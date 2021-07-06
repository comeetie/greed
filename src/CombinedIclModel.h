#ifndef COMBINEDICLMODEL
#define COMBINEDICLMODEL

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MergeMat.h"
#include "Partition.h"
#include "IclModel.h"
#include "IclModelEmission.h"

using namespace Rcpp;

class CombinedIclModel : public IclModel
{
public:
  CombinedIclModel(std::vector<IclModelEmission*> IclModelsi, S4 model,arma::vec cli,bool verb=false);
  void set_cl(arma::vec cli);
  double icl_emiss(const List & obs_stats);
  double icl_emiss(const List & obs_stats,int oldcl,int newcl);
  double icl_prop(arma::vec counts);
  double icl_prop(arma::vec counts,int oldcl,int newcl);
  arma::mat delta_prop_swap(int i,arma::uvec iclust);
  arma::mat delta_swap(int i,arma::uvec iclust);
  void swap_update(int i, int newcl);
  double delta_merge(int k, int l);
  double delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats);
  void merge_update(int k, int l);
  
  
  // main method for greedy swaping
  void greedy_swap(int nbpassmax, arma::vec workingset,arma::uvec iclust);
  // main method for greedy swaping with move constraints KxK sparse matrix
  void greedy_swap(int nbpassmax, arma::vec workingset,arma::sp_mat & move_mat);
  // methods  to compute merge matrix deltas 
  MergeMat delta_merge();
  // update version
  MergeMat delta_merge(arma::mat delta, int obk, int obl,const List & old_stats);
  // method  to compute merge matrix deltas with constraints on possible merge
  SpMergeMat nasty_delta_merge(const arma::sp_mat & merge_graph);
  SpMergeMat delta_merge(const arma::sp_mat & merge_graph);
  // method  to compute merge matrix deltas with constraints on possible merge (update version)
  SpMergeMat delta_merge(arma::sp_mat & merge_graph, int obk, int obl,const List & old_stats);
  // main method for greedy swaping
  void greedy_merge();
  arma::sp_mat greedy_merge(const arma::sp_mat & merge_graph);
  arma::sp_mat batch_greedy_merge(const arma::sp_mat & merge_graph,int nb_try,double reduction_factor);
  // main method for greedy merge
  List greedy_merge_path();
  // get posterior probs p(Zi|X,Z-i)
  arma::mat get_probs();
  // accessors
  List get_obs_stats();
  virtual List get_obs_stats_cst(){return List::create();};
  arma::vec get_cl(){return cl;};
  arma::vec get_counts(){return counts;};
  S4 get_model(){return model;};
  int get_K(){return K;};
  virtual ~CombinedIclModel(){};
private:
  std::vector<IclModelEmission*> IclModels;
  Partition clp;
};

#endif

