#ifndef ICLMODEL
#define ICLMODEL

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;




class IclModel
{
public:
  double icl(const List & obs_stats);
  double icl(const List & obs_stats,int oldcl,int newcl);
  virtual double icl_emiss(const List & obs_stats){};
  virtual double icl_emiss(const List & obs_stats,int oldcl,int newcl){};
  virtual arma::mat delta_swap(const int i){};
  virtual void swap_update(const int i,const int newcl){};
  void greedy_swap(int nbpassmax);
  virtual MergeMat delta_merge(){};
  virtual MergeMat delta_merge(arma::mat delta, int obk, int obl){};
  virtual void merge_update(const int k,const int l){};
  void greedy_merge();
  List greedy_merge_path();
  virtual List get_obs_stats(){};
  arma::vec get_cl(){return cl;};
  arma::vec get_counts(){return counts;};
  arma::mat get_x_counts(){return x_counts;};
protected:
  double alpha;
  int K;
  int N;
  arma::vec cl;
  arma::vec counts;
  arma::mat x_counts;
  double icl_value;
};

#endif

