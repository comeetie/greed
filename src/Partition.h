#ifndef PARTITION
#define PARTITION

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;



class Partition 
{
public:
  Partition(){};
  Partition(arma::vec & cl,int K);
  arma::vec get_cl();
  void swap(int i,int newcl);
  void merge(int k, int l);
  void erase(int k);
  int get(int i);
private:
  std::vector<int *> cl_pointers;
  std::vector<int *> cl_addr;
  std::vector<int> cl_labels;
  int K;
  int N;
};

#endif
