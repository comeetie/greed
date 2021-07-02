// [[Rcpp::depends(RcppArmadillo)]]
#include "Partition.h"
using namespace Rcpp;



Partition::Partition(arma::vec & cl, int Ki){
  K =Ki;
  cl_labels = std::vector<int>(K);
  cl_addr = std::vector<int*>(K);;
  N = cl.n_elem;
  cl_pointers = std::vector<int*>(N);
  
  for(int k=0;k<K;k++){
    cl_labels[k]=k;
    cl_addr[k]= &cl_labels[k]; 
  }
  arma::vec::iterator it;
  int i=0;
  for (it=cl.begin();it!=cl.end();++it){
    cl_pointers[i]=&cl_labels[*it];
    i++;
  }
  
}

arma::vec Partition::get_cl(){
  std::vector<int *>::iterator it;
  arma::vec cl(N);
  int i=0;
  for (it=cl_pointers.begin();it!=cl_pointers.end();++it){
    cl(i)=*(*it);
    i++;
  }
  return cl;
}


void Partition::swap(int i,int newcl){
  cl_pointers[i]=cl_addr[newcl];
}


void Partition::erase(int k){

  cl_addr.erase(cl_addr.begin()+k);
  for(int i=0;i<K;i++){
    if(cl_labels[i]==k){
      cl_labels[i]=-1;
    }
    if(cl_labels[i]>k){
      cl_labels[i]=cl_labels[i]-1;
    }
  }
}

void Partition::merge(int k, int l){
  for(int i=0;i<K;i++){
    if(cl_labels[i]==k){
      cl_labels[i]=l;
    }
  }
  erase(k);
}
  
int Partition::get(int i){
  return *(cl_pointers[i]);
}
