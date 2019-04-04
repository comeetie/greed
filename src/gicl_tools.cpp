// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat submatcross(int oldcl,int newcl,int K){
  arma::mat result(4*(K-1),2);
  int nbr = 0;
  result.fill(0);
  for(int i = 0; i < K; ++i) {
    result(i,0)=oldcl;
    result(i,1)=i;
    result(i+K,0)=newcl;
    result(i+K,1)=i;
    if(i==oldcl | i==newcl){
      nbr = nbr + 1;
    }else{
      result(i+2*K-nbr,1)=oldcl;
      result(i+2*K-nbr,0)=i;
      result(i+(3*K-2)-nbr,1)=newcl;
      result(i+(3*K-2)-nbr,0)=i;
    }
  }

  return result;
}

arma::vec count(arma::vec cl,int K){
  arma::vec result(K);
  result.fill(0);
  for(int i = 0; i < cl.n_elem; ++i) {
    result(cl(i),0)+=1;
  }
  return result;
}

arma::mat gsum_mat(arma::vec cl,const arma::sp_mat& x,int K) {
  arma::mat result(K,K);
  result.fill(0);
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    result(cl(i.row()),cl(i.col())) += *i;
  }
  return result;
}


arma::mat gsum_mm(arma::vec cl,const arma::sp_mat& x,int K) {
  arma::mat result(K,x.n_cols);
  result.fill(0);
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    result(cl(i.row()),i.col()) += *i;
  }
  return result;
}

arma::mat gsum_col(arma::vec cl,const arma::sp_mat& x,int i, int K) {
  arma::sp_mat ccol = x.col(i);
  arma::mat result(K,1);
  result.fill(0);
  for (arma::sp_mat::const_iterator i = ccol.begin(); i != ccol.end(); ++i) {
    result(cl(i.row()),0) += *i;
  }
  return result;
}


arma::mat update_count(arma::vec counts,int oldcl,int newcl) {
  counts(oldcl)=counts(oldcl)-1;
  counts(newcl)=counts(newcl)+1;
  return counts;
}


arma::vec sum_cols(const arma::sp_mat& x){
  arma::vec tots(x.n_cols);
  tots.fill(0);
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i)  {
    tots(i.col()) += *i;
  }
  return tots;
}

arma::vec sum_cols(const arma::mat x){
  int nbcols = x.n_cols;
  int nbrows = x.n_rows;
  arma::vec tots(x.n_cols);
  tots.fill(0);
  for (int j=0;j<nbcols;++j)  {
    for (int i=0;i<nbrows;++i)  {
      tots(j) += x(i,j);
    }
  }
  return tots;
}

arma::vec sum_rows(const arma::mat x){
  int nbcols = x.n_cols;
  int nbrows = x.n_rows;
  arma::vec tots(x.n_cols);
  tots.fill(0);
  for (int i=0;i<nbrows;++i) {
    for (int j=0;j<nbcols;++j) {
      tots(i) += x(i,j);
    }
  }
  return tots;
}

double sum_lfact(const arma::sp_mat & x){
  arma::sp_mat::const_iterator start = x.begin();
  arma::sp_mat::const_iterator end   = x.end();
  double cst = 0;
  for(arma::sp_mat::const_iterator it = start; it != end; ++it)
  {
     cst += 1; 
  }
  return cst;
}


