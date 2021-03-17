// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "MultSbmUndirected.h"
using namespace Rcpp;



double MultSbmUndirected::icl_emiss(const List & obs_stats){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::cube edges_counts =as<arma::cube>(obs_stats["x_counts"]);
  // lets go
  double icl_emiss = 0;
  for (int k=0;k<K;++k){
    for (int l=0;l<=k;++l){
      if(k==l){
        if(counts(k)>1){
          arma::vec klcounts = edges_counts.tube(k,l)/2;
          icl_emiss+=lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
        }
      }else{
        arma::vec klcounts = edges_counts.tube(k,l);
        icl_emiss+=lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
      }
      
    }
  }
  return icl_emiss+cst;
}

double MultSbmUndirected::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::cube edges_counts =as<arma::cube>(obs_stats["x_counts"]);
  arma::mat si = submatcross(oldcl,newcl,counts.n_rows);
  double icl_emiss = 0;
  int k = 0;
  int l = 0;
  int cc = 0;
  for (arma::uword i = 0;i<si.n_rows;++i){
    k=si(i,0);
    l=si(i,1);
    if(counts(k)*counts(l)!=0){
      if(k<l){
        arma::vec klcounts = edges_counts.tube(k,l);
        icl_emiss+=lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
      }
      if((k==l) && (counts(k)>1)){
        arma::vec klcounts = edges_counts.tube(k,l)/2;
        icl_emiss+=lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
      }
    }
  }
  return icl_emiss;
}







double MultSbmUndirected::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
  // here old refers to the stats before the fusion between obk and obl
  //Rcout << "Je calculs des corrections !!" << std::endl;
  //Rcout << obk << "---- " << obl << std::endl;
  int a,b,ao,bo,lo;
  double icl_cor = 0;
  int cc, cc_old;
  double oxc,xc;
  arma::vec old_counts =as<arma::vec>(old_stats["counts"]);
  arma::cube old_x_counts =as<arma::cube>(old_stats["x_counts"]);
  arma::vec klcounts;
  cc = counts(k)*counts(l);
  arma::uvec kl;
  kl << k << l << arma::endr;
  arma::uvec mkl;
  mkl << obk << obl << arma::endr;
  if(l>=obk){
    lo=l+1;
  }else{
    lo=l;
  }
  for(int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      
      
      a = kl(i);
      b = mkl(j);
      if(b>=obk){
        b=b-1;
      }
      
      
      
      // new stats no fusion k/l
      if(j==1){
        
        klcounts = x_counts.tube(a,b);
        icl_cor -= lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
        
      }
      
      // new stats fusion k/l
      if((j==1) & (i==0)){
        klcounts  = x_counts.tube(k,b)+x_counts.tube(l,b);
        icl_cor += lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
      }
      
      // handling matrix sizes differences
      ao = kl(i);
      bo = mkl(j);
      if(ao>=obk){
        ao=ao+1;
      }
      
      // old stats no fusion k/l
      klcounts  = old_x_counts.tube(ao,bo);
      icl_cor += lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
      // old stats fusion k/l
      if(i==0){
        klcounts  = old_x_counts.tube(ao,bo)+old_x_counts.tube(lo,bo);
        icl_cor -= lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
        
      }
      
      
    }
  }
  
  return icl_cor;
}