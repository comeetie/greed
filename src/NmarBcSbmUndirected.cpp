// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "NmarBcSbmUndirected.h"
using namespace Rcpp;

NmarBcSbmUndirected::NmarBcSbmUndirected(arma::sp_mat& xp,arma::sp_mat& xpobs,S4 modeli,arma::vec& clt,bool verb) : NmarBcSbm(xp,xpobs,modeli,clt,verb){
  node_obs= arma::zeros<arma::vec>(N);
  arma::sp_mat node_deg = arma::sum(xtobs,0);
  arma::rowvec test(node_deg);
  int N_obs = arma::accu(test==(N-1));
  node_obs.elem(arma::find(test==(N-1))) = arma::ones<arma::vec>(N_obs);
  set_cl(clt);
}

double NmarBcSbmUndirected::icl_emiss(const List & obs_stats){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::vec counts_obs =as<arma::vec>(obs_stats["counts_obs"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::mat obs_counts =as<arma::mat>(obs_stats["x_counts_obs"]);
  arma::mat cmat = lgamma(a0+edges_counts)+lgamma(obs_counts-edges_counts+b0)+lgamma(a0+b0);
  cmat = cmat - lgamma(a0) - lgamma(b0) - lgamma(obs_counts+a0+b0);
  
  arma::vec d_counts = obs_counts.diag()/2;
  arma::vec ec = edges_counts.diag()/2;
  cmat.diag() =  lgamma(a0+ec)+lgamma(d_counts-ec+b0)+lgamma(a0+b0) - lgamma(a0) - lgamma(b0) - lgamma(d_counts+a0+b0);
  arma::uvec lonely = arma::find(counts==1);
  for (arma::uword i = 0;i<lonely.n_elem;++i){
    cmat(lonely(i),lonely(i))=0;
  }
  arma::vec correction_bc = lgamma(a0obs+b0obs)+lgamma(a0obs+counts_obs)+lgamma(b0obs+counts-counts_obs)-lgamma(a0obs+b0obs+counts)-lgamma(a0obs)-lgamma(b0obs);
  double icl_emiss=accu(arma::trimatl(cmat))+accu(correction_bc);
  return icl_emiss;
}





double NmarBcSbmUndirected::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::mat obs_counts =as<arma::mat>(obs_stats["x_counts_obs"]);
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
        cc = obs_counts(k,l);
        icl_emiss += lgamma(a0+edges_counts(k,l))+lgamma(b0+cc-edges_counts(k,l))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);  
      }
      if((k==l) && (counts(k)>1)){
        cc = obs_counts(k,k)/2;
        icl_emiss += lgamma(a0+edges_counts(k,l)/2)+lgamma(b0+cc-edges_counts(k,l)/2)+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);  
      }
    }
    
  }
  if(counts(k)>0){
    icl_emiss += lgamma(a0obs+b0obs)+lgamma(a0obs+counts_obs(k))+lgamma(b0obs+counts(k)-counts_obs(k))-lgamma(a0obs+b0obs+counts(k))-lgamma(a0obs)-lgamma(b0obs);
  }
  
  if(counts(l)>0){
    icl_emiss += lgamma(a0obs+b0obs)+lgamma(a0obs+counts_obs(l))+lgamma(b0obs+counts(l)-counts_obs(l))-lgamma(a0obs+b0obs+counts(l))-lgamma(a0obs)-lgamma(b0obs);
  }
  return icl_emiss;
}



double NmarBcSbmUndirected::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
  // here old refers to the stats before the fusion between obk and obl
  //Rcout << "Je calculs des corrections !!" << std::endl;
  //Rcout << obk << "---- " << obl << std::endl;
  int a,b,ao,bo,lo;
  double icl_cor = 0;
  int cc, cc_old;
  double oxc,xc;
  arma::vec old_counts =as<arma::vec>(old_stats["counts"]);
  arma::mat old_x_counts =as<arma::mat>(old_stats["x_counts"]);
  arma::mat old_obs_counts =as<arma::mat>(old_stats["x_counts_obs"]);
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
        cc = x_counts_obs(a,b);
        icl_cor -= lgamma(a0+x_counts(a,b))+lgamma(b0+cc-x_counts(a,b))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);
      }
      
      // new stats fusion k/l
      if((j==1) & (i==0)){
        
        cc = x_counts_obs(k,b)+x_counts_obs(l,b);
        xc    = x_counts(k,b)+x_counts(l,b);
        icl_cor += lgamma(a0+xc)+lgamma(b0+cc-xc)+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);
      }
      
      // handling matrix sizes differences
      ao = kl(i);
      bo = mkl(j);
      if(ao>=obk){
        ao=ao+1;
      }
      
      // old stats no fusion k/l
      cc_old = old_obs_counts(ao,bo);
      icl_cor += lgamma(a0+old_x_counts(ao,bo))+lgamma(b0+cc_old-old_x_counts(ao,bo))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc_old);
      // old stats fusion k/l
      if(i==0){
        
        cc_old = old_obs_counts(ao,bo)+old_obs_counts(lo,bo);
        oxc    = old_x_counts(ao,bo)+old_x_counts(lo,bo);
        icl_cor -= lgamma(a0+oxc)+lgamma(b0+cc_old-oxc)+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc_old);
      }
      
      
    }
  }
  
  return icl_cor;
  
  
}



