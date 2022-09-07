// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "SbmUndirected.h"
using namespace Rcpp;


double SbmUndirected::icl_emiss(const List & obs_stats){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::mat matcount = counts*counts.t();
  arma::mat cmat = lgamma(a0+edges_counts)+lgamma(matcount-edges_counts+b0)+lgamma(a0+b0);
  cmat = cmat - lgamma(a0) - lgamma(b0) - lgamma(matcount+a0+b0);
  arma::vec d_counts = (matcount.diag()-counts)/2;
  arma::vec ec = edges_counts.diag()/2;
  cmat.diag() =  lgamma(a0+ec)+lgamma(d_counts-ec+b0)+lgamma(a0+b0) - lgamma(a0) - lgamma(b0) - lgamma(d_counts+a0+b0);
  arma::uvec lonely = arma::find(counts==1);
  for (arma::uword i = 0;i<lonely.n_elem;++i){
    cmat(lonely(i),lonely(i))=0;
  }
  double icl_emiss=accu(arma::trimatl(cmat));
  return icl_emiss;
}

double SbmUndirected::icl_emiss(const List & obs_stats,int oldcl,int newcl,bool dead_cluster){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::umat si = submatcross(oldcl,newcl,counts.n_rows);
  double icl_emiss = 0;
  int k = 0;
  int l = 0;
  int cc = 0;
  for (arma::uword i = 0;i<si.n_rows;++i){
    k=si(i,0);
    l=si(i,1);
    if(counts(k)*counts(l)!=0){
      if(k<l){
        cc = counts(k)*counts(l);
        icl_emiss += lgamma(a0+edges_counts(k,l))+lgamma(b0+cc-edges_counts(k,l))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);  
      }
      if((k==l) && (counts(k)>1)){
        cc = counts(k)*(counts(k)-1)/2;
        icl_emiss += lgamma(a0+edges_counts(k,l)/2)+lgamma(b0+cc-edges_counts(k,l)/2)+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);  
      }
    }
    
  }
  return icl_emiss;
}

double SbmUndirected::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
  // here old refers to the stats before the fusion between obk and obl
  //Rcout << "Je calculs des corrections !!" << std::endl;
  //Rcout << obk << "---- " << obl << std::endl;
  int a,b,ao,bo,lo;
  double icl_cor = 0;
  int cc, cc_old;
  double oxc,xc;
  arma::vec old_counts =as<arma::vec>(old_stats["counts"]);
  arma::mat old_x_counts =as<arma::mat>(old_stats["x_counts"]);
  cc = counts(k)*counts(l);
  arma::ivec kl({k, l});
  arma::ivec mkl({obk, obl});
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
        cc = counts(a)*counts(b);
        icl_cor -= lgamma(a0+x_counts(a,b))+lgamma(b0+cc-x_counts(a,b))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);
      }
      
      // new stats fusion k/l
      if((j==1) & (i==0)){
        cc = (counts(k)+counts(l))*counts(b);
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
      cc_old = old_counts(ao)*old_counts(bo);
      icl_cor += lgamma(a0+old_x_counts(ao,bo))+lgamma(b0+cc_old-old_x_counts(ao,bo))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc_old);
      // old stats fusion k/l
      if(i==0){
        cc_old = (old_counts(ao)+old_counts(lo))*old_counts(bo);
        oxc    = old_x_counts(ao,bo)+old_x_counts(lo,bo);
        icl_cor -= lgamma(a0+oxc)+lgamma(b0+cc_old-oxc)+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc_old);
      }
      
      
    }
  }
  
  return icl_cor;
  
  
}
