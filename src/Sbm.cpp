// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModelEmission.h"
#include "Sbm.h"
using namespace Rcpp;

Sbm::Sbm(const arma::sp_mat  & xp,S4 modeli,bool verb){
  model = modeli;
  x  = xp;
  xt = xp.t();
  N  = x.n_rows;
  a0 = model.slot("a0");
  b0 = model.slot("b0");
  x  = xp;
  xt = xp.t();
  verbose=verb;
}



void Sbm::set_cl(arma::uvec clt){
  K = arma::max(clt)+1;
  x_counts = gsum_mat(clt,x,K);
  counts = count(clt,K);
}

double Sbm::icl_emiss(const List & obs_stats){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::mat matcount = counts*counts.t();
  arma::mat cmat = lgamma(a0+edges_counts)+lgamma(matcount-edges_counts+b0)+lgamma(a0+b0);
  cmat = cmat - lgamma(a0) - lgamma(b0) - lgamma(matcount+a0+b0);
  double icl_emiss=accu(cmat);
  return icl_emiss;
}

double Sbm::icl_emiss(const List & obs_stats,int oldcl,int newcl, bool dead_cluster){
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
      cc = counts(k)*counts(l);
      icl_emiss += lgamma(a0+edges_counts(k,l))+lgamma(b0+cc-edges_counts(k,l))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);  
    }
    
  }
  return icl_emiss;
  
  
}



List Sbm::get_obs_stats(){
  return List::create(Named("counts", counts), Named("x_counts", x_counts));
}



arma::vec Sbm::delta_swap(int i,arma::uvec & cl, bool almost_dead_cluster, arma::uvec iclust, int K){
  int self=x(i,i);
  int oldcl = cl(i);
  arma::sp_mat col_sum = gsum_col(cl,x,i,K);
  col_sum(oldcl)=col_sum(oldcl)-self;
  arma::sp_mat row_sum = gsum_col(cl,xt,i,K);
  row_sum(oldcl)=row_sum(oldcl)-self;
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  int k = 0;
  // for each possible move
  for(arma::uword j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    if(k!=oldcl){
      arma::mat new_ec = x_counts;
      new_ec.col(k) = new_ec.col(k)+col_sum;
      new_ec.row(k) = new_ec.row(k)+row_sum.t();
      new_ec.col(oldcl) = new_ec.col(oldcl)-col_sum;
      new_ec.row(oldcl) = new_ec.row(oldcl)-row_sum.t();
      new_ec(k,k)=new_ec(k,k)+self;
      new_ec(oldcl,oldcl)=new_ec(oldcl,oldcl)-self;
      arma::vec new_counts = update_count(counts,oldcl,k);
      List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
      delta(k)=icl_emiss(new_stats,oldcl,k,almost_dead_cluster)-icl_emiss(old_stats,oldcl,k,false);
    }
  }
  return delta;
  
}


void Sbm::swap_update(const int i,arma::uvec &  cl,bool dead_cluster,const int newcl){
  int self=x(i,i);
  int oldcl = cl(i);
  arma::sp_mat col_sum = gsum_col(cl,x,i,K);
  col_sum(oldcl)=col_sum(oldcl)-self;
  arma::sp_mat row_sum = gsum_col(cl,xt,i,K);
  row_sum(oldcl)=row_sum(oldcl)-self;
  
  x_counts.col(newcl) = x_counts.col(newcl)+col_sum;
  x_counts.row(newcl) = x_counts.row(newcl)+row_sum.t();
  x_counts.col(oldcl) = x_counts.col(oldcl)-col_sum;
  x_counts.row(oldcl) = x_counts.row(oldcl)-row_sum.t();
  x_counts(newcl,newcl)=x_counts(newcl,newcl)+self;
  x_counts(oldcl,oldcl)=x_counts(oldcl,oldcl)-self;
  counts = update_count(counts,oldcl,newcl);
  if(dead_cluster){
    counts.shed_row(oldcl);
    x_counts = x_counts(arma::find(arma::linspace(0,K-1,K)!=oldcl),arma::find(arma::linspace(0,K-1,K)!=oldcl));
    --K;
  }
  
  
}


double Sbm::delta_merge(int k, int l){
  
  arma::mat new_ec = x_counts;
  arma::mat new_counts = counts;
  new_counts(l)=new_counts(k)+new_counts(l);
  new_counts(k)=0;
  new_ec.col(l) = new_ec.col(k)+new_ec.col(l);
  new_ec.row(l) = new_ec.row(k)+new_ec.row(l);
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
  double delta = icl_emiss(new_stats,k,l,true)-icl_emiss(old_stats,k,l,false);
  return delta;
  
}


double Sbm::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
  // here old refers to the stats before the fusion between obk and obl

  int a,b,ao,bo,lo;
  double icl_cor = 0;
  int cc, cc_old;
  double oxc,xc;
  arma::vec old_counts =as<arma::vec>(old_stats["counts"]);
  arma::mat old_x_counts =as<arma::mat>(old_stats["x_counts"]);
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
        cc = counts(a)*counts(b);
        icl_cor -= lgamma(a0+x_counts(a,b))+lgamma(b0+cc-x_counts(a,b))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);
        icl_cor -= lgamma(a0+x_counts(b,a))+lgamma(b0+cc-x_counts(b,a))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);
      }
      
      // new stats fusion k/l
      if((j==1) & (i==0)){
        cc = (counts(k)+counts(l))*counts(b);
        xc    = x_counts(k,b)+x_counts(l,b);
        icl_cor += lgamma(a0+xc)+lgamma(b0+cc-xc)+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);
        xc    = x_counts(b,k)+x_counts(b,l);
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
      icl_cor += lgamma(a0+old_x_counts(bo,ao))+lgamma(b0+cc_old-old_x_counts(bo,ao))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc_old);
      // old stats fusion k/l
      if(i==0){
        cc_old = (old_counts(ao)+old_counts(lo))*old_counts(bo);
        oxc    = old_x_counts(ao,bo)+old_x_counts(lo,bo);
        icl_cor -= lgamma(a0+oxc)+lgamma(b0+cc_old-oxc)+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc_old);
        oxc    = old_x_counts(bo,ao)+old_x_counts(bo,lo);
        icl_cor -= lgamma(a0+oxc)+lgamma(b0+cc_old-oxc)+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc_old);
      }
      
      
    }
  }
  return icl_cor;
}



void Sbm::merge_update(int k,int l){
  counts(l) = counts(k)+counts(l);
  counts.shed_row(k);
  
  x_counts.col(l) = x_counts.col(k)+x_counts.col(l);
  x_counts.row(l) = x_counts.row(k)+x_counts.row(l);
  x_counts = x_counts(arma::find(arma::linspace(0,K-1,K)!=k),arma::find(arma::linspace(0,K-1,K)!=k));
  
  --K;
}



