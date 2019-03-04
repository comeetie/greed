// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "Sbm.h"
using namespace Rcpp;


Sbm::Sbm(arma::sp_mat& xp,int Ki,double alphai,double a0i,double b0i){
  alpha = alphai;
  a0 = a0i;
  b0 = b0i;
  x  = xp;
  xt = xp.t();
  N  = x.n_rows;
  K  = Ki;
  cl = as<arma::vec>(sample(K,N,true));
  cl = cl -1;
  x_counts = gsum_mat(cl,x,K);
  counts = count(cl,K);
}

Sbm::Sbm(arma::sp_mat& xp,int Ki,double alphai,double a0i,double b0i,arma::vec& clt){
  alpha = alphai;
  a0 = a0i;
  b0 = b0i;
  x  = xp;
  xt = xp.t();
  N  = x.n_rows;
  K  = Ki;
  cl = clt;
  x_counts = gsum_mat(cl,x,K);
  counts = count(cl,K);
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

double Sbm::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::mat si = submatcross(oldcl,newcl,counts.n_rows);
  double icl_emiss = 0;
  int k = 0;
  int l = 0;
  int cc = 0;
  for (int i = 0;i<si.n_rows;++i){
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

arma::mat Sbm::delta_swap(int i){
  int self=x(i,i);
  int oldcl = cl(i);
  arma::mat col_sum = gsum_col(cl,x,i,K);
  col_sum(oldcl)=col_sum(oldcl)-self;
  arma::mat row_sum = gsum_col(cl,xt,i,K);
  row_sum(oldcl)=row_sum(oldcl)-self;
  arma::mat delta(K,1);
  delta.fill(0);
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  for(int k = 0; k < K; ++k) {
    if(k!=oldcl){
      arma::mat new_ec = x_counts;
      new_ec.col(k) = new_ec.col(k)+col_sum;
      new_ec.row(k) = new_ec.row(k)+row_sum.t();
      new_ec.col(oldcl) = new_ec.col(oldcl)-col_sum;
      new_ec.row(oldcl) = new_ec.row(oldcl)-row_sum.t();
      new_ec(k,k)=new_ec(k,k)+self;
      new_ec(oldcl,oldcl)=new_ec(oldcl,oldcl)-self;
      arma::mat new_counts = update_count(counts,oldcl,k);
      List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
      delta(k,0)=icl(new_stats,oldcl,k)-icl(old_stats,oldcl,k);
    }
  }
  return delta;
}


void Sbm::swap_update(const int i,const int newcl){
  int self=x(i,i);
  int oldcl = cl(i);
  arma::mat col_sum = gsum_col(cl,x,i,K);
  col_sum(oldcl)=col_sum(oldcl)-self;
  arma::mat row_sum = gsum_col(cl,xt,i,K);
  row_sum(oldcl)=row_sum(oldcl)-self;
  arma::mat new_ec = x_counts;
  new_ec.col(newcl) = new_ec.col(newcl)+col_sum;
  new_ec.row(newcl) = new_ec.row(newcl)+row_sum.t();
  new_ec.col(oldcl) = new_ec.col(oldcl)-col_sum;
  new_ec.row(oldcl) = new_ec.row(oldcl)-row_sum.t();
  new_ec(newcl,newcl)=new_ec(newcl,newcl)+self;
  new_ec(oldcl,oldcl)=new_ec(oldcl,oldcl)-self;
  arma::mat new_counts = update_count(counts,oldcl,newcl);
  cl(i)=newcl;
  if(new_counts(oldcl)==0){
    counts = new_counts.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    x_counts = new_ec(arma::find(arma::linspace(0,K-1,K)!=oldcl),arma::find(arma::linspace(0,K-1,K)!=oldcl));
    cl.elem(arma::find(cl>oldcl))=cl.elem(arma::find(cl>oldcl))-1;
    --K;
  }else{
    counts=new_counts;
    x_counts=new_ec;
  }

}


MergeMat Sbm::delta_merge(){
  arma::mat delta(K,K);
  delta.fill(0);
  int bk = 0;
  int bl = 0;
  double bv = -std::numeric_limits<double>::infinity();
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  for(int k = 1 ; k < K ; ++k) {
    for (int l = 0 ; l<k ; ++l){
      arma::mat new_ec = x_counts;
      arma::mat new_counts = counts;
      new_counts(l)=new_counts(k)+new_counts(l);
      new_counts(k)=0;
      new_ec.col(l) = new_ec.col(k)+new_ec.col(l);
      new_ec.row(l) = new_ec.row(k)+new_ec.row(l);
      List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
      delta(k,l)=icl(new_stats,k,l)-icl(old_stats,k,l);
      delta(l,k)=delta(k,l);
      if(delta(k,l)>bv){
        bk=k;
        bl=l;
        bv=delta(k,l);
      }
    }
  }
  return MergeMat(bk,bl,bv,delta);
}

MergeMat Sbm::delta_merge(arma::mat delta, int obk, int obl){
  delta = delta(arma::find(arma::linspace(0,K,K+1)!=obk),arma::find(arma::linspace(0,K,K+1)!=obk));
  int bk = 0;
  int bl = 0;
  double bv = -std::numeric_limits<double>::infinity();
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  for(int k = 1; k < K; ++k) {
    for (int l = 0;l<k;++l){
      if(k == obl | l == obl){
        arma::mat new_ec = x_counts;
        arma::mat new_counts = counts;
        new_counts(l)=new_counts(k)+new_counts(l);
        new_counts(k)=0;
        new_ec.col(l) = new_ec.col(k)+new_ec.col(l);
        new_ec.row(l) = new_ec.row(k)+new_ec.row(l);
        List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
        delta(k,l)=icl(new_stats,k,l)-icl(old_stats,k,l);
        delta(l,k)=delta(k,l);
      }
      if(delta(k,l)>bv){
        bk=k;
        bl=l;
        bv=delta(k,l);
      }
    }
  }
  return MergeMat(bk,bl,bv,delta);
}

void Sbm::merge_update(int k,int l){
  cl(arma::find(cl==k))=arma::ones(counts(k),1)*l;
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  counts(l) = counts(k)+counts(l);
  counts    = counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));

  x_counts.col(l) = x_counts.col(k)+x_counts.col(l);
  x_counts.row(l) = x_counts.row(k)+x_counts.row(l);
  x_counts = x_counts(arma::find(arma::linspace(0,K-1,K)!=k),arma::find(arma::linspace(0,K-1,K)!=k));
  --K;
}
