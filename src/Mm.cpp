// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "Mm.h"
using namespace Rcpp;


Mm::Mm(arma::sp_mat& xp,int Ki,double alphai,double betai){
  alpha = alphai;
  beta = betai;
  x  = xp;
  xt = x.t();
  N  = x.n_rows;
  K  = Ki;
  cl = as<arma::vec>(sample(K,N,true));
  cl = cl -1;
  x_counts = gsum_mm(cl,x,K);
  counts = count(cl,K);
}

Mm::Mm(arma::sp_mat& xp,int Ki,double alphai,double betai,arma::vec& clt){
  alpha = alphai;
  beta = betai;
  x  = xp;
  xt = x.t();
  N  = x.n_rows;
  K  = Ki;
  cl = clt;
  x_counts = gsum_mm(cl,x,K);
  counts = count(cl,K);
}


List Mm::get_obs_stats(){
  return List::create(Named("counts", counts), Named("x_counts", x_counts));
}

double Mm::icl_emiss(const List & obs_stats){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::mat x_counts =as<arma::vec>(obs_stats["x_counts"]);
  double icl_emiss=K*lgamma(x_counts.n_rows*beta)-K*x_counts.n_rows*lgamma(beta)+accu(lgamma(x_counts+beta))-accu(lgamma(counts+x_counts.n_rows*beta));
  return icl_emiss;
}

double Mm::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::mat x_counts =as<arma::vec>(obs_stats["x_counts"]);
  double icl_emiss = 0;
  for (int d = 0;d<x_counts.n_rows;++d){
      icl_emiss += lgamma(x_counts.n_rows*beta)-x_counts.n_rows*lgamma(beta)+accu(lgamma(x_counts.row(newcl)+beta))-lgamma(counts(newcl)+x_counts.n_rows*beta);
    if(counts(oldcl)!=0){
      icl_emiss += lgamma(x_counts.n_rows*beta)-x_counts.n_rows*lgamma(beta)+accu(lgamma(x_counts.row(oldcl)+beta))-lgamma(counts(oldcl)+x_counts.n_rows*beta);
    }
  }
  return icl_emiss;
}


arma::mat Mm::delta_swap(int i){
  // Rcout << "--" << i << "--" << std::endl;;
  int oldcl = cl(i);
  arma::sp_mat crow = xt.col(i).t();
  // Rcout << crow << std::endl;;
  arma::mat delta(K,1);
  delta.fill(0);
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  for(int k = 0; k < K; ++k) {
    if(k!=oldcl){
      arma::mat new_ec = x_counts;
      new_ec.row(k) = new_ec.row(k)+crow;
      new_ec.row(oldcl) = new_ec.row(oldcl)-crow;
      arma::mat new_counts = update_count(counts,oldcl,k);
      List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
      delta(k,0)=icl(old_stats,oldcl,k)-icl(new_stats,oldcl,k);
    }
  }
  return delta;
}


void Mm::swap_update(int i,int newcl){
  // Rcout << "swap !!" << i << newcl << std::endl;;
  int oldcl = cl(i);
  arma::sp_mat crow = xt.col(i).t();
  arma::mat new_ec = x_counts;
  
  new_ec.row(newcl) = new_ec.row(newcl)+crow;
  new_ec.row(oldcl) = new_ec.row(oldcl)-crow;
  
  arma::mat new_counts = update_count(counts,oldcl,newcl);
  cl(i)=newcl;
  
  if(new_counts(oldcl)==0){
    counts = new_counts.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    x_counts = new_ec(arma::find(arma::linspace(0,K-1,K)!=oldcl),arma::find(arma::linspace(0,K-1,K)<0));
    cl.elem(arma::find(cl>oldcl))=cl.elem(arma::find(cl>oldcl))-1;
    --K;
  }else{
    counts=new_counts;
    x_counts=new_ec;
  }

}


MergeMat Mm::delta_merge(){
  arma::mat delta(K,K);
  delta.fill(0);
  int bk = 0;
  int bl = 0;
  double bv = -std::numeric_limits<double>::infinity();
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  for(int k = 1; k < K; ++k) {
    for (int l = 0;l<k;++l){
      arma::mat new_ec = x_counts;
      arma::mat new_counts = counts;
      new_counts(l)=new_counts(k)+new_counts(l);
      new_counts(k)=0;
      new_ec.row(k) = new_ec.row(k)-new_ec.row(k);
      new_ec.row(l) = new_ec.row(l)+new_ec.row(k);
      List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
      delta(k,l)=icl(old_stats,k,l)-icl(new_stats,k,l);
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

MergeMat Mm::delta_merge(arma::mat delta, int obk, int obl){
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
        new_ec.row(k) = new_ec.row(k)-new_ec.row(k);
        new_ec.row(l) = new_ec.row(l)+new_ec.row(k);
        List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
        delta(k,l)=icl(old_stats,k,l)-icl(new_stats,k,l);
        delta(l,k)=delta(k,l);
        if(delta(k,l)>bv){
          bk=k;
          bl=l;
          bv=delta(k,l);
        }
      }
    }
  }
  return MergeMat(bk,bl,bv,delta);
}

void Mm::merge_update(int k,int l){
  cl(arma::find(cl==k))=arma::ones(counts(k),1)*l;
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  
  counts(l) = counts(k)+counts(l);
  
  

  counts    = counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));
  
  x_counts.row(l) = x_counts.row(k)+x_counts.row(l);
  
  x_counts = x_counts.rows(arma::find(arma::linspace(0,K-1,K)!=k));

  --K;
}
