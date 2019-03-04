// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "Mmvp.h"
using namespace Rcpp;


Mmvp::Mmvp(arma::sp_mat& xp,int Ki,double alphai){
  alpha = alphai;
  x  = xp;
  xt = x.t();
  N  = x.n_rows;
  K  = Ki;
  cl = as<arma::vec>(sample(K,N,true));
  cl = cl -1;
  p  = sum_cols(x);
  x_counts = gsum_mm(cl,x,K);
  Rcout << "xcounts : " << x_counts << std::endl;
  counts = count(cl,K);
  Rcout << "counts : " << counts << std::endl;
}

Mmvp::Mmvp(arma::sp_mat& xp,int Ki,double alphai,arma::vec& clt){
  alpha = alphai;
  x  = xp;
  xt = x.t();
  N  = x.n_rows;
  K  = Ki;
  cl = clt;
  p  = sum_cols(x);
  x_counts = gsum_mm(cl,x,K);
  counts = count(cl,K);
}


List Mmvp::get_obs_stats(){
  return List::create(Named("counts", counts), Named("x_counts", x_counts));
}

double Mmvp::icl_emiss(const List & obs_stats){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::mat x_counts =as<arma::mat>(obs_stats["x_counts"]);
  int K = x_counts.n_rows;
  int d = x_counts.n_cols;
  double icl_emiss = 0;
  for (int k = 0;k<K;++k){
    icl_emiss = arma::accu(lgamma(x_counts.row(k)+1))-arma::accu((x_counts.row(k)+1)*log(counts+p));
  }
  icl_emiss = icl_emiss + K*arma::accu(log(p));
  return icl_emiss;
}

double Mmvp::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  double icl_emiss = 0;
  return icl_emiss;
}


arma::mat Mmvp::delta_swap(int i){
  // Rcout << "--" << i << "--" << std::endl;;
  int oldcl = cl(i);
  arma::sp_mat crow = xt.col(i).t();


  arma::mat delta(K,1);
  delta.fill(0);
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  for(int k = 0; k < K; ++k) {
    if(k!=oldcl){
      arma::mat new_ec = x_counts;
      new_ec.row(k) = new_ec.row(k)+crow;
      new_ec.row(oldcl) = new_ec.row(oldcl)-crow;
      arma::mat new_counts = update_count(counts,oldcl,k);
      if(new_counts(oldcl)==0){
        new_counts=new_counts(arma::find(arma::linspace(0,K-1,K)!=oldcl));
        new_ec=new_ec(arma::find(arma::linspace(0,K-1,K)!=oldcl),arma::find(arma::linspace(0,K-1,K)>=0));
        Rcout << "death " << new_ec << std::endl;
      }
      List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
      delta(k,0)=icl(old_stats)-icl(new_stats);
    }
  }
  if(oldcl==2){
    Rcout << "CL 3" << std::endl;
    Rcout << crow << std::endl;
    Rcout << delta << std::endl;
  }
  return delta;
}


void Mmvp::swap_update(int i,int newcl){
 //  Rcout << "swap !!" << i << newcl << std::endl;;
  // Rcout << counts << std::endl;
  int oldcl = cl(i);
  arma::sp_mat crow = xt.col(i).t();
  arma::mat new_ec = x_counts;
  
  new_ec.row(newcl) = new_ec.row(newcl)+crow;
  new_ec.row(oldcl) = new_ec.row(oldcl)-crow;
  
  arma::mat new_counts = update_count(counts,oldcl,newcl);
  cl(i)=newcl;

  if(new_counts(oldcl)==0){
    counts = new_counts.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    x_counts = new_ec(arma::find(arma::linspace(0,K-1,K)!=oldcl),arma::find(arma::linspace(0,K-1,K)>=0));
    cl.elem(arma::find(cl>oldcl))=cl.elem(arma::find(cl>oldcl))-1;
    --K;
  }else{
    counts=new_counts;
    x_counts=new_ec;
  }
  //Rcout << counts << std::endl;

}


MergeMat Mmvp::delta_merge(){
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

MergeMat Mmvp::delta_merge(arma::mat delta, int obk, int obl){
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

void Mmvp::merge_update(int k,int l){
  cl(arma::find(cl==k))=arma::ones(counts(k),1)*l;
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  
  counts(l) = counts(k)+counts(l);
  
  

  counts    = counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));
  
  x_counts.row(l) = x_counts.row(k)+x_counts.row(l);
  
  x_counts = x_counts.rows(arma::find(arma::linspace(0,K-1,K)!=k));

  --K;
}
