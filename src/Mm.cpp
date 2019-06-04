// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "Mm.h"
using namespace Rcpp;


Mm::Mm(arma::sp_mat& xp,int Ki,double alphai,double betai,bool verb){
  // dirichlet prior parameter on proportion
  alpha = alphai;
  // dirichlet prior parameter on proportion
  beta = betai;
  // data
  x  = xp;
  // store also transpose for fast colum acces
  xt = x.t();
  // Number of individuals
  N  = x.n_rows;
  // First value for K
  K  = Ki;
  // sample Z
  cl = as<arma::vec>(sample(K,N,true));
  cl = cl -1;
  // construct oberved stats 
  // x_counts : col sums for each cluster
  x_counts = gsum_mm(cl,x,K);
  // counts : number of row in each cluster
  counts = count(cl,K);
  // TODO : add a filed to store the icl const ?
  verbose=verb;
}

Mm::Mm(arma::sp_mat& xp,int Ki,double alphai,double betai,arma::vec& clt,bool verb){
  // dirichlet prior parameter on proportion
  alpha = alphai;
  // dirichlet prior parameter on proportion
  beta = betai;
  // data
  x  = xp;
  // store also transpose for fast colum acces
  xt = x.t();
  // Number of individuals
  N  = x.n_rows;
  // First value for K
  K  = Ki;
  // init Z
  cl = clt;
  // construct oberved stats 
  // x_counts : col sums for each cluster
  x_counts = gsum_mm(cl,x,K);
  // counts : number of row in each cluster
  counts = count(cl,K);
  // TODO : add a filed to store the icl const ?
  verbose=verb;
}


List Mm::get_obs_stats(){
  // return observed stats
  return List::create(Named("counts", counts), Named("x_counts", x_counts));
}

double Mm::icl_emiss(const List & obs_stats){
  // compute log(p(X|Z))
  arma::mat x_counts =as<arma::mat>(obs_stats["x_counts"]);
  int K = x_counts.n_rows;
  int d = x_counts.n_cols;
  double icl_emiss = 0;
  for (int k = 0;k<K;++k){
    // B(X+beta)/B(beta)
    icl_emiss += lgamma(d*beta)+arma::accu(lgamma(x_counts.row(k)+beta))-d*lgamma(beta)-lgamma(arma::accu(x_counts.row(k)+beta));
  }
  return icl_emiss;
}

double Mm::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  // compute log(p(X|Z)) but only for the 2 classes which haved changed (oldcl,newcl)
  arma::mat x_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  double icl_emiss = 0;
  int d = x_counts.n_cols;
  for (int k = 0;k<x_counts.n_rows;++k){
    // only for the 2 classes which haved changed (oldcl,newcl) and not empty
    if((k==oldcl || k ==newcl) && counts(k)>0){
      // compute log(p(X|Z))
      icl_emiss += lgamma(d*beta)+arma::accu(lgamma(x_counts.row(k)+beta))-d*lgamma(beta)-lgamma(arma::accu(x_counts.row(k)+beta));
    }
  }
  return icl_emiss;
}


arma::mat Mm::delta_swap(int i){
  
  // old cluster
  int oldcl = cl(i);
  
  // extract current row 
  arma::sp_mat crow = xt.col(i).t();

  // initialize vecor of delta ICL
  arma::vec delta(K);
  delta.fill(0);
  // old stats
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  
  // for each possible move
  for(int k = 0; k < K; ++k) {
    if(k!=oldcl){
      // construct new stats
      arma::mat new_ec = x_counts;
      // siwtch current row
      new_ec.row(k) = new_ec.row(k)+crow;
      new_ec.row(oldcl) = new_ec.row(oldcl)-crow;
      // update cluster counts
      arma::vec new_counts = update_count(counts,oldcl,k);
      // new stats and delta
      List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
      delta(k)=icl(new_stats,oldcl,k)-icl(old_stats,oldcl,k);
    }
  }
  
  return delta;
}


void Mm::swap_update(int i,int newcl){
  // a swap is done !
  // old cluster
  int oldcl = cl(i);
  // current row
  arma::sp_mat crow = xt.col(i).t();
  // update x_counts
  arma::mat new_ec = x_counts;
  new_ec.row(newcl) = new_ec.row(newcl)+crow;
  new_ec.row(oldcl) = new_ec.row(oldcl)-crow;
  // update counts
  arma::mat new_counts = update_count(counts,oldcl,newcl);
  // update cl
  cl(i)=newcl;
  // if a cluster is dead
  if(new_counts(oldcl)==0){
    // remove from counts
    counts = new_counts.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    // remove from x_counts
    x_counts = new_ec.rows(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    // update cl to take into account de dead cluster
    cl.elem(arma::find(cl>oldcl))=cl.elem(arma::find(cl>oldcl))-1;
    // upate K
    --K;
  }else{
    // just update
    counts=new_counts;
    x_counts=new_ec;
  }
  

}


MergeMat Mm::delta_merge(){
  // inititalize delta merge matrix
  arma::mat delta(K,K);
  delta.fill(0);
  // index to store current best merge
  int bk = 0;
  int bl = 0;
  // initialize bv found to -infty
  double bv = -std::numeric_limits<double>::infinity();
  // store cuurent stats
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  for(int k = 1; k < K; ++k) {
    for (int l = 0;l<k;++l){
      // for each possible merge
      arma::mat new_ec = x_counts;
      arma::mat new_counts = counts;
      // counts after merge
      new_counts(l) = new_counts(k)+new_counts(l);
      new_counts(k) = 0;
      // x_counts after merge on l
      // row/col k will not be taken into account since counts(k)==0
      new_ec.row(l) = new_ec.row(l)+new_ec.row(k);
      List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
      // delta
      delta(k,l)=icl(new_stats,k,l)-icl(old_stats,k,l);
      delta(l,k)=delta(k,l);
      // best merge ?
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
  // optimized version to compute only new values of the merge mat
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
        new_counts(l) = new_counts(k)+new_counts(l);
        new_counts(k) = 0;
        new_ec.row(l) = new_ec.row(l)+new_ec.row(k);
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

// Merge update between k and l
void Mm::merge_update(int k,int l){
  // update cl
  cl(arma::find(cl==k))=arma::ones(counts(k),1)*l;
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  // update counts
  counts(l) = counts(k)+counts(l);
  counts    = counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));
  // update x_counts
  x_counts.row(l) = x_counts.row(k)+x_counts.row(l);
  x_counts = x_counts.rows(arma::find(arma::linspace(0,K-1,K)!=k));
  // update K
  --K;
}
