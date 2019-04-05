// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "Mreg.h"
using namespace Rcpp;


Mreg::Mreg(const arma::mat & Xi,const arma::colvec & yi, int Ki,double alphai,double regi, double a0i, double b0i,bool verb){
  // dirichlet prior parameter on proportion
  alpha = alphai;
  // dirichlet prior parameter on proportion
  reg = regi;
  a0 = a0i;
  b0 = b0i;
  // data
  X  = Xi;
  y  = yi;
  // Number of individuals
  N  = X.n_rows;
  // First value for K
  K  = Ki;
  // sample Z
  cl = as<arma::vec>(sample(K,N,true));
  cl = cl -1;
  // construct oberved stats 
  for(int k=0;k<K;k++){
    regs.push_back(lm_post(X.rows(arma::find(cl==k)),y.elem(arma::find(cl==k)),reg,a0,b0));
  }
  // counts : number of row in each cluster
  counts = count(cl,K);
  // TODO : add a filed to store the icl const ?
  verbose=verb;
}

Mreg::Mreg(const arma::mat & Xi,const arma::colvec & yi, int Ki,double alphai,double regi, double a0i, double b0i,arma::vec& cli,bool verb){
  // dirichlet prior parameter on proportion
  alpha = alphai;
  // dirichlet prior parameter on proportion
  reg = regi;
  a0 = a0i;
  b0 = b0i;
  // data
  X  = Xi;
  y  = yi;
  // Number of individuals
  N  = X.n_rows;
  // First value for K
  K  = Ki;
  cl= cli;
  // construct oberved stats 
  for(int k=0;k<K;k++){
    regs.push_back(lm_post(X.rows(arma::find(cl==k)),y.elem(arma::find(cl==k)),reg,a0,b0));
  }
  // counts : number of row in each cluster
  counts = count(cl,K);
  // TODO : add a filed to store the icl const ?
  verbose=verb;
}


List Mreg::get_obs_stats(){
  // return observed stats
  return List::create(Named("counts", counts), Named("regs", regs));
}

double Mreg::icl_emiss(const List & obs_stats){
  // compute log(p(X|Z))
  List regs =as<List>(obs_stats["regs"]);
  int K = regs.size();
  double icl_emiss = 0;
  for (int k = 0;k<K;++k){
    List rk = regs[k];
    double cle = rk["log_evidence"];
    icl_emiss += cle;
  }
  return icl_emiss;
}

double Mreg::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  // compute log(p(X|Z)) but only for the 2 classes which haved changed (oldcl,newcl)
  List regs =as<List>(obs_stats["regs"]);
  int K = regs.size();
  arma::vec counts = as<arma::vec>(obs_stats["counts"]);
  double icl_emiss = 0;
  for (int k = 0;k<K;++k){
    // only for the 2 classes which haved changed (oldcl,newcl) and not empty
    if((k==oldcl || k ==newcl) && counts(k)>0){
      // compute log(p(X|Z))
      List rk = regs[k];
      double cle = rk["log_evidence"];
      icl_emiss += cle;
    }
  }
  return icl_emiss;
}


arma::mat Mreg::delta_swap(int i){
  
  // old cluster
  int oldcl = cl(i);
  
  // extract current row 
  arma::mat xc = X.row(i);

  arma::colvec yc(1); 
  yc = yc.ones()*y(i);
  // initialize vecor of delta ICL
  arma::vec delta(K);
  delta.fill(0);
  // old stats
  List old_stats = List::create(Named("counts", counts), Named("regs", regs));
  

  // for each possible move
  for(int k = 0; k < K; ++k) {
    if(k!=oldcl){
      // construct new stats
      List new_regs = clone(regs);
      // switch row x from cluster oldcl to k
      
      
      List regk = new_regs[k];
      
      new_regs[k]=lm_post_add(regk,xc,yc,reg,a0,b0);
      
      List regold = new_regs[oldcl];
      
      new_regs[oldcl]=lm_post_del(regold,xc,yc,reg,a0,b0);
      // update cluster counts
      
 
      arma::vec new_counts = update_count(counts,oldcl,k);
      // new stats and delta
      

      List new_stats = List::create(Named("counts", new_counts), Named("regs", new_regs));
      delta(k)=icl(new_stats,oldcl,k)-icl(old_stats,oldcl,k);
    }
  }
  
  return delta;
}


void Mreg::swap_update(int i,int newcl){
  // a swap is done !
  // old cluster
  int oldcl = cl(i);
  // current row
  

  arma::mat x = X.row(i);
  arma::colvec yc(1); 
  yc(0) = y(i);
  // update regs
  List new_regs = clone(regs);
  // switch row x from cluster oldcl to k
  new_regs[newcl]=lm_post_add(new_regs[newcl],x,yc,reg,a0,b0);
  new_regs[oldcl]=lm_post_del(new_regs[oldcl],x,yc,reg,a0,b0);
  // update counts
  arma::mat new_counts = update_count(counts,oldcl,newcl);
  // update cl
  cl(i)=newcl;
  // if a cluster is dead
  if(new_counts(oldcl)==0){
    // remove from counts
    counts = new_counts.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    // remove from regs
    IntegerVector idx = seq_len(regs.length()) - 1;
    regs= new_regs[idx!=oldcl];
    // update cl to take into account de dead cluster
    cl.elem(arma::find(cl>oldcl))=cl.elem(arma::find(cl>oldcl))-1;
    // upate K
    --K;
  }else{
    // just update
    counts=new_counts;
    regs=new_regs;
  }
  

}


MergeMat Mreg::delta_merge(){
  // inititalize delta merge matrix
  arma::mat delta(K,K);
  delta.fill(0);
  // index to store current best merge
  int bk = 0;
  int bl = 0;
  // initialize bv found to -infty
  double bv = -std::numeric_limits<double>::infinity();
  // store cuurent stats
  List old_stats = List::create(Named("counts", counts), Named("regs", regs));
  for(int k = 1; k < K; ++k) {
    for (int l = 0;l<k;++l){
      // for each possible merge
      List new_regs = clone(regs);
      arma::mat new_counts = counts;
      // counts after merge
      new_counts(l) = new_counts(k)+new_counts(l);
      new_counts(k) = 0;
      // x_counts after merge on l
      // row/col k will not be taken into account since counts(k)==0
      new_regs[l] = lm_post_merge(regs[k],regs[l],reg,a0,b0);
      
      List new_stats = List::create(Named("counts", new_counts), Named("regs", new_regs));
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

MergeMat Mreg::delta_merge(arma::mat delta, int obk, int obl){
  // optimized version to compute only new values of the merge mat
  delta = delta(arma::find(arma::linspace(0,K,K+1)!=obk),arma::find(arma::linspace(0,K,K+1)!=obk));
  int bk = 0;
  int bl = 0;
  double bv = -std::numeric_limits<double>::infinity();
  List old_stats = List::create(Named("counts", counts), Named("regs", x_counts));
  for(int k = 1; k < K; ++k) {
    for (int l = 0;l<k;++l){
      if(k == obl | l == obl){
        List new_regs = clone(regs);
        arma::mat new_counts = counts;
        new_counts(l) = new_counts(k)+new_counts(l);
        new_counts(k) = 0;
        new_regs[l] = lm_post_merge(new_regs[k],new_regs[l],reg,a0,b0);
        List new_stats = List::create(Named("counts", new_counts), Named("regs", new_regs));
        delta(k,l)=icl(new_stats,k,l)-icl(old_stats,k,l);
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

// Merge update between k and l
void Mreg::merge_update(int k,int l){
  // update cl
  cl(arma::find(cl==k))=arma::ones(counts(k),1)*l;
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  // update counts
  counts(l) = counts(k)+counts(l);
  counts    = counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));
  // update x_counts
  regs[l] = lm_post_merge(regs[k],regs[l],reg,a0,b0);
  IntegerVector idx = seq_len(regs.length()) - 1;
  regs = regs[idx!=k];
  // update K
  --K;
}
