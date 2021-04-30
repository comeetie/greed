// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "Mvmregcomp.h"
using namespace Rcpp;



Mvmregcomp::Mvmregcomp(const arma::mat & Xi,const arma::mat & Yi, S4 modeli,arma::vec& cli,bool verb){
  // dirichlet prior parameter on proportion
  model = modeli;
  alpha = model.slot("alpha");

  // data
  X  = Xi;
  Y  = Yi;
  // Number of individuals
  N  = X.n_rows;


  M  = inv_sympd(X.t()*X)*X.t()*Y;
  double beta = model.slot("tau");
  Kp  =  beta*X.t()*X/N;
  

  arma::mat R  = Y-X*M;
  arma::mat RR = cov(R);
  S0 = as<arma::mat>(model.slot("epsilon"));
  if(S0.has_nan()){
    S0 = RR;
    S0.zeros();
    S0.diag() = 0.01*RR.diag();
    model.slot("epsilon") = S0;
  }
  
  if(Rcpp::traits::is_nan<REALSXP>(model.slot("N0"))){
    N0 = Y.n_cols; 
    model.slot("N0")=N0;
  }else{
    N0 = model.slot("N0");
  }
  
  set_cl(cli);

  // TODO : add a filed to store the icl const ?
  verbose=verb;
}


void Mvmregcomp::set_cl(arma::vec cli){
  cl = cli;
  K = arma::max(cl)+1;
  // construct oberved stats 
  for(int k=0;k<K;k++){
    regs.push_back(mvlm_post_comp(X.rows(arma::find(cl==k)),Y.rows(arma::find(cl==k)),M,Kp,S0,N0));
  }
  
  // counts : number of row in each cluster
  counts = count(cl,K);
}

List Mvmregcomp::get_obs_stats(){
  // return observed stats
  return List::create(Named("counts", counts), Named("regs", regs));
}

double Mvmregcomp::icl_emiss(const List & obs_stats){
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

double Mvmregcomp::icl_emiss(const List & obs_stats,int oldcl,int newcl){
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


arma::mat Mvmregcomp::delta_swap(int i,arma::uvec iclust){
  
  // old cluster
  int oldcl = cl(i);
  
  // extract current row 
  arma::rowvec xc = X.row(i);

  arma::rowvec yc = Y.row(i);
  // initialize vecor of delta ICL
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  // old stats
  List old_stats = List::create(Named("counts", counts), Named("regs", regs));
  
  int k = 0;
  // for each possible move
  for(arma::uword j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    if(k!=oldcl){
      // construct new stats
      List new_regs = clone(regs);
      // switch row x from cluster oldcl to k
      
      
      List regk = new_regs[k];
      
      new_regs[k]=mvlm_post_add1_comp(regk,xc,yc,M,Kp,S0,N0);
      
      List regold = new_regs[oldcl];
      
      new_regs[oldcl]=mvlm_post_del1_comp(regold,xc,yc,M,Kp,S0,N0);
      // update cluster counts
      
 
      arma::vec new_counts = update_count(counts,oldcl,k);
      // new stats and delta
      

      List new_stats = List::create(Named("counts", new_counts), Named("regs", new_regs));
      delta(k)=icl(new_stats,oldcl,k)-icl(old_stats,oldcl,k);
    }
  }
  
  return delta;
}


void Mvmregcomp::swap_update(int i,int newcl){
  // a swap is done !
  // old cluster
  int oldcl = cl(i);
  // current row
  

  arma::rowvec xc = X.row(i);
  arma::rowvec yc = Y.row(i);
  // update regs
  List new_regs = clone(regs);
  // switch row x from cluster oldcl to k
  new_regs[newcl]=mvlm_post_add1_comp(new_regs[newcl],xc,yc,M,Kp,S0,N0);
  new_regs[oldcl]=mvlm_post_del1_comp(new_regs[oldcl],xc,yc,M,Kp,S0,N0);
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


double Mvmregcomp::delta_merge(int k, int l){
  
  List old_stats = List::create(Named("counts", counts), Named("regs", regs));
  // for each possible merge
  List new_regs = clone(regs);
  arma::mat new_counts = counts;
  // counts after merge
  new_counts(l) = new_counts(k)+new_counts(l);
  new_counts(k) = 0;
  // x_counts after merge on l
  // row/col k will not be taken into account since counts(k)==0
  new_regs[l] = mvlm_post_merge_comp(regs[k],regs[l],M,Kp,S0,N0);
  
  List new_stats = List::create(Named("counts", new_counts), Named("regs", new_regs));
  // delta
  double delta=icl(new_stats,k,l)-icl(old_stats,k,l);
  return delta;
}



// Merge update between k and l
void Mvmregcomp::merge_update(int k,int l){
  // update cl
  cl(arma::find(cl==k))=arma::ones(counts(k),1)*l;
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  // update counts
  counts(l) = counts(k)+counts(l);
  counts    = counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));
  // update x_counts
  regs = clone(regs);
  regs[l] = mvlm_post_merge_comp(regs[k],regs[l],M,Kp,S0,N0);
  IntegerVector idx = seq_len(regs.length()) - 1;
  regs = regs[idx!=k];
  // update K
  --K;
}
