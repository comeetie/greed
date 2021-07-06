// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModelEmission.h"
#include "GmmE.h"
using namespace Rcpp;



GmmE::GmmE(const arma::mat & Xi,S4 modeli,bool verb){
  model= modeli;
  // dirichlet prior parameter on proportion

  // data
  X  = Xi;
  // Number of individuals
  N  = X.n_rows;

  tau = model.slot("tau");
  if(Rcpp::traits::is_nan<REALSXP>(model.slot("N0"))){
    N0 = X.n_cols; 
    model.slot("N0")=N0;
  }else{
    N0 = model.slot("N0");
  }
  mu = as<arma::rowvec>(model.slot("mu"));
  if(mu.has_nan()){
    mu = arma::mean(X,0);
    model.slot("mu")=mu;
  }
  epsilon = as<arma::mat>(model.slot("epsilon"));
  if(epsilon.has_nan()){
    epsilon = 0.01*arma::diagmat(cov(X));
    model.slot("epsilon") = epsilon;
  }
  
  // TODO : add a filed to store the icl const ?
  verbose=verb;
}

void GmmE::set_cl(arma::vec cl){
  // construct oberved stats 
  K = arma::max(cl)+1;
  for(int k=0;k<K;k++){
    regs.push_back(gmm_marginal(X.rows(arma::find(cl==k)),tau,N0,epsilon,mu));
  }
  // counts : number of row in each cluster
  counts = count(cl,K);
}



List GmmE::get_obs_stats(){
  // return observed stats
  return List::create(Named("counts", counts), Named("regs", regs));
}

double GmmE::icl_emiss(const List & obs_stats){
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

double GmmE::icl_emiss(const List & obs_stats,int oldcl,int newcl){
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


arma::mat GmmE::delta_swap(int i,int K,Partition clp, arma::uvec iclust){
  
  int oldcl = clp.get(i);
  // extract current row 
  arma::rowvec xc = X.row(i);

  // initialize vecor of delta ICL
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  // old stats
  List old_stats = List::create(Named("counts", counts), Named("regs", regs));
  List new_regs = List(K);
  int k = 0;
  // for each possible move
  for(arma::uword j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    if(k!=oldcl){
      // construct new stats
      new_regs[k]=gmm_marginal_add1(regs[k],xc,tau,N0,epsilon,mu);
      new_regs[oldcl]=gmm_marginal_del1(regs[oldcl],xc,tau,N0,epsilon,mu);
      // update cluster counts
      arma::vec new_counts = update_count(counts,oldcl,k);
      // new stats and delta
      List new_stats = List::create(Named("counts", new_counts), Named("regs", new_regs));
      delta(k)=icl_emiss(new_stats,oldcl,k)-icl_emiss(old_stats,oldcl,k);
    }
  }

  return delta;
}


void GmmE::swap_update(int i,Partition clp,bool dead_cluster,const int newcl){
  // a swap is done !
  int oldcl = clp.get(i);
  // current row
  arma::rowvec xc = X.row(i);
  // update regs
  // switch row x from cluster oldcl to k
  regs[newcl]=gmm_marginal_add1(regs[newcl],xc,tau,N0,epsilon,mu);
  regs[oldcl]=gmm_marginal_del1(regs[oldcl],xc,tau,N0,epsilon,mu);
  // update counts
  counts = update_count(counts,oldcl,newcl);

  // if a cluster is dead
  if(counts(oldcl)==0){
    // remove from counts
    counts.shed_row(oldcl);
    // remove from regs
    regs.erase(oldcl);
    // upate K
    --K;
  }
  

}


double GmmE::delta_merge(int k, int l){
  
  List old_stats = List::create(Named("counts", counts), Named("regs", regs));
  // for each possible merge
  List new_regs = List(K);
  arma::mat new_counts = counts;
  // counts after merge
  new_counts(l) = new_counts(k)+new_counts(l);
  new_counts(k) = 0;
  // x_counts after merge on l
  // row/col k will not be taken into account since counts(k)==0
  new_regs[l] = gmm_marginal_merge(regs[k],regs[l],tau,N0,epsilon,mu);
  
  List new_stats = List::create(Named("counts", new_counts), Named("regs", new_regs));
  // delta
  double delta=icl_emiss(new_stats,k,l)-icl_emiss(old_stats,k,l);
  return delta;
}



// Merge update between k and l
void GmmE::merge_update(int k,int l){

  // update counts
  counts(l) = counts(k)+counts(l);
  counts.shed_row(k);
  // update x_counts
  regs[l] = gmm_marginal_merge(regs[k],regs[l],tau,N0,epsilon,mu);
  regs.erase(k);
  // update K
  --K;
}
