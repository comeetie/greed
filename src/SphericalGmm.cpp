// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "SphericalGmm.h"
using namespace Rcpp;



SphericalGmm::SphericalGmm(const arma::mat & Xi,S4 modeli, arma::vec& cli,bool verb){
  // dirichlet prior parameter on proportion
  model= modeli;
  alpha = model.slot("alpha");

  // data
  X  = Xi;
  // Number of individuals
  N  = X.n_rows;

  tau = model.slot("tau");
  kappa = model.slot("kappa");

  
  beta= model.slot("beta");
  if(Rcpp::traits::is_nan<REALSXP>(beta)){
    beta = 0.1*arma::mean(arma::var(X)); 
    model.slot("beta")=beta;
  }
  
  mu = as<arma::rowvec>(model.slot("mu"));
  if(mu.has_nan()){
    mu = arma::mean(X,0);
    model.slot("mu")=mu;
  }
  set_cl(cli);

  verbose=verb;
}

void SphericalGmm::set_cl(arma::vec cli){
  // construct oberved stats 
  cl = cli;
  K = arma::max(cl)+1;
  for(int k=0;k<K;k++){
    regs.push_back(gmm_marginal_spherical(X.rows(arma::find(cl==k)),kappa,tau,beta,mu));
  }
  // counts : number of row in each cluster
  counts = count(cl,K);
}



List SphericalGmm::get_obs_stats(){
  // return observed stats
  return List::create(Named("counts", counts), Named("regs", regs));
}

double SphericalGmm::icl_emiss(const List & obs_stats){
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

double SphericalGmm::icl_emiss(const List & obs_stats,int oldcl,int newcl){
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


arma::mat SphericalGmm::delta_swap(int i,arma::uvec iclust){
  
  // old cluster
  int oldcl = cl(i);
  
  // extract current row 
  arma::rowvec xc = X.row(i);

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
      
      new_regs[k]=gmm_marginal_spherical_add1(regk,xc,kappa,tau,beta,mu);
      
      List regold = new_regs[oldcl];
      
      new_regs[oldcl]=gmm_marginal_spherical_del1(regold,xc,kappa,tau,beta,mu);
      // update cluster counts
      
 
      arma::vec new_counts = update_count(counts,oldcl,k);
      // new stats and delta
      

      List new_stats = List::create(Named("counts", new_counts), Named("regs", new_regs));
      delta(k)=icl(new_stats,oldcl,k)-icl(old_stats,oldcl,k);
    }
  }
  
  return delta;
}


void SphericalGmm::swap_update(int i,int newcl){
  // a swap is done !
  // old cluster
  int oldcl = cl(i);
  // current row
  

  arma::rowvec xc = X.row(i);
  // update regs
  List new_regs = clone(regs);
  // switch row x from cluster oldcl to k
  new_regs[newcl]=gmm_marginal_spherical_add1(new_regs[newcl],xc,kappa,tau,beta,mu);
  new_regs[oldcl]=gmm_marginal_spherical_del1(new_regs[oldcl],xc,kappa,tau,beta,mu);
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


double SphericalGmm::delta_merge(int k, int l){
  
  List old_stats = List::create(Named("counts", counts), Named("regs", regs));
  // for each possible merge
  List new_regs = clone(regs);
  arma::mat new_counts = counts;
  // counts after merge
  new_counts(l) = new_counts(k)+new_counts(l);
  new_counts(k) = 0;
  // x_counts after merge on l
  // row/col k will not be taken into account since counts(k)==0
  new_regs[l] = gmm_marginal_spherical_merge(regs[k],regs[l],kappa,tau,beta,mu);
  
  List new_stats = List::create(Named("counts", new_counts), Named("regs", new_regs));
  // delta
  double delta=icl(new_stats,k,l)-icl(old_stats,k,l);
  return delta;
}



// Merge update between k and l
void SphericalGmm::merge_update(int k,int l){
  // update cl
  cl(arma::find(cl==k))=arma::ones(counts(k),1)*l;
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  // update counts
  counts(l) = counts(k)+counts(l);
  counts    = counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));
  // update x_counts
  regs = clone(regs);
  regs[l] = gmm_marginal_spherical_merge(regs[k],regs[l],kappa,tau,beta,mu);
  IntegerVector idx = seq_len(regs.length()) - 1;
  regs = regs[idx!=k];
  // update K
  --K;
}
