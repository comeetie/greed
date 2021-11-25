// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModelEmission.h"
#include "DiagGmm.h"
using namespace Rcpp;



DiagGmm::DiagGmm(const arma::mat & Xi,S4 modeli,bool verb){
  
  if(Rcpp::traits::is_nan<REALSXP>(modeli.slot("beta"))||  as<arma::rowvec>(modeli.slot("mu")).has_nan()){
    model=clone(modeli);
  }else{
    model= modeli;
  }
  
  
  // data
  X  = Xi;
  
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
  verbose=verb;
}

void DiagGmm::set_cl(arma::uvec cl){
  // construct oberved stats 
  K = arma::max(cl)+1;
  for(int k=0;k<K;k++){
    regs.push_back(gmm_marginal_spherical(X.rows(arma::find(cl==k)),kappa,tau,beta,mu));
  }
}



List DiagGmm::get_obs_stats(){
  // return observed stats
  return clone(regs);
}

double DiagGmm::icl_emiss(const List & regs){
  // compute log(p(X|Z))
  int K = regs.size();
  double icl_emiss = 0;
  for (int k = 0;k<K;++k){
    List rk = regs[k];
    double cle = rk["log_evidence"];
    icl_emiss += cle;
  }
  return icl_emiss;
}

double DiagGmm::icl_emiss(const List & regs,int oldcl,int newcl,bool dead_cluster){
  // compute log(p(X|Z)) but only for the 2 classes which haved changed (oldcl,newcl)
  List reg_newcl = regs[newcl];
  double icl_emiss = reg_newcl["log_evidence"];
  if(!dead_cluster){
    List reg_oldcl = regs[oldcl];
    icl_emiss += as<double>(reg_oldcl["log_evidence"]);
  }
  return icl_emiss;
}


arma::vec DiagGmm::delta_swap(const int i,arma::uvec &  cl,bool almost_dead_cluster,arma::uvec iclust,int K){
  
  int oldcl = cl(i);
  // extract current row 
  arma::rowvec xc = X.row(i);

  // initialize vecor of delta ICL
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  // old stats
  List new_regs = List(K);
  int k = 0;
  // for each possible move
  for(arma::uword j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    if(k!=oldcl){
      // construct new stats
      new_regs[k]=gmm_marginal_spherical_add1(regs[k],xc,kappa,tau,beta,mu);
      if(!almost_dead_cluster){
        new_regs[oldcl]=gmm_marginal_spherical_del1(regs[oldcl],xc,kappa,tau,beta,mu);
      }
      delta(k)=icl_emiss(new_regs,oldcl,k,almost_dead_cluster)-icl_emiss(regs,oldcl,k,false);
    }
  }

  return delta;
}


void DiagGmm::swap_update(const int i,arma::uvec & cl,bool dead_cluster,const int newcl){
  // a swap is done !
  int oldcl = cl(i);
  // current row
  arma::rowvec xc = X.row(i);
  // update regs
  // switch row x from cluster oldcl to k
  regs[newcl]=gmm_marginal_spherical_add1(regs[newcl],xc,kappa,tau,beta,mu);
  if(!dead_cluster){
    regs[oldcl]=gmm_marginal_spherical_del1(regs[oldcl],xc,kappa,tau,beta,mu);
  }else{
    // remove from regs
    regs.erase(oldcl);
    // upate K
    --K;
  }
  

}


double DiagGmm::delta_merge(int k, int l){
  
  // for each possible merge
  List new_regs = List(K);
  new_regs[l] = gmm_marginal_spherical_merge(regs[k],regs[l],kappa,tau,beta,mu);
  // delta
  double delta=icl_emiss(new_regs,k,l,true)-icl_emiss(regs,k,l,false);
  return delta;
}



// Merge update between k and l
void DiagGmm::merge_update(int k,int l){
  // update regs
  regs[l] = gmm_marginal_spherical_merge(regs[k],regs[l],kappa,tau,beta,mu);
  regs.erase(k);
  // update K
  --K;
}
