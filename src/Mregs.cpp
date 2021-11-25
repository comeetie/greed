// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "Mregs.h"
using namespace Rcpp;

Mregs::Mregs(const arma::mat & Xi,const arma::mat & Yi,S4 modeli,bool verb){
  // dirichlet prior parameter on proportion
  
  S0 = as<arma::mat>(modeli.slot("epsilon"));
  if(Rcpp::traits::is_nan<REALSXP>(modeli.slot("N0")) ||  S0.has_nan()){
    model = clone(modeli);
  }else{
    model = modeli;
  }
  
  // data
  X  = Xi;
  Y  = Yi;
  // Number of individuals
  N  = X.n_rows;
  
  
  M  = inv_sympd(X.t()*X)*X.t()*Y;
  M.zeros();
  double beta = model.slot("tau");
  Kp  =  beta*X.t()*X/N;
  
  
  arma::mat R  = Y-X*M;
  arma::mat RR = cov(R);

  if(S0.has_nan()){
    S0 = RR;
    S0.zeros();
    S0.diag() = 0.1*RR.diag();
    model.slot("epsilon") = S0;
  }
  
  if(Rcpp::traits::is_nan<REALSXP>(model.slot("N0"))){
    N0 = Y.n_cols; 
    model.slot("N0")=N0;
  }else{
    N0 = model.slot("N0");
  }
  
  
  // TODO : add a filed to store the icl const ?
  verbose=verb;
}




void Mregs::set_cl(arma::uvec cl){
  K = arma::max(cl)+1;
  // construct oberved stats 
  for(int k=0;k<K;k++){
    regs.push_back(mvlm_post_comp(X.rows(arma::find(cl==k)),Y.rows(arma::find(cl==k)),M,Kp,S0,N0));
  }

}

List Mregs::get_obs_stats(){
  // return observed stats
  return clone(regs);
}

double Mregs::icl_emiss(const List & regs){
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

double Mregs::icl_emiss(const List & regs,int oldcl,int newcl, bool dead_cluster){
  // compute log(p(X|Z)) but only for the 2 classes which haved changed (oldcl,newcl)
  List rk = regs[newcl];
  double icl_emiss = rk["log_evidence"];
  if(!dead_cluster){
    rk = regs[oldcl];
    icl_emiss += as<double>(rk["log_evidence"]);
  }
  
  return icl_emiss;
}


arma::vec Mregs::delta_swap(const int i,arma::uvec & cl, bool almost_dead_cluster, arma::uvec iclust, int K){
  int oldcl = cl(i);
  // extract current row 
  arma::rowvec xc = X.row(i);
  
  arma::rowvec yc = Y.row(i);
  // initialize vecor of delta ICL
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  // old stats
  List old_stats = regs;
  
  int k = 0;
  // for each possible move
  for(arma::uword j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    if(k!=oldcl){
      // construct new stats
      List new_regs = List(K);
      // switch row x from cluster oldcl to k
      
      new_regs[k]=mvlm_post_add1_comp(regs[k],xc,yc,M,Kp,S0,N0);
      
      new_regs[oldcl]=mvlm_post_del1_comp(regs[oldcl],xc,yc,M,Kp,S0,N0);
    
      // new stats and delta
      
      
      delta(k)=icl_emiss(new_regs,oldcl,k,almost_dead_cluster)-icl_emiss(old_stats,oldcl,k,false);
    }
  }
  
  return delta;
}


void Mregs::swap_update(const int i,arma::uvec &  cl,bool dead_cluster,const int newcl){
  // a swap is done !
  // old cluster
  int oldcl = cl(i);
  // current row
  
  
  arma::rowvec xc = X.row(i);
  arma::rowvec yc = Y.row(i);
  // switch row x from cluster oldcl to k
  regs[newcl]=mvlm_post_add1_comp(regs[newcl],xc,yc,M,Kp,S0,N0);
  regs[oldcl]=mvlm_post_del1_comp(regs[oldcl],xc,yc,M,Kp,S0,N0);
  if(dead_cluster){
    regs.erase(oldcl);
    // upate K
    --K;
  }
  
}


double Mregs::delta_merge(int k, int l){
  
  List new_regs = List(K);
  new_regs[l] = mvlm_post_merge_comp(regs[k],regs[l],M,Kp,S0,N0);
  
  double delta=icl_emiss(new_regs,k,l,true)-icl_emiss(regs,k,l,false);
  return delta;
}



// Merge update between k and l
void Mregs::merge_update(int k,int l){
  regs[l] = mvlm_post_merge_comp(regs[k],regs[l],M,Kp,S0,N0);
  regs.erase(k);
  // update K
  --K;
}