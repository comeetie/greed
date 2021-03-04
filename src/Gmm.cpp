// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "Gmm.h"
using namespace Rcpp;



Gmm::Gmm(const arma::mat & Xi,double alphai,double taui,int N0i, const arma::mat epsiloni, const arma::rowvec mui, arma::vec& cli,bool verb){
  // dirichlet prior parameter on proportion
  alpha = alphai;

  // data
  X  = Xi;
  // Number of individuals
  N  = X.n_rows;

  tau = taui;
  N0 = N0i;
  epsilon = epsiloni;
  mu = mui;
  set_cl(cli);
  S = cov(X);
  //normfact = -N/2*log(det(S));
  List normref = as<List>(gmm_marginal_eb(X,tau,N0,epsilon,mu));
  normfact = normref["log_evidence"];

  // TODO : add a filed to store the icl const ?
  verbose=verb;
}

void Gmm::set_cl(arma::vec cli){
  // construct oberved stats 
  cl = cli;
  K = arma::max(cl)+1;
  for(int k=0;k<K;k++){
    regs.push_back(gmm_marginal_eb(X.rows(arma::find(cl==k)),tau,N0,epsilon,mu));
  }
  // counts : number of row in each cluster
  counts = count(cl,K);
}



List Gmm::get_obs_stats(){
  // return observed stats
  return List::create(Named("counts", counts), Named("regs", regs));
}

double Gmm::icl_emiss(const List & obs_stats){
  // compute log(p(X|Z))
  List regs =as<List>(obs_stats["regs"]);
  int K = regs.size();
  double icl_emiss = 0;
  for (int k = 0;k<K;++k){
    List rk = regs[k];
    double cle = rk["log_evidence"];
    icl_emiss += cle;
  }
  return icl_emiss-normfact;
}

double Gmm::icl_emiss(const List & obs_stats,int oldcl,int newcl){
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
  return icl_emiss-normfact;
}


arma::mat Gmm::delta_swap(int i,arma::uvec iclust){
  
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
      
      new_regs[k]=gmm_marginal_add1_eb(regk,xc,tau,N0,epsilon,mu);
      
      List regold = new_regs[oldcl];
      
      new_regs[oldcl]=gmm_marginal_del1_eb(regold,xc,tau,N0,epsilon,mu);
      // update cluster counts
      
 
      arma::vec new_counts = update_count(counts,oldcl,k);
      // new stats and delta
      

      List new_stats = List::create(Named("counts", new_counts), Named("regs", new_regs));
      delta(k)=icl(new_stats,oldcl,k)-icl(old_stats,oldcl,k);
    }
  }
  
  return delta;
}


void Gmm::swap_update(int i,int newcl){
  // a swap is done !
  // old cluster
  int oldcl = cl(i);
  // current row
  

  arma::rowvec xc = X.row(i);
  // update regs
  List new_regs = clone(regs);
  // switch row x from cluster oldcl to k
  new_regs[newcl]=gmm_marginal_add1_eb(new_regs[newcl],xc,tau,N0,epsilon,mu);
  new_regs[oldcl]=gmm_marginal_del1_eb(new_regs[oldcl],xc,tau,N0,epsilon,mu);
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


double Gmm::delta_merge(int k, int l){
  
  List old_stats = List::create(Named("counts", counts), Named("regs", regs));
  // for each possible merge
  List new_regs = clone(regs);
  arma::mat new_counts = counts;
  // counts after merge
  new_counts(l) = new_counts(k)+new_counts(l);
  new_counts(k) = 0;
  // x_counts after merge on l
  // row/col k will not be taken into account since counts(k)==0
  new_regs[l] = gmm_marginal_merge_eb(regs[k],regs[l],tau,N0,epsilon,mu);
  
  List new_stats = List::create(Named("counts", new_counts), Named("regs", new_regs));
  // delta
  double delta=icl(new_stats,k,l)-icl(old_stats,k,l);
  return delta;
}



// Merge update between k and l
void Gmm::merge_update(int k,int l){
  // update cl
  cl(arma::find(cl==k))=arma::ones(counts(k),1)*l;
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  // update counts
  counts(l) = counts(k)+counts(l);
  counts    = counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));
  // update x_counts
  regs = clone(regs);
  regs[l] = gmm_marginal_merge_eb(regs[k],regs[l],tau,N0,epsilon,mu);
  IntegerVector idx = seq_len(regs.length()) - 1;
  regs = regs[idx!=k];
  // update K
  --K;
}
