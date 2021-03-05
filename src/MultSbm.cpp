// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "MultSbm.h"
using namespace Rcpp;



MultSbm::MultSbm(const arma::cube& xp,double alphai,arma::vec& clt,bool verb){
  alpha = alphai;
  x  = xp;
  N  = x.n_rows;
  M = x.n_slices;
  set_cl(clt);
  verbose=verb;
  // cst to add 
  cst = 0;
}

void MultSbm::set_cl(arma::vec clt){
  cl = clt;
  K = arma::max(cl)+1;
  x_counts = gsum_cube(cl,x,K);
  counts = count(cl,K);
}

double MultSbm::icl_emiss(const List & obs_stats){

  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::mat edges_counts =as<arma::cube>(obs_stats["x_counts"]);
  arma::mat matcount = counts*counts.t();
  // lets go


  double icl_emiss = 0;
  return icl_emiss+cst;
}

double MultSbm::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  
  return 0;
}





List MultSbm::get_obs_stats(){
  return List::create(Named("counts", counts), Named("x_counts", x_counts));
}

arma::mat MultSbm::delta_swap(int i,arma::uvec iclust){
  arma::vec self=x.tube(i,i);
  int oldcl = cl(i);
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  int k = 0;
  // for each possible move
  for(int j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    if(k!=oldcl){
      arma::mat new_ec = x_counts;
      arma::vec new_counts = update_count(counts,oldcl,k);
      
      // todo
      List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
      delta(k)=icl(new_stats,oldcl,k)-icl(old_stats,oldcl,k);
    }
    
  }
  return delta;
}


void MultSbm::swap_update(const int i,const int newcl){
 
//todo

}


double MultSbm::delta_merge(int k, int l){
  
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  
  arma::mat new_ec = x_counts;
  arma::mat new_counts = counts;
  new_counts(l)=new_counts(k)+new_counts(l);
  new_counts(k)=0;
  // todo
  List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
  double delta=icl(new_stats,k,l)-icl(old_stats,k,l);
  
  return delta;
}


double MultSbm::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
  return 0;
}
