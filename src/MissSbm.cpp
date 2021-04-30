// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "MissSbm.h"
using namespace Rcpp;



MissSbm::MissSbm(arma::sp_mat& xp,arma::sp_mat& xpobs,S4 modeli,arma::vec& clt,bool verb){
  model = modeli;
  alpha = model.slot("alpha");
  a0 = model.slot("a0");
  b0 = model.slot("b0");
  List sampling_priors = as<List>(model.slot("sampling_priors"));
  a0obs = sampling_priors["a0obs"];
  b0obs = sampling_priors["b0obs"];
  x  = xp;
  xt = xp.t();
  xobs  = xpobs;
  xtobs = xpobs.t();
  N  = x.n_rows;
  set_cl(clt);
  verbose=verb;
}

void MissSbm::set_cl(arma::vec clt){
  cl = clt;
  K = arma::max(cl)+1;
  x_counts = gsum_mat(cl,x,K);
  x_counts_obs = gsum_mat(cl,xobs,K);
  counts = count(cl,K);
}


double MissSbm::icl_emiss(const List & obs_stats){
  return 0;
}

double MissSbm::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  return 0;
}


List MissSbm::get_obs_stats(){
  return List::create(Named("counts", counts), Named("x_counts", x_counts),Named("x_counts_obs", x_counts_obs));
}

arma::mat MissSbm::delta_swap(int i,arma::uvec iclust){
  int self=x(i,i);
  int selfobs=xobs(i,i);
  int oldcl = cl(i);
  arma::sp_mat col_sum = gsum_col(cl,x,i,K);
  col_sum(oldcl)=col_sum(oldcl)-self;
  
  arma::sp_mat col_sum_obs = gsum_col(cl,xobs,i,K);
  col_sum_obs(oldcl)=col_sum_obs(oldcl)-selfobs;
  
  arma::sp_mat row_sum = gsum_col(cl,xt,i,K);
  row_sum(oldcl)=row_sum(oldcl)-self;
  
  arma::sp_mat row_sum_obs = gsum_col(cl,xtobs,i,K);
  row_sum_obs(oldcl)=row_sum_obs(oldcl)-selfobs;
  
  
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts), Named("x_counts_obs", x_counts_obs));
  int k = 0;
  // for each possible move
  for(arma::uword j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    if(k!=oldcl){
      
      arma::mat new_ec = x_counts;
      new_ec.col(k) = new_ec.col(k)+col_sum;
      new_ec.row(k) = new_ec.row(k)+row_sum.t();
      new_ec.col(oldcl) = new_ec.col(oldcl)-col_sum;
      new_ec.row(oldcl) = new_ec.row(oldcl)-row_sum.t();
      new_ec(k,k)=new_ec(k,k)+self;
      new_ec(oldcl,oldcl)=new_ec(oldcl,oldcl)-self;
      
      
      arma::mat new_ec_obs = x_counts_obs;
      new_ec_obs.col(k) = new_ec_obs.col(k)+col_sum_obs;
      new_ec_obs.row(k) = new_ec_obs.row(k)+row_sum_obs.t();
      new_ec_obs.col(oldcl) = new_ec_obs.col(oldcl)-col_sum_obs;
      new_ec_obs.row(oldcl) = new_ec_obs.row(oldcl)-row_sum_obs.t();
      new_ec_obs(k,k)=new_ec_obs(k,k)+selfobs;
      new_ec_obs(oldcl,oldcl)=new_ec_obs(oldcl,oldcl)-selfobs;
      
      arma::vec new_counts = update_count(counts,oldcl,k);
      List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec), Named("x_counts_obs", new_ec_obs));
      delta(k)=icl(new_stats,oldcl,k)-icl(old_stats,oldcl,k);
    }
  }
  //Rcout << delta << std::endl;
  return delta;
}


void MissSbm::swap_update(const int i,const int newcl){
  int self=x(i,i);
  int selfobs=xobs(i,i);
  int oldcl = cl(i);
  arma::sp_mat col_sum = gsum_col(cl,x,i,K);
  col_sum(oldcl)=col_sum(oldcl)-self;
  
  arma::sp_mat col_sum_obs = gsum_col(cl,xobs,i,K);
  col_sum_obs(oldcl)=col_sum_obs(oldcl)-selfobs;
  
  arma::sp_mat row_sum = gsum_col(cl,xt,i,K);
  row_sum(oldcl)=row_sum(oldcl)-self;
  
  arma::sp_mat row_sum_obs = gsum_col(cl,xtobs,i,K);
  row_sum_obs(oldcl)=row_sum_obs(oldcl)-selfobs;
  
  arma::mat new_ec = x_counts;
  new_ec.col(newcl) = new_ec.col(newcl)+col_sum;
  new_ec.row(newcl) = new_ec.row(newcl)+row_sum.t();
  new_ec.col(oldcl) = new_ec.col(oldcl)-col_sum;
  new_ec.row(oldcl) = new_ec.row(oldcl)-row_sum.t();
  new_ec(newcl,newcl)=new_ec(newcl,newcl)+self;
  new_ec(oldcl,oldcl)=new_ec(oldcl,oldcl)-self;
  
  
  arma::mat new_ec_obs = x_counts_obs;
  new_ec_obs.col(newcl) = new_ec_obs.col(newcl)+col_sum_obs;
  new_ec_obs.row(newcl) = new_ec_obs.row(newcl)+row_sum_obs.t();
  new_ec_obs.col(oldcl) = new_ec_obs.col(oldcl)-col_sum_obs;
  new_ec_obs.row(oldcl) = new_ec_obs.row(oldcl)-row_sum_obs.t();
  new_ec_obs(newcl,newcl)=new_ec_obs(newcl,newcl)+selfobs;
  new_ec_obs(oldcl,oldcl)=new_ec_obs(oldcl,oldcl)-selfobs;
  arma::mat new_counts = update_count(counts,oldcl,newcl);
  cl(i)=newcl;
  if(new_counts(oldcl)==0){
    counts = new_counts.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    x_counts = new_ec(arma::find(arma::linspace(0,K-1,K)!=oldcl),arma::find(arma::linspace(0,K-1,K)!=oldcl));
    x_counts_obs = new_ec_obs(arma::find(arma::linspace(0,K-1,K)!=oldcl),arma::find(arma::linspace(0,K-1,K)!=oldcl));
    cl.elem(arma::find(cl>oldcl))=cl.elem(arma::find(cl>oldcl))-1;
    --K;
  }else{
    counts=new_counts;
    x_counts=new_ec;
    x_counts_obs=new_ec_obs;
  }

}


double MissSbm::delta_merge(int k, int l){
  
  arma::mat new_ec = x_counts;
  arma::mat new_ec_obs = x_counts_obs;
  arma::mat new_counts = counts;
  new_counts(l)=new_counts(k)+new_counts(l);
  new_counts(k)=0;
  new_ec.col(l) = new_ec.col(k)+new_ec.col(l);
  new_ec.row(l) = new_ec.row(k)+new_ec.row(l);
  new_ec_obs.col(l) = new_ec_obs.col(k)+new_ec_obs.col(l);
  new_ec_obs.row(l) = new_ec_obs.row(k)+new_ec_obs.row(l);
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts), Named("x_counts_obs", x_counts_obs));
  List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec), Named("x_counts_obs", new_ec_obs));
  double delta = icl(new_stats,k,l)-icl(old_stats,k,l);
  
  
  return delta;
}


void MissSbm::merge_update(int k,int l){
  cl(arma::find(cl==k))=arma::ones(counts(k),1)*l;
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  counts(l) = counts(k)+counts(l);
  counts    = counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));

  x_counts.col(l) = x_counts.col(k)+x_counts.col(l);
  x_counts.row(l) = x_counts.row(k)+x_counts.row(l);
  x_counts = x_counts(arma::find(arma::linspace(0,K-1,K)!=k),arma::find(arma::linspace(0,K-1,K)!=k));
  x_counts_obs.col(l) = x_counts_obs.col(k)+x_counts_obs.col(l);
  x_counts_obs.row(l) = x_counts_obs.row(k)+x_counts_obs.row(l);
  x_counts_obs = x_counts_obs(arma::find(arma::linspace(0,K-1,K)!=k),arma::find(arma::linspace(0,K-1,K)!=k));
  --K;
}

double MissSbm::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){

  return 0;
  
}



