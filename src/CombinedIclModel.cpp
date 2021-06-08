// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "IclModel.h"
#include "CombinedIclModel.h"
using namespace Rcpp;



CombinedIclModel::CombinedIclModel(std::vector<IclModelEmission*> IclModelsi,S4 modeli, arma::vec cli, bool verb){

  model=modeli;
  alpha = model.slot("alpha");
  IclModels = IclModelsi;
  set_cl(cli);
  verbose=verb;
}

void CombinedIclModel::set_cl(arma::vec cli){

  cl = cli;
  N = cl.n_elem; 
  K = arma::max(cl)+1;
  counts = count(cl,K);
  for(int i=0;i<IclModels.size();i++){
    IclModelEmission * Mp = IclModels[i];
    Mp->set_cl(cl);
  }
}



List CombinedIclModel::get_obs_stats(){
  List obs_stats = List::create(Named("counts", counts));
  for(int i=0;i<IclModels.size();i++){
    IclModelEmission * Mp = IclModels[i];
    obs_stats.push_back(Mp->get_obs_stats());
  }
  //return observed stats
  return obs_stats;
}



double CombinedIclModel::icl_emiss(const List & obs_stats){
  double icl = 0;
  for(int i=0;i<IclModels.size();i++){
    IclModelEmission * Mp = IclModels[i];
    List cobs = as<List>(obs_stats[i+1]);
    icl += Mp->icl_emiss(cobs);
  }
  return icl;
}





double CombinedIclModel::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  double icl = 0;
  for(int i=0;i<IclModels.size();i++){
    IclModelEmission * Mp = IclModels[i];
    List cobs = as<List>(obs_stats[i+1]);
    icl += Mp->icl_emiss(cobs,oldcl,newcl);
  }
  return icl;
}



// main function to compute ICL from observed stats
double CombinedIclModel::icl_prop(arma::vec counts){
  // number of cluster
  int K = counts.n_elem;
  // log(p(Z))
  double icl_prop = lgamma(K*alpha)+arma::accu(lgamma(alpha+counts))-K*lgamma(alpha)-lgamma(arma::accu(counts+alpha));
  return icl_prop;
}

// main function to compute ICL from observed stats optimized version for computing delta which only invlove change in 2 clusters
double CombinedIclModel::icl_prop(arma::vec counts,int oldcl,int newcl){
  // compute the first part p(Z) from clusters counts
  int K = counts.n_elem;
  double icl_prop = 0;
  if(counts(oldcl)!=0){
    // both clusters are healthy
    icl_prop = lgamma(K*alpha)+lgamma(alpha+counts(oldcl))+lgamma(alpha+counts(newcl))-K*lgamma(alpha)-lgamma(K*alpha+N);
  }else{
    // cluster oldclass is dead, count(oldcl)==0 chnage of dimension
    icl_prop = lgamma((K-1)*alpha)+lgamma(alpha+counts(newcl,0))-(K-1)*lgamma(alpha)-lgamma((K-1)*alpha+N);
  }
  // complete with log(p(X|X)) from derived class
  return icl_prop;
}

arma::mat CombinedIclModel::delta_prop_swap(int i,arma::uvec iclust){
  int oldcl = cl(i);
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  int k = 0;
  // for each possible move
  for(arma::uword j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    if(k!=oldcl){
      arma::vec new_counts = update_count(counts,oldcl,k);
      delta(k)=icl_prop(new_counts,oldcl,k)-icl_prop(counts,oldcl,k);
    }
  }
  return delta;
}

arma::mat CombinedIclModel::delta_swap(int i,arma::uvec iclust){
  arma::vec delta(K);
  delta.fill(0);
  for(int m=0;m<IclModels.size();m++){
    IclModelEmission * Mp = IclModels[m];
    delta += Mp->delta_swap(i,K,cl,iclust);
  }
  delta += delta_prop_swap(i,iclust);
  return delta;
}


void CombinedIclModel::swap_update(int i,int newcl){
  bool dead_cluster = false;
  if(counts(cl(i))==1){
    dead_cluster=true;
  }
  for(int m=0;m<IclModels.size();m++){
    IclModelEmission * Mp = IclModels[m];
    Mp->swap_update(i,cl,dead_cluster,newcl);
  }
  int oldcl = cl(i);
  arma::mat new_counts = update_count(counts,oldcl,newcl);
  cl(i)=newcl;
  if(new_counts(oldcl)==0){
    counts = new_counts.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    cl.elem(arma::find(cl>oldcl))=cl.elem(arma::find(cl>oldcl))-1;
    --K;
  }else{
    counts=new_counts;
  }
}


double CombinedIclModel::delta_merge(int k,int l){
  double delta = 0;
  for(int i=0;i<IclModels.size();i++){
    IclModelEmission * Mp = IclModels[i];
    delta += Mp->delta_merge(k,l);
  }

  arma::mat new_counts = counts;
  new_counts(l)=new_counts(k)+new_counts(l);
  new_counts(k)=0;
  delta += icl_prop(new_counts,k,l)-icl_prop(counts,k,l);
  return delta;
}



// Merge update between k and l
void CombinedIclModel::merge_update(int k,int l){
  for(int i=0;i<IclModels.size();i++){
    IclModelEmission * Mp = IclModels[i];
    Mp->merge_update(k,l);
  }
  //update cl
  cl(arma::find(cl==k))=arma::ones(counts(k),1)*l;
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  // update counts
  counts(l) = counts(k)+counts(l);
  counts.shed_row(k);
  --K;
}





double CombinedIclModel::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
  double correction = 0;
  for(int i=0;i<IclModels.size();i++){
    IclModelEmission * Mp = IclModels[i];
    correction += Mp->delta_merge_correction(k,l,obk,obl,old_stats[i+1]);
  }
  return correction;
}




