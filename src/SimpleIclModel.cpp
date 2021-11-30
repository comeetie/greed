// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "IclModel.h"
#include "SimpleIclModel.h"
#include "IclModelEmission.h"
using namespace Rcpp;



SimpleIclModel::SimpleIclModel(IclModelEmission * emission_modeli,S4 modeli, arma::uvec cli, bool verb){

  model=modeli;
  alpha = model.slot("alpha");
  emission_model = emission_modeli;
  set_cl(cli);
  verbose=verb;
}

void SimpleIclModel::set_cl(arma::uvec cli){
  N = cli.n_elem; 
  K = arma::max(cli)+1;
  cl=cli;
  counts = count(cli,K);
  emission_model->set_cl(cli);
}

// main function to compute ICL from observed stats
double SimpleIclModel::icl(const List & obs_stats){
  // compute the first part p(Z) from clusters counts
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  // number of cluster
  int K = counts.n_elem;
  // log(p(Z))
  double icl_prop = lgamma(K*alpha)+arma::accu(lgamma(alpha+counts))-K*lgamma(alpha)-lgamma(arma::accu(counts+alpha));
  // complete with log(p(X|X)) from derived class
  double icl_e = this->icl_emiss(obs_stats);
  double icl = icl_prop+icl_e;
  return icl;
}

// main function to compute ICL from observed stats optimized version for computing delta which only invlove change in 2 clusters
double SimpleIclModel::icl(const List & obs_stats,int oldcl,int newcl){
  // compute the first part p(Z) from clusters counts
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  int K = counts.n_elem;
  double icl_prop = 0;
  if(counts(oldcl)!=0){
    // both clusters are healthy
    icl_prop = lgamma(K*alpha)+lgamma(alpha+counts(oldcl))+lgamma(alpha+counts(newcl))-K*lgamma(alpha)-lgamma(K*alpha+N);
  }else{
    // cluster oldclass is dead, count(oldcl)==0 chnage of dimension
    icl_prop = lgamma((K-1)*alpha)+lgamma(alpha+counts(newcl))-(K-1)*lgamma(alpha)-lgamma((K-1)*alpha+N);
  }
  // complete with log(p(X|X)) from derived class
  double icl_e = this->icl_emiss(obs_stats,oldcl,newcl);
  double icl = icl_prop+icl_e;
  return icl;
}



List SimpleIclModel::get_obs_stats(){
  List obs_stats = List::create(counts);
  obs_stats.push_back(emission_model->get_obs_stats());
  CharacterVector compn = CharacterVector(2);
  compn[0]="counts";
  compn[1]=as<std::string>(model.attr("class"));
  obs_stats.names() = compn;
  //return observed stats
  return obs_stats;
}



double SimpleIclModel::icl_emiss(const List & obs_stats){
  return(emission_model->icl_emiss(obs_stats[1]));
}





double SimpleIclModel::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  return(emission_model->icl_emiss(obs_stats[1],oldcl,newcl,false));
}



// main function to compute ICL from observed stats
double SimpleIclModel::icl_prop(arma::vec counts){
  // number of cluster
  int K = counts.n_elem;
  // log(p(Z))
  double icl_prop = lgamma(K*alpha)+arma::accu(lgamma(alpha+counts))-K*lgamma(alpha)-lgamma(arma::accu(counts+alpha));
  return icl_prop;
}

// main function to compute ICL from observed stats optimized version for computing delta which only invlove change in 2 clusters
double SimpleIclModel::icl_prop(arma::vec counts,int oldcl,int newcl){
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
  return icl_prop;
}

arma::vec SimpleIclModel::delta_prop_swap(int i,arma::uvec iclust){
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

arma::vec SimpleIclModel::delta_swap(int i,arma::uvec iclust){
  arma::vec delta(K);
  delta = delta_prop_swap(i,iclust);
  bool dead_cluster=false;
  if(counts(cl(i))==1){
    dead_cluster=true;
  }
  delta += emission_model->delta_swap(i,cl,dead_cluster,iclust,K);
  return delta;
}


void SimpleIclModel::swap_update(int i,int newcl){
  bool almost_dead_cluster = false;
  if(counts(cl(i))==1){
    almost_dead_cluster=true;
  }
  emission_model->swap_update(i,cl,almost_dead_cluster,newcl);

  int oldcl = cl(i);
  counts = update_count(counts,oldcl,newcl);
  cl(i)=newcl;
  if(counts(oldcl)==0){
    counts.shed_row(oldcl);
    cl.elem(arma::find(cl>oldcl))=cl.elem(arma::find(cl>oldcl))-1;
    --K;
  }
}


double SimpleIclModel::delta_merge(int k,int l){
  double delta = emission_model->delta_merge(k,l);
  arma::vec new_counts = counts;
  new_counts(l)=new_counts(k)+new_counts(l);
  new_counts(k)=0;
  delta += icl_prop(new_counts,k,l)-icl_prop(counts,k,l);
  return delta;
}



// Merge update between k and l
void SimpleIclModel::merge_update(int k,int l){
  emission_model->merge_update(k,l);
  
  //update cl
  cl(arma::find(cl==k)).fill(l);
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  // update counts
  counts(l) = counts(k)+counts(l);
  counts.shed_row(k);
  --K;
}


double SimpleIclModel::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
  return(emission_model->delta_merge_correction(k,l,obk,obl,old_stats[1]));
}


double SimpleIclModel::delta_merge_correction_prop(int k,int l,int obk,int obl,const List & old_stats){
  int Kold = K+1;
  double cor_prop = lgamma((Kold-2)*alpha)-2*lgamma((Kold-1)*alpha)+lgamma(Kold*alpha)+2*lgamma((Kold-1)*alpha+N)-lgamma((Kold-2)*alpha+N)-lgamma(Kold*alpha+N);
  return cor_prop;
}

