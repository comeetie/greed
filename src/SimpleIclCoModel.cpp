// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "IclModel.h"
#include "SimpleIclCoModel.h"
#include "IclModelEmission.h"
using namespace Rcpp;



SimpleIclCoModel::SimpleIclCoModel(IclModelEmission * emission_modeli,S4 modeli, arma::uvec cli,int Nri, int Nci, bool verb){
  N  = Nri+Nci;
  Nc=Nci;
  Nr=Nri;
  model=modeli;
  alpha = model.slot("alpha");
  emission_model = emission_modeli;
  set_cl(cli);
  verbose=verb;

}

void SimpleIclCoModel::set_cl(arma::uvec cli){
  cl = cli;
  K =arma::max(cl)+1;
  arma::uvec clr = cl.subvec(0,Nr-1);
  row_clusts = arma::unique(clr);
  Kr = row_clusts.n_elem;
  arma::uvec clc = cl.subvec(Nr,N-1);
  col_clusts = arma::unique(clc);
  Kc = col_clusts.n_elem;

  counts = count(cl,K);
  arma::vec types(K);
  types.fill(0);
  for (arma::uword i = 0;i<row_clusts.n_elem;++i){
    types(row_clusts(i))=1; 
  }
  for (arma::uword i = 0;i<col_clusts.n_elem;++i){
    if(types(col_clusts(i))!=0){
      Rcpp::stop("Invalid partition in co-clustering");      
    }
    types(col_clusts(i))=2;
  }
  clusttypes=types;
  emission_model->set_cl(cli);
}

// main function to compute ICL from observed stats
double SimpleIclCoModel::icl(const List & obs_stats){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  double iprop = this->icl_prop(counts);
  double iemiss = this->icl_emiss(obs_stats);
  return iprop+iemiss;
}

// main function to compute ICL from observed stats optimized version for computing delta which only invlove change in 2 clusters
double SimpleIclCoModel::icl(const List & obs_stats,int oldcl,int newcl){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  double iprop = this->icl_prop(counts,oldcl,newcl);
  double iemiss = this->icl_emiss(obs_stats,oldcl,newcl);
  return iprop+iemiss;
}



List SimpleIclCoModel::get_obs_stats(){
  List obs_stats = List::create(counts);
  obs_stats.push_back(emission_model->get_obs_stats());
  CharacterVector compn = CharacterVector(2);
  compn[0]="counts";
  compn[1]=as<std::string>(model.attr("class"));
  obs_stats.names() = compn;
  //return observed stats
  return obs_stats;
}



double SimpleIclCoModel::icl_emiss(const List & obs_stats){
  return(emission_model->icl_emiss(obs_stats[1]));
}





double SimpleIclCoModel::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  return(emission_model->icl_emiss(obs_stats[1],oldcl,newcl,false));
}



// main function to compute ICL from observed stats
double SimpleIclCoModel::icl_prop(arma::vec counts){
  // log(p(Z))
  double icl_prop = lgamma(Kr*alpha)+lgamma(Kc*alpha)+arma::accu(lgamma(alpha+counts))-K*lgamma(alpha)-lgamma(Nr+alpha*Kr)-lgamma(Nc+alpha*Kc);
  return icl_prop;
}

// main function to compute ICL from observed stats optimized version for computing delta which only invlove change in 2 clusters
double SimpleIclCoModel::icl_prop(arma::vec counts,int oldcl,int newcl){
  double icl_prop = 0;
  if(counts(oldcl)!=0){
    // both clusters are healthy
    icl_prop = lgamma(Kr*alpha)+lgamma(Kc*alpha)+lgamma(alpha+counts(oldcl))+lgamma(alpha+counts(newcl))-K*lgamma(alpha)-lgamma(Nr+alpha*Kr)-lgamma(Nc+alpha*Kc);
  }else{
    // cluster oldclass is dead, count(oldcl)==0 chnage of dimension
    if(clusttypes(oldcl)==1){
      icl_prop = lgamma((Kr-1)*alpha)+lgamma(Kc*alpha)+lgamma(alpha+counts(newcl))-(K-1)*lgamma(alpha)-lgamma(Nr+alpha*(Kr-1))-lgamma(Nc+alpha*Kc);
    }
    if(clusttypes(oldcl)==2){
      icl_prop = lgamma(Kr*alpha)+lgamma((Kc-1)*alpha)+lgamma(alpha+counts(newcl))-(K-1)*lgamma(alpha)-lgamma(Nr+alpha*Kr)-lgamma(Nc+alpha*(Kc-1));
    }
  }
  
  return icl_prop;
}

arma::vec SimpleIclCoModel::delta_prop_swap(int i,arma::uvec iclust){
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

arma::vec SimpleIclCoModel::delta_swap(int i,arma::uvec iclust){
  arma::vec delta(K);
  delta = delta_prop_swap(i,iclust);
  bool dead_cluster=false;
  if(counts(cl(i))==1){
    dead_cluster=true;
  }
  delta += emission_model->delta_swap(i,cl,dead_cluster,iclust,K);
  return delta;
}


void SimpleIclCoModel::swap_update(int i,int newcl){
  bool almost_dead_cluster = false;
  if(counts(cl(i))==1){
    almost_dead_cluster=true;
  }
  emission_model->swap_update(i,cl,almost_dead_cluster,newcl);

  int oldcl = cl(i);
  int cd = 0;
  counts = update_count(counts,oldcl,newcl);
  cl(i)=newcl;
  if(counts(oldcl)==0){
    counts.shed_row(oldcl);
    cl.elem(arma::find(cl>oldcl))=cl.elem(arma::find(cl>oldcl))-1;
    if(clusttypes(oldcl)==1){
      Kr=Kr-1;
    }
    if(clusttypes(oldcl)==2){
      Kc=Kc-1;
    }
    clusttypes.shed_row(oldcl);
    row_clusts = arma::find(clusttypes==1);
    col_clusts = arma::find(clusttypes!=1);
    --K;
  }
}


double SimpleIclCoModel::delta_merge(int k,int l){
  double delta = emission_model->delta_merge(k,l);
  arma::vec new_counts = counts;
  new_counts(l)=new_counts(k)+new_counts(l);
  new_counts(k)=0;
  delta += icl_prop(new_counts,k,l)-icl_prop(counts,k,l);
  return delta;
}



// Merge update between k and l
void SimpleIclCoModel::merge_update(int k,int l){
  emission_model->merge_update(k,l);
  
  //update cl
  cl(arma::find(cl==k)).fill(l);
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  // update counts
  counts(l) = counts(k)+counts(l);
  counts.shed_row(k);
  --K;
  
  
  if(clusttypes(k)==1){
    Kr=Kr-1;
  }else{
    Kc=Kc-1;
  }
  
  clusttypes.shed_row(k);
  row_clusts = arma::find(clusttypes==1);
  col_clusts = arma::find(clusttypes!=1);
}


double SimpleIclCoModel::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
  return(emission_model->delta_merge_correction(k,l,obk,obl,old_stats[1]));
}

double SimpleIclCoModel::delta_merge_correction_prop(int k,int l,int obk,int obl,const List & old_stats){
  int Kold;
  double cor_prop = 0;
  if((clusttypes(k)==1) &&  (clusttypes(l)==1) && (clusttypes(obl)==1)){
    Kold = Kr+1;
    cor_prop = lgamma((Kold-2)*alpha)-2*lgamma((Kold-1)*alpha)+lgamma(Kold*alpha)+2*lgamma((Kold-1)*alpha+N)-lgamma((Kold-2)*alpha+N)-lgamma(Kold*alpha+N);
  }
  if((clusttypes(k)==2) &&  (clusttypes(l)==2) && (clusttypes(obl)==2)){
    Kold = Kc+1;
    cor_prop = lgamma((Kold-2)*alpha)-2*lgamma((Kold-1)*alpha)+lgamma(Kold*alpha)+2*lgamma((Kold-1)*alpha+N)-lgamma((Kold-2)*alpha+N)-lgamma(Kold*alpha+N);
  }
  return cor_prop;
}
