// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "LcaE.h"
using namespace Rcpp;



LcaE::LcaE(arma::mat xi,S4 modeli,bool verb){
  model = modeli;
  beta = model.slot("beta");
  x  = xi;
  verbose=verb;

}

void LcaE::set_cl(arma::vec cl){

  K = arma::max(cl)+1;
  int nbvar = x.n_cols;
  counts = count(cl,K);
  nbmod = arma::zeros(nbvar);
  for(int i=0;i<nbvar;i++){
    nbmod(i)=arma::max(x.col(i))+1;
    arma::mat ct =table_count(cl,x.col(i),K,nbmod(i));
    x_counts.push_back(ct);
  }

}


double LcaE::icl_emiss(const List & obs_stats){
  arma::vec counts_cur = as<arma::vec>(obs_stats["counts"]);
  List x_counts_cur = as<List>(obs_stats["x_counts"]);
  double icl_emiss=0;
  
  for(int j=0;j<x.n_cols;j++){
    arma::mat temp_mat =  as<arma::mat>(x_counts_cur[j]);
    icl_emiss+=arma::accu(lgamma(temp_mat+beta));
    icl_emiss+= K*lgamma(beta*nbmod(j));
    icl_emiss-= K*nbmod(j)*lgamma(beta);
    icl_emiss-= arma::accu(lgamma(counts_cur+beta*nbmod(j)));
  }
  
  return icl_emiss;
}

double LcaE::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  arma::vec counts_cur = as<arma::vec>(obs_stats["counts"]);
  List x_counts_cur = as<List>(obs_stats["x_counts"]);
  double icl_emiss=0;
  int Kc =K;
  if(counts_cur(oldcl)==0){
    Kc--;
  }
  for(int j=0;j<x.n_cols;j++){
    arma::mat temp_mat =  as<arma::mat>(x_counts_cur[j]);
    icl_emiss+=arma::accu(lgamma(temp_mat.row(newcl)+beta));
    if(counts_cur(oldcl)!=0){
      icl_emiss+=arma::accu(lgamma(temp_mat.row(oldcl)+beta));
    }
    icl_emiss+=arma::accu(lgamma(temp_mat+beta));
    icl_emiss+= Kc*lgamma(beta*nbmod(j));
    icl_emiss-= Kc*nbmod(j)*lgamma(beta);
    icl_emiss-= lgamma(counts_cur(newcl)+beta*nbmod(j));
    if(counts_cur(oldcl)!=0){
      icl_emiss-= lgamma(counts_cur(oldcl)+beta*nbmod(j));
    }
    
  }
  
  return icl_emiss;
}


List LcaE::get_obs_stats(){
  return List::create(Named("x_counts", clone(x_counts)),Named("counts", counts));
}


arma::mat LcaE::delta_swap(int i,int K, Partition clp,arma::uvec iclust){
  

  int oldcl = clp.get(i);
 
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  
  List old_stats = List::create(Named("x_counts", x_counts),Named("counts", counts));
  arma::vec new_counts = counts;
  int k = 0;
  // for each possible move
  for(arma::uword j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    
    if(k!=oldcl){
      List new_x_counts;
      arma::vec new_counts = update_count(counts,oldcl,k);
      for(int v=0;v<x.n_cols;v++){

        arma::mat temp_mat =  as<arma::mat>(x_counts[v]);
        temp_mat(oldcl,x(i,v))=temp_mat(oldcl,x(i,v))-1;
        temp_mat(k,x(i,v))=temp_mat(k,x(i,v))+1;
        new_x_counts.push_back(temp_mat);
      }

      List new_stats = List::create(Named("x_counts", new_x_counts),Named("counts", new_counts));
      delta(k)=icl_emiss(new_stats,oldcl,k)-icl_emiss(old_stats,oldcl,k);

    }
  }
  return delta;
}



void LcaE::swap_update(const int i,Partition clp,bool dead_cluster, const int newcl){
  

  
  int oldcl = clp.get(i); 

  for(int v=0;v<x.n_cols;v++){

    arma::mat temp_mat =  as<arma::mat>(x_counts[v]);
    temp_mat(oldcl,x(i,v))=temp_mat(oldcl,x(i,v))-1;
    temp_mat(newcl,x(i,v))=temp_mat(newcl,x(i,v))+1;
    
    if(dead_cluster){
      temp_mat.shed_row(oldcl);
    }
    x_counts[v] = temp_mat;
  }
  
  counts = update_count(counts,oldcl,newcl);
  if(dead_cluster){
    counts = counts.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    K--;
  }
  
  

}


double LcaE::delta_merge(int k, int l){
  
  // k,l -> l and k is removed
  
  List old_stats = List::create(Named("x_counts", x_counts),Named("counts", counts));
  List new_x_counts;
  for(int v=0;v<x.n_cols;v++){
    arma::mat temp_mat =  as<arma::mat>(x_counts[v]);
    temp_mat.row(l) =temp_mat.row(l)+temp_mat.row(k);
    temp_mat.row(k)-=temp_mat.row(k);
    new_x_counts.push_back(temp_mat);
  }
  arma::vec new_counts=counts;
  new_counts(l)+=new_counts(k);
  new_counts(k)=0;
  List new_stats = List::create(Named("x_counts", new_x_counts),Named("counts", new_counts));
  
  double delta = icl_emiss(new_stats,k,l)-icl_emiss(old_stats,k,l);

  return delta;
}


void LcaE::merge_update(int k,int l){
  
  for(int v=0;v<x.n_cols;v++){
    arma::mat temp_mat =  as<arma::mat>(x_counts[v]);
    temp_mat.row(l)=temp_mat.row(l)+temp_mat.row(k);
    temp_mat.shed_row(k);
    x_counts[v] = temp_mat;
  }
  counts(l)+=counts(k);
  counts = counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));
  --K;
}

