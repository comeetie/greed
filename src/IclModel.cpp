// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
using namespace Rcpp;

double IclModel::icl(const List & obs_stats){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  
  int N = arma::accu(counts);
  int K = x_counts.n_elem;
  double icl_prop = lgamma(K*alpha)+arma::accu(lgamma(alpha+counts))-K*lgamma(alpha)-lgamma(arma::accu(counts+alpha));
  double icl_e = this->icl_emiss(obs_stats);

  double icl = icl_prop+icl_e;
  return icl;
}

double IclModel::icl(const List & obs_stats,int oldcl,int newcl){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  double icl_prop = 0;
  if(counts(oldcl)!=0){
    icl_prop = lgamma(K*alpha)+lgamma(alpha+counts(oldcl,0))+lgamma(alpha+counts(newcl,0))-K*lgamma(alpha)-lgamma(K*alpha+N);
  }else{
    icl_prop = lgamma((K-1)*alpha)+lgamma(alpha+counts(newcl,0))-(K-1)*lgamma(alpha)-lgamma((K-1)*alpha+N);
  }
  double icl_e = this->icl_emiss(obs_stats,oldcl,newcl);
  double icl = icl_prop+icl_e;
  return icl;
}


void IclModel::greedy_swap(int nbpassmax){
  int nbpass = 0;
  int nbmove = 0;
  int cnode = 0;
  bool hasMoved = true;
  while (hasMoved && nbpass < nbpassmax ){
    // 
    arma::vec pass= as<arma::vec>(sample(N,N))-1;
    hasMoved=false;
    nbmove=0;
    for (int i=0;i<N ;++i){
      cnode=pass(i);
      
      arma::mat delta = this->delta_swap(cnode);
      
      int ncl = delta.index_max();
      if(ncl!=cl(cnode)){
        //Rcout << delta << std::endl;
        this->swap_update(cnode,ncl);
        hasMoved=true;
        ++nbmove;
      }
    }
    ++nbpass;
    icl_value = icl(this->get_obs_stats());
    /* Rcout << "##################################"<< std::endl;
    Rcout << "Pass N°"<< nbpass << " completed with " << nbmove << " moves, icl :" << icl_value << std::endl;
    Rcout << "##################################"<< std::endl; */
  }
}


void IclModel::greedy_merge(){
  
  MergeMat merge_mat = this->delta_merge();
  int nbmerge = 0;
  
  while(merge_mat.getValue()>0){
    ++nbmerge;
    /* Rcout << "##################################"<< std::endl;
    Rcout << "Merge N°"<< nbmerge << " with delta :" << merge_mat.getValue() << std::endl;
    Rcout << "##################################"<< std::endl;  */
    
    this->merge_update(merge_mat.getK(),merge_mat.getL());
  
    merge_mat = this->delta_merge(merge_mat.getMergeMat(),merge_mat.getK(),merge_mat.getL());

  }
  icl_value = icl(this->get_obs_stats());
  /* Rcout << "##################################"<< std::endl;
  Rcout << "Final icl : "<< icl_value << std::endl;
  Rcout << "##################################"<< std::endl; */
}

List IclModel::greedy_merge_path(){
  MergeMat merge_mat = this->delta_merge();
  List path = List();
  int nbmerge = 0;
  while(K>1){
    ++nbmerge;
    /* Rcout << "##################################"<< std::endl;
    Rcout << "Merge N°"<< nbmerge << " with delta :" << merge_mat.getValue() << std::endl;
    Rcout << "##################################"<< std::endl; */
    int k = merge_mat.getK();
    int l = merge_mat.getL();
    this->merge_update(k,l);
    double icl = this->icl(this->get_obs_stats());
    path.push_back(List::create(Named("counts") = counts, 
                                Named("x_counts") = x_counts, 
                                Named("cl") = cl+1, 
                                Named("K") = K,
                                Named("icl")=icl,
                                Named("k")=k+1,
                                Named("l")=l+1));
    merge_mat = this->delta_merge(merge_mat.getMergeMat(),merge_mat.getK(),merge_mat.getL());
  }
  return path;
}
