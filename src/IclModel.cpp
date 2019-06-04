// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
using namespace Rcpp;


// main function to compute ICL from observed stats
double IclModel::icl(const List & obs_stats){
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
double IclModel::icl(const List & obs_stats,int oldcl,int newcl){
  // compute the first part p(Z) from clusters counts
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
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
  double icl_e = this->icl_emiss(obs_stats,oldcl,newcl);
  double icl = icl_prop+icl_e;
  return icl;
}

// main function for greedy swaping
void IclModel::greedy_swap(int nbpassmax){
  // number of pass over data
  int nbpass = 0;
  // number of move during the current pass
  int nbmove = 0;
  // current node to swap
  int cnode = 0;
  // boolean to test if a move occurs
  bool hasMoved = true;
  // while their are moves
  arma::vec workingset(N);
  workingset.ones();
  while (hasMoved && nbpass < nbpassmax){
    // suffle the index 
    arma::vec pass= as<arma::vec>(sample(N,N))-1;
    // reinit move counter
    hasMoved=false;
    nbmove=0;
    // perform a pass
    for (int i=0;i<N ;++i){
      // current node
      cnode=pass(i);
      //Rcout << cnode << "K:" << K << std::endl;
      if (workingset(cnode)==1 && K>1){
        // compute delta swap
        arma::vec delta = this->delta_swap(cnode);
        //Rcout << delta << std::endl;
        // best swap
        int ncl = delta.index_max();
        if(ncl!=cl(cnode)){
  
          // if best swap corresponds to a move
          // update the stats
          this->swap_update(cnode,ncl);
  
          // update the move counters
          hasMoved=true;
          ++nbmove;
          // workingset(cnode)=1;
        }else{
          arma::vec deltaneg = delta.elem(arma::find(delta<0));
          int bmn= deltaneg.index_max();
          if(deltaneg(bmn) <  -10 ){
            workingset(cnode) = 0;
            //Rcout << "BMN :"<< bmn << "val" << deltaneg(bmn) << std::endl;
          }
        }
      }
    }
    // update the pass counter
    ++nbpass;
    // compute icl after the pass
    icl_value = icl(this->get_obs_stats());

    if(verbose){
      Rcout << "##################################"<< std::endl;
      Rcout << "Pass N°"<< nbpass << " completed with " << nbmove << " moves, icl :" << icl_value << "K :" << K << ", working set size :" << arma::accu(workingset)  << std::endl;
      Rcout << "##################################"<< std::endl; 
    } 
  }
}


// get p(z_i|X,z-i)
arma::mat IclModel::get_probs(){

  // perform a pass
   arma::mat probs(N,K);
   
  for (int i=0;i<N ;++i){
    // compute delta swap
    arma::vec delta = this->delta_swap(i);
    // transform to probabilities
    arma::rowvec pr =  (exp(delta)/arma::accu(exp(delta))).t();
    probs.row(i) = pr;
  } 
  return probs;
}


// main function for greedy merging
void IclModel::greedy_merge(){
  
  // init the merge matrix(K,K) with the delta icl of each merge 
  MergeMat merge_mat = this->delta_merge();
  // init merge counter
  int nbmerge = 0;
  
  // while a positive merge exists
  while(merge_mat.getValue()>0){
    
    // increment 
    ++nbmerge;
    
      
    
    // perform the merge and update the stats
    this->merge_update(merge_mat.getK(),merge_mat.getL());
    if(verbose){
      Rcout << "##################################"<< std::endl;
      Rcout << "Merge icl : "<< icl(this->get_obs_stats()) << std::endl;
      Rcout << "##################################"<< std::endl; 
    }
    // update the merge matrix
    merge_mat = this->delta_merge(merge_mat.getMergeMat(),merge_mat.getK(),merge_mat.getL());

  }
  // compute final icl value
  icl_value = icl(this->get_obs_stats());
  if(verbose){
    Rcout << "##################################"<< std::endl;
    Rcout << "Final icl : "<< icl_value << std::endl;
    Rcout << "##################################"<< std::endl; 
  }
}


// main function for greedy merge path perform all merge until its remain only one cluster and store the successive models
List IclModel::greedy_merge_path(){
  // set alpha to 1 in order to compute icl and log(alpha)
  alpha = 1;
  // init the merge matrix(K,K) with the delta icl of each merge 
  MergeMat merge_mat = this->delta_merge();
  // init list to store merge path
  List path = List();
  // init merge counter
  int nbmerge = 0;
  double icl = this->icl(this->get_obs_stats());
  double iclold = icl;
  while(K>1){
    // arma::mat mm = this->delta_merge().getMergeMat();
    // Rcout << mm << std::endl;
    ++nbmerge;
    // get best merge
    int k = merge_mat.getK();
    int l = merge_mat.getL();
    // update stats
    // Rcout << k<< l << ":  "<< mm(k,l) << std::endl;
    this->merge_update(k,l);
    // compute new icl
    iclold = icl;
    icl = this->icl(this->get_obs_stats());
    // store current solution
    path.push_back(List::create(Named("obs_stats") = this->get_obs_stats(), 
                                Named("cl") = cl+1, 
                                Named("K") = K,
                                Named("icl1")=icl,
                                Named("logalpha")=icl-iclold,
                                Named("k")=k+1,
                                Named("l")=l+1));
    // update merge matrix
    merge_mat = this->delta_merge(merge_mat.getMergeMat(),merge_mat.getK(),merge_mat.getL());
  }
  return path;
}



