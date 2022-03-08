// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "SpMergeMat.h"
#include "IclModel.h"

using namespace Rcpp;



// main function for greedy swaping
void IclModel::greedy_swap(int nbpassmax, arma::vec workingset,arma::uvec iclust){

  // number of pass over data
  int nbpass = 0;
  // number of move during the current pass
  int nbmove = 0;
  // current node to swap
  int cnode = 0;
  // boolean to test if a move occurs
  bool hasMoved = true;
  // while their are moves

  while (hasMoved && nbpass < nbpassmax && K>1 && iclust.n_elem>1){
    

    // suffle the index 
    arma::uvec pass= as<arma::uvec>(sample(N,N))-1;
    // reinit move counter
    hasMoved=false;
    nbmove=0;
    // perform a pass
    for (int i=0;i<N ;++i){
      // current node
      cnode=pass(i);

      if (workingset(cnode)==1){
        // compute delta swap

        arma::vec delta = this->delta_swap(cnode,iclust);

        // best swap
        int ncl=cl(cnode);
        try{
          ncl = delta.index_max();
        }catch(std::exception &ex){
          warning("Undefined move cost, please check the priors parameters.");
          ncl=cl(cnode);
        }

        if(ncl!=cl(cnode)){
  
          //if best swap corresponds to a move
          //update the stats
          if(counts(cl(cnode))==1){
            //Rcout << iclust << std::endl;
            iclust = iclust.elem(arma::find(iclust!=cl(cnode)));
            iclust.elem(arma::find(iclust>cl(cnode)))=iclust.elem(arma::find(iclust>cl(cnode)))-1;
            //Rcout << iclust << std::endl;
          }
          this->swap_update(cnode,ncl);
  
          // update the move counters
          hasMoved=true;
          ++nbmove;
          // workingset(cnode)=1;
          if(K==1){
            break;
          }
        }else{
          arma::vec deltaneg = delta.elem(arma::find(delta<0));
          int bmn= deltaneg.index_max();
          if(deltaneg(bmn) <  -4 ){
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
      Rcout << "Swap convergence in " << nbpass << " epochs with " << nbmove << " moves, icl :" << icl_value << "K :" << K << ", working set size :" << arma::accu(workingset)  << std::endl;
      //Rcout << "Swap convergence, with an ICL of "<< icl_value << " and " << K << " clusters." << std::endl;
      Rcout << "##################################"<< std::endl; 
    } 

  }

}


// main function for greedy swaping with move constraints given as a sparse KxK matrix 
void IclModel::greedy_swap(int nbpassmax, arma::vec workingset,arma::sp_mat & moves_mat){
  // number of pass over data
  int nbpass = 0;
  // number of move during the current pass
  int nbmove = 0;
  // current node to swap
  int cnode = 0;
  // boolean to test if a move occurs
  bool hasMoved = true;
  // while their are moves
  while (hasMoved && nbpass < nbpassmax && K>1){
    // suffle the index 
    arma::uvec pass= as<arma::uvec>(sample(N,N))-1;
    // reinit move counter
    hasMoved=false;
    nbmove=0;
    // perform a pass
    for (int i=0;i<N ;++i){
      // current node
      cnode=pass(i);
      // Rcout << cnode << "K:" << K << std::endl;
      if (workingset(cnode)==1){
        // compute delta swap
        arma::uvec iclust = possible_moves(cl(cnode),moves_mat);
        arma::vec delta = this->delta_swap(cnode,iclust);
        //Rcout << delta << std::endl;
        // best swap
        int ncl=cl(cnode);
        try{
          ncl = delta.index_max();
        }catch(std::exception &ex){
          warning("Undefined move cost, please check the input data and the priors parameters.");
          ncl=cl(cnode);
        }
        
        if(ncl!=cl(cnode)){
          
          // if best swap corresponds to a move
          // update the stats and deal with cluster death
          if(counts(cl(cnode))==1){
              // remove the cluster from the move matrix
              moves_mat=delrowcol_copy(moves_mat,cl(cnode));

          }
          this->swap_update(cnode,ncl);
          
          // update the move counters
          hasMoved=true;
          ++nbmove;
          // workingset(cnode)=1;
          if(K==1){
            break;
          }
        }else{
          try{
            arma::vec deltaneg = delta.elem(arma::find(delta<0));
            int bmn= deltaneg.index_max();
            if(deltaneg(bmn) <  -2 ){
              workingset(cnode) = 0;
              //Rcout << "BMN :"<< bmn << "val" << deltaneg(bmn) << std::endl;
            }
          }catch(std::exception &ex){
            warning("Undefined move cost, please check the input data and the priors parameters.");
            workingset(cnode) = 0;
          }
        }
      }
    }
    // update the pass counter
    ++nbpass;
    // compute icl after the pass
    icl_value = icl(this->get_obs_stats());
    
  }
  if(verbose){
    Rcout << "##################################"<< std::endl;
    //plaRcout << "Swap convergence in " << nbpass << " epochs with " << nbmove << " moves, icl :" << icl_value << "K :" << K << ", working set size :" << arma::accu(workingset)  << std::endl;
    Rcout << "swap convergence, with an ICL of "<< icl_value << " and " << K << " clusters." << std::endl;
    Rcout << "##################################"<< std::endl; 
  } 
}


// get p(z_i|X,z-i)
arma::mat IclModel::get_probs(){

  // perform a pass
  arma::mat probs(N,K);
  arma::uvec iclust = arma::find(arma::ones(K));
  for (int i=0;i<N ;++i){
    // compute delta swap
    arma::vec delta = this->delta_swap(i,iclust);
    // transform to probabilities
    arma::rowvec pr =  (exp(delta)/arma::accu(exp(delta))).t();
    probs.row(i) = pr;
  } 
  return probs;
}


// init merge matrix
MergeMat IclModel::delta_merge(){
  // inititalize delta merge matrix
  arma::mat delta(K,K);
  delta.fill(0);
  // index to store current best merge
  int bk = 0;
  int bl = 0;
  // initialize bv found to -infty
  double bv = -std::numeric_limits<double>::infinity();
  for(int k = 1; k < K; ++k) {
    for (int l = 0;l<k;++l){
      // for each possible merge
      delta(k,l)=this->delta_merge(k,l);
      // best merge ?
      if(delta(k,l)>bv){
        bk=k;
        bl=l;
        bv=delta(k,l);
      }
    }
  }
  return MergeMat(bk,bl,bv,delta);
}

// update merge matrix after merge of obk/obl
// obl < obk so didn't change when removing row/col obk
MergeMat IclModel::delta_merge(arma::mat delta, int obk, int obl, const List & old_stats){
  // optimized version to compute only new values of the merge mat
  delta = delta(arma::find(arma::linspace(0,K,K+1)!=obk),arma::find(arma::linspace(0,K,K+1)!=obk));
  int bk = 0;
  int bl = 0;
  double bv = -std::numeric_limits<double>::infinity();
  for(int k = 1; k < K; ++k) {
    for (int l = 0;l<k;++l){
      if((k == obl) | (l == obl)){
        delta(k,l)=this->delta_merge(k,l);
      }else{
        //Rcout << k <<" : " <<l << " - " << obk <<" : " << obl <<std::endl;
        //Rcout << delta(k,l) << std::endl;
        delta(k,l)=delta(k,l)+this->delta_merge_correction_prop(k,l,obk,obl,old_stats)+this->delta_merge_correction(k,l,obk,obl,old_stats);
        //Rcout << delta(k,l) << std::endl;
      }
      if(delta(k,l)>bv){
        bk=k;
        bl=l;
        bv=delta(k,l);
      }
      
    }
  }
  return MergeMat(bk,bl,bv,delta);
}

// init merge matrix sparse
SpMergeMat IclModel::delta_merge(const arma::sp_mat & merge_graph){
  
  // inititalize delta merge matrix
  arma::sp_mat delta = merge_graph;
  // index to store current best merge
  int bk = 0;
  int bl = 0;
  // initialize bv found to -infty
  double bv = -std::numeric_limits<double>::infinity();
  // store cuurent stats
  for (arma::sp_mat::iterator i = delta.begin(); i != delta.end(); ++i) {
    if(i.col()<i.row()){
      delta(i.row(),i.col())=this->delta_merge(i.row(),i.col());
      delta(i.col(),i.row())=delta(i.row(),i.col());
      // best merge ?
      if(delta(i.row(),i.col())>bv){
        bk=i.row();
        bl=i.col();
        bv=delta(i.row(),i.col());
      }
    }
   
  }

  return SpMergeMat(bk,bl,bv,delta);
}





// update merge matrix after merge of obk/obl sparse
SpMergeMat IclModel::delta_merge(arma::sp_mat & merge_graph, int obk, int obl,const List & old_stats){
  // optimized version to compute only new values of the merge mat
  //delta = delta(arma::find(arma::linspace(0,K,K+1)!=obk),arma::find(arma::linspace(0,K,K+1)!=obk));
  
  
  arma::sp_mat delta = delrowcol_copy(merge_graph,obk);
  int bk = 0;
  int bl = 0;
  int k,l;
  double bv = -std::numeric_limits<double>::infinity();
  arma::sp_mat::iterator i = delta.begin();
  arma::sp_mat::iterator end = delta.end();
  for (; i != end; ++i) {
     if(i.col()<i.row()){
      k=i.row();
      l=i.col();
      if((k == obl) | (l == obl)){
        delta(k,l)=this->delta_merge(k,l);
      }else{
        delta(k,l)=delta(k,l)+this->delta_merge_correction_prop(k,l,obk,obl,old_stats)+this->delta_merge_correction(k,l,obk,obl,old_stats);
      }
      delta(l,k)=delta(k,l);

      if(delta(k,l)>bv){
        bk=k;
        bl=l;
        bv=delta(k,l);
      }
    }
  }
  
  return SpMergeMat(bk,bl,bv,delta);
}



// main function for greedy merging with prior merge graph
arma::sp_mat IclModel::greedy_merge(const arma::sp_mat & merge_graph){
  // init the merge matrix(K,K) with the delta icl of each merge 
  SpMergeMat merge_mat = this->delta_merge(merge_graph);
  arma::sp_mat delta = merge_mat.getMergeMat();
  // init merge counter
  int nbmerge = 0;
  double cicl = this->icl(this->get_obs_stats());
  double bicl= cicl;
  arma::sp_mat best_merge_mat = delta;
  arma::uvec bcl = cl;
  // while their are merge to explore
  while(delta.n_nonzero>0 && K>1){
    
    // increment
    ++nbmerge;
    List old_stats = this->get_obs_stats();
    // perform the merge and update the stats
    this->merge_update(merge_mat.getK(),merge_mat.getL());
    if(verbose){
       Rcout << "##################################"<< std::endl;
       Rcout << "Merge icl : "<< icl(this->get_obs_stats()) << std::endl;
       Rcout << "##################################"<< std::endl;
    }
    // update the merge matrix
    //merge_mat = this->delta_merge(merge_mat.getMergeMat(),merge_mat.getK(),merge_mat.getL(),old_stats);

    //Rcout << delta.n_nonzero << std::endl;
    
    cicl = cicl + merge_mat.getValue();

    // Rcout << merge_mat.getMergeMat() << std::endl;
    // int Ko = merge_mat.getK();
    merge_mat = this->delta_merge(delta,merge_mat.getK(),merge_mat.getL(),old_stats);
    
    delta = merge_mat.getMergeMat();
    if(cicl > bicl){
      bicl=cicl;
      bcl = cl;
      best_merge_mat=delta;
    }
    
    
    //check test for merge mat correction
    // delrowcol(delta,Ko);
    // Rcout << "--- check correction ---" << std::endl;
    // Rcout << K << std::endl;
    // SpMergeMat merge_mat_comp = this->delta_merge(delta);
    // Rcout << arma::max(arma::max(arma::abs(merge_mat.getMergeMat()-merge_mat_comp.getMergeMat()))) << std::endl;
    // Rcout << merge_mat.getMergeMat() << std::endl;
    // Rcout << merge_mat_comp.getMergeMat() << std::endl;
  }
  //compute final icl value
  this->set_cl(bcl);
  icl_value = icl(this->get_obs_stats());
  
  if(verbose){
    Rcout << "##################################"<< std::endl;
    Rcout << "Merge convergence, with an ICL of "<< icl_value << " and " << K << " clusters." << std::endl;
    Rcout << "##################################"<< std::endl; 
  }
  return best_merge_mat;
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
    
    //Rcout << merge_mat.getValue()<< std::endl;
    List old_stats = this->get_obs_stats();
    // perform the merge and update the stats
    this->merge_update(merge_mat.getK(),merge_mat.getL());
    // if(verbose){
    //   Rcout << "##################################"<< std::endl;
    //   Rcout << "Merge icl : "<< icl(this->get_obs_stats()) << std::endl;
    //   Rcout << "##################################"<< std::endl; 
    // }
    // update the merge matrix
    merge_mat = this->delta_merge(merge_mat.getMergeMat(),merge_mat.getK(),merge_mat.getL(),old_stats);
    //Rcout << merge_mat.getValue()<< std::endl;
    // check test for merge mat correction
    // MergeMat merge_mat_comp = this->delta_merge();
    // Rcout << arma::max(merge_mat.getMergeMat()-merge_mat_comp.getMergeMat()) << std::endl;
    // Rcout << arma::min(merge_mat.getMergeMat()-merge_mat_comp.getMergeMat()) << std::endl;
    
  }
  // compute final icl value
  icl_value = icl(this->get_obs_stats());
  if(verbose){
    Rcout << "##################################"<< std::endl;
    Rcout << "Merge convergence, with an ICL of "<< icl_value << " and " << K << " clusters." << std::endl;
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
  List obs_stats = this->get_obs_stats();
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  double icl = this->icl_emiss(obs_stats)-log(static_cast<double>(K))+arma::accu(lgamma(counts))-lgamma(N);
  double iclold = icl;
  int k=0;
  int l=0;
  while(K>1){
    // arma::mat mm = this->delta_merge().getMergeMat();
    // Rcout << mm << std::endl;
    ++nbmerge;
    // get best merge
    if(K==2){
      k=1;
      l=0;
    }else{
      k = merge_mat.getK();
      l = merge_mat.getL();
    }
    // update stats
    // Rcout << k<< l << ":  "<< mm(k,l) << std::endl;
    List old_stats = this->get_obs_stats();
    this->merge_update(k,l);
    // compute new icl
    iclold = icl;
    obs_stats = this->get_obs_stats();
    if(merge_mat.getValue()>-std::numeric_limits<double>::infinity()){
      counts =as<arma::vec>(obs_stats["counts"]);
      icl = this->icl_emiss(obs_stats)-log(static_cast<double>(K))+arma::accu(lgamma(counts))-lgamma(N);
    }else{
      icl = -std::numeric_limits<double>::infinity();
    }
    // icl = icl+merge_mat.getValue();
    // this->icl(this->get_obs_stats());
    // store current solution

    path.push_back(List::create(Named("obs_stats") = this->get_obs_stats(), 
                                Named("K") = K,
                                Named("icl1")=icl,
                                Named("logalpha")=icl-iclold,
                                Named("k")=k+1,
                                Named("l")=l+1,
                                // useless? : 
                                Named("merge_mat") = arma::trimatl(merge_mat.getMergeMat())));
                                
    // update merge matrix
    merge_mat = this->delta_merge(merge_mat.getMergeMat(),merge_mat.getK(),merge_mat.getL(),old_stats);
  }
  return path;
}



