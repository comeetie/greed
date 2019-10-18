// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "SpMergeMat.h"
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

// main function to compute ICL from observed stats optimized version for computing delta which only invlove change in 2 clusters ans sparse vector
double IclModel::icl(const List & obs_stats,const List & up_stats,int oldcl,int newcl){
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
  double icl_e = this->icl_emiss(obs_stats,up_stats,oldcl,newcl);
  double icl = icl_prop+icl_e;
  return icl;
}


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
    arma::vec pass= as<arma::vec>(sample(N,N))-1;
    // reinit move counter
    hasMoved=false;
    nbmove=0;
    // perform a pass
    for (int i=0;i<N ;++i){
      // current node
      cnode=pass(i);
      //Rcout << cnode << "K:" << K << std::endl;
      if (workingset(cnode)==1){
        // compute delta swap
        arma::vec delta = this->delta_swap(cnode,iclust);
        //Rcout << delta << std::endl;
        // best swap
        int ncl = delta.index_max();
        if(ncl!=cl(cnode)){
  
          // if best swap corresponds to a move
          // update the stats
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


  }
  if(verbose){
    Rcout << "##################################"<< std::endl;
    Rcout << "Swap convergence in " << nbpass << " epochs with " << nbmove << " moves, icl :" << icl_value << "K :" << K << ", working set size :" << arma::accu(workingset)  << std::endl;
    //Rcout << "Swap convergence, with an ICL of "<< icl_value << " and " << K << " clusters." << std::endl;
    Rcout << "##################################"<< std::endl; 
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
    arma::vec pass= as<arma::vec>(sample(N,N))-1;
    // reinit move counter
    hasMoved=false;
    nbmove=0;
    // perform a pass
    for (int i=0;i<N ;++i){
      // current node
      cnode=pass(i);
      //Rcout << cnode << "K:" << K << std::endl;
      if (workingset(cnode)==1){
        // compute delta swap
        arma::uvec iclust = possible_moves(cl(cnode),moves_mat);
        arma::vec delta = this->delta_swap(cnode,iclust);
        //Rcout << delta << std::endl;
        // best swap
        int ncl = delta.index_max();
        if(ncl!=cl(cnode)){
          
          // if best swap corresponds to a move
          // update the stats and deal with cluster death
          if(counts(cl(cnode))==1){
              // remove the cluster from the move matrix
              delrowcol(moves_mat,cl(cnode));

          }
          this->swap_update(cnode,ncl);
          
          // update the move counters
          hasMoved=true;
          ++nbmove;
          // workingset(cnode)=1;
        }else{
          arma::vec deltaneg = delta.elem(arma::find(delta<0));
          int bmn= deltaneg.index_max();
          if(deltaneg(bmn) <  -5 ){
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
      if(k == obl | l == obl){
        delta(k,l)=this->delta_merge(k,l);
      }else{
        //Rcout << k <<" : " <<l << " - " << obk <<" : " << obl <<std::endl;
        //Rcout << delta(k,l) << std::endl;
        delta(k,l)=delta(k,l)+this->delta_merge_correction(k,l,obk,obl,old_stats);
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


// init merge matrix sparse
SpMergeMat IclModel::nasty_delta_merge(const arma::sp_mat & merge_graph){
  // inititalize delta merge matrix
  arma::sp_mat delta = merge_graph;
  arma::sp_mat delta_pattern = merge_graph;
  // index to store current best merge
  int bk = 0;
  int bl = 0;
  // initialize bv found to -infty
  double bv = -std::numeric_limits<double>::infinity();
  // store cuurent stats
  arma::sp_mat::iterator i= delta.begin();
  bool move = true;
  while(move){
    move = false;
    for (int d=0;d<delta.n_cols;++d ){
      bv = -std::numeric_limits<double>::infinity();
      for (arma::sp_mat::iterator i = delta.begin_col(d); i != delta.end_col(d); ++i) {
        delta(i.row(),i.col())=this->delta_merge(i.row(),i.col());
        delta(i.col(),i.row())=delta(i.row(),i.col());
        // best merge ?
        if(delta(i.row(),i.col())>bv){
          bv=delta(i.row(),i.col());
          bk=i.row();
          bl=i.col();
        }
      }
      if(bv>0){
        this->merge_update(bk,bl);
        delrowcol(delta,bk);
        move=true;
      }
    }

    
  }
  
  return SpMergeMat(bk,bl,bv,delta);
}




// update merge matrix after merge of obk/obl sparse
SpMergeMat IclModel::delta_merge(const arma::sp_mat & merge_graph, int obk, int obl,const List & old_stats){
  // optimized version to compute only new values of the merge mat
  //delta = delta(arma::find(arma::linspace(0,K,K+1)!=obk),arma::find(arma::linspace(0,K,K+1)!=obk));
  arma::sp_mat deltaO = merge_graph;
  delrowcol(deltaO,obk);
  arma::sp_mat delta = deltaO;
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
      if(k == obl | l == obl){
        delta(k,l)=this->delta_merge(k,l);
      }else{
        //delta(k,l)=delta(k,l)+this->delta_merge_correction(k,l,obk,obl,old_stats);
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

// init merge matrix sparse
arma::sp_mat IclModel::batch_greedy_merge(const arma::sp_mat & merge_graph,int nb_try,double reduction_factor){
  // inititalize delta merge matrix
  arma::sp_mat delta = merge_graph;
  arma::uvec col0;
  col0<< 0 << arma::endr;
  arma::uvec col1;
  col1<< 1 << arma::endr;
  arma::uvec col01;
  col01 << 0 << 1 << arma::endr;
  // index to store current best merge
  int nb_sols=K*reduction_factor;
  int pos,k,l;
  bool move = true;
  while(move){
    move = false;
    // initialize bv found to -infty
    arma::umat best_kls(nb_sols,2);
    arma::vec best_kls_values(nb_sols);
    
    best_kls_values.fill(-std::numeric_limits<double>::infinity());
    
    for (int l=0;l<delta.n_cols;++l ){
      int j = 0;
      arma::vec Ks(delta.col(l).n_nonzero);
      for (arma::sp_mat::iterator i = delta.begin_col(l); i != delta.end_col(l); ++i) {
        Ks(j)=i.row();
        j++;
      }
      arma::vec pass;
      if(nb_try< delta.col(l).n_nonzero){
        pass = as<arma::vec>(sample(delta.col(l).n_nonzero,nb_try))-1;
      }else{
        pass = arma::linspace(0,delta.col(l).n_nonzero-1, delta.col(l).n_nonzero);
      }

      for (int i=0;i<pass.n_elem;i++) {
        int k=Ks(pass(i));
        delta(k,l) = this->delta_merge(k,l);
        // best merge ?
        pos = nb_sols-1;
        while(pos >=0 && delta(k,l)>best_kls_values(pos) ){
          pos--;
        }
        if(pos<(nb_sols-1)){
          //insertion

            if(pos<(nb_sols-3)){
              best_kls_values.subvec(pos+2,best_kls.n_rows-1)=best_kls_values.subvec(pos+1,best_kls.n_rows-2);
              best_kls.submat(pos+2,0,best_kls.n_rows-1,1)=best_kls.submat(pos+1,0,best_kls.n_rows-2,1);
            }

            best_kls(pos+1,0)=k;
            best_kls(pos+1,1)=l;
            best_kls_values(pos+1) = delta(k,l);

      
        }
      }
    }
    pos = 0;

    while(pos<best_kls_values.n_elem && best_kls_values(pos)>0){
      int bk = best_kls(pos,0);
      int bl = best_kls(pos,1);
      k = std::max(bk,bl);
      l = std::min(bk,bl);
      if(k!=l){
        if(delta(k,l)!=best_kls_values(pos)){
          Rcout << k  <<" ; " << l <<" ; " << delta(k,l)<< ";" << best_kls_values(pos) << std::endl;
        }

        this->merge_update(k,l);
        delta.shed_col(k);
        delta.shed_row(k);
        best_kls_values=best_kls_values(arma::find(best_kls.col(0)!=k));
        best_kls=best_kls.submat(arma::find(best_kls.col(0)!=k),col01);
        best_kls_values=best_kls_values(arma::find(best_kls.col(1)!=k));
        best_kls=best_kls.submat(arma::find(best_kls.col(1)!=k),col01);
        //best_kls.submat(arma::find(best_kls.col(0)==k),col0).fill(l);
        //best_kls.submat(arma::find(best_kls.col(1)==k),col1).fill(l);
        best_kls.submat(arma::find(best_kls.col(0)>k),col0) = best_kls.submat(arma::find(best_kls.col(0)>k),col0)-1;
        best_kls.submat(arma::find(best_kls.col(1)>k),col1) = best_kls.submat(arma::find(best_kls.col(1)>k),col1)-1;
        //Rcout << best_kls << std::endl;
        move=true;
      }
      pos++;
    }
    if(move){
       nb_sols = K*reduction_factor;
       nb_sols = std::max(nb_sols,1);
    }
  }
  return delta;
}


// main function for greedy merging with prior merge graph
arma::sp_mat IclModel::greedy_merge(const arma::sp_mat & merge_graph){
  // init the merge matrix(K,K) with the delta icl of each merge 
  SpMergeMat merge_mat = this->delta_merge(merge_graph);
  // init merge counter
  int nbmerge = 0;
  //while a positive merge exists
  while(merge_mat.getValue()>0){

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
    arma::sp_mat delta = merge_mat.getMergeMat();
    merge_mat = this->delta_merge(delta,merge_mat.getK(),merge_mat.getL(),old_stats);

    //delrowcol(delta,merge_mat.getK());
    //merge_mat = this->delta_merge(delta);
  }
  //compute final icl value
  icl_value = icl(this->get_obs_stats());
  if(verbose){
    Rcout << "##################################"<< std::endl;
    Rcout << "Merge convergence, with an ICL of "<< icl_value << " and " << K << " clusters." << std::endl;
    Rcout << "##################################"<< std::endl; 
  }
  return merge_mat.getMergeMat();
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
    // check test
    //MergeMat merge_mat_comp = this->delta_merge();
    //Rcout << arma::max(merge_mat.getMergeMat()-merge_mat_comp.getMergeMat()) << std::endl;
    //Rcout << arma::min(merge_mat.getMergeMat()-merge_mat_comp.getMergeMat()) << std::endl;
    
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
    List old_stats = this->get_obs_stats();
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
                                Named("l")=l+1,
                                Named("merge_mat") = arma::trimatl(merge_mat.getMergeMat())));
    // update merge matrix
    merge_mat = this->delta_merge(merge_mat.getMergeMat(),merge_mat.getK(),merge_mat.getL(),old_stats);
  }
  return path;
}



