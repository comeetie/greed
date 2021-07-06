// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "IclModel.h"
#include "Partition.h"
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


  N = cli.n_elem; 
  K = arma::max(cli)+1;
  clp = Partition(cli,K);
  counts = count(cli,K);
  for(int i=0;i<IclModels.size();i++){
    IclModelEmission * Mp = IclModels[i];
    Mp->set_cl(cli);
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
  return icl_prop;
}

arma::mat CombinedIclModel::delta_prop_swap(int i,arma::uvec iclust){
  int oldcl = clp.get(i);
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
    delta += Mp->delta_swap(i,K,clp,iclust);
  }

  delta += delta_prop_swap(i,iclust);

  return delta;
}


void CombinedIclModel::swap_update(int i,int newcl){
  bool dead_cluster = false;
  if(counts(clp.get(i))==1){
    dead_cluster=true;
  }
  for(int m=0;m<IclModels.size();m++){
    IclModelEmission * Mp = IclModels[m];
    Mp->swap_update(i,clp,dead_cluster,newcl);
  }
  int oldcl = clp.get(i);
  counts = update_count(counts,oldcl,newcl);
  clp.swap(i,newcl);
  if(counts(oldcl)==0){
    counts.shed_row(oldcl);
    clp.erase(oldcl);
    --K;
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
  clp.swap(k,l);
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



// main function for greedy swaping
void CombinedIclModel::greedy_swap(int nbpassmax, arma::vec workingset,arma::uvec iclust){

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

        // best swap
        int ncl = delta.index_max();
        int currentcl = clp.get(cnode); 
        // if best swap corresponds to a move
        // update the stats
        if(ncl!=currentcl){
          // one cluster will die
          if(counts(currentcl)==1){
            iclust = iclust.elem(arma::find(iclust!=currentcl));
            iclust.elem(arma::find(iclust>currentcl))=iclust.elem(arma::find(iclust>currentcl))-1;
          }
          this->swap_update(cnode,ncl);
          
          // update the move counters
          hasMoved=true;
          ++nbmove;
          if(K==1){
            break;
          }

        }else{
          arma::vec deltaneg = delta.elem(arma::find(delta<0));
          int bmn= deltaneg.index_max();
          if(deltaneg(bmn) <  -4 ){
            workingset(cnode) = 0;
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
void CombinedIclModel::greedy_swap(int nbpassmax, arma::vec workingset,arma::sp_mat & moves_mat){
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
      // Rcout << cnode << "K:" << K << std::endl;
      if (workingset(cnode)==1){
        // compute delta swap
        arma::uvec iclust = possible_moves(clp.get(cnode),moves_mat);
        Rcout << "0" << std::endl;
        arma::vec delta = this->delta_swap(cnode,iclust);
        Rcout << "3" << std::endl;
        Rcout << delta << std::endl;
        //Rcout << delta << std::endl;
        // best swap
        int ncl = delta.index_max();
        // if best swap corresponds to a move
        // update the stats and deal with cluster death
        if(ncl!=clp.get(cnode)){
          // one cluster will die
          if(counts(clp.get(cnode))==1){
            // remove the cluster from the move matrix
            moves_mat=delrowcol_copy(moves_mat,clp.get(cnode));
          }
          this->swap_update(cnode,ncl);
          
          // update the move counters
          hasMoved=true;
          ++nbmove;
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
    
  }
  if(verbose){
    Rcout << "##################################"<< std::endl;
    //plaRcout << "Swap convergence in " << nbpass << " epochs with " << nbmove << " moves, icl :" << icl_value << "K :" << K << ", working set size :" << arma::accu(workingset)  << std::endl;
    Rcout << "swap convergence, with an ICL of "<< icl_value << " and " << K << " clusters." << std::endl;
    Rcout << "##################################"<< std::endl; 
  } 
}


// get p(z_i|X,z-i)
arma::mat CombinedIclModel::get_probs(){
  
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
MergeMat CombinedIclModel::delta_merge(){
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
MergeMat CombinedIclModel::delta_merge(arma::mat delta, int obk, int obl, const List & old_stats){
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
SpMergeMat CombinedIclModel::delta_merge(const arma::sp_mat & merge_graph){
  
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
SpMergeMat CombinedIclModel::delta_merge(arma::sp_mat & merge_graph, int obk, int obl,const List & old_stats){
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
        delta(k,l)=delta(k,l)+this->delta_merge_correction(k,l,obk,obl,old_stats);
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
arma::sp_mat CombinedIclModel::greedy_merge(const arma::sp_mat & merge_graph){
  // init the merge matrix(K,K) with the delta icl of each merge 
  SpMergeMat merge_mat = this->delta_merge(merge_graph);
  arma::sp_mat delta = merge_mat.getMergeMat();
  // init merge counter
  int nbmerge = 0;
  double cicl = this->icl(this->get_obs_stats());
  double bicl= cicl;
  arma::sp_mat best_merge_mat = delta;
  arma::vec bcl = cl;
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
void CombinedIclModel::greedy_merge(){
  
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
List CombinedIclModel::greedy_merge_path(){
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
                                Named("merge_mat") = arma::trimatl(merge_mat.getMergeMat())));
    // update merge matrix
    merge_mat = this->delta_merge(merge_mat.getMergeMat(),merge_mat.getK(),merge_mat.getL(),old_stats);
  }
  return path;
}


