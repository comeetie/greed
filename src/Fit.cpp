// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include "IclModel.h"
#include "MergeMat.h"
#include "Sbm.h"
#include "SbmUndirected.h"
#include "MarSbm.h"
#include "MarSbmUndirected.h"
#include "NmarBdSbm.h"
#include "NmarBdSbmUndirected.h"
#include "DcSbm.h"
#include "DcSbmUndirected.h"
#include "MultSbm.h"
#include "MultSbmUndirected.h"
#include "CoDcSbm.h"
#include "Mm.h"
#include "Gmm.h"
#include "SphericalGmm.h"
#include "Mvmregcomp.h"
using namespace Rcpp;


IclModel * init(S4 model,List data, arma::vec clt, bool verbose) {
  
  IclModel * M;
  int N = clt.n_elem;
  clt = clt-arma::ones(N);
  S4 sol;
  try{
    if((strcmp(model.slot("name"),"sbm")!=0) && 
       (strcmp(model.slot("name"),"misssbm")!=0) && 
       (strcmp(model.slot("name"),"dcsbm")!=0) &&
       (strcmp(model.slot("name"),"multsbm")!=0) &&
       (strcmp(model.slot("name"),"co_dcsbm")!=0) && 
       (strcmp(model.slot("name"),"mm")!=0) &&
       (strcmp(model.slot("name"),"mvmreg")!=0) &&
       (strcmp(model.slot("name"),"diaggmm")!=0) &&
       (strcmp(model.slot("name"),"gmm")!=0)){
       stop("Unsuported model");
    }
    if(strcmp(model.slot("name"),"sbm")==0){
      if((strcmp(model.slot("type"),"directed")!=0) && (strcmp(model.slot("type"),"undirected")!=0)){
        stop("Unsuported model type only directed / undirected are allowed");
      } 
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      if(strcmp(model.slot("type"),"directed")==0){
        M = new Sbm(xp,model,clt,verbose);
      }
      if(strcmp(model.slot("type"),"undirected")==0){
        M = new SbmUndirected(xp,model,clt,verbose);
      }
    }
    if(strcmp(model.slot("name"),"misssbm")==0){
      if((strcmp(model.slot("type"),"directed")!=0) && (strcmp(model.slot("type"),"undirected")!=0)){
        stop("Unsuported model type only directed / undirected are allowed");
      } 
      if((strcmp(model.slot("sampling"),"dyad")!=0) && (strcmp(model.slot("sampling"),"block-dyad")!=0)){
        stop("Unsuported sampling scheme only  'dyad' / 'block-dyad' are allowed");
      } 
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      arma::sp_mat xpobs = as<arma::sp_mat>(data["Xobs"]);
      if(strcmp(model.slot("sampling"),"dyad")==0){
        if(strcmp(model.slot("type"),"directed")==0){
          M = new MarSbm(xp,xpobs,model,clt,verbose);
        }
        if(strcmp(model.slot("type"),"undirected")==0){
          M = new MarSbmUndirected(xp,xpobs,model,clt,verbose);
        }
      }
      if(strcmp(model.slot("sampling"),"block-dyad")==0){
        if(strcmp(model.slot("type"),"directed")==0){
          M = new NmarBdSbm(xp,xpobs,model,clt,verbose);
        }
        if(strcmp(model.slot("type"),"undirected")==0){
          M = new NmarBdSbmUndirected(xp,xpobs,model,clt,verbose);
        }
      }
    }
    if(strcmp(model.slot("name"),"dcsbm")==0){
      if((strcmp(model.slot("type"),"directed")!=0) && (strcmp(model.slot("type"),"undirected")!=0)){
        stop("Unsuported model type only directed / undirected are allowed");
      } 
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      if(strcmp(model.slot("type"),"directed")==0){
        M = new DcSbm(xp,model,clt,verbose);
      }
      if(strcmp(model.slot("type"),"undirected")==0){
        M = new DcSbmUndirected(xp,model,clt,verbose);
      }
    }
    if(strcmp(model.slot("name"),"multsbm")==0){
      if((strcmp(model.slot("type"),"directed")!=0) && (strcmp(model.slot("type"),"undirected")!=0)){
        stop("Unsuported model type only directed / undirected are allowed");
      } 
      arma::cube xp = as<arma::cube>(data["X"]);
      if(strcmp(model.slot("type"),"directed")==0){
        M = new MultSbm(xp,model,clt,verbose);
      }
      if(strcmp(model.slot("type"),"undirected")==0){
        M = new MultSbmUndirected(xp,model,clt,verbose);
      }
    }
    if(strcmp(model.slot("name"),"co_dcsbm")==0){
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      int Nr = static_cast<int>(data["Nrows"]);
      int Nc = static_cast<int>(data["Ncols"]);
      M = new CoDcSbm(xp,Nr,Nc,model,clt,verbose);
    }
    if(strcmp(model.slot("name"),"mm")==0){
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      M = new Mm(xp,model,clt,verbose);
    }
    if(strcmp(model.slot("name"),"gmm")==0){
      arma::mat X = as<arma::mat>(data["X"]);
      M = new Gmm(X,model,clt,verbose);
    }
    
    if(strcmp(model.slot("name"),"diaggmm")==0){
      arma::mat X = as<arma::mat>(data["X"]);
      M = new SphericalGmm(X,model,clt,verbose);
    }
    if(strcmp(model.slot("name"),"mvmreg")==0 ){
      arma::mat X = as<arma::mat>(data["X"]);
      arma::mat Y = as<arma::mat>(data["Y"]);
      M = new Mvmregcomp(X,Y,model,clt,verbose);
    }
    
    
    return(M);
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }
  return NULL;
}

S4 init_sol(S4 model,String type="fit") {
  String mname = model.slot("name");
  mname+="_";
  mname+=type;
  S4 sol(mname);
  sol.slot("name")=mname;
  return(sol);
}
// post_probs
// @param model icl_model
// @param data list with clustering data (fields depend on model type)
// @param clt cluster labels in 1,..,K //' @export
// [[Rcpp::export]]
arma::mat post_probs(S4 model,List data,  arma::vec& clt) {
  IclModel * M = init(model,data,clt,false);
  arma::mat probs = M->get_probs();
  delete M;
  return(probs);
}

// fit_greed_cstr
// @param model icl_model
// @param data list with clustering data (fileds depend on model type)
// @param clt cluster labels {0,...,K-1}
// @param workingset 
// @param iclust 
// @param type: merge, swap, none, or both (default)  
// @param nb_max_pass maximum number of pass for greedy swap
// @param verbose boolean for verbose mode default to false
// @return a model_fit object  //' @export
// [[Rcpp::export]]
S4 fit_greed_cstr(S4 model,List data,  arma::vec& clt,arma::vec workingset,arma::uvec iclust, std::string type="both", int nb_max_pass = 50,bool verbose=false) {
  IclModel * M = init(model,data,clt,verbose);
  S4 sol = init_sol(model);
  if(type!="merge" && type!="swap" && type!="both" && type!="none"){
    stop("Unsuported algorithm");
  }
  if(type=="swap" || type=="both"){
    M->greedy_swap(nb_max_pass,workingset,iclust-1);
  }
  if(type=="merge" || type=="both"){
    M->greedy_merge();
  }

  List obs_stats = M->get_obs_stats();
  List obs_stats_cst = M->get_obs_stats_cst();
  sol.slot("obs_stats_cst") = obs_stats_cst;
  sol.slot("model") = M->get_model();
  sol.slot("obs_stats") = obs_stats;
  sol.slot("cl") = M->get_cl()+1 ;
  sol.slot("icl") = M->icl(obs_stats);
  sol.slot("K") = M->get_K();
  delete M;
  return(sol);
}



// [[Rcpp::export]]
S4 merge_cstr(S4 model,List data,  arma::vec& clt,arma::sp_mat & merge_graph,bool verbose=false) {
  IclModel * M = init(model,data,clt,verbose);
  S4 sol = init_sol(model);
  arma::sp_mat move_mat = M->greedy_merge(merge_graph);
  List obs_stats = M->get_obs_stats();
  
  List obs_stats_cst = M->get_obs_stats_cst();
  sol.slot("obs_stats_cst") = obs_stats_cst;
  double bicl = M->icl(obs_stats);
  sol.slot("model") = M->get_model();
  sol.slot("obs_stats") = obs_stats;
  sol.slot("cl") = M->get_cl()+1 ;
  sol.slot("icl") = bicl;
  sol.slot("move_mat") = move_mat;
  sol.slot("K") = M->get_K();

  delete M;
  return(sol);
}




// [[Rcpp::export]]
S4 swap_cstr(S4 model,List data,  arma::vec& clt,arma::sp_mat & move_mat, int nb_max_pass = 50, bool verbose=false) {
  
  
  IclModel * M = init(model,data,clt,verbose);
  S4 sol = init_sol(model);
  int N = clt.n_elem;
  arma::vec workingset = arma::ones(N);
  M->greedy_swap(nb_max_pass,workingset,move_mat);
  
  List obs_stats = M->get_obs_stats();
  
  List obs_stats_cst = M->get_obs_stats_cst();
  sol.slot("obs_stats_cst") = obs_stats_cst;
  double bicl = M->icl(obs_stats);
  sol.slot("model") = M->get_model();
  sol.slot("obs_stats") = obs_stats;
  sol.slot("cl") = M->get_cl()+1 ;
  sol.slot("icl") = bicl;
  sol.slot("K") = M->get_K();
  
  delete M;
  return(sol);
}


// fit_greed
// @param model icl_model
// @param data list with clustering data (fileds depend on model type)
// @param clt cluster labels {0,...,K-1}
// @param type merge, swap, none, or both (default)  
// @param nb_max_pass maximum number of pass for greedy swap
// @param verbose boolean for verbose mode default to false
// @return a model_fit object  //' @export
// [[Rcpp::export]]
S4 fit_greed(S4 model,List data,  arma::vec& clt,std::string type="both", int nb_max_pass = 50,bool verbose=false) {
  int N = clt.n_elem;
  arma::vec workingset = arma::ones(N);
  int Ki = arma::max(clt);
  arma::uvec iclust = arma::find(arma::ones(Ki))+1;
  S4 sol = fit_greed_cstr(model,data,clt,workingset,iclust,type,nb_max_pass,verbose);
  return(sol);
}

// fit_greed_path
// @param data list with clustering data depnds on model type
// @param init_fit initial fit object
// @return a model_path object //' @export
// [[Rcpp::export]]
S4 fit_greed_path(List data, S4 init_fit) {
  S4 model = init_fit.slot("model");
  arma::vec clt = init_fit.slot("cl");
  IclModel * M = init(model,data,clt,false);
  S4 sol = init_sol(model,"path");
  sol.slot("model")=model;
  List obs_stats = M->get_obs_stats();
  sol.slot("obs_stats") = obs_stats;
  List obs_stats_cst = M->get_obs_stats_cst();
  sol.slot("obs_stats_cst") = obs_stats_cst;
  sol.slot("cl") = M->get_cl()+1 ;
  sol.slot("icl") = M->icl(obs_stats);
  sol.slot("K") = M->get_K();
  sol.slot("path") = M->greedy_merge_path();
  delete M;
  return(sol);
  
}

// merge_mat
// @param data list with clustering data depnds on model type
// @param init_fit initial fit object
// @return a cost merge matrix //' @export
// [[Rcpp::export]]
arma::mat merge_mat(List data, S4 init_fit) {
  S4 model = init_fit.slot("model");
  arma::vec clt = init_fit.slot("cl");
  IclModel * M = init(model,data,clt,false);
  MergeMat merge_mat = M->delta_merge();
  arma::mat mm = merge_mat.getMergeMat();
  delete M;
  return(mm);
}



bool test_swap(List data, S4 model,arma::vec& clt) {
  //IclModel * M = init(model,data,clt,false);
  //List obs_stats = M->get_obs_stats();
  //double icl_init =  M->icl(obs_stats);
  return(true);
}
