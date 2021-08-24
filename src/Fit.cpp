// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>


// #include "MarSbm.h"
// #include "MarSbmUndirected.h"
// #include "NmarBdSbm.h"
// #include "NmarBdSbmUndirected.h"
// #include "NmarBcSbm.h"
// #include "NmarBcSbmUndirected.h"
// #include "MarDyadSbmE.h"
// #include "MissSbmE.h"







#include "MergeMat.h"
#include "gicl_tools.h"
#include "IclModel.h"
#include "IclModelEmission.h"
#include "CombinedIclModel.h"
#include "SimpleIclModel.h"
#include "Gmm.h"
#include "DiagGmm.h"
#include "Lca.h"
#include "DcSbm.h"
#include "DcSbmUndirected.h"
#include "Sbm.h"
#include "SbmUndirected.h"
#include "MultSbm.h"
#include "MultSbmUndirected.h"
#include "Mm.h"
#include "Mregs.h"
#include "CoDcSbm.h"

using namespace Rcpp;


IclModel * init(S4 model,List data, arma::uvec clr, bool verbose) {
  
  IclModel * M;
  int N = clr.n_elem;
  arma::uvec clt = to_zero_based(clr);
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
       (strcmp(model.slot("name"),"gmm")!=0) && 
       (strcmp(model.slot("name"),"lca")!=0) &&
       (strcmp(model.slot("name"),"mmm")!=0)){
       stop("Unsuported model");
    }
    if(strcmp(model.slot("name"),"sbm")==0){
      if((strcmp(model.slot("type"),"directed")!=0) && (strcmp(model.slot("type"),"undirected")!=0) && (strcmp(model.slot("type"),"guess")!=0)){
        stop("Unsuported model type only directed / undirected are allowed");
      } 
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      if((strcmp(model.slot("type"),"guess")==0)){
        if(arma::accu(abs(xp-xp.t()))==0){
          model.slot("type")="undirected";
        }else{
          model.slot("type")="directed";
        }
      }

      if(strcmp(model.slot("type"),"directed")==0){
        IclModelEmission * sbm = new Sbm(xp,model,verbose);
        M = new SimpleIclModel(sbm,model,clt,verbose);
      }
      if(strcmp(model.slot("type"),"undirected")==0){
        IclModelEmission * sbm = new SbmUndirected(xp,model,verbose);
        M = new SimpleIclModel(sbm,model,clt,verbose);
      }
    }
    if(strcmp(model.slot("name"),"misssbm")==0){
      if((strcmp(model.slot("type"),"directed")!=0) && (strcmp(model.slot("type"),"undirected")!=0) && (strcmp(model.slot("type"),"guess")!=0)){
        stop("Unsuported model type only directed / undirected are allowed");
      } 
      
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      arma::sp_mat xpobs = as<arma::sp_mat>(data["Xobs"]);
      if((strcmp(model.slot("type"),"guess")==0)){
        if(arma::accu(abs(xp-xp.t()))==0 && arma::accu(abs(xpobs-xpobs.t()))==0){
          model.slot("type")="undirected";
        }else{
          model.slot("type")="directed";
        }
      } 
      if((strcmp(model.slot("sampling"),"dyad")!=0) && (strcmp(model.slot("sampling"),"block-dyad")!=0) && (strcmp(model.slot("sampling"),"block-node")!=0)){
        stop("Unsuported sampling scheme only  'dyad' / 'block-dyad / 'block-node' are allowed");
      } 
      if(strcmp(model.slot("sampling"),"dyad")==0){
        if(strcmp(model.slot("type"),"directed")==0){
  //        IclModelEmission * Mobs = new MissSbmE(xp,xpobs,model,verbose);
  //        IclModelEmission * Msampling = new MarDyadSbmE(xpobs,model,verbose);
  //        std::vector<IclModelEmission*> IclModels;
  //        IclModels.push_back(Mobs);
  //        IclModels.push_back(Msampling);
  //        M = new CombinedIclModel(IclModels,model,clt,verbose);
          //M = new MarSbm(xp,xpobs,model,clt,verbose);
          
        }
        if(strcmp(model.slot("type"),"undirected")==0){
    //      M = new MarSbmUndirected(xp,xpobs,model,clt,verbose);
        }
      }
      if(strcmp(model.slot("sampling"),"block-dyad")==0){

        if(strcmp(model.slot("type"),"directed")==0){
      //    M = new NmarBdSbm(xp,xpobs,model,clt,verbose);
        }
        if(strcmp(model.slot("type"),"undirected")==0){
      //    M = new NmarBdSbmUndirected(xp,xpobs,model,clt,verbose);
        }
      }
      if(strcmp(model.slot("sampling"),"block-node")==0){
        if(strcmp(model.slot("type"),"directed")==0){
        //  M = new NmarBcSbm(xp,xpobs,model,clt,verbose);
        }
        if(strcmp(model.slot("type"),"undirected")==0){
      //    M = new NmarBcSbmUndirected(xp,xpobs,model,clt,verbose);
        }
      }
    }
    if(strcmp(model.slot("name"),"dcsbm")==0){
      if((strcmp(model.slot("type"),"directed")!=0) && (strcmp(model.slot("type"),"undirected")!=0) && (strcmp(model.slot("type"),"guess")!=0)){
        stop("Unsuported model type only directed / undirected are allowed");
      } 
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      if((strcmp(model.slot("type"),"guess")==0)){
        if(arma::accu(abs(xp-xp.t()))==0){
          model.slot("type")="undirected";
        }else{
          model.slot("type")="directed";
        }
      }
      if(strcmp(model.slot("type"),"directed")==0){
        IclModelEmission * dcsbm = new DcSbm(xp,model,verbose);
        M = new SimpleIclModel(dcsbm,model,clt,verbose);
      }
      if(strcmp(model.slot("type"),"undirected")==0){
        IclModelEmission * dcsbm = new DcSbmUndirected(xp,model,verbose);
        M = new SimpleIclModel(dcsbm,model,clt,verbose);
      }
    }
    if(strcmp(model.slot("name"),"multsbm")==0){
      if((strcmp(model.slot("type"),"directed")!=0) && (strcmp(model.slot("type"),"undirected")!=0) && (strcmp(model.slot("type"),"guess")!=0)){
        stop("Unsuported model type only directed / undirected are allowed");
      } 
      arma::cube xp = as<arma::cube>(data["X"]);
      
      if((strcmp(model.slot("type"),"guess")==0)){
        int diff = 0;
        for (arma::uword s=0;s<xp.n_slices;s++){
          diff+=arma::accu(abs(xp.slice(s)-xp.slice(s).t()));
        }
        if(diff==0){
          model.slot("type")="undirected";
        }else{
          model.slot("type")="directed";
        }
      }
      
      if(strcmp(model.slot("type"),"directed")==0){
        IclModelEmission * mult = new MultSbm(xp,model,verbose);
        M = new SimpleIclModel(mult,model,clt,verbose);
      }
      if(strcmp(model.slot("type"),"undirected")==0){
        IclModelEmission * mult = new MultSbmUndirected(xp,model,verbose);
        M = new SimpleIclModel(mult,model,clt,verbose);
      }
    }
    if(strcmp(model.slot("name"),"co_dcsbm")==0){
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      int Nr = static_cast<int>(data["Nrows"]);
      int Nc = static_cast<int>(data["Ncols"]);
      IclModelEmission * codcsbm = new CoDcSbm(xp,Nr,Nc,model,verbose);
      M = new SimpleIclModel(codcsbm,model,clt,verbose);
    }
    if(strcmp(model.slot("name"),"mm")==0){
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      IclModelEmission * mm = new Mm(xp,model,verbose);
      M = new SimpleIclModel(mm,model,clt,verbose);
    }
    if(strcmp(model.slot("name"),"gmm")==0){
      arma::mat X = as<arma::mat>(data["X"]);
      IclModelEmission * gmm = new Gmm(X,model,verbose);
      M = new SimpleIclModel(gmm,model,clt,verbose);
    }
    
    if(strcmp(model.slot("name"),"diaggmm")==0){
      arma::mat X = as<arma::mat>(data["X"]);
      IclModelEmission * dgmm = new DiagGmm(X,model,verbose);
      M = new SimpleIclModel(dgmm,model,clt,verbose);
    }
    if(strcmp(model.slot("name"),"mvmreg")==0 ){
      arma::mat X = as<arma::mat>(data["X"]);
      arma::mat Y = as<arma::mat>(data["Y"]);
      IclModelEmission * mregs = new Mregs(X,Y,model,verbose);
      M = new SimpleIclModel(mregs,model,clt,verbose);
    }
    
    if(strcmp(model.slot("name"),"lca")==0 ){
      arma::umat X = as<arma::umat>(data["X"]);
      IclModelEmission * lca = new Lca(X,model,verbose);
      M = new SimpleIclModel(lca,model,clt,verbose);
    }
    if(strcmp(model.slot("name"),"mmm")==0 ){
      
    std::vector<IclModelEmission*> icl_models;
    arma::umat Xcat = as<arma::umat>(data["Xcat"]);
    IclModelEmission * lca = new Lca(Xcat,model,verbose);
    icl_models.push_back(lca);
    arma::mat Xnum = as<arma::mat>(data["Xnum"]);
    IclModelEmission * gmm = new Gmm(Xnum,model,verbose);
    icl_models.push_back(gmm);
    M = new CombinedIclModel(icl_models,model,clt,verbose);
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
arma::mat post_probs(S4 model,List data,  arma::uvec& clt) {
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
S4 fit_greed_cstr(S4 model,List data,  arma::uvec& clt,arma::vec workingset,arma::uvec iclust, std::string type="both", int nb_max_pass = 50,bool verbose=false) {
  IclModel * M = init(model,data,clt,verbose);
  S4 sol = init_sol(model);
  if(type!="merge" && type!="swap" && type!="both" && type!="none"){
    stop("Unsuported algorithm");
  }
  if(type=="swap" || type=="both"){
    //dynamic_cast<SimpleIclModel*>(M)->greedy_swap(nb_max_pass,workingset,iclust-1);

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
S4 merge_cstr(S4 model,List data,  arma::uvec& clt,arma::sp_mat & merge_graph,bool verbose=false) {
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
S4 swap_cstr(S4 model,List data,  arma::uvec& clt,arma::sp_mat & move_mat, int nb_max_pass = 50, bool verbose=false) {
  
  
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
S4 fit_greed(S4 model,List data,  arma::uvec& clt,std::string type="both", int nb_max_pass = 50,bool verbose=false) {
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
  arma::uvec clt = init_fit.slot("cl");
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
  arma::uvec clt = init_fit.slot("cl");
  IclModel * M = init(model,data,clt,false);
  MergeMat merge_mat = M->delta_merge();
  arma::mat mm = merge_mat.getMergeMat();
  delete M;
  return(mm);
}


