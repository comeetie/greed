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
#include "SimpleIclCoModel.h"
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
#include "DcLbm.h"
using namespace Rcpp;




IclModelEmission * init_emission_model(S4 model,List data, arma::uvec clt, bool verbose) {

  IclModelEmission * Memission;
  S4 sol;
  try{
    // SBM
    if(model.is("SbmPrior")){
      if((strcmp(model.slot("type"),"directed")!=0) && (strcmp(model.slot("type"),"undirected")!=0) && (strcmp(model.slot("type"),"guess")!=0)){
        stop("Unsuported model type only directed / undirected are allowed");
      } 
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      if((strcmp(model.slot("type"),"guess")==0)){
        model=clone(model);
        if(arma::accu(abs(xp-xp.t()))==0){
          model.slot("type")="undirected";
        }else{
          model.slot("type")="directed";
        }
      }

      if(strcmp(model.slot("type"),"directed")==0){
        Memission = new Sbm(xp,model,verbose);
      }
      if(strcmp(model.slot("type"),"undirected")==0){
        Memission = new SbmUndirected(xp,model,verbose);

      }
    }
  
  
    // DcSbm
    if(model.is("DcSbmPrior")){
      if((strcmp(model.slot("type"),"directed")!=0) && (strcmp(model.slot("type"),"undirected")!=0) && (strcmp(model.slot("type"),"guess")!=0)){
        stop("Unsuported model type only directed / undirected are allowed");
      } 
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      if((strcmp(model.slot("type"),"guess")==0)){
        model=clone(model);
        if(arma::accu(abs(xp-xp.t()))==0){
          model.slot("type")="undirected";
        }else{
          model.slot("type")="directed";
        }
      }
      if(strcmp(model.slot("type"),"directed")==0){
        Memission = new DcSbm(xp,model,verbose);
      }
      if(strcmp(model.slot("type"),"undirected")==0){
        Memission = new DcSbmUndirected(xp,model,verbose);
      }
    }
    
    
    
    // MultSbm
    if(model.is("MultSbmPrior")){
      if((strcmp(model.slot("type"),"directed")!=0) && (strcmp(model.slot("type"),"undirected")!=0) && (strcmp(model.slot("type"),"guess")!=0)){
        stop("Unsuported model type only directed / undirected are allowed");
      } 
      arma::cube xp = as<arma::cube>(data["X"]);
      
      if((strcmp(model.slot("type"),"guess")==0)){
        int diff = 0;
        model=clone(model);
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
        Memission = new MultSbm(xp,model,verbose);
      }
      if(strcmp(model.slot("type"),"undirected")==0){
        Memission = new MultSbmUndirected(xp,model,verbose);
      }
    }
    
    // DcLbm
    if(model.is("DcLbmPrior")){
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      int Nr = static_cast<int>(data["Nrows"]);
      int Nc = static_cast<int>(data["Ncols"]);
      Memission = new DcLbm(xp,Nr,Nc,model,verbose);
    }
    
    // MoM
    if(model.is("MoMPrior")){
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      Memission = new Mm(xp,model,verbose);
    }
    
    // Gmm
    if(model.is("GmmPrior")){
      arma::mat X = as<arma::mat>(data["X"]);
      Memission = new Gmm(X,model,verbose);
    }
    
    // DiagGmm
    if(model.is("DiagGmmPrior")){
      arma::mat X = as<arma::mat>(data["X"]);
      Memission = new DiagGmm(X,model,verbose);
    }
    
    // MoR
    if(model.is("MoRPrior")){
      arma::mat X = as<arma::mat>(data["X"]);
      arma::mat Y = as<arma::mat>(data["Y"]);
      Memission = new Mregs(X,Y,model,verbose);

    }
    
    // Lca
    if(model.is("LcaPrior")){
      arma::umat X = as<arma::umat>(data["X"]);
      Memission = new Lca(X,model,verbose);
    }
    

    return(Memission);
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }
  return NULL;
}

IclModel * init(S4 model,List data, arma::uvec clr, bool verbose){
  
  int N = clr.n_elem;
  arma::uvec clt = to_zero_based(clr);
  
  IclModel * M;
  if(!model.is("DlvmPrior")){
    stop("Unsuported model");
  }
  if(model.is("DlvmCoPrior")){
    int Nr = static_cast<int>(data["Nrows"]);
    int Nc = static_cast<int>(data["Ncols"]);
    IclModelEmission * m = init_emission_model(model,data, clt, verbose);
    M = new SimpleIclCoModel(m,model,clt,Nr,Nc,verbose);
  }else{
    if(model.is("DlvmPrior")){
      if(model.is("CombinedModels")){
        List models = as<List>(model.slot("models"));
        std::vector<IclModelEmission*> icl_models;
        CharacterVector models_names = models.names();
        for( int m=0;m<models.length();m++){
          String name = models_names[m];
          IclModelEmission * cm = init_emission_model(models[name],data[name], clt, verbose);
          icl_models.push_back(cm);
        }
        M = new CombinedIclModel(icl_models,model,clt,verbose);
      }else{
        IclModelEmission * m = init_emission_model(model,data, clt, verbose);
        M = new SimpleIclModel(m,model,clt,verbose);
      }
    }
  }

  return M;
}




S4 init_sol(S4 model,String type="Fit") {
  String mname = model.attr("class");
  mname+=type;
  S4 sol(mname);
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
  S4 sol = init_sol(model,"Path");
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
// test_swap
// @param model icl_model
// @param data list with clustering data depends on model type
// @param cl initital clustering
// @return a DeltaICL(cl,cl_swap)-DeltaICLOpt(cl_swap)
// [[Rcpp::export]]
double test_swap(S4 model,List data,arma::uvec& cl,int i, int newcl){
  int Ki = arma::max(cl);
  IclModel * M = init(model,data,cl,FALSE);
  double icli = M->icl(M->get_obs_stats());

  arma::uvec iclust = arma::find(arma::ones(Ki));
  arma::vec delta_swap = M->delta_swap(i-1,iclust);

  arma::uvec & cl_swap = cl;
  cl_swap(i-1)=newcl;
  IclModel * Mswap = init(model,data,cl,FALSE);
  double icls = Mswap->icl(Mswap->get_obs_stats());

  double delta = icls-icli;

  return std::abs(delta-delta_swap[newcl-1]);
}

// test_merge
// @param model icl_model
// @param data list with clustering data depends on model type
// @param cl initital clustering
// @return a DeltaICL(cl,cl_swap)-DeltaICLOpt(cl_swap)
// [[Rcpp::export]]
double test_merge(S4 model,List data,arma::uvec& cl,int k, int l){
  if(k<l){
    int temp = k;
    k=l;
    l=temp;
  }
  IclModel * M = init(model,data,cl,FALSE);
  arma::uvec & clm = cl;
  clm(arma::find(clm==k)).fill(l);
  clm.elem(arma::find(clm>k))=cl.elem(arma::find(clm>k))-1;
  k=k-1;
  l=l-1;

  double dmerge = M->delta_merge(k,l);

  IclModel * Mmerge = init(model,data,clm,FALSE);
  
  double delta = Mmerge->icl(Mmerge->get_obs_stats())-M->icl(M->get_obs_stats());

  return std::abs(dmerge-delta);
}


// test_merge_correction
// @param model icl_model
// @param data list with clustering data depends on model type
// @param cl initital clustering
// @return a 
// [[Rcpp::export]]
arma::mat test_merge_correction(S4 model,List data,arma::uvec& cl,int k, int l){
  if(k<l){
    int temp = k;
    k=l;
    l=temp;
  }
  IclModel * M = init(model,data,cl,FALSE);
  
  
  MergeMat init_merge_mat = M->delta_merge();
  List old_stats = M->get_obs_stats();
  M->merge_update(k-1,l-1);
  MergeMat fusopt_merge_mat = M->delta_merge(init_merge_mat.getMergeMat(),k-1,l-1,old_stats);
  
  
  arma::uvec & clm = cl;
  clm(arma::find(clm==k)).fill(l);
  clm.elem(arma::find(clm>k))=cl.elem(arma::find(clm>k))-1;

  
  IclModel * Mmerge = init(model,data,clm,FALSE);
  
  MergeMat fusnoopt_merge_mat = Mmerge->delta_merge(); 
  
  return fusopt_merge_mat.getMergeMat()-fusnoopt_merge_mat.getMergeMat(); 
}

