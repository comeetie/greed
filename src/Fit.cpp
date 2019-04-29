// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include "IclModel.h"
#include "Sbm.h"
#include "DcSbm.h"
#include "Mm.h"
#include "Mreg.h"
using namespace Rcpp;


IclModel * init(S4 model,List data, arma::vec& clt) {
  
  IclModel * M;
  int Ki = arma::max(clt);
  int N = clt.n_elem;
  clt = clt-arma::ones(N);
  S4 sol;
  try{
    if(strcmp(model.slot("name"),"sbm")!=0 && 
       strcmp(model.slot("name"),"dcsbm")!=0 && 
       strcmp(model.slot("name"),"mm")!=0 && 
       strcmp(model.slot("name"),"mreg")!=0){
      stop("Unsuported model");
    }
    if(strcmp(model.slot("name"),"sbm")==0){
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      M = new Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"),clt,false);
    }
    if(strcmp(model.slot("name"),"dcsbm")==0){
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      M = new DcSbm(xp,Ki,model.slot("alpha"),clt,false);
    }
    if(strcmp(model.slot("name"),"mm")==0){
      arma::sp_mat xp = as<arma::sp_mat>(data["X"]);
      M = new Mm(xp,Ki,model.slot("alpha"),model.slot("beta"),clt,false);
    }
    if(strcmp(model.slot("name"),"mreg")==0){
      arma::mat X = as<arma::mat>(data["X"]);
      arma::colvec y = as<arma::colvec>(data["y"]);
      M = new Mreg(X,y,Ki,model.slot("alpha"),model.slot("reg"),model.slot("a0"),model.slot("b0"),clt);
    }
    
    
    return(M);
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }
}

S4 init_sol(S4 model,String type="fit") {
  String mname = model.slot("name");
  mname+="_";
  mname+=type;
  S4 sol(mname);
  sol.slot("name")=mname;
  return(sol);
}

//' post_probs
//' @param model icl_model
//' @param xp sparseMatrix
//' @param clt cluster labels {0,...,K-1}
//' @export
// [[Rcpp::export]]
arma::mat post_probs(S4 model,List data,  arma::vec& clt) {
  IclModel * M = init(model,data,clt);
  return(M->get_probs());
}

//' fit_greed_init
//' @param model icl_model
//' @param xp sparseMatrix
//' @param clt cluster labels {0,...,K-1}
//' @param type : merge, swap, none, or both (default)  
//' @param nb_max_pass : maximum number of pass for greedy swap 
//' @export
// [[Rcpp::export]]
S4 fit_greed(S4 model,List data,  arma::vec& clt, std::string type="both", int nb_max_pass = 50,bool verbose=false) {
  IclModel * M = init(model,data,clt);
  S4 sol = init_sol(model);
  if(type!="merge" && type!="swap" && type!="both" && type!="none"){
    stop("Unsuported algorithm");
  }
  if(type=="swap" || type=="both"){
    M->greedy_swap(nb_max_pass);
  }
  if(type=="merge" || type=="both"){
    M->greedy_merge();
  }

  List obs_stats = M->get_obs_stats();
  sol.slot("model") = model;
  sol.slot("obs_stats") = obs_stats;
  sol.slot("cl") = M->get_cl()+1 ;
  sol.slot("icl") = M->icl(obs_stats);
  sol.slot("K") = M->get_K();
  return(sol);
}



//' fit_greed 
//' @param xp sparseMatrix
//' @param init initial fit
//' @export
// [[Rcpp::export]]
S4 fit_greed_path(List data, S4 init_fit) {
  S4 model = init_fit.slot("model");
  arma::vec clt = init_fit.slot("cl");
  IclModel * M = init(model,data,clt);
  S4 sol = init_sol(model,"path");
  List obs_stats = M->get_obs_stats();
  sol.slot("obs_stats") = obs_stats;
  sol.slot("cl") = M->get_cl()+1 ;
  sol.slot("icl") = M->icl(obs_stats);
  sol.slot("K") = M->get_K();
  sol.slot("path") = M->greedy_merge_path();
  return(sol);
  
}
