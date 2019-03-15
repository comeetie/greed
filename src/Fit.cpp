// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include "Sbm.h"
#include "DcSbm.h"
#include "Mm.h"
using namespace Rcpp;

//' fit_greed
//' @param model icl_model
//' @param xp sparseMatrix
//' @param Ki initia guess for K
//' @param type : merge, swap or both (default)  
//' @param nb_max_pass : maximum number of pass for greedy swap 
//' @export
// [[Rcpp::export]]
S4 fit_greed(S4 model,arma::sp_mat& xp, int Ki, std::string type="both", int nb_max_pass = 50,bool verbose = false) {
  
  IclModel * M;
  S4 sol("sbm_fit");
  try{
    if(strcmp(model.slot("name"),"sbm")!=0 && strcmp(model.slot("name"),"dcsbm")!=0 && strcmp(model.slot("name"),"mm")!=0){
      stop("Unsuported model");
    }
    if(strcmp(model.slot("name"),"sbm")==0){
      M = new Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"),verbose);
      S4 solt("sbm_fit");
      solt.slot("name") = "sbm_fit";
      sol = solt;
    }
    if(strcmp(model.slot("name"),"dcsbm")==0){
      M = new DcSbm(xp,Ki,model.slot("alpha"),verbose);
      S4 solt("dcsbm_fit");
      solt.slot("name") = "dcsbm_fit";
      sol = solt;
    }
    if(strcmp(model.slot("name"),"mm")==0){
      M = new Mm(xp,Ki,model.slot("alpha"),model.slot("beta"),verbose);
      S4 solt("mm_fit");
      solt.slot("name") = "mm_fit";
      sol = solt;
    }
    if(type!="merge" && type!="swap" && type!="both"){
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
    return(sol);
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }
}



//' fit_greed_init_type
//' @param model icl_model
//' @param xp sparseMatrix
//' @param clt cluster labels {0,...,K-1}
//' @param type : merge, swap or both (default)  
//' @param nb_max_pass : maximum number of pass for greedy swap 
//' @export
// [[Rcpp::export]]
S4 fit_greed_init(S4 model,arma::sp_mat& xp,  arma::vec& clt, std::string type="both", int nb_max_pass = 50,bool verbose=false) {
  
  IclModel * M;
  int Ki = arma::max(clt);
  int N = clt.n_elem;
  clt = clt-arma::ones(N);
  S4 sol("sbm_fit");
  try{
    if(strcmp(model.slot("name"),"sbm")!=0 && strcmp(model.slot("name"),"dcsbm")!=0 && strcmp(model.slot("name"),"mm")!=0){
      stop("Unsuported model");
    }
    if(strcmp(model.slot("name"),"sbm")==0){
      M = new Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"),clt,verbose);
      S4 solt("sbm_fit");
      solt.slot("name") = "sbm_fit";
      sol = solt;
    }
    if(strcmp(model.slot("name"),"dcsbm")==0){
      M = new DcSbm(xp,Ki,model.slot("alpha"),clt,verbose);
      S4 solt("dcsbm_fit");
      solt.slot("name") = "dcsbm_fit";
      sol = solt;
    }
    if(strcmp(model.slot("name"),"mm")==0){
      M = new Mm(xp,Ki,model.slot("alpha"),model.slot("beta"),clt,verbose);
      S4 solt("mm_fit");
      solt.slot("name") = "mm_fit";
      sol = solt;
    }
    if(type!="merge" && type!="swap" && type!="both"){
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
    return(sol);
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }
}



//' fit_greed 
//' @param xp sparseMatrix
//' @param init initial fit
//' @export
// [[Rcpp::export]]
S4 fit_greed_path(arma::sp_mat& xp, S4 init) {
  S4 model = init.slot("model");
  int Ki = arma::max(as<arma::vec>(init.slot("cl")));
  int N  = as<arma::vec>(init.slot("cl")).n_elem;
  arma::vec clt = as<arma::vec>(init.slot("cl"))-arma::ones(N);
 
  IclModel * M;
  S4 sol("sbm_path");
   try{
     if(strcmp(model.slot("name"),"sbm")!=0 && strcmp(model.slot("name"),"dcsbm")!=0 && strcmp(model.slot("name"),"mm")!=0){
       stop("Unsuported model");
     }
     if(strcmp(model.slot("name"),"sbm")==0){
       M = new Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"),clt);
       S4 solt("sbm_path");
       solt.slot("name") = "sbm_path";
       sol = solt;
     }
     if(strcmp(model.slot("name"),"dcsbm")==0){
       M = new DcSbm(xp,Ki,model.slot("alpha"),clt);
       S4 solt("dcsbm_path");
       solt.slot("name") = "dcsbm_path";
       sol = solt;
     }
     if(strcmp(model.slot("name"),"mm")==0){
       M = new Mm(xp,Ki,model.slot("alpha"),model.slot("beta"),clt);
       S4 solt("mm_path");
       solt.slot("name") = "mm_path";
       sol = solt;
     }
     List obs_stats = M->get_obs_stats();
     sol.slot("obs_stats") = obs_stats;
     sol.slot("cl") = M->get_cl()+1 ;
     sol.slot("icl") = M->icl(obs_stats);
     sol.slot("path") = M->greedy_merge_path();
     return(sol);
   }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }
  
}
