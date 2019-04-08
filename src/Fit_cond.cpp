// [[Rcpp::depends(RcppArmadillo)]]
#include "Mreg.h"
using namespace Rcpp;


//' init
//' @param model icl_model
//' @param xp sparseMatrix
//' @param clt cluster labels {0,...,K-1}
//' @export
// [[Rcpp::export]]
S4 init_cond(S4 model,arma::mat& X, arma::colvec& y,  arma::vec& clt) {
  
  IclModel * M;
  int Ki = arma::max(clt);
  int N = clt.n_elem;
  clt = clt-arma::ones(N);
  S4 sol("mreg_fit");
  try{
    if(strcmp(model.slot("name"),"mreg")!=0){
      stop("Unsuported model");
    }
    if(strcmp(model.slot("name"),"mreg")==0){
      M = new Mreg(X,y,Ki,model.slot("alpha"),model.slot("reg"),model.slot("a0"),model.slot("b0"),clt,false);
      S4 solt("mreg_fit");
      solt.slot("name") = "mreg_fit";
      sol = solt;
    }

    
    List obs_stats = M->get_obs_stats();
    sol.slot("model") = model;
    sol.slot("obs_stats") = obs_stats;
    sol.slot("cl") = M->get_cl()+1 ;
    sol.slot("icl") = M->icl(obs_stats);
    sol.slot("K") = M->get_K();
    return(sol);
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }
}

//' fit_greed
//' @param model icl_model
//' @param xp sparseMatrix
//' @param Ki initia guess for K
//' @param type : merge, swap or both (default)  
//' @param nb_max_pass : maximum number of pass for greedy swap 
//' @export
// [[Rcpp::export]]
S4 fit_greed_cond(S4 model,arma::mat& X, arma::colvec& y,int Ki, std::string type="both", int nb_max_pass = 50,bool verbose = false) {
  
  IclModel * M;
  S4 sol("mreg_fit");
  try{
    if(strcmp(model.slot("name"),"mreg")!=0){
      stop("Unsuported model");
    }
    if(strcmp(model.slot("name"),"mreg")==0){
      M = new Mreg(X,y,Ki,model.slot("alpha"),model.slot("reg"),model.slot("a0"),model.slot("b0"),verbose);
      S4 solt("mreg_fit");
      solt.slot("name") = "mreg_fit";
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
    sol.slot("K") = M->get_K();
    return(sol);
  }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }
}



//' fit_greed_init
//' @param model icl_model
//' @param xp sparseMatrix
//' @param clt cluster labels {0,...,K-1}
//' @param type : merge, swap or both (default)  
//' @param nb_max_pass : maximum number of pass for greedy swap 
//' @export
// [[Rcpp::export]]
S4 fit_greed_init_cond(S4 model,arma::mat& X, arma::colvec& y, arma::vec& clt, std::string type="both", int nb_max_pass = 50,bool verbose=false) {
  
  IclModel * M;
  int Ki = arma::max(clt);
  int N = clt.n_elem;
  clt = clt-arma::ones(N);
  S4 sol("mreg_fit");
  try{
    if(strcmp(model.slot("name"),"mreg")!=0){
      stop("Unsuported model");
    }
    if(strcmp(model.slot("name"),"mreg")==0){
      M = new Mreg(X,y,Ki,model.slot("alpha"),model.slot("reg"),model.slot("a0"),model.slot("b0"),clt,verbose);
      S4 solt("mreg_fit");
      solt.slot("name") = "mreg_fit";
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
    sol.slot("K") = M->get_K();
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
S4 fit_greed_path_cond(arma::mat& X,arma::colvec & y, S4 init) {
  S4 model = init.slot("model");
  int Ki = arma::max(as<arma::vec>(init.slot("cl")));
  int N  = as<arma::vec>(init.slot("cl")).n_elem;
  arma::vec clt = as<arma::vec>(init.slot("cl"))-arma::ones(N);
 
  IclModel * M;
  S4 sol("mreg_path");
   try{
     if(strcmp(model.slot("name"),"mreg")!=0){
       stop("Unsuported model");
     }
     if(strcmp(model.slot("name"),"mreg")==0){
       M = new Mreg(X,y,Ki,model.slot("alpha"),model.slot("reg"),model.slot("a0"),model.slot("b0"),clt);
       S4 solt("mreg_path");
       solt.slot("name") = "mreg_path";
       sol = solt;
     }
     List obs_stats = M->get_obs_stats();
     sol.slot("obs_stats") = obs_stats;
     sol.slot("cl") = M->get_cl()+1 ;
     sol.slot("icl") = M->icl(obs_stats);
     sol.slot("K") = M->get_K();
     sol.slot("path") = M->greedy_merge_path();
     return(sol);
   }catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  }
  
}
