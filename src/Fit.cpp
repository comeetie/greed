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
//' @export
// [[Rcpp::export]]
S4 fit_greed(S4 model,arma::sp_mat& xp, int Ki) {

  if(strcmp(model.slot("name"),"sbm")==0){

    Sbm alg = Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"));
    alg.greedy_swap(25);
    alg.greedy_merge();
    List obs_stats = alg.get_obs_stats();
    S4 sol("sbm_fit");
    sol.slot("model") = model;
    sol.slot("name") = "sbm_fit";
    sol.slot("obs_stats") = obs_stats;
    sol.slot("cl") = alg.get_cl()+1 ;
    sol.slot("icl") = alg.icl(obs_stats);
    Rcout << "Run :" << alg.icl(obs_stats) << std::endl;
    return(sol);
  }
  if(strcmp(model.slot("name"),"dcsbm")==0){
    
    DcSbm alg = DcSbm(xp,Ki,model.slot("alpha"));
    alg.greedy_swap(25);
    alg.greedy_merge();
    List obs_stats = alg.get_obs_stats();
    S4 sol("dcsbm_fit");
    sol.slot("model") = model;
    sol.slot("name") = "dcsbm_fit";
    sol.slot("obs_stats") = obs_stats;
    sol.slot("cl") = alg.get_cl()+1 ;
    sol.slot("icl") = alg.icl(obs_stats);
    Rcout << "Run :" << alg.icl(obs_stats) << std::endl;
    return(sol);
  }
  if(strcmp(model.slot("name"),"mm")==0){

    Mm alg = Mm(xp,Ki,model.slot("alpha"),model.slot("beta"));
    alg.greedy_swap(25);
    alg.greedy_merge();
    List obs_stats = alg.get_obs_stats();
    S4 sol("mm_fit");
    sol.slot("model") = model;
    sol.slot("name") = "mm_fit";
    sol.slot("obs_stats") = obs_stats;
    sol.slot("cl") = alg.get_cl()+1 ;
    sol.slot("icl") = alg.icl(obs_stats);
    Rcout << "Run :" << alg.icl(obs_stats) << std::endl;
    return(sol);
  }

}



//' fit_greed_init
//' @param model icl_model
//' @param xp sparseMatrix
//' @param Ki initia guess for K
//' @param clt luster labels 
//' @export
// [[Rcpp::export]]
S4 fit_greed_init(S4 model,arma::sp_mat& xp, int Ki, arma::vec& clt) {
  
  if(strcmp(model.slot("name"),"sbm")==0){
    
    Sbm alg = Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"),clt);
    //alg.greedy_swap(100);
    alg.greedy_merge();
    List obs_stats = alg.get_obs_stats();
    S4 sol("sbm_fit");
    sol.slot("model") = model;
    sol.slot("name") = "sbm_fit";
    sol.slot("obs_stats") = obs_stats;
    sol.slot("cl") = alg.get_cl()+1 ;
    sol.slot("icl") = alg.icl(obs_stats);
    Rcout << "Run :" << alg.icl(obs_stats) << std::endl;
    return(sol);
  }
  if(strcmp(model.slot("name"),"dcsbm")==0){
    
    DcSbm alg = DcSbm(xp,Ki,model.slot("alpha"),clt);
    //alg.greedy_swap(100);
    alg.greedy_merge();
    List obs_stats = alg.get_obs_stats();
    S4 sol("dcsbm_fit");
    sol.slot("model") = model;
    sol.slot("name") = "dcsbm_fit";
    sol.slot("obs_stats") = obs_stats;
    sol.slot("cl") = alg.get_cl()+1 ;
    sol.slot("icl") = alg.icl(obs_stats);
    Rcout << "Run :" << alg.icl(obs_stats) << std::endl;
    return(sol);
  }
  if(strcmp(model.slot("name"),"mm")==0){
    
    Mm alg = Mm(xp,Ki,model.slot("alpha"),model.slot("beta"),clt);
    //alg.greedy_swap(100);
    alg.greedy_merge();
    List obs_stats = alg.get_obs_stats();
    S4 sol("mm_fit");
    sol.slot("name") = "sbm_fit";
    sol.slot("model") = model;
    sol.slot("obs_stats") = obs_stats;
    sol.slot("cl") = alg.get_cl()+1 ;
    sol.slot("icl") = alg.icl(obs_stats);
    Rcout << "Run :" << alg.icl(obs_stats) << std::endl;
    return(sol);
  }
}



//' fit_greed 
//' @param xp sparseMatrix
//' @param init initial fit
//' @export
// [[Rcpp::export]]
S4 fit_greed_path(arma::sp_mat& xp, S4 init) {
  
  if(strcmp(init.slot("name"),"sbm_fit")==0){
    S4 model = init.slot("model");
    int Ki = arma::max(as<arma::vec>(init.slot("cl")));
    int N  = as<arma::vec>(init.slot("cl")).n_elem;
    arma::vec clt = as<arma::vec>(init.slot("cl"))-arma::ones(N);
    Sbm alg = Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"),clt);
    List obs_stats = alg.get_obs_stats();
    S4 sol("sbm_path");
    sol.slot("name") = "sbm_path";
    sol.slot("obs_stats") = obs_stats;
    sol.slot("cl") = alg.get_cl()+1 ;
    sol.slot("icl") = alg.icl(obs_stats);
    sol.slot("path") = alg.greedy_merge_path();
    return(sol);
  }
  if(strcmp(init.slot("name"),"dcsbm_fit")==0){
    S4 model = init.slot("model");
    int Ki = arma::max(as<arma::vec>(init.slot("cl")));
    int N  = as<arma::vec>(init.slot("cl")).n_elem;
    arma::vec clt = as<arma::vec>(init.slot("cl"))-arma::ones(N);
    DcSbm alg = DcSbm(xp,Ki,model.slot("alpha"),clt);
    List obs_stats = alg.get_obs_stats();
    S4 sol("dcsbm_path");
    sol.slot("name") = "dcsbm_path";
    sol.slot("obs_stats") = obs_stats;
    sol.slot("cl") = alg.get_cl()+1 ;
    sol.slot("icl") = alg.icl(obs_stats);
    sol.slot("path") = alg.greedy_merge_path();
    return(sol);
  }
  if(strcmp(init.slot("name"),"mm_fit")==0){
    S4 model = init.slot("model");
    int Ki = arma::max(as<arma::vec>(init.slot("cl")));
    int N  = as<arma::vec>(init.slot("cl")).n_elem;
    arma::vec clt = as<arma::vec>(init.slot("cl"))-arma::ones(N);
    Mm alg = Mm(xp,Ki,model.slot("alpha"),model.slot("beta"),clt);
    List obs_stats = alg.get_obs_stats();
    S4 sol("mm_path");
    sol.slot("name") = "mm_path";
    sol.slot("obs_stats") = obs_stats;
    sol.slot("cl") = alg.get_cl()+1 ;
    sol.slot("icl") = alg.icl(obs_stats);
    sol.slot("path") = alg.greedy_merge_path();
    return(sol);
  }
  
}

//' fit_greed_init_type
//' @param model icl_model
//' @param xp sparseMatrix
//' @param Ki initia guess for K
//' @param clt luster labels 
//' @export
// [[Rcpp::export]]
S4 fit_greed_init_type(S4 model,arma::sp_mat& xp, int Ki, arma::vec& clt, std::string type) {
  
  IclModel * M;
  S4 sol("sbm_fit");
  try{
    if(strcmp(model.slot("name"),"sbm")!=0 && strcmp(model.slot("name"),"dcsbm")!=0 && strcmp(model.slot("name"),"mm")!=0){
      stop("Unsuported model");
    }
    if(strcmp(model.slot("name"),"sbm")==0){
      M = new Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"),clt);
      S4 solt("sbm_fit");
      solt.slot("name") = "sbm_fit";
      sol = solt;
    }
    if(strcmp(model.slot("name"),"dcsbm")==0){
      M = new DcSbm(xp,Ki,model.slot("alpha"),clt);
      S4 solt("dcsbm_fit");
      solt.slot("name") = "dcsbm_fit";
      sol = solt;
    }
    if(strcmp(model.slot("name"),"mm")==0){
      M = new Mm(xp,Ki,model.slot("alpha"),model.slot("beta"),clt);
      S4 solt("mm_fit");
      solt.slot("name") = "mm_fit";
      sol = solt;
    }
    Rcout << type << std::endl;
    if(type!="merge" && type!="swap" && type!="both"){
      stop("Unsuported algorithm");
    }
    if(type=="swap" || type=="both"){
      M->greedy_swap(100);
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
