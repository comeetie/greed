// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include "Sbm.h"
#include "Mm.h"
using namespace Rcpp;


// [[Rcpp::export]]
List fit_greed(S4 model,arma::sp_mat& xp, int Ki) {

  if(strcmp(model.slot("name"),"sbm")==0){

    Sbm alg = Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"));
    alg.greedy_swap(100);
    alg.greedy_merge();
    List obs_stats = alg.get_obs_stats();
    return List::create(Named("counts", obs_stats["counts"]), 
                        Named("x_counts", obs_stats["x_counts"]), 
                        Named("cl",alg.get_cl()+1),
                        Named("icl", alg.icl(obs_stats)));
  }
  if(strcmp(model.slot("name"),"mm")==0){

    Mm alg = Mm(xp,Ki,model.slot("alpha"),model.slot("beta"));
    alg.greedy_swap(100);
    alg.greedy_merge();
    List obs_stats = alg.get_obs_stats();
    return List::create(Named("counts", obs_stats["counts"]), 
                        Named("x_counts", obs_stats["x_counts"]), 
                        Named("cl",alg.get_cl()+1),
                        Named("icl", alg.icl(obs_stats)));
  }

  return List::create();
}



// [[Rcpp::export]]
List fit_greed_init(S4 model,arma::sp_mat& xp, int Ki, arma::vec& clt) {
  
  if(strcmp(model.slot("name"),"sbm")==0){
    
    Sbm alg = Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"),clt);
    alg.greedy_swap(100);
    alg.greedy_merge();
    List obs_stats = alg.get_obs_stats();
    return List::create(Named("counts", obs_stats["counts"]), 
                        Named("x_counts", obs_stats["x_counts"]), 
                        Named("cl",alg.get_cl()),
                        Named("icl", alg.icl(obs_stats)));
  }
  if(strcmp(model.slot("name"),"mm")==0){
    
    Mm alg = Mm(xp,Ki,model.slot("alpha"),model.slot("beta"),clt);
    alg.greedy_swap(100);
    alg.greedy_merge();
    List obs_stats = alg.get_obs_stats();
    return List::create(Named("counts", obs_stats["counts"]), 
                        Named("x_counts", obs_stats["x_counts"]), 
                        Named("cl",alg.get_cl()),
                        Named("icl", alg.icl(obs_stats)));
  }
  
  return List::create();
}


