// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include "Sbm.h"
#include "Mm.h"
using namespace Rcpp;


// [[Rcpp::export]]
List fit_icl(S4 model,arma::sp_mat& xp, int Ki) {
  Sbm alg = Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"));
 // alg.greedy_swap(100);
  alg.greedy_merge();
  List obs_stats = alg.get_obs_stats();
  return List::create(Named("counts", obs_stats["counts"]), 
                      Named("x_counts", obs_stats["x_counts"]), 
                      Named("cl",alg.get_cl()),
                      Named("icl", alg.icl(obs_stats)));
}



// [[Rcpp::export]]
List fit_icl_init(S4 model,arma::sp_mat& xp, int Ki, arma::vec& clt) {
  Sbm alg = Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"),clt);
  alg.greedy_swap(100);
  alg.greedy_merge();
  List obs_stats = alg.get_obs_stats();
  return List::create(Named("counts", obs_stats["counts"]), 
                      Named("x_counts", obs_stats["x_counts"]), 
                      Named("cl",alg.get_cl()),
                      Named("icl", alg.icl(obs_stats)));
}



// [[Rcpp::export]]
List fit_greed_sbm(S4 model,arma::sp_mat& xp, int Ki) {
  Sbm alg = Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"));
  alg.greedy_swap(100);
  alg.greedy_merge();
  List obs_stats = alg.get_obs_stats();
  return List::create(Named("counts", obs_stats["counts"]), 
                      Named("x_counts", obs_stats["x_counts"]), 
                      Named("cl",alg.get_cl()),
                      Named("icl", alg.icl(obs_stats)));
}

// [[Rcpp::export]]
List fit_greed_mm(S4 model,arma::sp_mat& xp, int Ki) {
  Mm alg = Mm(xp,Ki,model.slot("alpha"),model.slot("beta"));
  Rcout << "Lets Go" << std::endl;
  Rcout << "##################################"<< std::endl;
  alg.greedy_swap(100);
  alg.greedy_merge(); 
  List obs_stats = alg.get_obs_stats();
  return List::create(Named("counts", obs_stats["counts"]), 
                      Named("x_counts", obs_stats["x_counts"]), 
                      Named("cl",alg.get_cl()),
                      Named("icl", alg.icl(obs_stats)));
}