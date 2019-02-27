
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include "Sbm.h"
using namespace Rcpp;


// [[Rcpp::export]]
List fit_icl(S4 model,arma::sp_mat& xp, int Ki) {
  Sbm alg = Sbm(xp,Ki,model.slot("alpha"),model.slot("a0"),model.slot("b0"));
  alg.greedy_swap(100);
  alg.greedy_merge();
  List obs_stats = alg.get_obs_stats();
  return List::create(Named("counts", obs_stats["counts"]), Named("x_counts", obs_stats["x_counts"]), Named("cl",alg.get_cl()));;
}
