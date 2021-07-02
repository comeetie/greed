// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "Partition.h" 
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::vec test_partition(arma::vec cl, int K){
  Partition Part(cl,K);
  Part.swap(1,0);
  Part.swap(0,0);
  Part.swap(2,0);
  Part.erase(5);
  Part.swap(12,0);
  Part.swap(10,0);
  Part.swap(1,5);
  return Part.get_cl();
}


