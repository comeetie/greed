// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat submatcross(int oldcl,int newcl,int K){
  arma::mat result(4*(K-1),2);
  int nbr = 0;
  result.fill(0);
  for(int i = 0; i < K; ++i) {
    result(i,0)=oldcl;
    result(i,1)=i;
    result(i+K,0)=newcl;
    result(i+K,1)=i;
    if(i==oldcl | i==newcl){
      nbr = nbr + 1;
    }else{
      result(i+2*K-nbr,1)=oldcl;
      result(i+2*K-nbr,0)=i;
      result(i+(3*K-2)-nbr,1)=newcl;
      result(i+(3*K-2)-nbr,0)=i;
    }
  }

  return result;
}

arma::vec count(arma::vec cl,int K){
  arma::vec result(K);
  result.fill(0);
  for(int i = 0; i < cl.n_elem; ++i) {
    result(cl(i),0)+=1;
  }
  return result;
}

arma::mat gsum_mat(arma::vec cl,const arma::sp_mat& x,int K) {
  arma::mat result(K,K);
  result.fill(0);
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    result(cl(i.row()),cl(i.col())) += *i;
  }
  return result;
}


arma::mat gsum_mm(arma::vec cl,const arma::sp_mat& x,int K) {
  arma::mat result(K,x.n_cols);
  result.fill(0);
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    result(cl(i.row()),i.col()) += *i;
  }
  return result;
}

arma::mat gsum_col(arma::vec cl,const arma::sp_mat& x,int i, int K) {
  arma::sp_mat ccol = x.col(i);
  arma::mat result(K,1);
  result.fill(0);
  for (arma::sp_mat::const_iterator i = ccol.begin(); i != ccol.end(); ++i) {
    result(cl(i.row()),0) += *i;
  }
  return result;
}


arma::mat update_count(arma::vec counts,int oldcl,int newcl) {
  counts(oldcl)=counts(oldcl)-1;
  counts(newcl)=counts(newcl)+1;
  return counts;
}


arma::vec sum_cols(const arma::sp_mat& x){
  arma::vec tots(x.n_cols);
  tots.fill(0);
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i)  {
    tots(i.col()) += *i;
  }
  return tots;
}

arma::vec sum_cols(const arma::mat x){
  int nbcols = x.n_cols;
  int nbrows = x.n_rows;
  arma::vec tots(x.n_cols);
  tots.fill(0);
  for (int j=0;j<nbcols;++j)  {
    for (int i=0;i<nbrows;++i)  {
      tots(j) += x(i,j);
    }
  }
  return tots;
}

arma::vec sum_rows(const arma::mat x){
  int nbcols = x.n_cols;
  int nbrows = x.n_rows;
  arma::vec tots(x.n_cols);
  tots.fill(0);
  for (int i=0;i<nbrows;++i) {
    for (int j=0;j<nbcols;++j) {
      tots(i) += x(i,j);
    }
  }
  return tots;
}

//' lm_post
//' @param X
//' @param y
//' @param regu
//' @param a0
//' @param b0
//' @export
// [[Rcpp::export]]
List lm_post(const arma::mat X,const arma::colvec& y,double regu, double a0, double b0) {
  int n = X.n_rows, d = X.n_cols;
  arma::mat Sprior(d,d);
  Sprior.zeros();
  Sprior.diag() = arma::ones<arma::vec>(d)*regu;
    
  arma::mat S = X.t()*X+Sprior;
  arma::colvec Xty = X.t()*y;
  arma::colvec mu =  inv_sympd(S)*Xty;
  
  double a = a0+n/2;
  double yty = arma::as_scalar(y.t()*y);
  double b = b0+0.5*(yty-arma::as_scalar(mu.t()*S*mu));
  
  double log_evidence = -n/2*log(2*M_PI)-0.5*d*log(regu)+0.5*det(S)+a0*log(b0)-a*log(b)+lgamma(a)-lgamma(a0);
  return List::create(Named("S")  = S,
                      Named("mu") = mu,
                      Named("a")  = a,
                      Named("b")  = b,
                      Named("n")  = n,
                      Named("yty")  = yty,
                      Named("Xty")  = Xty,
                      Named("log_evidence")=log_evidence);
}
  
  
  //' lm_post_add
  //' @param X
  //' @param y
  //' @param regu
  //' @param a0
  //' @param b0
  //' @export
  // [[Rcpp::export]]
  List lm_post_add(List current, const arma::mat X,const arma::colvec& y,double regu, double a0, double b0) {
    int n = as<int>(current["n"])+X.n_rows;
    int d = as<arma::mat>(current["S"]).n_cols;
    arma::mat S = as<arma::mat>(current["S"])+X.t()*X;
    arma::colvec Xty = as<arma::colvec>(current["Xty"])+X.t()*y;
    arma::colvec mu =  inv_sympd(S)*Xty;
    
    double yty = as<double>(current["yty"])+arma::as_scalar(y.t()*y);
    double a = a0+n/2;
    double b = b0+0.5*(yty-arma::as_scalar(mu.t()*S*mu));
    double log_evidence = -n/2*log(2*M_PI)-0.5*d*log(regu)+0.5*det(S)+a0*log(b0)-a*log(b)+lgamma(a)-lgamma(a0);
    return List::create(Named("S")  = S,
                        Named("mu") = mu,
                        Named("a")  = a,
                        Named("b")  = b,
                        Named("n")  = n,
                        Named("yty")  = yty,
                        Named("Xty")  = Xty,
                        Named("log_evidence")=log_evidence);
  }



//' lm_post_add
//' @param X
//' @param y
//' @param regu
//' @param a0
//' @param b0
//' @export
// [[Rcpp::export]]
List lm_post_del(List current, const arma::mat X,const arma::colvec& y,double regu, double a0, double b0) {
  int n = as<int>(current["n"])-X.n_rows;
  int d = as<arma::mat>(current["S"]).n_cols;
  arma::mat S = as<arma::mat>(current["S"])-X.t()*X;
  arma::colvec Xty = as<arma::colvec>(current["Xty"])-X.t()*y;
  arma::colvec mu =  inv_sympd(S)*Xty;
  
  double yty = as<double>(current["yty"])-arma::as_scalar(y.t()*y);
  double a = a0+n/2;
  double b = b0+0.5*(yty-arma::as_scalar(mu.t()*S*mu));
  double log_evidence = -n/2*log(2*M_PI)-0.5*d*log(regu)+0.5*det(S)+a0*log(b0)-a*log(b)+lgamma(a)-lgamma(a0);
  return List::create(Named("S")  = S,
                      Named("mu") = mu,
                      Named("a")  = a,
                      Named("b")  = b,
                      Named("n")  = n,
                      Named("yty")  = yty,
                      Named("Xty")  = Xty,
                      Named("log_evidence")=log_evidence);
}



