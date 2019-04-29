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

double sum_lfact(const arma::sp_mat & x){
  arma::sp_mat::const_iterator start = x.begin();
  arma::sp_mat::const_iterator end   = x.end();
  double cst = 0;
  for(arma::sp_mat::const_iterator it = start; it != end; ++it)
  {
     cst += 1; 
  }
  return cst;
}

//' lm_post
//' @param X data matrix of covariates Nxd
//' @param y target Nx1
//' @param regu prior precision parameter
//' @param a0 prior parameter
//' @param b0 prior parameter
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
  
  double log_evidence = -n/2*log(2*M_PI)-0.5*d*log(regu)+0.5*log(det(S))+a0*log(b0)-a*log(b)+lgamma(a)-lgamma(a0);
  return List::create(Named("S")  = S,
                      Named("mu") = mu,
                      Named("a")  = a,
                      Named("b")  = b,
                      Named("n")  = n,
                      Named("yty")  = yty,
                      Named("Xty")  = Xty,
                      Named("log_evidence")=log_evidence,
                      Named("detS")=det(S),
                      Named("iS")=inv_sympd(S));
}
  
  
  //' lm_post_add1
  //' @param current gaussian linear model to update
  //' @param X data matrix of covariates 1xd
  //' @param y target 1x1
  //' @param regu prior precision parameter
  //' @param a0 prior parameter
  //' @param b0 prior parameter
  //' @export
  // [[Rcpp::export]]
  List lm_post_add1(List current, const arma::rowvec X,double y,double regu, double a0, double b0) {
    

    
    int n = as<int>(current["n"])+X.n_rows;
    int d = as<arma::mat>(current["S"]).n_cols;
    arma::mat S = as<arma::mat>(current["S"])+X.t()*X;
    arma::mat iSo = as<arma::mat>(current["iS"]);
    // algo de mise a jour sequentiel pour 1 point
    // https://en.wikipedia.org/wiki/Matrix_determinant_lemma
    double detS = as_scalar(1+X*iSo*X.t())*as<double>(current["detS"]);
    
    // https://en.wikipedia.org/wiki/Woodbury_matrix_identity
    // https://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula
    arma::mat iS = iSo-(iSo*X.t()*X*iSo)/as_scalar(1+X*iSo*X.t()); 
    arma::colvec Xty = as<arma::colvec>(current["Xty"])+X.t()*y;
    arma::colvec mu =  iS*Xty;
    
    
    
    
    double yty = as<double>(current["yty"])+y*y;
    double a = a0+n/2;
    double b = b0+0.5*(yty-arma::as_scalar(mu.t()*S*mu));
    double log_evidence = -n/2*log(2*M_PI)-0.5*d*log(regu)+0.5*log(detS)+a0*log(b0)-a*log(b)+lgamma(a)-lgamma(a0);
    return List::create(Named("S")  = S,
                        Named("mu") = mu,
                        Named("a")  = a,
                        Named("b")  = b,
                        Named("n")  = n,
                        Named("yty")  = yty,
                        Named("Xty")  = Xty,
                        Named("log_evidence")=log_evidence,
                        Named("detS")=detS,
                        Named("iS")=iS);
  }



//' lm_post_del1
//' @param current gaussian linear model to update
//' @param X data matrix of covariates 1xd
//' @param y target 1x1
//' @param regu prior precision parameter
//' @param a0 prior parameter
//' @param b0 prior parameter
//' @export
// [[Rcpp::export]]
List lm_post_del1(List current, const arma::rowvec X,double y,double regu, double a0, double b0) {
  int n = as<int>(current["n"])-X.n_rows;
  int d = as<arma::mat>(current["S"]).n_cols;
  arma::mat S = as<arma::mat>(current["S"])-X.t()*X;
  arma::mat iSo = as<arma::mat>(current["iS"]);
  // algo de mise a jour sequentiel pour 1 point
  // https://en.wikipedia.org/wiki/Matrix_determinant_lemma
  double detS = as_scalar(1-X*iSo*X.t())*as<double>(current["detS"]);
  
  // https://en.wikipedia.org/wiki/Woodbury_matrix_identity
  // https://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula
  arma::mat Xn = -X;
  arma::mat iS = iSo-(iSo*X.t()*Xn*iSo)/as_scalar(1+Xn*iSo*X.t()); 
  arma::colvec Xty = as<arma::colvec>(current["Xty"])-X.t()*y;
  arma::colvec mu =  iS*Xty;
  
  double yty = as<double>(current["yty"])-y*y;
  double a = a0+n/2;
  double b = b0+0.5*(yty-arma::as_scalar(mu.t()*S*mu));
  double log_evidence = -n/2*log(2*M_PI)-0.5*d*log(regu)+0.5*log(detS)+a0*log(b0)-a*log(b)+lgamma(a)-lgamma(a0);
  return List::create(Named("S")  = S,
                      Named("mu") = mu,
                      Named("a")  = a,
                      Named("b")  = b,
                      Named("n")  = n,
                      Named("yty")  = yty,
                      Named("Xty")  = Xty,
                      Named("log_evidence")=log_evidence,
                      Named("detS")=det(S),
                      Named("iS")=inv_sympd(S));
}

//' lm_post_merge
//' @param current_k gaussian linear model to merge
//' @param current_l gaussian linear model to merge
//' @param regu prior precision parameter
//' @param a0 prior parameter
//' @param b0 prior parameter
//' @export
// [[Rcpp::export]]
List lm_post_merge(List current_k,List current_l,double regu, double a0, double b0) {
  int n = as<int>(current_k["n"])+as<int>(current_l["n"]);
  int d = as<arma::mat>(current_k["S"]).n_cols;
  arma::mat Sprior(d,d);
  Sprior.zeros();
  Sprior.diag() = arma::ones<arma::vec>(d)*regu;
  arma::mat S = as<arma::mat>(current_k["S"])+as<arma::mat>(current_l["S"])-Sprior;
  arma::colvec Xty = as<arma::colvec>(current_k["Xty"])+as<arma::colvec>(current_l["Xty"]);
  arma::colvec mu =  inv_sympd(S)*Xty;
  
  double yty = as<double>(current_k["yty"])+as<double>(current_l["yty"]);
  double a = a0+n/2;
  double b = b0+0.5*(yty-arma::as_scalar(mu.t()*S*mu));
  double log_evidence = -n/2*log(2*M_PI)-0.5*d*log(regu)+0.5*log(det(S))+a0*log(b0)-a*log(b)+lgamma(a)-lgamma(a0);
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
//' @param current gaussian linear model to update
//' @param X data matrix of covariates Ntxd
//' @param y target Ntx1
//' @param regu prior precision parameter
//' @param a0 prior parameter
//' @param b0 prior parameter
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
  double log_evidence = -n/2*log(2*M_PI)-0.5*d*log(regu)+0.5*log(det(S))+a0*log(b0)-a*log(b)+lgamma(a)-lgamma(a0);
  return List::create(Named("S")  = S,
                      Named("mu") = mu,
                      Named("a")  = a,
                      Named("b")  = b,
                      Named("n")  = n,
                      Named("yty")  = yty,
                      Named("Xty")  = Xty,
                      Named("log_evidence")=log_evidence);
}



//' lm_post_del
//' @param current gaussian linear model to update
//' @param X data matrix of covariates Ntxd
//' @param y target Ntx1
//' @param regu prior precision parameter
//' @param a0 prior parameter
//' @param b0 prior parameter
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
  double log_evidence = -n/2*log(2*M_PI)-0.5*d*log(regu)+0.5*log(det(S))+a0*log(b0)-a*log(b)+lgamma(a)-lgamma(a0);
  return List::create(Named("S")  = S,
                      Named("mu") = mu,
                      Named("a")  = a,
                      Named("b")  = b,
                      Named("n")  = n,
                      Named("yty")  = yty,
                      Named("Xty")  = Xty,
                      Named("log_evidence")=log_evidence);
}
