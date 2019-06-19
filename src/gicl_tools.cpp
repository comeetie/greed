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

// [[Rcpp::export]]
arma::sp_mat add_sppat(const arma::sp_mat & a, const arma::sp_mat & b){
  arma::sp_mat result(a.n_rows,1);
  for (arma::sp_mat::const_iterator i = b.begin(); i != b.end(); ++i) {
    result(i.row(),0) = a(i.row(),0)+*i;
  }
  return result;
}

// [[Rcpp::export]]
arma::sp_mat delcol(const arma::sp_mat & a, int ci){
  arma::sp_mat result(a.n_rows,a.n_cols-1);
  for (arma::sp_mat::const_iterator i = a.begin(); i != a.end(); ++i) {
    if(i.col()<ci){
      result(i.row(),i.col()) = a(i.row(),i.col());
    }
    if(i.col()>ci){
      result(i.row(),i.col()-1) = a(i.row(),i.col());
    }    

  }
  return result;
}

// [[Rcpp::export]]
arma::sp_mat delrowcol(const arma::sp_mat & a, int ci){
  arma::sp_mat result(a.n_rows-1,a.n_cols-1);
  int k=0;
  int l=0;
  for (arma::sp_mat::const_iterator i = a.begin(); i != a.end(); ++i) {
    k = i.row();
    if(i.row()>ci){
      k --; 
    }
    l=i.col();
    if(i.col()>ci){
      l--;
    }
    if(i.row()!=ci & i.col()!=ci){
      result(k,l) = a(i.row(),i.col());
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

// [[Rcpp::export]]
arma::sp_mat gsum_mm(arma::vec cl,const arma::sp_mat& x,int K) {
  arma::sp_mat result(x.n_rows,K);
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    result(i.row(),cl(i.col())) += *i;
  }
  return result;
}

arma::sp_mat gsum_col(arma::vec cl,const arma::sp_mat& x,int i, int K) {
  arma::sp_mat ccol = x.col(i);
  arma::sp_mat result(K,1);
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


// [[Rcpp::export]]
List mvlm_post(const arma::mat X,const arma::mat Y,double alpha, double N0) {
  int n = X.n_rows, m = X.n_cols, d=Y.n_cols;
  
  
  arma::mat Sprior(m,m);
  Sprior.zeros();
  Sprior.diag() = arma::ones<arma::vec>(m)*alpha;
  
  
  arma::mat SMprior(d,d);
  SMprior.zeros();
  SMprior.diag() = arma::ones<arma::vec>(d)*N0;
  
  arma::mat S = X.t()*X+Sprior;
  arma::mat Xty = X.t()*Y;
  arma::mat iS = inv_sympd(S);
  arma::mat mu =  iS*Xty;
  arma::mat Yty = Y.t()*Y;
  arma::mat Syx = Yty - Y.t()*X*mu;
  
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  double log_evidence = arma::accu(lgamma((n+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));
  log_evidence = log_evidence + m*d/2*log(alpha) - d/2*log(det(S)) - n*d/2*log(M_PI);
  log_evidence = log_evidence + N0*d/2*log(N0)-(n+N0)/2*log(det(Syx+SMprior));
   // log_evidence = arma::accu(di);
    //+ m*d/2*log(alpha) - d/2*log(det(S)) - n*d/2*log(M_PI)+N0*d/2*log(N0)-(n+N0)/2*log(det(Syx+SMprior));
  return List::create(Named("S")  = S,
                      Named("mu") = mu,
                      Named("n")  = n,
                      Named("Xty")  = Xty,
                      Named("Yty")  = Yty,
                      Named("Syx")=Syx,
                      Named("iS")=iS,
                      Named("log_evidence")=log_evidence);
}


// [[Rcpp::export]]
List mvlm_post_add1(List current, const arma::rowvec X,const arma::rowvec Y,double alpha, double N0) {
  int m = X.n_cols, d=Y.n_cols;
  int n = as<int>(current["n"])+1;
  
  arma::mat S = as<arma::mat>(current["S"])+X.t()*X;
  arma::mat Xty = as<arma::mat>(current["Xty"])+X.t()*Y;
  arma::mat Yty = as<arma::mat>(current["Yty"])+Y.t()*Y;
  
  
  arma::mat SMprior(d,d);
  SMprior.zeros();
  SMprior.diag() = arma::ones<arma::vec>(d)*N0;
  
  
  arma::mat iSo = as<arma::mat>(current["iS"]);
  // algo de mise a jour sequentiel pour 1 point
  // https://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula
  arma::mat iS = iSo-(iSo*X.t()*X*iSo)/as_scalar(1+X*iSo*X.t()); 
  arma::mat mu =  iS*Xty;
  
  arma::mat Syx = Yty - Xty.t()*mu;
  
  
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  double log_evidence = arma::accu(lgamma((n+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));
  log_evidence = log_evidence + m*d/2*log(alpha) - d/2*log(det(S)) - n*d/2*log(M_PI);
  log_evidence = log_evidence + N0*d/2*log(N0)-(n+N0)/2*log(det(Syx+SMprior));
  return List::create(Named("S")  = S,
                      Named("mu") = mu,
                      Named("n")  = n,
                      Named("Xty")  = Xty,
                      Named("Yty")  = Yty,
                      Named("Syx")=Syx,
                      Named("iS")=iS,
                      Named("log_evidence")=log_evidence);
}

// [[Rcpp::export]]
List mvlm_post_del1(List current, const arma::rowvec X,const arma::rowvec Y,double alpha, double N0) {
  int m = X.n_cols, d=Y.n_cols;
  int n = as<int>(current["n"])-1;
  
  arma::mat S = as<arma::mat>(current["S"])-X.t()*X;
  arma::mat Xty = as<arma::mat>(current["Xty"])-X.t()*Y;
  arma::mat Yty = as<arma::mat>(current["Yty"])-Y.t()*Y;
  
  
  arma::mat SMprior(d,d);
  SMprior.zeros();
  SMprior.diag() = arma::ones<arma::vec>(d)*N0;
  
  
  arma::mat iSo = as<arma::mat>(current["iS"]);
  // algo de mise a jour sequentiel pour 1 point
  // https://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula
  arma::mat Xn = -X;
  arma::mat iS = iSo-(iSo*X.t()*Xn*iSo)/as_scalar(1+Xn*iSo*X.t());
  arma::mat mu =  iS*Xty;
  
  arma::mat Syx = Yty - Xty.t()*mu;
  
  
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  double log_evidence = arma::accu(lgamma((n+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));
  log_evidence = log_evidence + m*d/2*log(alpha) - d/2*log(det(S)) - n*d/2*log(M_PI);
  log_evidence = log_evidence + N0*d/2*log(N0)-(n+N0)/2*log(det(Syx+SMprior));
  return List::create(Named("S")  = S,
                      Named("mu") = mu,
                      Named("n")  = n,
                      Named("Xty")  = Xty,
                      Named("Yty")  = Yty,
                      Named("Syx")=Syx,
                      Named("iS")=iS,
                      Named("log_evidence")=log_evidence);
}

// [[Rcpp::export]]
List mvlm_post_merge(List current1, List current2,double alpha, double N0) {
  int m = as<arma::mat>(current1["S"]).n_cols, d=as<arma::mat>(current1["Yty"]).n_cols;
  int n = as<int>(current1["n"])+as<int>(current2["n"]);
  
  arma::mat Sprior(m,m);
  Sprior.zeros();
  Sprior.diag() = arma::ones<arma::vec>(m)*alpha;
  
  
  arma::mat S = as<arma::mat>(current1["S"])+as<arma::mat>(current2["S"])-Sprior;
  arma::mat Xty = as<arma::mat>(current1["Xty"])+as<arma::mat>(current2["Xty"]);
  arma::mat Yty = as<arma::mat>(current1["Yty"])+as<arma::mat>(current2["Yty"]);
  
  
  arma::mat SMprior(d,d);
  SMprior.zeros();
  SMprior.diag() = arma::ones<arma::vec>(d)*N0;
  
  arma::mat iS = inv_sympd(S);
  arma::mat mu =  iS*Xty;
  
  arma::mat Syx = Yty - Xty.t()*mu;
  
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  double log_evidence = arma::accu(lgamma((n+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));
  log_evidence = log_evidence + m*d/2*log(alpha) - d/2*log(det(S)) - n*d/2*log(M_PI);
  log_evidence = log_evidence + N0*d/2*log(N0)-(n+N0)/2*log(det(Syx+SMprior));
  return List::create(Named("S")  = S,
                      Named("mu") = mu,
                      Named("n")  = n,
                      Named("Xty")  = Xty,
                      Named("Yty")  = Yty,
                      Named("Syx")=Syx,
                      Named("iS")=iS,
                      Named("log_evidence")=log_evidence);
}


