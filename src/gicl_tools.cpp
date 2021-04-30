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
    if((i==oldcl) || (i==newcl)){
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
arma::sp_mat sp_cross(arma::sp_mat colvec,arma::sp_mat rowvec,int self, int oldcl, int newcl, int K){
  arma::sp_mat result(K,K);
  result.col(oldcl)=result.col(oldcl)-colvec;
  result.col(newcl)=result.col(newcl)+colvec;
  result.row(oldcl)=result.row(oldcl)-rowvec;
  result.row(newcl)=result.row(newcl)+rowvec;
  result(oldcl,oldcl)=result(oldcl,oldcl)-self;
  result(newcl,newcl)=result(newcl,newcl)+self;
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
arma::sp_mat add_spmatpat(const arma::sp_mat & a, const arma::sp_mat & b){
  arma::sp_mat result(a.n_rows,a.n_cols);
  for (arma::sp_mat::const_iterator i = b.begin(); i != b.end(); ++i) {
    result(i.row(),i.col()) = a(i.row(),i.col())+*i;
  }
  return result;
}


// [[Rcpp::export]]
arma::sp_mat which_spmatpat(const arma::sp_mat & a, const arma::sp_mat & b){
  arma::sp_mat result(a.n_rows,a.n_cols);
  for (arma::sp_mat::const_iterator i = b.begin(); i != b.end(); ++i) {
    result(i.row(),i.col()) = a(i.row(),i.col());
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


arma::sp_mat delrowcol_copy(const arma::sp_mat & a, int ci){
//   //a.shed_row(ci);
//   //a.shed_col(ci);
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
    if((i.row()!=ci) & (i.col()!=ci)){
      result(k,l) = a(i.row(),i.col());
    }

  }
  return result;
}



// [[Rcpp::export]]
void delrowcol(arma::sp_mat & a, int ci){
  a.shed_row(ci);
  a.shed_col(ci);
}


// [[Rcpp::export]]
arma::cube gsum_cube(arma::vec cl,const arma::cube& x, int K){
  arma::cube res = arma::cube(K,K,x.n_slices);
  res.fill(0);
  for(int i = 0; i < cl.n_elem; ++i) {
    for(int j = 0; j < cl.n_elem; ++j) {
      res.tube( cl(i), cl(j) )=res.tube( cl(i), cl(j) )+ x.tube(i,j);
    }  
  }
  return res;
}

arma::vec count(arma::vec cl,int K){
  arma::vec result(K);
  result.fill(0);
  for(int i = 0; i < cl.n_elem; ++i) {
    result(cl(i),0)+=1;
  }
  return result;
}


// [[Rcpp::export]]
arma::mat gsum_mat(arma::vec cl,const arma::sp_mat& x,int K) {
  arma::mat result(K,K);
  result.fill(0);
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    result(cl(i.row()),cl(i.col())) += *i;
  }
  return result;
}

// [[Rcpp::export]]
arma::mat gsum_bimat(arma::vec clr,arma::vec clc, const arma::sp_mat& x,int K) {
  arma::mat result(K,K);
  result.fill(0);
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    result(clr(i.row()),clc(i.col())) += *i;
  }
  return result;
}



// [[Rcpp::export]]
arma::sp_mat gsum_mat_sp(arma::vec cl,const arma::sp_mat& x,int K) {
  arma::sp_mat result(K,K);
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
List mvlm_post_comp(const arma::mat X,const arma::mat Y,const arma::mat M,const arma::mat K, const arma::mat S0, double N0) {
  
  
  
  // https://tminka.github.io/papers/minka-linear.pdf
  // Bayesian linear regression
  int n = X.n_rows, m = X.n_cols, d=Y.n_cols;


  arma::mat S = X.t()*X+K;

  arma::mat Xty = X.t()*Y+(M.t()*K).t();

  arma::mat iS = inv_sympd(S);
  arma::mat mu =  iS*Xty;
 
  arma::mat Yty = Y.t()*Y+M.t()*K*M;

  arma::mat Syx = Yty - Y.t()*X*mu;
  
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  double log_evidence = arma::accu(lgamma((n+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));

  log_evidence = log_evidence + d/2*log(det(K))- d/2*log(det(S)) - n*d/2*log(M_PI);

  log_evidence = log_evidence + N0/2*log(det(S0))-(n+N0)/2*log(det(Syx+S0));

  
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
List mvlm_post_add1_comp(List current, const arma::rowvec X,const arma::rowvec Y,const arma::mat M,const arma::mat K, const arma::mat S0, double N0) {
  int m = X.n_cols, d=Y.n_cols;
  int n = as<int>(current["n"])+1;
  
  arma::mat S = as<arma::mat>(current["S"])+X.t()*X;
  arma::mat Xty = as<arma::mat>(current["Xty"])+X.t()*Y;
  arma::mat Yty = as<arma::mat>(current["Yty"])+Y.t()*Y;
  
  arma::mat iSo = as<arma::mat>(current["iS"]);
  // algo de mise a jour sequentiel pour 1 point
  // https://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula
  arma::mat iS = iSo-(iSo*X.t()*X*iSo)/as_scalar(1+X*iSo*X.t()); 
  arma::mat mu =  iS*Xty;
  
  arma::mat Syx = Yty - Xty.t()*mu;
  
  
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  double log_evidence = arma::accu(lgamma((n+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));
  log_evidence = log_evidence + d/2*log(det(K))- d/2*log(det(S)) - n*d/2*log(M_PI);
  log_evidence = log_evidence + N0/2*log(det(S0))-(n+N0)/2*log(det(Syx+S0));
  
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
List mvlm_post_del1_comp(List current, const arma::rowvec X,const arma::rowvec Y,const arma::mat M,const arma::mat K, const arma::mat S0, double N0) {

  int m = X.n_cols, d=Y.n_cols;
  int n = as<int>(current["n"])-1;
  
  arma::mat S = as<arma::mat>(current["S"])-X.t()*X;
  arma::mat Xty = as<arma::mat>(current["Xty"])-X.t()*Y;
  arma::mat Yty = as<arma::mat>(current["Yty"])-Y.t()*Y;
  

  arma::mat iSo = as<arma::mat>(current["iS"]);
  // algo de mise a jour sequentiel pour 1 point
  // https://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula
  arma::mat Xn = -X;
  arma::mat iS = iSo-(iSo*X.t()*Xn*iSo)/as_scalar(1+Xn*iSo*X.t());
  arma::mat mu =  iS*Xty;
  
  arma::mat Syx = Yty - Xty.t()*mu;
  
  
  
  
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  double log_evidence = arma::accu(lgamma((n+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));
  log_evidence = log_evidence + d/2*log(det(K))- d/2*log(det(S)) - n*d/2*log(M_PI);
  log_evidence = log_evidence + N0/2*log(det(S0))-(n+N0)/2*log(det(Syx+S0));
  
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
List mvlm_post_merge_comp(List current1, List current2,const arma::mat M,const arma::mat K, const arma::mat S0, double N0) {

  int m = as<arma::mat>(current1["S"]).n_cols, d=as<arma::mat>(current1["Yty"]).n_cols;
  int n = as<int>(current1["n"])+as<int>(current2["n"]);
  

  
  arma::mat S = as<arma::mat>(current1["S"])+as<arma::mat>(current2["S"])-K;
  arma::mat Xty = as<arma::mat>(current1["Xty"])+as<arma::mat>(current2["Xty"])-(M.t()*K).t();
  arma::mat Yty = as<arma::mat>(current1["Yty"])+as<arma::mat>(current2["Yty"])-M.t()*K*M;
  

  arma::mat iS = inv_sympd(S);
  arma::mat mu =  iS*Xty;
  
  arma::mat Syx = Yty - Xty.t()*mu;
  
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  double log_evidence = arma::accu(lgamma((n+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));
  log_evidence = log_evidence + d/2*log(det(K))- d/2*log(det(S)) - n*d/2*log(M_PI);
  log_evidence = log_evidence + N0/2*log(det(S0))-(n+N0)/2*log(det(Syx+S0));
  
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
List gmm_marginal(const arma::mat X,double tau,int N0i, const arma::mat epsilon, const arma::rowvec mu) {
  
  
  double ng = X.n_rows, d = X.n_cols;
  double N0 = N0i;
  arma::rowvec m = arma::mean(X,0);
  arma::mat M(ng,d);
  M.each_row() = m;
  arma::mat S = (X-M).t()*(X-M); 
  arma::mat Sp = epsilon+tau*ng/(tau+ng)*(m-mu).t()*(m-mu)+S;
  
  
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  double log_evidence = arma::accu(lgamma((ng+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));
  log_evidence = log_evidence - ng*d/2*log(M_PI) + d/2*log(tau) - d/2*log(tau+ng) ;
  log_evidence = log_evidence + N0/2*log(det(epsilon))-(ng+N0)/2*log(det(Sp));
  
  return List::create(Named("S")  = S,
                      Named("m") = m,
                      Named("ng")  = ng,
                      Named("Sp") = Sp,
                      Named("log_evidence")=log_evidence);
}


List gmm_marginal_spherical(const arma::mat X,double kappa,double tau,double beta, const arma::rowvec mu) {
  
  
  double ng = X.n_rows, d = X.n_cols;
  arma::rowvec m = arma::mean(X,0);
  arma::mat M(ng,d);
  M.each_row() = m;
  arma::rowvec S = arma::sum(arma::pow(X-M,2),0);
  arma::rowvec betan = beta +0.5*S + (tau*ng)/(2*(tau+ng))*arma::pow(m-mu,2);
  double taun = tau+ng;
  double kappan = kappa+(ng/2);
  double log_evidence = arma::accu(lgamma(kappan)-lgamma(kappa)+kappa*log(beta)-kappan*log(betan)+0.5*log(tau)-0.5*log(taun)-ng/2*log(2*M_PI));
    
  return List::create(Named("S")  = S,
                      Named("m") = m,
                      Named("ng")  = ng,
                      Named("log_evidence")=log_evidence);
}



double gauss_evidence(double d, double ng,double N0, double tau,double eps,arma::mat epsilon, arma::mat S, arma::rowvec m, arma::rowvec mu){
  arma::mat Sp = eps*epsilon+S;
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  //+tau*ng/(tau+ng)*(m-mu).t()*(m-mu)
  //+ d/2*log(tau) - d/2*log(tau+ng)
  double log_evidence = arma::accu(lgamma((ng+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));
  log_evidence = log_evidence - ng*d/2*log(M_PI)  ;
  log_evidence = log_evidence + N0/2*log(det(eps*epsilon))-(ng+N0)/2*log(det(Sp));
  return log_evidence; 
}



// [[Rcpp::export]]
List gmm_marginal_add1(List current, const arma::rowvec X,double tau,int N0i, const arma::mat epsilon, const arma::rowvec mu) {
  

  double d = X.n_cols;
  
  arma::mat mold = as<arma::mat>(current["m"]);
  double ngold = as<double>(current["ng"]);
  double ng = ngold +1;
  double N0 = N0i;
  arma::rowvec m = (mold*ngold+X)/ng;
  arma::mat S =  as<arma::mat>(current["S"])+(X-m).t()*(X-mold);
  arma::mat Sp = epsilon+tau*ng/(tau+ng)*(m-mu).t()*(m-mu)+S;
  
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  double log_evidence = arma::accu(lgamma((ng+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));
  log_evidence = log_evidence - ng*d/2*log(M_PI) + d/2*log(tau) - d/2*log(tau+ng) ;
  log_evidence = log_evidence + N0/2*log(det(epsilon))-(ng+N0)/2*log(det(Sp));
  
  return List::create(Named("S")  = S,
                      Named("m") = m,
                      Named("ng")  = ng,
                      Named("log_evidence")=log_evidence);
}



List gmm_marginal_spherical_add1(List current, const arma::rowvec X,double kappa,double tau,double beta, const arma::rowvec mu) {
  
  
  arma::mat mold = as<arma::mat>(current["m"]);
  double ngold = as<double>(current["ng"]);
  double ng = ngold +1;
  arma::rowvec m = (mold*ngold+X)/ng;
  arma::rowvec S = as<arma::rowvec>(current["S"])+((X-mold)%(X-m));
  arma::rowvec betan = beta +0.5*S + (tau*ng)/(2*(tau+ng))*arma::pow(m-mu,2);
  double taun = tau+ng;
  double kappan = kappa+(ng/2);
  double log_evidence = arma::accu(lgamma(kappan)-lgamma(kappa)+kappa*log(beta)-kappan*log(betan)+0.5*log(tau)-0.5*log(taun)-ng/2*log(2*M_PI));
  
  return List::create(Named("S")  = S,
                      Named("m") = m,
                      Named("ng")  = ng,
                      Named("log_evidence")=log_evidence);
}




// [[Rcpp::export]]
List gmm_marginal_del1(List current, const arma::rowvec X,double tau,int N0i, const arma::mat epsilon, const arma::rowvec mu) {
  
  
  double d = X.n_cols;
  
  arma::mat mold = as<arma::mat>(current["m"]);
  double ngold = as<int>(current["ng"]);
  double ng = ngold -1;
  double N0 = N0i;
  arma::rowvec m = (mold*ngold-X)/ng;
  arma::mat S =  as<arma::mat>(current["S"])-(X-m).t()*(X-mold);
  arma::mat Sp = epsilon+tau*ng/(tau+ng)*(m-mu).t()*(m-mu)+S;
  
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  double log_evidence = arma::accu(lgamma((ng+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));
  log_evidence = log_evidence - ng*d/2*log(M_PI) + d/2*log(tau) - d/2*log(tau+ng) ;
  log_evidence = log_evidence + N0/2*log(det(epsilon))-(ng+N0)/2*log(det(Sp));
  
  return List::create(Named("S")  = S,
                      Named("m") = m,
                      Named("ng")  = ng,
                      Named("log_evidence")=log_evidence);
}

List gmm_marginal_spherical_del1(List current, const arma::rowvec X,double kappa,double tau,double beta, const arma::rowvec mu) {
  
  
  arma::mat mold = as<arma::mat>(current["m"]);
  double ngold = as<double>(current["ng"]);
  double ng = ngold  - 1;
  arma::rowvec m = (mold*ngold-X)/ng;
  arma::rowvec S = as<arma::rowvec>(current["S"])-((X-mold)%(X-m));
  arma::rowvec betan = beta +0.5*S + (tau*ng)/(2*(tau+ng))*arma::pow(m-mu,2);
  double taun = tau+ng;
  double kappan = kappa+(ng/2);
  double log_evidence = arma::accu(lgamma(kappan)-lgamma(kappa)+kappa*log(beta)-kappan*log(betan)+0.5*log(tau)-0.5*log(taun)-ng/2*log(2*M_PI));

  return List::create(Named("S")  = S,
                      Named("m") = m,
                      Named("ng")  = ng,
                      Named("log_evidence")=log_evidence);
}








// [[Rcpp::export]]
List gmm_marginal_merge(List current1, List current2,double tau,int N0i, const arma::mat epsilon, const arma::rowvec mu) {
  
  
  double N0 = N0i;

  double ng1 = as<double>(current1["ng"]);
  double ng2 = as<double>(current2["ng"]);
  double ng  = ng1+ng2;
  
  arma::rowvec m1 = as<arma::rowvec>(current1["m"]);
  arma::rowvec m2 = as<arma::rowvec>(current2["m"]);
  arma::rowvec m = m1*(ng1/ng)+m2*(ng2/ng);

  double d = m1.n_cols;
  arma::mat S1 =  as<arma::mat>(current1["S"]);
  arma::mat S2 =  as<arma::mat>(current2["S"]);
  arma::mat S = S1+ng1*(m1-m).t()*(m1-m)+S2+ng2*(m2-m).t()*(m2-m);
  
  arma::mat Sp = epsilon+tau*ng/(tau+ng)*(m-mu).t()*(m-mu)+S;
  
  arma::vec di = arma::linspace<arma::vec>(1, d,d);
  double log_evidence = arma::accu(lgamma((ng+N0+1-di)/2)) - arma::accu(lgamma((N0+1-di)/2));
  log_evidence = log_evidence - ng*d/2*log(M_PI) + d/2*log(tau) - d/2*log(tau+ng) ;
  log_evidence = log_evidence + N0/2*log(det(epsilon))-(ng+N0)/2*log(det(Sp));
  
  return List::create(Named("S")  = S,
                      Named("m") = m,
                      Named("ng")  = ng,
                      Named("log_evidence")=log_evidence);
}


List gmm_marginal_spherical_merge(List current1, List current2,double kappa,double tau,double beta, const arma::rowvec mu) {
  
  double ng1 = as<double>(current1["ng"]);
  double ng2 = as<double>(current2["ng"]);
  double ng  = ng1+ng2;
  
  arma::rowvec m1 = as<arma::rowvec>(current1["m"]);
  arma::rowvec m2 = as<arma::rowvec>(current2["m"]);
  arma::rowvec m = m1*(ng1/ng)+m2*(ng2/ng);
  
  double d = m1.n_cols;
  arma::rowvec S1 =  as<arma::rowvec>(current1["S"]);
  arma::rowvec S2 =  as<arma::rowvec>(current2["S"]);
  arma::rowvec S = S1+ng1*arma::pow(m1-m,2)+S2+ng2*arma::pow(m2-m,2);
  
  
  arma::rowvec betan = beta +0.5*S + (tau*ng)/(2*(tau+ng))*arma::pow(m-mu,2);
  double taun = tau+ng;
  double kappan = kappa+(ng/2);
  double log_evidence = arma::accu(lgamma(kappan)-lgamma(kappa)+kappa*log(beta)-kappan*log(betan)+0.5*log(tau)-0.5*log(taun)-ng/2*log(2*M_PI));
  
  return List::create(Named("S")  = S,
                      Named("m") = m,
                      Named("ng")  = ng,
                      Named("log_evidence")=log_evidence);
}




// [[Rcpp::export]] 
arma::uvec possible_moves(int k,arma::sp_mat & move_mat){
  return arma::find(move_mat.col(k));
}




