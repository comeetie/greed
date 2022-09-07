// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "DcSbmUndirected.h"
using namespace Rcpp;



void DcSbmUndirected::set_cl(arma::uvec clt){
  K = arma::max(clt)+1;
  x_counts = gsum_mat(clt,x,K);
  counts = count(clt,K);
  din = sum(x_counts).t();
  dout = sum(x_counts.t()).t();
  double cst_denom = 0;
  int nlinks =0;
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    cst_denom+=lgamma(*i+1);
    nlinks += *i;
  }
  cst = arma::accu(lgamma(sum(x)+1))-cst_denom/2 + (nlinks/2)*log(p);
}

double DcSbmUndirected::icl_emiss(const List & obs_stats){

  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::vec din =as<arma::vec>(obs_stats["din"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::mat matcount = counts*counts.t();
  // lets go

  double icl_emiss = accu(lgamma(counts)-lgamma(counts+din)+din % log(counts));
  edges_counts.diag()=edges_counts.diag()/2;
  matcount.diag() = (matcount.diag()-counts)/2;
  arma::mat cmat = lgamma(edges_counts+1)-(edges_counts+1) % log(p*matcount+1);
  arma::uvec lonely = arma::find(counts==1);
  for (arma::uword i = 0;i<lonely.n_elem;++i){
    cmat(lonely(i),lonely(i))=0;
  }
  icl_emiss=icl_emiss + arma::accu(arma::trimatl(cmat));
  return icl_emiss + cst;
}

double DcSbmUndirected::icl_emiss(const List & obs_stats,int oldcl,int newcl, bool dead){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::vec din =as<arma::vec>(obs_stats["din"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::umat si = submatcross(oldcl,newcl,counts.n_rows);
  double icl_emiss = (lgamma(counts(newcl))-lgamma(counts(newcl)+din(newcl))+din(newcl)*log(counts(newcl)));
  if(!dead){
    icl_emiss += (lgamma(counts(oldcl))-lgamma(counts(oldcl)+din(oldcl))+din(oldcl)*log(counts(oldcl)));
  }
  int k = 0;
  int l = 0;
  int cc = 0;
  for (arma::uword i = 0;i<si.n_rows;++i){
    k=si(i,0);
    l=si(i,1);
    if(counts(k)*counts(l)!=0){
      if(k<l){
        cc = counts(k)*counts(l);
        // lets go
        icl_emiss += lgamma(edges_counts(k,l)+1)-(edges_counts(k,l)+1)*log(p*cc+1);
      }
      if((k==l) && (counts(k)>1)){
        cc = (counts(k)*(counts(k)-1))/2;
        // lets go
        icl_emiss += lgamma(edges_counts(k,l)/2+1)-(edges_counts(k,l)/2+1)*log(p*cc+1);
      }
      
    }
    
  }
  return icl_emiss+cst;
}



double DcSbmUndirected::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
  // here old refers to the stats before the fusion between obk and obl
  //Rcout << "Je calculs des corrections !!" << std::endl;
  //Rcout << obk << "---- " << obl << std::endl;
  int a,b,ao,bo,lo;
  double icl_cor = 0;
  int cc, cc_old;
  double oxc,xc;
  arma::vec old_counts =as<arma::vec>(old_stats["counts"]);
  arma::mat old_x_counts =as<arma::mat>(old_stats["x_counts"]);
  cc = counts(k)*counts(l);
  arma::ivec kl({k, l});
  arma::ivec mkl({obk, obl});
  if(l>=obk){
    lo=l+1;
  }else{
    lo=l;
  }
  for(int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      
      
      a = kl(i);
      b = mkl(j);
      if(b>=obk){
        b=b-1;
      }
      
      
      
      // new stats no fusion k/l
      if(j==1){
        cc = counts(a)*counts(b);
        icl_cor -= lgamma(x_counts(a,b)+1)-(x_counts(a,b)+1)*log(p*cc+1);
      }
      
      // new stats fusion k/l
      if((j==1) & (i==0)){
        cc = (counts(k)+counts(l))*counts(b);
        xc    = x_counts(k,b)+x_counts(l,b);
        icl_cor += lgamma(xc+1)-(xc+1)*log(p*cc+1);
      }
      
      // handling matrix sizes differences
      ao = kl(i);
      bo = mkl(j);
      if(ao>=obk){
        ao=ao+1;
      }
      // old stats no fusion k/l
      cc_old = old_counts(ao)*old_counts(bo);
      oxc = old_x_counts(ao,bo);
      icl_cor += lgamma(oxc+1)-(oxc+1)*log(p*cc_old+1);
      // old stats fusion k/l
      if(i==0){
        cc_old = (old_counts(ao)+old_counts(lo))*old_counts(bo);
        oxc    = old_x_counts(ao,bo)+old_x_counts(lo,bo);
        icl_cor -= lgamma(oxc+1)-(oxc+1)*log(p*cc_old+1);
      }
      
      
    }
  }

  return icl_cor;
  
  
}

