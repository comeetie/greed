// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "DcSbm.h"
using namespace Rcpp;


DcSbm::DcSbm(arma::sp_mat& xp,int Ki,double alphai,bool verb){
  alpha = alphai;
  x  = xp;
  xt = xp.t();
  N  = x.n_rows;
  K  = Ki;
  cl = as<arma::vec>(sample(K,N,true));
  cl = cl -1;
  x_counts = gsum_mat(cl,x,K);
  counts = count(cl,K);
  din = sum_cols(x_counts);
  dout = sum_rows(x_counts);
  p= arma::accu(x_counts)/(N*N);
  verbose=verb;
  
  cst = - sum_lfact(xp);
}

DcSbm::DcSbm(arma::sp_mat& xp,int Ki,double alphai,arma::vec& clt,bool verb){
  alpha = alphai;
  x  = xp;
  xt = xp.t();
  N  = x.n_rows;
  K  = Ki;
  cl = clt;
  x_counts = gsum_mat(cl,x,K);
  counts = count(cl,K);
  din = sum_cols(x_counts);
  dout = sum_rows(x_counts);
  p= arma::accu(x_counts)/(N*N);
  verbose=verb;
  cst = - sum_lfact(xp);
}


double DcSbm::icl_emiss(const List & obs_stats){

  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::vec din =as<arma::vec>(obs_stats["din"]);
  arma::vec dout =as<arma::vec>(obs_stats["dout"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::mat matcount = counts*counts.t();
  // lets go


  double icl_emiss = accu(lgamma(counts)-lgamma(counts+din)+din % log(counts))+accu(lgamma(counts)-lgamma(counts+dout)+dout % log(counts));

  icl_emiss=icl_emiss + arma::accu(lgamma(edges_counts+1)-(edges_counts+1) % log(p*matcount+1));

  return icl_emiss+cst;
}

double DcSbm::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::vec din =as<arma::vec>(obs_stats["din"]);
  arma::vec dout =as<arma::vec>(obs_stats["dout"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::mat si = submatcross(oldcl,newcl,counts.n_rows);
  double icl_emiss = lgamma(counts(newcl))-lgamma(counts(newcl)+din(newcl))+din(newcl)*log(counts(newcl));
  icl_emiss += lgamma(counts(newcl))-lgamma(counts(newcl)+dout(newcl))+dout(newcl)*log(counts(newcl));
  if(counts(oldcl)!=0){
    icl_emiss += lgamma(counts(oldcl))-lgamma(counts(oldcl)+dout(oldcl))+dout(oldcl)*log(counts(oldcl));
    icl_emiss += lgamma(counts(oldcl))-lgamma(counts(oldcl)+din(oldcl))+din(oldcl)*log(counts(oldcl));
  }
  int k = 0;
  int l = 0;
  int cc = 0;
  for (int i = 0;i<si.n_rows;++i){
    k=si(i,0);
    l=si(i,1);
    if(counts(k)*counts(l)!=0){
        cc = counts(k)*counts(l);
      // lets go
        icl_emiss += lgamma(edges_counts(k,l)+1)-(edges_counts(k,l)+1)*log(p*cc+1);
    }
    
  }
  return icl_emiss+cst;
}


List DcSbm::get_obs_stats(){
  return List::create(Named("counts", counts), Named("din", din),Named("dout", dout),Named("x_counts", x_counts));
}

arma::mat DcSbm::delta_swap(int i){
  int self=x(i,i);
  int oldcl = cl(i);
  arma::mat col_sum = gsum_col(cl,x,i,K);
  col_sum(oldcl)=col_sum(oldcl)-self;
  arma::mat row_sum = gsum_col(cl,xt,i,K);
  row_sum(oldcl)=row_sum(oldcl)-self;
  arma::vec delta(K);
  delta.fill(0);
  int cdin = arma::accu(col_sum)+self;
  int cdout = arma::accu(row_sum)+self;
  List old_stats = List::create(Named("counts", counts),Named("din", din),Named("dout", dout), Named("x_counts", x_counts));
  for(int k = 0; k < K; ++k) {
    if(k!=oldcl){
      arma::mat new_ec = x_counts;
      new_ec.col(k) = new_ec.col(k)+col_sum;
      new_ec.row(k) = new_ec.row(k)+row_sum.t();
      new_ec.col(oldcl) = new_ec.col(oldcl)-col_sum;
      new_ec.row(oldcl) = new_ec.row(oldcl)-row_sum.t();
      new_ec(k,k)=new_ec(k,k)+self;
      new_ec(oldcl,oldcl)=new_ec(oldcl,oldcl)-self;
      arma::vec new_counts = update_count(counts,oldcl,k);
      arma::vec new_din = din;
      new_din(oldcl) = new_din(oldcl)-cdin;
      new_din(k) = new_din(k)+cdin;
      arma::vec new_dout = dout;
      new_dout(oldcl) = new_dout(oldcl)-cdout;
      new_dout(k) = new_dout(k)+cdout;
      List new_stats = List::create(Named("counts", new_counts),Named("din", new_din),Named("dout", new_dout), Named("x_counts", new_ec));
      delta(k)=icl(new_stats,oldcl,k)-icl(old_stats,oldcl,k);
    }
  }
  return delta;
}


void DcSbm::swap_update(const int i,const int newcl){
  int self=x(i,i);
  int oldcl = cl(i);
  arma::mat col_sum = gsum_col(cl,x,i,K);
  col_sum(oldcl)=col_sum(oldcl)-self;
  arma::mat row_sum = gsum_col(cl,xt,i,K);
  row_sum(oldcl)=row_sum(oldcl)-self;
  arma::mat new_ec = x_counts;
  int cdin = arma::accu(col_sum)+self;
  int cdout = arma::accu(row_sum)+self;
  new_ec.col(newcl) = new_ec.col(newcl)+col_sum;
  new_ec.row(newcl) = new_ec.row(newcl)+row_sum.t();
  new_ec.col(oldcl) = new_ec.col(oldcl)-col_sum;
  new_ec.row(oldcl) = new_ec.row(oldcl)-row_sum.t();
  new_ec(newcl,newcl)=new_ec(newcl,newcl)+self;
  new_ec(oldcl,oldcl)=new_ec(oldcl,oldcl)-self;
  arma::mat new_counts = update_count(counts,oldcl,newcl);
  arma::vec new_din = din;
  new_din(oldcl) = new_din(oldcl)-cdin;
  new_din(newcl) = new_din(newcl)+cdin;
  arma::vec new_dout = dout;
  new_dout(oldcl) = new_dout(oldcl)-cdout;
  new_dout(newcl) = new_dout(newcl)+cdout;
  cl(i)=newcl;
  if(new_counts(oldcl)==0){
    counts = new_counts.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    din = new_din.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    dout = new_dout.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    x_counts = new_ec(arma::find(arma::linspace(0,K-1,K)!=oldcl),arma::find(arma::linspace(0,K-1,K)!=oldcl));
    cl.elem(arma::find(cl>oldcl))=cl.elem(arma::find(cl>oldcl))-1;
    --K;
  }else{
    counts=new_counts;
    x_counts=new_ec;
    din=new_din;
    dout=new_dout;
  }


}


MergeMat DcSbm::delta_merge(){
  arma::mat delta(K,K);
  delta.fill(0);
  int bk = 0;
  int bl = 0;
  double bv = -std::numeric_limits<double>::infinity();
  List old_stats = List::create(Named("counts", counts), Named("din", din),Named("dout", dout),Named("x_counts", x_counts));
  for(int k = 1 ; k < K ; ++k) {
    for (int l = 0 ; l<k ; ++l){
      arma::mat new_ec = x_counts;
      arma::mat new_counts = counts;
      new_counts(l)=new_counts(k)+new_counts(l);
      new_counts(k)=0;
      new_ec.col(l) = new_ec.col(k)+new_ec.col(l);
      new_ec.row(l) = new_ec.row(k)+new_ec.row(l);
      arma::vec new_din = din;
      new_din(l) = new_din(k)+new_din(l);
      arma::vec new_dout = dout;
      new_dout(l) = new_dout(k)+new_dout(l);
      List new_stats = List::create(Named("counts", new_counts),Named("din", new_din),Named("dout", new_dout), Named("x_counts", new_ec));
      delta(k,l)=icl(new_stats,k,l)-icl(old_stats,k,l);
      delta(l,k)=delta(k,l);
      if(delta(k,l)>bv){
        bk=k;
        bl=l;
        bv=delta(k,l);
      }
    }
  }
  return MergeMat(bk,bl,bv,delta);
}

MergeMat DcSbm::delta_merge(arma::mat delta, int obk, int obl){
  delta = delta(arma::find(arma::linspace(0,K,K+1)!=obk),arma::find(arma::linspace(0,K,K+1)!=obk));
  int bk = 0;
  int bl = 0;
  double bv = -std::numeric_limits<double>::infinity();
  List old_stats = List::create(Named("counts", counts), Named("din", din),Named("dout", dout),Named("x_counts", x_counts));
  for(int k = 1; k < K; ++k) {
    for (int l = 0;l<k;++l){
      if(k == obl | l == obl){
        arma::mat new_ec = x_counts;
        arma::mat new_counts = counts;
        new_counts(l)=new_counts(k)+new_counts(l);
        new_counts(k)=0;
        new_ec.col(l) = new_ec.col(k)+new_ec.col(l);
        new_ec.row(l) = new_ec.row(k)+new_ec.row(l);
        arma::vec new_din = din;
        new_din(l) = new_din(k)+new_din(l);
        arma::vec new_dout = dout;
        new_dout(l) = new_dout(k)+new_dout(l);
        List new_stats = List::create(Named("counts", new_counts),Named("din", new_din),Named("dout", new_dout), Named("x_counts", new_ec));
        delta(k,l)=icl(new_stats,k,l)-icl(old_stats,k,l);
        delta(l,k)=delta(k,l);
      }
      if(delta(k,l)>bv){
        bk=k;
        bl=l;
        bv=delta(k,l);
      }
    }
  }
  return MergeMat(bk,bl,bv,delta);
}

void DcSbm::merge_update(int k,int l){
  cl(arma::find(cl==k))=arma::ones(counts(k),1)*l;
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  counts(l) = counts(k)+counts(l);
  counts    = counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));

  x_counts.col(l) = x_counts.col(k)+x_counts.col(l);
  x_counts.row(l) = x_counts.row(k)+x_counts.row(l);
  x_counts = x_counts(arma::find(arma::linspace(0,K-1,K)!=k),arma::find(arma::linspace(0,K-1,K)!=k));
  
  din(l) = din(l)+din(k);
  din = din.elem(arma::find(arma::linspace(0,K-1,K)!=k));
  dout(l) = dout(l)+dout(k);
  dout = dout.elem(arma::find(arma::linspace(0,K-1,K)!=k));
  
  
  --K;
}
