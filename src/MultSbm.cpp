// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "MultSbm.h"
using namespace Rcpp;



MultSbm::MultSbm(const arma::cube& xp,double alphai,double betai,arma::vec& clt,bool verb){
  alpha = alphai;
  beta = betai;
  x  = xp;
  N  = x.n_rows;
  M = x.n_slices;
  set_cl(clt);
  verbose=verb;
  // cst to add if needed 
  cst = 0;
}

void MultSbm::set_cl(arma::vec clt){
  cl = clt;
  K = arma::max(cl)+1;
  x_counts = gsum_cube(cl,x,K);
  counts = count(cl,K);
}

double MultSbm::icl_emiss(const List & obs_stats){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::cube edges_counts =as<arma::cube>(obs_stats["x_counts"]);
  arma::mat matcount = counts*counts.t();
  // lets go
  double icl_emiss = 0;
  for (int k=0;k<K;++k){
    for (int l=0;l<K;++l){
      arma::vec klcounts = edges_counts.tube(k,l);
      icl_emiss+=lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
    }
  }
  return icl_emiss/2+cst;
}

double MultSbm::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::cube edges_counts =as<arma::cube>(obs_stats["x_counts"]);
  arma::mat si = submatcross(oldcl,newcl,counts.n_rows);
  double icl_emiss = 0;
  int k = 0;
  int l = 0;
  int cc = 0;
  for (arma::uword i = 0;i<si.n_rows;++i){
    k=si(i,0);
    l=si(i,1);
    if(counts(k)*counts(l)!=0){
      arma::vec klcounts = edges_counts.tube(k,l);
      icl_emiss+=lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
    }
    
  }
  return icl_emiss/2;
}





List MultSbm::get_obs_stats(){
  return List::create(Named("counts", counts), Named("x_counts", x_counts));
}

arma::mat MultSbm::delta_swap(int i,arma::uvec iclust){
  arma::vec self=x.tube(i,i);
  int oldcl = cl(i);
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  arma::mat cur_row(M,K);
  cur_row.fill(0);
  arma::mat cur_col(M,K);
  cur_col.fill(0);
  for(int j = 0; j < x.n_rows; ++j) {
    if(j!=i){
      arma::colvec curr = x.tube(i,j);
      cur_row.col(cl(j))=cur_row.col(cl(j))+curr;
      arma::colvec curc = x.tube(j,i);
      cur_col.col(cl(j))=cur_col.col(cl(j))+curc;
    }
  }
  
  
  int k = 0;
  // for each possible move
  for(int j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    if(k!=oldcl){
      arma::cube new_ec = x_counts;
      arma::vec new_counts = update_count(counts,oldcl,k);
      // todo
      
      for (int m=0;m<M;++m){
        new_ec(oldcl,oldcl,m)=new_ec(oldcl,oldcl,m)-self(m);
        new_ec(k,k,m)=new_ec(k,k,m)+self(m);
      }
      
      for (int l=0; l<K;++l){
        for (int m=0;m<M;++m){
          new_ec(oldcl,l,m)=new_ec(oldcl,l,m)-cur_row(m,l);
          new_ec(k,l,m)=new_ec(k,l,m)+cur_row(m,l);
          new_ec(l,oldcl,m)=new_ec(l,oldcl,m)-cur_col(m,l);
          new_ec(l,k,m)=new_ec(l,k,m)+cur_col(m,l);
        } 
      }
      
      List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
      delta(k)=icl(new_stats,oldcl,k)-icl(old_stats,oldcl,k);
    }
    
  }
  return delta;
}


void MultSbm::swap_update(const int i,const int newcl){
  arma::vec self=x.tube(i,i);
  int oldcl = cl(i);
  arma::mat cur_row(M,K);
  cur_row.fill(0);
  arma::mat cur_col(M,K);
  cur_col.fill(0);
  for(int j = 0; j < x.n_rows; ++j) {
    if(j!=i){
      arma::colvec curr = x.tube(i,j);
      cur_row.col(cl(j))=cur_row.col(cl(j))+curr;
      arma::colvec curc = x.tube(j,i);
      cur_col.col(cl(j))=cur_col.col(cl(j))+curc;
    }
  }
  counts = update_count(counts,oldcl,newcl);
  
  for (int m=0;m<M;m++){
    x_counts(oldcl,oldcl,m)=x_counts(oldcl,oldcl,m)-self(m);
    x_counts(newcl,newcl,m)=x_counts(newcl,newcl,m)+self(m);
  }
  for (int l=0; l<K;++l){
    for (int m=0;m<M;m++){
      x_counts(oldcl,l,m)=x_counts(oldcl,l,m)-cur_row(m,l);
      x_counts(newcl,l,m)=x_counts(newcl,l,m)+cur_row(m,l);
      x_counts(l,oldcl,m)=x_counts(l,oldcl,m)-cur_col(m,l);
      x_counts(l,newcl,m)=x_counts(l,newcl,m)+cur_col(m,l);
    } 
  }

  cl(i)=newcl;
  if(counts(oldcl)==0){
    counts = counts.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    x_counts.shed_row(oldcl);
    x_counts.shed_col(oldcl);
    cl.elem(arma::find(cl>oldcl))=cl.elem(arma::find(cl>oldcl))-1;
    --K;
  }
}


double MultSbm::delta_merge(int k, int l){
  
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts));
  
  arma::cube new_ec = x_counts;
  arma::mat new_counts = counts;
  new_counts(l)=new_counts(k)+new_counts(l);
  new_counts(k)=0;
  
  new_ec.tube(l,l)=new_ec.tube(l,l)+new_ec.tube(k,k);
  new_ec.tube(k,k)=new_ec.tube(k,k)-new_ec.tube(k,k);
  for (int h=0; h<K;++h){
    for (int m=0;m<M;m++){
      new_ec(l,h,m)=new_ec(l,h,m)+new_ec(k,h,m);
      new_ec(k,h,m)=new_ec(k,h,m)-new_ec(k,h,m);
      new_ec(h,l,m)=new_ec(h,l,m)+new_ec(h,k,m);
      new_ec(h,k,m)=new_ec(h,k,m)-new_ec(h,k,m);
    } 
  }
  
  List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec));
  double delta=icl(new_stats,k,l)-icl(old_stats,k,l);
  
  return delta;
}


double MultSbm::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
  // here old refers to the stats before the fusion between obk and obl
  //Rcout << "Je calculs des corrections !!" << std::endl;
  //Rcout << obk << "---- " << obl << std::endl;
  int a,b,ao,bo,lo;
  double icl_cor = 0;
  int cc, cc_old;
  double oxc,xc;
  arma::vec old_counts =as<arma::vec>(old_stats["counts"]);
  arma::cube old_x_counts =as<arma::cube>(old_stats["x_counts"]);
  arma::vec klcounts;
  cc = counts(k)*counts(l);
  arma::uvec kl;
  kl << k << l << arma::endr;
  arma::uvec mkl;
  mkl << obk << obl << arma::endr;
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
        
        klcounts = x_counts.tube(a,b);
        icl_cor -= lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
        klcounts = x_counts.tube(b,a);
        icl_cor -= lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
        
      }
      
      // new stats fusion k/l
      if((j==1) & (i==0)){
        klcounts  = x_counts.tube(k,b)+x_counts.tube(l,b);
        icl_cor += lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
        klcounts  = x_counts.tube(b,k)+x_counts.tube(b,l);
        icl_cor += lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
      }
      
      // handling matrix sizes differences
      ao = kl(i);
      bo = mkl(j);
      if(ao>=obk){
        ao=ao+1;
      }
      
      // old stats no fusion k/l
      klcounts  = old_x_counts.tube(ao,bo);
      icl_cor += lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
      klcounts  = old_x_counts.tube(bo,ao);
      icl_cor += lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
      // old stats fusion k/l
      if(i==0){
        klcounts  = old_x_counts.tube(ao,bo)+old_x_counts.tube(lo,bo);
        icl_cor -= lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
        klcounts  = old_x_counts.tube(bo,ao)+old_x_counts.tube(bo,lo);
        icl_cor -= lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
        
      }
      
      
    }
  }
  
  return icl_cor;
}


void MultSbm::merge_update(int k, int l){
  arma::cube new_ec = x_counts;
  arma::mat new_counts = counts;
  cl(arma::find(cl==k))=arma::ones(counts(k),1)*l;
  new_counts(l)=new_counts(k)+new_counts(l);
  new_counts(k)=0;
  new_ec.tube(l,l)=new_ec.tube(l,l)+new_ec.tube(k,k);
  new_ec.tube(k,k)=new_ec.tube(k,k)-new_ec.tube(k,k);
  for (int h=0; h<K;++h){
    for (int m=0;m<M;m++){
      new_ec(l,h,m)=new_ec(l,h,m)+new_ec(k,h,m);
      new_ec(k,h,m)=new_ec(k,h,m)-new_ec(k,h,m);
      new_ec(h,l,m)=new_ec(h,l,m)+new_ec(h,k,m);
      new_ec(h,k,m)=new_ec(h,k,m)-new_ec(h,k,m);
    } 
  }
  
  counts = new_counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));
  x_counts = new_ec;
  x_counts.shed_row(k);
  x_counts.shed_col(k);
  cl.elem(arma::find(cl>k))=cl.elem(arma::find(cl>k))-1;
  --K;
}