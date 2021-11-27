// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MultSbm.h"
using namespace Rcpp;

MultSbm::MultSbm(const arma::cube  & xp,S4 modeli,bool verb){
  model = modeli;
  beta = model.slot("beta");
  x  = xp;
  N  = x.n_rows;
  M = x.n_slices;
  verbose=verb;
  // cst to add if needed 
  cst = 0;
}



void MultSbm::set_cl(arma::uvec clt){
  K = arma::max(clt)+1;
  x_counts = gsum_cube(clt,x,K);
}

double MultSbm::icl_emiss(const List & obs_stats){
  arma::cube edges_counts =as<arma::cube>(obs_stats["x_counts"]);
  // lets go
  double icl_emiss = 0;
  for (int k=0;k<K;++k){
    for (int l=0;l<K;++l){
      arma::vec klcounts = edges_counts.tube(k,l);
      icl_emiss+=lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
    }
  }
  return icl_emiss+cst;
}

double MultSbm::icl_emiss(const List & obs_stats,int oldcl,int newcl, bool dead_cluster){
  arma::cube edges_counts =as<arma::cube>(obs_stats["x_counts"]);
  arma::umat si = submatcross(oldcl,newcl,K);
  double icl_emiss = 0;
  int k = 0;
  int l = 0;
  for (arma::uword i = 0;i<si.n_rows;++i){
    k=si(i,0);
    l=si(i,1);
    if(!dead_cluster){
      arma::vec klcounts = edges_counts.tube(k,l);
      icl_emiss+=lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
    }else{
      if((k!=oldcl) & (l!=oldcl)){
        arma::vec klcounts = edges_counts.tube(k,l);
        icl_emiss+=lgamma(M*beta)+arma::accu(lgamma(beta+klcounts))-M*lgamma(beta)-lgamma(arma::accu(klcounts+beta));
      }
    }
    
  }
  return icl_emiss;
  
}



List MultSbm::get_obs_stats(){
  return List::create(Named("x_counts", x_counts));
}



arma::vec MultSbm::delta_swap(int i,arma::uvec & cl, bool almost_dead_cluster, arma::uvec iclust, int K){
  arma::vec self=x.tube(i,i);
  int oldcl = cl(i);
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  List old_stats = List::create( Named("x_counts", x_counts));
  arma::mat cur_row(M,K);
  cur_row.fill(0);
  arma::mat cur_col(M,K);
  cur_col.fill(0);
  for(arma::uword j = 0; j < x.n_rows; ++j) {
    if(j!=i){
      arma::colvec curr = x.tube(i,j);
      cur_row.col(cl(j))=cur_row.col(cl(j))+curr;
      arma::colvec curc = x.tube(j,i);
      cur_col.col(cl(j))=cur_col.col(cl(j))+curc;
    }
  }
  
  
  int k = 0;
  // for each possible move
  for(arma::uword j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    if(k!=oldcl){
      arma::cube new_ec = x_counts;

      // todo
      
      for (arma::uword  m=0;m<M;++m){
        new_ec(oldcl,oldcl,m)=new_ec(oldcl,oldcl,m)-self(m);
        new_ec(k,k,m)=new_ec(k,k,m)+self(m);
      }
      
      for (arma::uword  l=0; l<K;++l){
        for (arma::uword  m=0;m<M;++m){
          new_ec(oldcl,l,m)=new_ec(oldcl,l,m)-cur_row(m,l);
          new_ec(k,l,m)=new_ec(k,l,m)+cur_row(m,l);
          new_ec(l,oldcl,m)=new_ec(l,oldcl,m)-cur_col(m,l);
          new_ec(l,k,m)=new_ec(l,k,m)+cur_col(m,l);
        } 
      }
      
      List new_stats = List::create(Named("x_counts", new_ec));
      delta(k)=icl_emiss(new_stats,oldcl,k,almost_dead_cluster)-icl_emiss(old_stats,oldcl,k,false);
    }
    
  }
  return delta;  
}


void MultSbm::swap_update(const int i,arma::uvec &  cl,bool dead_cluster,const int newcl){
  arma::vec self=x.tube(i,i);
  int oldcl = cl(i);
  arma::mat cur_row(M,K);
  cur_row.fill(0);
  arma::mat cur_col(M,K);
  cur_col.fill(0);
  for(arma::uword  j = 0; j < x.n_rows; ++j) {
    if(j!=i){
      arma::colvec curr = x.tube(i,j);
      cur_row.col(cl(j))=cur_row.col(cl(j))+curr;
      arma::colvec curc = x.tube(j,i);
      cur_col.col(cl(j))=cur_col.col(cl(j))+curc;
    }
  }

  
  for (arma::uword m=0;m<M;m++){
    x_counts(oldcl,oldcl,m)=x_counts(oldcl,oldcl,m)-self(m);
    x_counts(newcl,newcl,m)=x_counts(newcl,newcl,m)+self(m);
  }
  for (arma::uword  l=0; l<K;++l){
    for (arma::uword  m=0;m<M;m++){
      x_counts(oldcl,l,m)=x_counts(oldcl,l,m)-cur_row(m,l);
      x_counts(newcl,l,m)=x_counts(newcl,l,m)+cur_row(m,l);
      x_counts(l,oldcl,m)=x_counts(l,oldcl,m)-cur_col(m,l);
      x_counts(l,newcl,m)=x_counts(l,newcl,m)+cur_col(m,l);
    } 
  }
  

  if(dead_cluster){
    x_counts.shed_row(oldcl);
    x_counts.shed_col(oldcl);
    --K;
  }
  
}


double MultSbm::delta_merge(int k, int l){
  
  List old_stats = List::create(Named("x_counts", x_counts));
  
  arma::cube new_ec = x_counts;
  
  new_ec.tube(l,l)=new_ec.tube(l,l)+new_ec.tube(k,k);
  new_ec.tube(k,k)=new_ec.tube(k,k)-new_ec.tube(k,k);
  for (arma::uword  h=0; h<K;++h){
    for (arma::uword  m=0;m<M;m++){
      new_ec(l,h,m)=new_ec(l,h,m)+new_ec(k,h,m);
      new_ec(k,h,m)=new_ec(k,h,m)-new_ec(k,h,m);
      new_ec(h,l,m)=new_ec(h,l,m)+new_ec(h,k,m);
      new_ec(h,k,m)=new_ec(h,k,m)-new_ec(h,k,m);
    } 
  }
  
  List new_stats = List::create( Named("x_counts", new_ec));
  double delta=icl_emiss(new_stats,k,l,true)-icl_emiss(old_stats,k,l,false);
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
  arma::cube old_x_counts =as<arma::cube>(old_stats["x_counts"]);
  arma::vec klcounts;

  arma::uvec kl;
  kl << k << l << arma::endr;
  arma::uvec mkl;
  mkl << obk << obl << arma::endr;
  if(l>=obk){
    lo=l+1;
  }else{
    lo=l;
  }
  for(arma::uword  i=0;i<2;i++){
    for (arma::uword  j=0;j<2;j++){
      
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



void MultSbm::merge_update(int k,int l){
  x_counts.tube(l,l)=x_counts.tube(l,l)+x_counts.tube(k,k);
  x_counts.tube(k,k)=x_counts.tube(k,k)-x_counts.tube(k,k);
  for (arma::uword  h=0; h<K;++h){
    for (arma::uword  m=0;m<M;m++){
      x_counts(l,h,m)=x_counts(l,h,m)+x_counts(k,h,m);
      x_counts(k,h,m)=x_counts(k,h,m)-x_counts(k,h,m);
      x_counts(h,l,m)=x_counts(h,l,m)+x_counts(h,k,m);
      x_counts(h,k,m)=x_counts(h,k,m)-x_counts(h,k,m);
    } 
  }
  

  x_counts.shed_row(k);
  x_counts.shed_col(k);

  --K;
}



