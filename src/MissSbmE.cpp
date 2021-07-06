// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "MissSbmE.h"
using namespace Rcpp;



MissSbmE::MissSbmE(arma::sp_mat& xp,arma::sp_mat& xpobs,S4 modeli,bool verb){
  model = modeli;
  a0 = model.slot("a0");
  b0 = model.slot("b0");
  x  = xp;
  xt = xp.t();
  xobs  = xpobs;
  xtobs = xpobs.t();
  verbose=verb;

}

void MissSbmE::set_cl(arma::vec cl){
  K = arma::max(cl)+1;
  counts = count(cl,K);
  x_counts = gsum_mat(cl,x,K);
  x_counts_obs = gsum_mat(cl,xobs,K);
}


double MissSbmE::icl_emiss(const List & obs_stats){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::mat obs_counts =as<arma::mat>(obs_stats["x_counts_obs"]);
  arma::mat cmat = lgamma(a0+edges_counts)+lgamma(obs_counts-edges_counts+b0)+lgamma(a0+b0);
  cmat = cmat - lgamma(a0) - lgamma(b0) - lgamma(obs_counts+a0+b0);
  double icl_emiss=accu(cmat);
  return icl_emiss;
}

double MissSbmE::icl_emiss(const List & obs_stats,int oldcl,int newcl){
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::mat obs_counts =as<arma::mat>(obs_stats["x_counts_obs"]);
  arma::mat si = submatcross(oldcl,newcl,counts.n_rows);
  double icl_emiss = 0;
  int k = 0;
  int l = 0;
  int cc = 0;
  for (arma::uword i = 0;i<si.n_rows;++i){
    k=si(i,0);
    l=si(i,1);
    if(counts(k)*counts(l)!=0){
      cc = obs_counts(k,l);
      icl_emiss += lgamma(a0+edges_counts(k,l))+lgamma(b0+cc-edges_counts(k,l))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);  
    }
    
  }
  return icl_emiss;
}


List MissSbmE::get_obs_stats(){
  return List::create(Named("counts", counts),Named("x_counts", x_counts),Named("x_counts_obs", x_counts_obs));
}
arma::mat MissSbmE::delta_swap(int i,int K,Partition clp,arma::uvec iclust){
  
  int self=x(i,i);
  int selfobs=xobs(i,i);
  int oldcl = clp.get(i);
  

  
  arma::sp_mat col_sum = gsum_partition(clp,x,i,K);
  col_sum(oldcl)=col_sum(oldcl)-self;
  
  arma::sp_mat col_sum_obs = gsum_partition(clp,xobs,i,K);
  col_sum_obs(oldcl)=col_sum_obs(oldcl)-selfobs;
  
  arma::sp_mat row_sum = gsum_partition(clp,xt,i,K);
  row_sum(oldcl)=row_sum(oldcl)-self;
  
  arma::sp_mat row_sum_obs = gsum_partition(clp,xtobs,i,K);
  row_sum_obs(oldcl)=row_sum_obs(oldcl)-selfobs;
  
  
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts), Named("x_counts_obs", x_counts_obs));
  int k = 0;
  // for each possible move
  for(arma::uword j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    if(k!=oldcl){
      
      arma::mat new_ec = x_counts;
      new_ec.col(k) = new_ec.col(k)+col_sum;
      new_ec.row(k) = new_ec.row(k)+row_sum.t();
      new_ec.col(oldcl) = new_ec.col(oldcl)-col_sum;
      new_ec.row(oldcl) = new_ec.row(oldcl)-row_sum.t();
      new_ec(k,k)=new_ec(k,k)+self;
      new_ec(oldcl,oldcl)=new_ec(oldcl,oldcl)-self;
      
      
      arma::mat new_ec_obs = x_counts_obs;
      new_ec_obs.col(k) = new_ec_obs.col(k)+col_sum_obs;
      new_ec_obs.row(k) = new_ec_obs.row(k)+row_sum_obs.t();
      new_ec_obs.col(oldcl) = new_ec_obs.col(oldcl)-col_sum_obs;
      new_ec_obs.row(oldcl) = new_ec_obs.row(oldcl)-row_sum_obs.t();
      new_ec_obs(k,k)=new_ec_obs(k,k)+selfobs;
      new_ec_obs(oldcl,oldcl)=new_ec_obs(oldcl,oldcl)-selfobs;
      
      arma::vec new_counts = update_count(counts,oldcl,k);
      List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec), Named("x_counts_obs", new_ec_obs));
      delta(k)=icl_emiss(new_stats,oldcl,k)-icl_emiss(old_stats,oldcl,k);
    }
  }
  return delta;
}



void MissSbmE::swap_update(const int i,Partition clp,bool dead_cluster,  const int newcl){
  int self=x(i,i);
  int oldcl = clp.get(i);
  int selfobs=xobs(i,i);
  arma::sp_mat col_sum = gsum_partition(clp,x,i,K);
  col_sum(oldcl)=col_sum(oldcl)-self;
  
  arma::sp_mat col_sum_obs = gsum_partition(clp,xobs,i,K);
  col_sum_obs(oldcl)=col_sum_obs(oldcl)-selfobs;
  
  arma::sp_mat row_sum = gsum_partition(clp,xt,i,K);
  row_sum(oldcl)=row_sum(oldcl)-self;
  
  arma::sp_mat row_sum_obs = gsum_partition(clp,xtobs,i,K);
  row_sum_obs(oldcl)=row_sum_obs(oldcl)-selfobs;
  
  arma::mat new_ec = x_counts;
  new_ec.col(newcl) = new_ec.col(newcl)+col_sum;
  new_ec.row(newcl) = new_ec.row(newcl)+row_sum.t();
  new_ec.col(oldcl) = new_ec.col(oldcl)-col_sum;
  new_ec.row(oldcl) = new_ec.row(oldcl)-row_sum.t();
  new_ec(newcl,newcl)=new_ec(newcl,newcl)+self;
  new_ec(oldcl,oldcl)=new_ec(oldcl,oldcl)-self;
  
  
  arma::mat new_ec_obs = x_counts_obs;
  new_ec_obs.col(newcl) = new_ec_obs.col(newcl)+col_sum_obs;
  new_ec_obs.row(newcl) = new_ec_obs.row(newcl)+row_sum_obs.t();
  new_ec_obs.col(oldcl) = new_ec_obs.col(oldcl)-col_sum_obs;
  new_ec_obs.row(oldcl) = new_ec_obs.row(oldcl)-row_sum_obs.t();
  new_ec_obs(newcl,newcl)=new_ec_obs(newcl,newcl)+selfobs;
  new_ec_obs(oldcl,oldcl)=new_ec_obs(oldcl,oldcl)-selfobs;
  arma::mat new_counts = update_count(counts,oldcl,newcl);
  if(new_counts(oldcl)==0){
    counts = new_counts.elem(arma::find(arma::linspace(0,K-1,K)!=oldcl));
    x_counts = new_ec(arma::find(arma::linspace(0,K-1,K)!=oldcl),arma::find(arma::linspace(0,K-1,K)!=oldcl));
    x_counts_obs = new_ec_obs(arma::find(arma::linspace(0,K-1,K)!=oldcl),arma::find(arma::linspace(0,K-1,K)!=oldcl));
     --K;
  }else{
    counts=new_counts;
    x_counts=new_ec;
    x_counts_obs=new_ec_obs;
  }

}


double MissSbmE::delta_merge(int k, int l){
  
  arma::mat new_ec = x_counts;
  arma::mat new_ec_obs = x_counts_obs;
  arma::mat new_counts = counts;
  new_counts(l)=new_counts(k)+new_counts(l);
  new_counts(k)=0;
  new_ec.col(l) = new_ec.col(k)+new_ec.col(l);
  new_ec.row(l) = new_ec.row(k)+new_ec.row(l);
  new_ec_obs.col(l) = new_ec_obs.col(k)+new_ec_obs.col(l);
  new_ec_obs.row(l) = new_ec_obs.row(k)+new_ec_obs.row(l);
  List old_stats = List::create(Named("counts", counts), Named("x_counts", x_counts), Named("x_counts_obs", x_counts_obs));
  List new_stats = List::create(Named("counts", new_counts), Named("x_counts", new_ec), Named("x_counts_obs", new_ec_obs));
  double delta = icl_emiss(new_stats,k,l)-icl_emiss(old_stats,k,l);
  
  return delta;
}


void MissSbmE::merge_update(int k,int l){

  counts(l) = counts(k)+counts(l);
  counts    = counts.elem(arma::find(arma::linspace(0,K-1,K)!=k));

  x_counts.col(l) = x_counts.col(k)+x_counts.col(l);
  x_counts.row(l) = x_counts.row(k)+x_counts.row(l);
  x_counts = x_counts(arma::find(arma::linspace(0,K-1,K)!=k),arma::find(arma::linspace(0,K-1,K)!=k));
  x_counts_obs.col(l) = x_counts_obs.col(k)+x_counts_obs.col(l);
  x_counts_obs.row(l) = x_counts_obs.row(k)+x_counts_obs.row(l);
  x_counts_obs = x_counts_obs(arma::find(arma::linspace(0,K-1,K)!=k),arma::find(arma::linspace(0,K-1,K)!=k));
  --K;
}




double MissSbmE::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
  // here old refers to the stats before the fusion between obk and obl
  //Rcout << "Je calculs des corrections !!" << std::endl;
  //Rcout << obk << "---- " << obl << std::endl;
  int a,b,ao,bo,lo;
  double icl_cor = 0;
  int cc, cc_old;
  double oxc,xc;
  arma::vec old_counts =as<arma::vec>(old_stats["counts"]);
  arma::mat old_x_counts =as<arma::mat>(old_stats["x_counts"]);
  arma::mat old_obs_counts =as<arma::mat>(old_stats["x_counts_obs"]);
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
        cc = x_counts_obs(a,b);
        icl_cor -= lgamma(a0+x_counts(a,b))+lgamma(b0+cc-x_counts(a,b))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);
        cc = x_counts_obs(b,a);
        icl_cor -= lgamma(a0+x_counts(b,a))+lgamma(b0+cc-x_counts(b,a))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);
      }
      
      // new stats fusion k/l
      if((j==1) & (i==0)){
        
        cc = x_counts_obs(k,b)+x_counts_obs(l,b);
        xc    = x_counts(k,b)+x_counts(l,b);
        icl_cor += lgamma(a0+xc)+lgamma(b0+cc-xc)+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc);
        cc = x_counts_obs(b,k)+x_counts_obs(b,l);
        xc    = x_counts(b,k)+x_counts(b,l);
        icl_cor += lgamma(a0+xc)+lgamma(b0+cc-xc)+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc); 
      }
      
      // handling matrix sizes differences
      ao = kl(i);
      bo = mkl(j);
      if(ao>=obk){
        ao=ao+1;
      }
      
      // old stats no fusion k/l
      cc_old = old_obs_counts(ao,bo);
      icl_cor += lgamma(a0+old_x_counts(ao,bo))+lgamma(b0+cc_old-old_x_counts(ao,bo))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc_old);
      cc_old = old_obs_counts(bo,ao);
      icl_cor += lgamma(a0+old_x_counts(bo,ao))+lgamma(b0+cc_old-old_x_counts(bo,ao))+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc_old);
      // old stats fusion k/l
      if(i==0){
        
        cc_old = old_obs_counts(ao,bo)+old_obs_counts(lo,bo);
        oxc    = old_x_counts(ao,bo)+old_x_counts(lo,bo);
        icl_cor -= lgamma(a0+oxc)+lgamma(b0+cc_old-oxc)+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc_old);
        cc_old = old_obs_counts(bo,ao)+old_obs_counts(bo,lo);
        oxc    = old_x_counts(bo,ao)+old_x_counts(bo,lo);
        icl_cor -= lgamma(a0+oxc)+lgamma(b0+cc_old-oxc)+lgamma(a0+b0)-lgamma(a0)-lgamma(b0)-lgamma(a0+b0+cc_old);
      }
      
      
    }
  }
  
  return icl_cor;
  
  
}


