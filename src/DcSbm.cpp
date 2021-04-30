// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "MergeMat.h"
#include "IclModel.h"
#include "DcSbm.h"
using namespace Rcpp;



DcSbm::DcSbm(const arma::sp_mat& xp,S4 modeli,arma::vec& clt,bool verb){
  model = modeli;
  alpha = modeli.slot("alpha");
  x  = xp;
  xt = xp.t();
  N  = x.n_rows;
  set_cl(clt);
  if(Rcpp::traits::is_nan<REALSXP>(model.slot("p"))){
    p= arma::accu(x_counts)/(N*N);
    model.slot("p") = p;
  }else{
    p = model.slot("p");
  }
  verbose=verb;
  // cst to add 
  cst = 0;
}

void DcSbm::set_cl(arma::vec clt){
  cl = clt;
  K = arma::max(cl)+1;
  x_counts = gsum_mat(cl,x,K);
  counts = count(cl,K);
  din = sum(x_counts).t();
  dout = sum(x_counts.t()).t();
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
  for (arma::uword i = 0;i<si.n_rows;++i){
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

List DcSbm::get_obs_stats_cst(){
  arma::sp_mat din_node = sum(x).t();
  arma::sp_mat dout_node = sum(x.t()).t();
  return List::create(Named("din_node", din_node),Named("dout_node", dout_node));
}


arma::mat DcSbm::delta_swap(int i,arma::uvec iclust){
  int self=x(i,i);
  int oldcl = cl(i);
  arma::sp_mat col_sum = gsum_col(cl,x,i,K);
  col_sum(oldcl)=col_sum(oldcl)-self;
  arma::sp_mat row_sum = gsum_col(cl,xt,i,K);
  row_sum(oldcl)=row_sum(oldcl)-self;
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  int cdin = arma::accu(col_sum)+self;
  int cdout = arma::accu(row_sum)+self;
  List old_stats = List::create(Named("counts", counts),Named("din", din),Named("dout", dout), Named("x_counts", x_counts));
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
  arma::sp_mat col_sum = gsum_col(cl,x,i,K);
  col_sum(oldcl)=col_sum(oldcl)-self;
  arma::sp_mat row_sum = gsum_col(cl,xt,i,K);
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


double DcSbm::delta_merge(int k, int l){
  
  List old_stats = List::create(Named("counts", counts), Named("din", din),Named("dout", dout),Named("x_counts", x_counts));
  
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
  double delta=icl(new_stats,k,l)-icl(old_stats,k,l);
  
  return delta;
}


double DcSbm::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
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
        cc = counts(a)*counts(b);
        icl_cor -= lgamma(x_counts(a,b)+1)-(x_counts(a,b)+1)*log(p*cc+1);
        icl_cor -= lgamma(x_counts(b,a)+1)-(x_counts(b,a)+1)*log(p*cc+1);
      }
      
      // new stats fusion k/l
      if((j==1) & (i==0)){
        cc = (counts(k)+counts(l))*counts(b);
        xc    = x_counts(k,b)+x_counts(l,b);
        icl_cor += lgamma(xc+1)-(xc+1)*log(p*cc+1);
        xc    = x_counts(b,k)+x_counts(b,l);
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
      oxc = old_x_counts(bo,ao);
      icl_cor += lgamma(oxc+1)-(oxc+1)*log(p*cc_old+1);
      // old stats fusion k/l
      if(i==0){
        cc_old = (old_counts(ao)+old_counts(lo))*old_counts(bo);
        oxc    = old_x_counts(ao,bo)+old_x_counts(lo,bo);
        icl_cor -= lgamma(oxc+1)-(oxc+1)*log(p*cc_old+1);
        oxc    = old_x_counts(bo,ao)+old_x_counts(bo,lo);
        icl_cor -= lgamma(oxc+1)-(oxc+1)*log(p*cc_old+1);
      }
      
      
    }
  }
  if(l==obk){
    icl_cor=0;
  }
  return icl_cor;
  
  
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
