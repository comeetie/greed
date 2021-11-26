// [[Rcpp::depends(RcppArmadillo)]]
#include "gicl_tools.h"
#include "CoDcSbm.h"
using namespace Rcpp;

CoDcSbm::CoDcSbm(const arma::sp_mat  & xp,int Nri, int Nci,S4 modeli,bool verb){
  model=modeli;
  

  x  = xp;
  N  = Nri+Nci;
  Nc=Nci;
  Nr=Nri;
  alpha = model.slot("alpha");

  if(Rcpp::traits::is_nan<REALSXP>(model.slot("p"))){
    p = arma::accu(x)/(Nr*Nc);
    model.slot("p") = p;
  }else{
    p = model.slot("p");
  }
  
  
  verbose=verb;
  // cst to add 
  cst = 0;
}



void CoDcSbm::set_cl(arma::uvec cl){

  K =arma::max(cl)+1;
  arma::uvec clr = cl.subvec(0,Nr-1);
  row_clusts = arma::unique(clr);
  Kr = row_clusts.n_elem;
  arma::uvec clc = cl.subvec(Nr,N-1);
  col_clusts = arma::unique(clc);
  Kc = col_clusts.n_elem;
  x_counts = gsum_bimat(clr,clc,x,K);
  counts = count(cl,K);
  dr = sum(x_counts.t()).t();
  dc = sum(x_counts).t();
  arma::vec types(K);
  types.fill(0);
  for (arma::uword i = 0;i<row_clusts.n_elem;++i){
    types(row_clusts(i))=1; 
  }
  for (arma::uword i = 0;i<col_clusts.n_elem;++i){
    if(types(col_clusts(i))!=0){
      Rcpp::stop("Invalid partition in co-clustering");      
    }
    types(col_clusts(i))=2;
  }
  clusttypes=types;
}



// a checker !!
double CoDcSbm::icl_prop_cor(const List & obs_stats){
  // compute the first part p(Zr,Zc) from clusters counts
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  // log(p(Z))
  double icl_prop = lgamma(Kr*alpha)+lgamma(Kc*alpha)+arma::accu(lgamma(alpha+counts))-K*lgamma(alpha)-lgamma(Nr+alpha*Kr)-lgamma(Nc+alpha*Kc);
  icl_prop -= lgamma(K*alpha)+arma::accu(lgamma(alpha+counts))-K*lgamma(alpha)-lgamma(arma::accu(counts+alpha));
  // complete with log(p(X|X)) from derived class
  return icl_prop;
}


double CoDcSbm::icl_prop_cor(const List & obs_stats,int oldcl,int newcl, bool dead){

  // compute the first part p(Zr,Zc) from clusters counts
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  double icl_prop = 0;
  if(!dead){
    // both clusters are healthy
    icl_prop -= lgamma(K*alpha)+lgamma(alpha+counts(oldcl))+lgamma(alpha+counts(newcl))-K*lgamma(alpha)-lgamma(K*alpha+N);
  }else{
    // cluster oldclass is dead, count(oldcl)==0 chnage of dimension
    icl_prop -= lgamma((K-1)*alpha)+lgamma(alpha+counts(newcl,0))-(K-1)*lgamma(alpha)-lgamma((K-1)*alpha+N);
  }
  if(!dead){
    // both clusters are healthy
    icl_prop += lgamma(Kr*alpha)+lgamma(Kc*alpha)+lgamma(alpha+counts(oldcl))+lgamma(alpha+counts(newcl))-K*lgamma(alpha)-lgamma(Nr+alpha*Kr)-lgamma(Nc+alpha*Kc);
  }else{
    // cluster oldclass is dead, count(oldcl)==0 chnage of dimension
    if(clusttypes(oldcl)==1){
      icl_prop += lgamma((Kr-1)*alpha)+lgamma(Kc*alpha)+lgamma(alpha+counts(newcl))-(K-1)*lgamma(alpha)-lgamma(Nr+alpha*(Kr-1))-lgamma(Nc+alpha*Kc);
    }
    if(clusttypes(oldcl)==2){
      icl_prop += lgamma(Kr*alpha)+lgamma((Kc-1)*alpha)+lgamma(alpha+counts(newcl))-(K-1)*lgamma(alpha)-lgamma(Nr+alpha*Kr)-lgamma(Nc+alpha*(Kc-1));
    }
  }

  return icl_prop;
}



double CoDcSbm::icl_emiss(const List & obs_stats){
  
  double icl_prop = icl_prop_cor(obs_stats);
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::vec dr =as<arma::vec>(obs_stats["dr"]);
  arma::vec dc =as<arma::vec>(obs_stats["dc"]);
  arma::vec d = dr+dc;
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  arma::mat sub_edges = edges_counts.submat(arma::find(clusttypes==1),arma::find(clusttypes==2)); 
  arma::mat matcount = counts*counts.t();
  
  arma::mat sub_mc = matcount.submat(arma::find(clusttypes==1),arma::find(clusttypes==2)); 
  
  // deggree correction
  double icl_emiss = arma::accu(lgamma(counts)-lgamma(counts+d)+d % log(counts));

  // gamma poisson
  icl_emiss=icl_emiss + arma::accu(lgamma(sub_edges+1)-(sub_edges+1) % log(p*sub_mc+1));

  return icl_emiss+icl_prop;
}

double CoDcSbm::icl_emiss(const List & obs_stats,int oldcl,int newcl, bool dead){
  if(clusttypes(oldcl)!=clusttypes(newcl)){
    return -std::numeric_limits<double>::infinity();
  }
  double icl_prop = icl_prop_cor(obs_stats,oldcl,newcl,dead);
  arma::vec counts =as<arma::vec>(obs_stats["counts"]);
  arma::vec ncounts =as<arma::vec>(obs_stats["counts"]);  
  
  
  arma::vec dr =as<arma::vec>(obs_stats["dr"]);
  arma::vec dc =as<arma::vec>(obs_stats["dc"]);
  arma::vec d = dr+dc;
  arma::mat edges_counts =as<arma::mat>(obs_stats["x_counts"]);
  
  
  // degree correction
  double icl_emiss = lgamma(counts(newcl))-lgamma(counts(newcl)+d(newcl))+d(newcl)*log(counts(newcl));
  if(!dead){
    icl_emiss += lgamma(counts(oldcl))-lgamma(counts(oldcl)+d(oldcl))+d(oldcl)*log(counts(oldcl));
  }
  
  // gamma poisson
  arma::uvec opp_clusts;
  if(clusttypes(oldcl)==1){
    opp_clusts = col_clusts;
  }else{
    opp_clusts = row_clusts;
  }
  for (arma::uword i = 0;i<opp_clusts.n_elem;++i){
    int l=opp_clusts(i);
    double cc=counts(newcl)*counts(l);
    
    if(clusttypes(oldcl)==1){
      icl_emiss += lgamma(edges_counts(newcl,l)+1)-(edges_counts(newcl,l)+1)*log(p*cc+1);
    }else{
      icl_emiss += lgamma(edges_counts(l,newcl)+1)-(edges_counts(l,newcl)+1)*log(p*cc+1);
    }
    if(!dead){
      cc = counts(oldcl)*counts(l);
      if(clusttypes(oldcl)==1){
        icl_emiss += lgamma(edges_counts(oldcl,l)+1)-(edges_counts(oldcl,l)+1)*log(p*cc+1);
      }else{
        icl_emiss += lgamma(edges_counts(l,oldcl)+1)-(edges_counts(l,oldcl)+1)*log(p*cc+1);
      }
    }
  }
  
  return icl_emiss+icl_prop;
}





List CoDcSbm::get_obs_stats(){
  return List::create(Named("counts", counts), Named("dr", dr),Named("dc", dc),Named("x_counts", x_counts));
}

List CoDcSbm::get_obs_stats_cst(){
  arma::sp_mat din_node = sum(x).t();
  arma::sp_mat dout_node = sum(x.t()).t();
  return List::create(Named("dcol", din_node),Named("drow", dout_node));
}

arma::vec CoDcSbm::delta_swap(const int i,arma::uvec & cl, bool almost_dead_cluster, arma::uvec iclust, int K){
  int oldcl = cl(i);
  int cd = 0;
  arma::sp_mat delta_counts;
  
  if(clusttypes(oldcl)==1){
    delta_counts =  gsum_col(cl.subvec(arma::span(Nr,N-1)),x.t(),i,K);
    cd = arma::accu(delta_counts);
  }else{
    delta_counts =  gsum_col(cl.subvec(arma::span(0,Nr-1)),x,i-Nr,K);
    cd = arma::accu(delta_counts);
  }
  
  
  arma::vec delta(K);
  delta.fill(-std::numeric_limits<double>::infinity());
  delta(oldcl)=0;
  
  List old_stats = List::create(Named("counts", counts),Named("dr", dr),Named("dc", dc), Named("x_counts", x_counts));
  int k = 0;
  // for each possible move
  for(arma::uword j = 0; j < iclust.n_elem; ++j) {
    k=iclust(j);
    if((k!=oldcl) & (clusttypes(k)==clusttypes(oldcl))){
      arma::mat new_ec = x_counts;
      arma::vec new_counts = update_count(counts,oldcl,k);
      arma::vec new_dr = dr;
      arma::vec new_dc = dc;
      if(clusttypes(oldcl)==1){
        new_dr(oldcl)=new_dr(oldcl)-cd;
        new_dr(k)=new_dr(k)+cd;
        new_ec.row(oldcl)=new_ec.row(oldcl)-delta_counts.t();
        new_ec.row(k)=new_ec.row(k)+delta_counts.t();
      }else{
        new_dc(oldcl)=new_dc(oldcl)-cd;
        new_dc(k)=new_dc(k)+cd;
        new_ec.col(oldcl)=new_ec.col(oldcl)-delta_counts;
        new_ec.col(k)=new_ec.col(k)+delta_counts;
      }
      List new_stats = List::create(Named("counts", new_counts),Named("dr", new_dr),Named("dc", new_dc), Named("x_counts", new_ec));
      delta(k)=icl_emiss(new_stats,oldcl,k,almost_dead_cluster)-icl_emiss(old_stats,oldcl,k,false);
      
    }
    
  }
  
  return delta;
}


void CoDcSbm::swap_update(const int i,arma::uvec &  cl,bool dead_cluster,const int newcl){
  
  int oldcl = cl(i);
  int cd = 0;
  counts = update_count(counts,oldcl,newcl);
  arma::sp_mat delta_counts;
  if(clusttypes(oldcl)==1){
    delta_counts =  gsum_col(cl.subvec(arma::span(Nr,N-1)),x.t(),i,K);
    cd = arma::accu(delta_counts);
    dr(oldcl)=dr(oldcl)-cd;
    dr(newcl)=dr(newcl)+cd;
    x_counts.row(oldcl)=x_counts.row(oldcl)-delta_counts.t();
    x_counts.row(newcl)=x_counts.row(newcl)+delta_counts.t();
  }else{
    delta_counts =  gsum_col(cl.subvec(arma::span(0,Nr-1)),x,i-Nr,K);
    cd = arma::accu(delta_counts);
    dc(oldcl)=dc(oldcl)-cd;
    dc(newcl)=dc(newcl)+cd;
    x_counts.col(oldcl)=x_counts.col(oldcl)-delta_counts;
    x_counts.col(newcl)=x_counts.col(newcl)+delta_counts;
  }
  
  
  if(dead_cluster){
    counts.shed_row(oldcl);
    dr.shed_row(oldcl);
    dc.shed_row(oldcl);
    x_counts = x_counts(arma::find(arma::linspace(0,K-1,K)!=oldcl),arma::find(arma::linspace(0,K-1,K)!=oldcl));
    
    if(clusttypes(oldcl)==1){
      Kr=Kr-1;
    }
    if(clusttypes(oldcl)==2){
      Kc=Kc-1;
    }
    clusttypes.shed_row(oldcl);
    row_clusts = arma::find(clusttypes==1);
    col_clusts = arma::find(clusttypes!=1);
    --K;
  }
}


double CoDcSbm::delta_merge(int k, int l){
  if(clusttypes(k)!=clusttypes(l)){
    return -std::numeric_limits<double>::infinity();
  }
  
  List old_stats = List::create(Named("counts", counts), Named("dr", dr),Named("dc", dc),Named("x_counts", x_counts));
  
  arma::mat new_ec = x_counts;
  arma::mat new_counts = counts;
  new_counts(l)=new_counts(k)+new_counts(l);
  new_counts(k)=0;
  
  
  arma::vec new_dr = dr;
  arma::vec new_dc = dc;
  if(clusttypes(k)==1){
    new_dr(l)=new_dr(k)+new_dr(l);
    new_ec.row(l)=new_ec.row(k)+new_ec.row(l);
    
  }else{
    new_dc(l)=new_dc(k)+new_dc(l);
    new_ec.col(l)=new_ec.col(k)+new_ec.col(l);
    
  }
  
  
  List new_stats = List::create(Named("counts", new_counts),Named("dr", new_dr),Named("dc", new_dc), Named("x_counts", new_ec));
  
  double delta=icl_emiss(new_stats,k,l,true)-icl_emiss(old_stats,k,l,false);
  
  return delta;
}


double CoDcSbm::delta_merge_correction(int k,int l,int obk,int obl,const List & old_stats){
  // here old refers to the stats before the fusion between obk and obl
  //Rcout << "Je calculs des corrections !!" << std::endl;
  //Rcout << obk << "---- " << obl << std::endl;
  arma::vec old_counts =as<arma::vec>(old_stats["counts"]);
  arma::mat edges_counts =as<arma::mat>(old_stats["x_counts"]);
  double correction = 0;
  int lo,ko;
  double cc;
  
  if((clusttypes(k)==clusttypes(l)) & (clusttypes(k)!=clusttypes(obl))){
    // k,l position in old_stats
    //Rcout << "Corr" << std::endl;
    if(l>=obk){
      lo=l+1;
    }else{
      lo=l;
    }
    if(k>=obk){
      ko=k+1;
    }else{
      ko=k;
    }
    if(clusttypes(k)==1){
      // remove old values
      cc = old_counts(ko)*old_counts(obk);
      correction += lgamma(edges_counts(ko,obk)+1)-(edges_counts(ko,obk)+1)*log(p*cc+1);
      cc = old_counts(ko)*old_counts(obl);
      correction += lgamma(edges_counts(ko,obl)+1)-(edges_counts(ko,obl)+1)*log(p*cc+1);
      cc = old_counts(lo)*old_counts(obk);
      correction += lgamma(edges_counts(lo,obk)+1)-(edges_counts(lo,obk)+1)*log(p*cc+1);
      cc = old_counts(lo)*old_counts(obl);
      correction += lgamma(edges_counts(lo,obl)+1)-(edges_counts(lo,obl)+1)*log(p*cc+1);
      //
      cc = (old_counts(ko)+old_counts(lo))*old_counts(obk);
      correction -= lgamma(edges_counts(lo,obk)+edges_counts(ko,obk)+1)-(edges_counts(lo,obk)+edges_counts(ko,obk)+1)*log(p*cc+1);
      cc = (old_counts(ko)+old_counts(lo))*old_counts(obl);
      correction -= lgamma(edges_counts(lo,obl)+edges_counts(ko,obl)+1)-(edges_counts(lo,obl)+edges_counts(ko,obl)+1)*log(p*cc+1);
      // add new values
      cc = counts(k)*counts(obl);
      correction -= lgamma(x_counts(k,obl)+1)-(x_counts(k,obl)+1)*log(p*cc+1);
      cc = counts(l)*counts(obl);
      correction -= lgamma(x_counts(l,obl)+1)-(x_counts(l,obl)+1)*log(p*cc+1);
      // merge k,l et merge obk,obl
      cc = (counts(k)+counts(l))*counts(obl);
      correction += lgamma(x_counts(l,obl)+x_counts(k,obl)+1)-(x_counts(l,obl)+x_counts(k,obl)+1)*log(p*cc+1);
      
    }else{
      // remove old values no fusion k,l no fusion obk,obl
      cc = old_counts(ko)*old_counts(obk);
      correction += lgamma(edges_counts(obk,ko)+1)-(edges_counts(obk,ko)+1)*log(p*cc+1);
      cc = old_counts(ko)*old_counts(obl);
      correction += lgamma(edges_counts(obl,ko)+1)-(edges_counts(obl,ko)+1)*log(p*cc+1);
      cc = old_counts(lo)*old_counts(obk);
      correction += lgamma(edges_counts(obk,lo)+1)-(edges_counts(obk,lo)+1)*log(p*cc+1);
      cc = old_counts(lo)*old_counts(obl);
      correction += lgamma(edges_counts(obl,lo)+1)-(edges_counts(obl,lo)+1)*log(p*cc+1);
      // remove old values fusion k,l, no fusion obk,obl
      cc = (old_counts(ko)+old_counts(lo))*old_counts(obk);
      correction -= lgamma(edges_counts(obk,lo)+edges_counts(obk,ko)+1)-(edges_counts(obk,lo)+edges_counts(obk,ko)+1)*log(p*cc+1);
      cc = (old_counts(ko)+old_counts(lo))*old_counts(obl);
      correction -= lgamma(edges_counts(obl,lo)+edges_counts(obl,ko)+1)-(edges_counts(obl,lo)+edges_counts(obl,ko)+1)*log(p*cc+1);
      // add new values no fusion k,l fusion obk,obl
      cc = counts(k)*counts(obl);
      correction -= lgamma(x_counts(obl,k)+1)-(x_counts(obl,k)+1)*log(p*cc+1);
      cc = counts(l)*counts(obl);
      correction -= lgamma(x_counts(obl,l)+1)-(x_counts(obl,l)+1)*log(p*cc+1);
      // add new values fusion k,l fusion obk,obl
      cc = (counts(k)+counts(l))*counts(obl);
      correction += lgamma(x_counts(obl,l)+x_counts(obl,k)+1)-(x_counts(obl,l)+x_counts(obl,k)+1)*log(p*cc+1);
    }
    
    
    
  }
  
  
  return correction;
  
}



void CoDcSbm::merge_update(int k,int l){
  

  counts(l) = counts(k)+counts(l);
  counts.shed_row(k);
  if(clusttypes(k)==1){
    Kr=Kr-1;
    dr(l)=dr(k)+dr(l);
    x_counts.row(l) = x_counts.row(k)+x_counts.row(l);
  }else{
    Kc=Kc-1;
    dc(l)=dc(k)+dc(l);
    x_counts.col(l) = x_counts.col(k)+x_counts.col(l);
  }
  x_counts = x_counts(arma::find(arma::linspace(0,K-1,K)!=k),arma::find(arma::linspace(0,K-1,K)!=k));
  dr.shed_row(k);
  dc.shed_row(k);
  
  clusttypes.shed_row(k);
  row_clusts = arma::find(clusttypes==1);
  col_clusts = arma::find(clusttypes!=1);
  --K;
}