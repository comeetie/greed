#' @importFrom graphics plot
#' @include models_classes.R fit_classes.R





#' @title print an icl_path object
#' 
#' @description
#' Print an \code{\link{icl_path-class}} object, model type and number of found clusters are provided.
#' @param x \code{\link{icl_path-class}} object to print
#' @return None (invisible NULL). No return value, called for side effects.
#' @export
setMethod(f = "print", 
          signature = signature("icl_path"),
          definition = function(x){
            cat(paste0("ICL clustering with a ",toupper(x@model@name)," model, ",length(x@obs_stats$counts), " clusters and an icl of ", round(x@icl),"."))
          })



pprint =function(x,M,l){
  K = length(x@obs_stats$counts)
  na = colnames(M)
  D=Matrix::rowSums(M)
  for (k in 1:K){
    ii=which(x@cl==k)
    topk=order(D[ii],decreasing = TRUE)[1:l]
    print(na[ii[topk]])
  }  
}

spy = function(x){
  ij=Matrix::which(x!=0,arr.ind = TRUE)
  gg=data.frame(i=ij[,1],j=ij[,2],v=x[ij])
  ggplot2::ggplot(gg)+ggplot2::geom_point(ggplot2::aes_(y=~-i,x=~j,size=~v))+
    ggplot2::scale_x_continuous("",c())+
    ggplot2::scale_y_continuous("",c())+
    ggplot2::scale_size_area(max_size=1,guide='none')
}

groupspy = function(x,clust){
  x=x[order(clust),order(clust)]
  lims = c(1,cumsum(table(clust)))
  ij=Matrix::which(x!=0,arr.ind = TRUE)
  gg=data.frame(i=ij[,1],j=ij[,2],v=x[ij])
  ggplot2::ggplot(gg)+ggplot2::geom_point(ggplot2::aes_(y=~-i,x=~j,size=~v),alpha=0.5)+
    ggplot2::scale_x_continuous("",breaks = lims,labels = c(),minor_breaks = c(),)+
    ggplot2::scale_y_continuous("",breaks = -lims,labels = c(),minor_breaks = c())+
    ggplot2::scale_size_area(max_size=1,guide='none')
}


plot_front = function(sol){
  if(sol@K<3){
    message("The fit contains only less than 3 clusters, an empty plot was produced.")
    return(ggplot2::ggplot())
  }
  icl = c(sol@icl,sapply(sol@path,function(v){v$icl1}))
  ggicl = data.frame(icl = icl[length(icl):1], K  = 1:length(icl))
  #ggfront= sol@ggtree %>% mutate(x=-H) %>% select(x,K) %>% arrange(x) %>% head(sol@K) %>% left_join(ggicl) %>% mutate(xp=lag(x)) 
  ggfront = merge(sol@ggtree[,c("H","K")],ggicl)
  ggfront$x = -ggfront$H
  ggfront = ggfront[order(ggfront$x)[1:sol@K],]
  ggfront$xp = c(min(ggfront$x)-0.05*diff(range(ggfront$x)), ggfront$x[1:(nrow(ggfront)-1)])
  ggfront = ggfront[ggfront$x!=ggfront$xp,]
  ggplot2::ggplot()+ggplot2::geom_abline(data=ggicl,ggplot2::aes_(intercept=~icl,slope=~K-1),alpha=0.2)+
    ggplot2::geom_point(data=ggfront,ggplot2::aes_(x=~x,y=~icl+x*(K-1)))+
    ggplot2::geom_segment(data=ggfront,ggplot2::aes_(x=~x,y=~icl+x*(K-1),xend=~xp,yend=~icl+xp*(K-1)))+
    ggplot2::scale_x_continuous(expression(paste("log(",alpha,")")),limits = c(min(ggfront$xp),0))+
    ggplot2::ylab("ICL")+
    ggplot2::ggtitle(paste0(toupper(sol@model@name)," model with : ",max(sol@cl)," clusters."))+
    ggplot2::theme_bw()
}

lapath = function(x){
  if(x@K<3){
    message("The fit contains less than 3 clusters, an empty plot was produced.")
    return(ggplot2::ggplot())
  }
  gg = data.frame(k=sapply(x@path,function(p){p$K}),logalpha=sapply(x@path,function(p){p$logalpha}))
  gg = rbind(gg,data.frame(k=length(x@obs_stats$counts),logalpha=x@logalpha)) 
  ggplot2::ggplot(data=gg)+ggplot2::geom_line(ggplot2::aes_(x=~k,y=~-logalpha))+
    ggplot2::geom_point(ggplot2::aes_(x=~k,y=~-logalpha))+
    ggplot2::ylab(expression(paste("-log(",alpha,")")))+
    ggplot2::ggtitle(paste0(toupper(x@model@name)," ",length(x@obs_stats$counts)," clusters"))+
    ggplot2::theme_bw()
}


iclpath = function(x){
  if(x@K==1){
    message("The fit contains only one cluster, an empty plot was produced.")
    return(ggplot2::ggplot())
  }
  gg = data.frame(k=sapply(x@path,function(p){length(p$counts)}),icl=sapply(x@path,function(p){p$icl}))
  gg = rbind(gg,data.frame(k=length(x@obs_stats$counts),icl=x@icl)) 
  ggplot2::ggplot(data=gg)+ggplot2::geom_line(ggplot2::aes_(x=~k,y=~icl))+
    ggplot2::geom_point(ggplot2::aes_(x=~k,y=~icl))+
    ggplot2::ylab(expression(paste("ICL")))+
    ggplot2::ggtitle(paste0(toupper(x@model@name)," ",length(x@obs_stats$counts)," clusters"))+
    ggplot2::theme_bw()
}


# Node link visualisations

nodelink = function(sol){
  ij = Matrix::which(sol@obs_stats$x_counts>0,arr.ind = TRUE)
  ld = sol@obs_stats$x_counts
  #/(sol@obs_stats$counts%*%t(sol@obs_stats$counts))
  ij = ij[ij[,1]!=ij[,2],]
  gglink = data.frame(from=ij[,1],to=ij[,2],p=ld[ij])
  gglink$y=ifelse(gglink$from<gglink$to,-0.3,0.3)
  ggnode = data.frame(i=1:length(sol@obs_stats$counts),pi=diag(sol@obs_stats$x_counts))
  gl = ggplot2::guide_legend()
  ggplot2::ggplot()+ggplot2::geom_curve(data=gglink,ggplot2::aes_(x=~from,xend=~to,y=~y,yend=~y,size=~p,alpha=~p),arrow=grid::arrow(length = grid::unit(2,"mm")),curvature = 0.7)+
    ggplot2::scale_x_continuous("",c())+
    ggplot2::scale_y_continuous("",c(),limits = c(-5,5))+
    ggplot2::scale_alpha("Link density:",limits=c(0,max(gglink$p)),guide="none")+
    ggplot2::scale_size_area("Clusters size:",limits=c(0,max(ggnode$pi)),max_size = 4,guide="none")+
    ggplot2::geom_point(data=ggnode,ggplot2::aes_(x=~i,y=~0,size=~pi))+
    ggplot2::ggtitle(paste0(toupper(sol@model@name)," model with : ",max(sol@cl)," clusters."))+
    ggplot2::theme_minimal()
}

#' nodelinklab
#' @param sol \code{\link{mm_path-class}} object to be plot
#' @param labels a vector of cluster labels
#' @param s threshold for links
#' @return a ggplot2 graph
#' @export 
nodelinklab = function(sol,labels,s=0){
  ij = Matrix::which(sol@obs_stats$x_counts>0,arr.ind = TRUE)
  ld = sol@obs_stats$x_counts
  #/(sol@obs_stats$counts%*%t(sol@obs_stats$counts))
  ij = ij[ij[,1]!=ij[,2],]
  gglink = data.frame(from=ij[,1],to=ij[,2],p=ld[ij])
  gglink$y=ifelse(gglink$from<gglink$to,-0.3,0.3)
  ggnode = data.frame(i=1:length(sol@obs_stats$counts),pi=diag(sol@obs_stats$x_counts))
  gl = ggplot2::guide_legend()
  ggplot2::ggplot()+ggplot2::geom_curve(data=gglink[gglink$p>s,],ggplot2::aes_(x=~from,xend=~to,y=~y,yend=~y,size=~p,alpha=~p),arrow=grid::arrow(length = grid::unit(2,"mm")),curvature = 0.7)+
    ggplot2::scale_x_continuous("",c())+
    ggplot2::scale_y_continuous("",c(),limits = c(-6,6))+
    ggplot2::scale_alpha("Link density:",limits=c(0,max(gglink$p)),guide="none")+
    ggplot2::scale_size_area("Clusters size:",limits=c(0,max(c(ggnode$pi,gglink$p))),max_size = 4,guide="none")+
    ggplot2::geom_point(data=ggnode,ggplot2::aes_(x=~i,y=~-0.1,size=~pi))+
    ggplot2::ggtitle(paste0(toupper(sol@model@name)," model with : ",max(sol@cl)," clusters."))+
    ggplot2::geom_text(data=data.frame(x=1:length(labels),label=labels),ggplot2::aes_(x=~x,label=~label,y=~0.05))+
    ggplot2::theme_minimal()
}

co_nodelink = function(sol){
  ij = Matrix::which(sol@obs_stats$co_x_counts>0,arr.ind = TRUE)
  ld = sol@obs_stats$co_x_counts
  cccol = as.numeric(table(sol@clcol))
  ccrow = as.numeric(table(sol@clrow))
  
  gglink = data.frame(from=ij[,1],to=ij[,2],p=ld[ij]) 
  gglink$nrow = ccrow[gglink$from]
  gglink$ncol = cccol[gglink$to]
  gglink$pn=gglink$p/(gglink$nrow*gglink$ncol)
  gglink=gglink[order(gglink$pn),]
  ggnode = rbind(data.frame(type="col",n=cccol,i=1:sol@Kcol),data.frame(type="row",n=ccrow,i=1:sol@Krow) )
  
  
  ggplot2::ggplot()+ggplot2::geom_curve(data=gglink,ggplot2::aes_(y=~from,yend=~max(from)+1,x=-5,xend=~to,size=~p,alpha=~p),curvature = 0.35)+
    ggplot2::geom_point(data=ggnode,ggplot2::aes_(x=~ifelse(type=="row",-5,i),y=~ifelse(type=="row",i,max(gglink$from)+1),size=~n^2))+
    ggplot2::theme_minimal()+
    ggplot2::scale_size_area("Clusters size:",limits=c(0,max(ggnode$n)^2),max_size = 7,guide="none")+
    ggplot2::scale_y_continuous("",c())+
    ggplot2::scale_x_continuous("",c())+
    ggplot2::scale_alpha("Link density:",limits=c(0,max(gglink$p)),guide="none")+
    ggplot2::ggtitle("Poisson, co-clustering")
}

# dendogram visualisation

dendo = function(x){
  if(x@K<3){
    message("The fit contains less than 3 clusters, an empty plot was produced.")
    return(ggplot2::ggplot())
  }
  ggtree = x@ggtree
  tree=ggplot2::ggplot()+ggplot2::geom_segment(data=ggtree[ggtree$node %in% ggtree$tree,],ggplot2::aes_(x=~xmin,y=~H,xend=~xmax,yend=~H))+
    ggplot2::geom_segment(data=ggtree[-1,],ggplot2::aes_(x=~x,y=~H,xend=~x,yend=~Hend))+
    ggplot2::scale_x_continuous("",breaks=c())+
    ggplot2::ylab(expression(paste("-log(",alpha,")")))+
    ggplot2::ggtitle(paste0(toupper(x@model@name)," ",length(x@obs_stats$counts)," clusters, dendogram"))+
    ggplot2::theme_bw()
  if(x@K<max(x@ggtree$K)){
    hc = x@ggtree$H[x@ggtree$K==x@K]
    tree=tree+ggplot2::geom_hline(ggplot2::aes(yintercept=hc),color='red',size=0.8,linetype="dashed",alpha=0.8)
  }
  tree
}

co_dendo = function(x){
  ggtree = x@ggtreerow
  rowt = ggplot2::ggplot()+ggplot2::geom_segment(data=ggtree[ggtree$node %in% ggtree$tree,],ggplot2::aes_(x=~xmin,y=~H,xend=~xmax,yend=~H))+
    ggplot2::geom_segment(data=ggtree,ggplot2::aes_(x=~x,y=~H,xend=~x,yend=~Hend))+
    ggplot2::scale_x_continuous("",breaks=c())+
    ggplot2::ylab(expression(paste("-log(",alpha,")")))+
    ggplot2::ggtitle(paste0(x@Krow," row clusters"))+
    ggplot2::theme_bw()
  if(x@K<max(x@ggtree$K)){
    hc = x@ggtree$H[x@ggtree$K==x@K]
    rowt=rowt+ggplot2::geom_hline(ggplot2::aes(yintercept=hc),color='red',size=0.8,linetype="dashed",alpha=0.8)
  }
  
  ggtree = x@ggtreecol
  colt = ggplot2::ggplot()+ggplot2::geom_segment(data=ggtree[ggtree$node %in% ggtree$tree,],ggplot2::aes_(x=~xmin,y=~H,xend=~xmax,yend=~H))+
    ggplot2::geom_segment(data=ggtree,ggplot2::aes_(x=~x,y=~H,xend=~x,yend=~Hend))+
    ggplot2::scale_x_continuous("",breaks=c())+
    ggplot2::ylab(expression(paste("-log(",alpha,")")))+
    ggplot2::ggtitle(paste0(x@Kcol," column clusters"))+
    ggplot2::theme_bw()
  if(x@K<max(x@ggtree$K)){
    hc = x@ggtree$H[x@ggtree$K==x@K]
    colt=colt+ggplot2::geom_hline(ggplot2::aes(yintercept=hc),color='red',size=0.8,linetype="dashed",alpha=0.8)
  }
  ggpubr::ggarrange(rowt,colt)
}


# matrice blocks visualisation

graph_blocks = function(x){
  K = length(x@obs_stats$counts)
  gg=data.frame(kc=rep(cumsum(x@obs_stats$counts),each=K),
                lc=rep(cumsum(x@obs_stats$counts),K),
                sizek = rep(x@obs_stats$counts,each=K),
                sizel = rep(x@obs_stats$counts,K), 
                count=as.vector(x@obs_stats$x_counts))
  ggplot2::ggplot(gg[gg$count>0,])+ggplot2::geom_tile(ggplot2::aes_(x=~kc-sizek/2,y=~lc-sizel/2,width=~sizek,height=~sizel,fill=~count/(sizek*sizel),alpha=~count/(sizek*sizel)))+
    ggplot2::scale_fill_distiller("Link density",palette="YlOrRd",direction = 1,guide = ggplot2::guide_legend(),limits=c(0,max(gg$count/(gg$sizek*gg$sizel))))+
    ggplot2::scale_alpha("Link density",range=c(0,1),limits=c(0,max(gg$count/(gg$sizek*gg$sizel))))+
    ggplot2::ggtitle(paste0(toupper(x@model@name)," model with : ",max(x@cl)," clusters."))+
    ggplot2::scale_x_continuous("",breaks=cumsum(x@obs_stats$counts),
                                labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),
                                minor_breaks = NULL,expand = ggplot2::expansion(mult = 0, add = 0))+
    ggplot2::scale_y_continuous("",breaks=cumsum(x@obs_stats$counts),
                                labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),
                                minor_breaks = NULL,expand = ggplot2::expansion(mult = 0, add = 0))+
    ggplot2::coord_fixed()+ggplot2::theme_bw()
}


co_blocks = function(x){
  K = x@Krow
  D = x@Kcol
  ccrow = cumsum(table(x@clrow))
  cccol = cumsum(table(x@clcol))
  gg=data.frame(kc=rep(ccrow,D),
                lc=rep(cccol,each=K),
                sizek = rep(table(x@clrow),D),
                sizel = rep(table(x@clcol),each=K), 
                count=as.vector(x@obs_stats$co_x_counts))
  ylab  = round(100*table(x@clcol)/sum(table(x@clcol)))
  xlab = round(100*table(x@clrow)/sum(table(x@clrow)))
  
  ggplot2::ggplot(gg)+ggplot2::geom_tile(ggplot2::aes_(y=~kc-sizek/2,x=~lc-sizel/2,height=~sizek,width=~sizel,fill=~count/(sizek*sizel),alpha=~count/(sizek*sizel)))+
    ggplot2::scale_fill_distiller("E[X]",palette="YlOrRd",direction = 1,guide = ggplot2::guide_legend(),limits=c(0,max(gg$count/(gg$sizek*gg$sizel))))+
    ggplot2::scale_alpha("E[X]",range=c(0,1),limits=c(0,max(gg$count/(gg$sizek*gg$sizel))))+
    ggplot2::ggtitle(paste0("Co-clustering with : ",max(x@cl)," clusters."))+
    ggplot2::scale_x_continuous("Col clusters",breaks=cccol,labels=ifelse(ylab>5,paste0(ylab,"%"),""),minor_breaks = NULL,expand = ggplot2::expansion(mult = 0, add = 0))+
    ggplot2::scale_y_continuous("Row clusters",breaks=ccrow,labels =ifelse(xlab>5,paste0(xlab,"%"),""),minor_breaks = NULL,expand = ggplot2::expansion(mult = 0, add = 0))+
    ggplot2::theme_bw()      
}

mat_blocks = function(x){
  K = length(x@obs_stats$counts)
  D = dim(x@obs_stats$x_counts)[1]
  gg=data.frame(kc=rep(cumsum(x@obs_stats$counts),D),
                lc=rep(1:D,each=K),
                sizek = rep(x@obs_stats$counts,D),
                sizel = rep(1,K*D), 
                count=as.vector(Matrix::t(x@obs_stats$x_counts)/Matrix::rowSums(Matrix::t(x@obs_stats$x_counts))))
  
  ggplot2::ggplot(gg)+ggplot2::geom_tile(ggplot2::aes_(y=~kc-sizek/2,x=~lc-sizel/2,height=~sizek,width=~sizel,fill=~log(count),alpha=~count))+
    ggplot2::scale_fill_distiller("E[X]",palette="YlOrRd",direction = 1,guide = ggplot2::guide_legend(),limits=c(1,log(max(gg$count))))+
    ggplot2::scale_alpha("E[X]",range=c(0,1),limits=c(0,max(gg$count)))+
    ggplot2::ggtitle(paste0("MM Model with : ",max(x@cl)," clusters."))+
    ggplot2::scale_x_continuous("Features",breaks=1:D,labels=rep("",D),minor_breaks = NULL,expand = ggplot2::expansion(mult = 0, add = 0))+
    ggplot2::scale_y_continuous("Clusters",breaks=cumsum(x@obs_stats$counts),labels = paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),minor_breaks = NULL,expand = ggplot2::expansion(mult = 0, add = 0))+
    ggplot2::theme_bw()      
}



#' graph_balance
#' @param x a \code{\link{sbm_fit-class}} object to be plot
#' @return a ggplot2 graph
#' @export 
graph_balance = function(x){
  K = length(x@obs_stats$counts)
  B=x@obs_stats$x_counts-t(x@obs_stats$x_counts)
  gg=data.frame(kc=rep(cumsum(x@obs_stats$counts),each=K),
                lc=rep(cumsum(x@obs_stats$counts),K),
                sizek = rep(x@obs_stats$counts,each=K),
                sizel = rep(x@obs_stats$counts,K),
                dk = as.vector(x@obs_stats$dout%*%t(x@obs_stats$din)),
                count=as.vector(B))
  vm =max(abs(gg$count))
  ggplot2::ggplot(gg)+ggplot2::geom_tile(ggplot2::aes_(x=~kc-sizek/2,y=~lc-sizel/2,width=~sizek,height=~sizel,fill=~count),alpha=0.7)+
    ggplot2::scale_fill_distiller("Balance :",direction = 1,palette="RdBu",limits=c(-vm,vm))+
    ggplot2::ggtitle(paste0(toupper(x@model@name)," model with : ",max(x@cl)," clusters."))+
    ggplot2::scale_x_continuous("",breaks=cumsum(x@obs_stats$counts),labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),minor_breaks = NULL)+
    ggplot2::scale_y_continuous("",breaks=cumsum(x@obs_stats$counts),labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),minor_breaks = NULL)+
    ggplot2::coord_fixed()+ggplot2::theme_bw()
}




mat_reg = function(x){
  K = length(x@obs_stats$counts)
  D = length(x@obs_stats$regs[[1]]$mu)
  gg=data.frame(kc=rep(cumsum(x@obs_stats$counts),D),
                lc=rep(1:D,each=K),
                sizek = rep(x@obs_stats$counts,D),
                sizel = rep(1,K*D), 
                count=as.vector(sapply(x@obs_stats$regs,function(reg){reg$mu})))
  ggplot2::ggplot(gg)+ggplot2::geom_tile(ggplot2::aes_(y=~kc-sizek/2,x=~lc-sizel/2,height=~sizek,width=~sizel,fill=~count,alpha=~count))+
    ggplot2::scale_fill_distiller(expression(paste(" ",beta," ")),type="seq",direction = 1,palette = 2)+
    ggplot2::scale_alpha(expression(paste(" ",beta," ")))+
    ggplot2::ggtitle(paste0("Mixture of Regression Model with : ",max(x@cl)," clusters."))+
    ggplot2::scale_x_continuous("Features",breaks=1:D,labels=rep("",D),minor_breaks = NULL)+
    ggplot2::scale_y_continuous("Clusters",breaks=cumsum(x@obs_stats$counts),labels = paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),minor_breaks = NULL)+
    ggplot2::theme_bw()      
}




mat_reg_line = function(x,Xd,yd){
  K = length(x@obs_stats$counts)
  D = length(x@obs_stats$regs[[1]]$mu)
  ggp= data.frame(x=Xd[,1],y=yd,K=factor(x@cl,levels=1:K))
  gg=data.frame(y=as.vector(cbind(seq(min(ggp$x),max(ggp$x),length.out = 20),rep(1,20))%*%sapply(x@obs_stats$regs,function(reg){reg$mu})),
                x=rep(seq(min(ggp$x),max(ggp$x),length.out = 20),K),K=rep(1:K,each=20))
 ggplot2::ggplot()+
    ggplot2::geom_point(data=ggp[,1:2],ggplot2::aes_(x=~x,y=~y),alpha=0.05)+
    ggplot2::geom_path(data=gg,ggplot2::aes_(x=~x,y=~y,group=~K))+
    ggplot2::geom_point(data=ggp,ggplot2::aes_(x=~x,y=~y,col=~K))+
    ggplot2::ggtitle(paste0("Mixture of Regression Model with : ",max(x@cl)," clusters."))+
    ggplot2::facet_wrap(~K)+
    ggplot2::theme_bw()

}







bi_plot = function(x){
  
  K = x@Krow
  D = x@Kcol
  ccrow = cumsum(table(x@clrow))
  cccol = cumsum(table(x@clcol))
  gg=data.frame(kc=rep(ccrow,D),
                lc=rep(cccol,each=K),
                sizek = rep(table(x@clrow),D),
                sizel = rep(table(x@clcol),each=K), 
                count=as.vector(x@obs_stats$co_x_counts))
  ylab  = round(100*table(x@clcol)/sum(table(x@clcol)))
  xlab = round(100*table(x@clrow)/sum(table(x@clrow)))
  
  
  theme_mm = ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                             axis.text.x=ggplot2::element_blank(),
                             axis.ticks.x=ggplot2::element_blank(),
                             axis.title.y=ggplot2::element_blank(),
                             axis.text.y=ggplot2::element_blank(),
                             axis.ticks.y=ggplot2::element_blank(),
                            panel.border = ggplot2::element_blank())      
  mat = ggplot2::ggplot(gg)+ggplot2::geom_tile(ggplot2::aes_(y=~kc-sizek/2,x=~lc-sizel/2,height=~sizek,width=~sizel,fill=~count/(sizek*sizel),alpha=~count/(sizek*sizel)))+
    ggplot2::scale_fill_distiller("E[X]",type="seq",direction = 1)+
    ggplot2::scale_alpha("E[X]")+
    ggplot2::scale_x_continuous("Col clusters",breaks=cccol,labels=ifelse(ylab>5,paste0(ylab,"%"),""),minor_breaks = NULL)+
    ggplot2::scale_y_continuous("Row clusters",breaks=ccrow,labels =ifelse(xlab>5,paste0(xlab,"%"),""),minor_breaks = NULL)+
    ggplot2::theme_bw()+theme_mm
  
  ggtree = x@ggtreerow
  rowt = ggplot2::ggplot()+ggplot2::geom_segment(data=ggtree[ggtree$node %in% ggtree$tree,],ggplot2::aes_(y=~xmin,x=~H,yend=~xmax,xend=~H))+
    ggplot2::geom_segment(data=ggtree[ggtree$Hend!=-1,],ggplot2::aes_(x=~H,y=~x,yend=~x,xend=~Hend))+
    ggplot2::scale_y_continuous("",breaks=c())+
    ggplot2::scale_x_reverse("",breaks=c())+
    ggplot2::theme_bw()+theme_mm
  if(x@K<max(x@ggtree$K)){
    hc = x@ggtree$H[x@ggtree$K==x@K]
    rowt=rowt+ggplot2::geom_hline(ggplot2::aes(yintercept=hc),color='red',size=0.8,linetype="dashed",alpha=0.8)
  }
  
  ggtree = x@ggtreecol
  colt = ggplot2::ggplot()+ggplot2::geom_segment(data=ggtree[ggtree$node %in% ggtree$tree,],ggplot2::aes_(x=~xmin,y=~H,xend=~xmax,yend=~H))+
    ggplot2::geom_segment(data=ggtree[ggtree$Hend!=-1,],ggplot2::aes_(x=~x,y=~H,xend=~x,yend=~Hend))+
    ggplot2::scale_x_continuous("",breaks=c())+
    ggplot2::ylab(expression(paste("-log(",alpha,")")))+
    ggplot2::ggtitle(paste0(x@Kcol," column clusters"))+
    ggplot2::theme_bw()
  if(x@K<max(x@ggtree$K)){
    hc = x@ggtree$H[x@ggtree$K==x@K]
    colt=colt+ggplot2::geom_hline(ggplot2::aes(yintercept=hc),color='red',size=0.8,linetype="dashed",alpha=0.8)
  }
  ggpubr::ggarrange(rowt,mat,common.legend = TRUE,widths = c(1,2))
}



graph_blocks_cube = function(x){
  K = length(x@obs_stats$counts)
  M=dim(x@obs_stats$x_counts)[3]
  ggl = lapply(1:M,function(m){data.frame(kc=rep(cumsum(x@obs_stats$counts),each=K),
                                     lc=rep(cumsum(x@obs_stats$counts),K),
                                     sizek = rep(x@obs_stats$counts,each=K),
                                     sizel = rep(x@obs_stats$counts,K), 
                                     count=as.vector(x@obs_stats$x_counts[,,m]),m=paste("Slice ",m),stringsAsFactors = FALSE)})
  
  gg=do.call(rbind,ggl);
  
  ggplot2::ggplot(gg[gg$count>0,])+ggplot2::geom_tile(ggplot2::aes_(x=~kc-sizek/2,y=~lc-sizel/2,width=~sizek,height=~sizel,fill=~count/(sizek*sizel),alpha=~count/(sizek*sizel)))+
    ggplot2::facet_wrap(~m)+
    ggplot2::scale_fill_distiller("Link density",palette="YlOrRd",direction = 1,guide = ggplot2::guide_legend(),limits=c(0,max(gg$count/(gg$sizek*gg$sizel))))+
    ggplot2::scale_alpha("Link density",range=c(0,1),limits=c(0,max(gg$count/(gg$sizek*gg$sizel))))+
    ggplot2::ggtitle(paste0(toupper(x@model@name)," model with : ",max(x@cl)," clusters."))+
    ggplot2::scale_x_continuous("",breaks=cumsum(x@obs_stats$counts),
                                labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),
                                minor_breaks = NULL,expand = ggplot2::expansion(mult = 0, add = 0))+
    ggplot2::scale_y_continuous("",breaks=cumsum(x@obs_stats$counts),
                                labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),
                                minor_breaks = NULL,expand = ggplot2::expansion(mult = 0, add = 0))+
    ggplot2::coord_fixed()+ggplot2::theme_bw()
}


nodelink_cube = function(sol){
  M=dim(sol@obs_stats$x_counts)[3]
  gglink_list = lapply(1:M,function(m){
  ij = which(sol@obs_stats$x_counts[,,m]>0,arr.ind = TRUE)
  ld = sol@obs_stats$x_counts[,,m]
  #/(sol@obs_stats$counts%*%t(sol@obs_stats$counts))
  ij = ij[ij[,1]!=ij[,2],]
  gglink = data.frame(from=ij[,1],to=ij[,2],p=ld[ij])
  gglink$y=ifelse(gglink$from<gglink$to,-0.3,0.3)
  gglink$m=paste("Slice",m)
  gglink
  })
  gglink =do.call(rbind,gglink_list)
  ggnode_list=lapply(1:M, function(m){data.frame(i=1:length(sol@obs_stats$counts),pi=diag(sol@obs_stats$x_counts[,,m]),m=paste("Slice",m),stringsAsFactors = FALSE)})
  ggnode =do.call(rbind,ggnode_list)
  gl = ggplot2::guide_legend()
  ggplot2::ggplot()+ggplot2::geom_curve(data=gglink,ggplot2::aes_(x=~from,xend=~to,y=~y,yend=~y,size=~p,alpha=~p),arrow=grid::arrow(length = grid::unit(2,"mm")),curvature = 0.7)+
    ggplot2::facet_wrap(~m)+
    ggplot2::scale_x_continuous("",c())+
    ggplot2::scale_y_continuous("",c(),limits = c(-5,5))+
    ggplot2::scale_alpha("Link density:",limits=c(0,max(gglink$p)),guide="none")+
    ggplot2::scale_size_area("Clusters size:",limits=c(0,max(ggnode$pi)),max_size = 4,guide="none")+
    ggplot2::geom_point(data=ggnode,ggplot2::aes_(x=~i,y=~0,size=~pi))+
    ggplot2::ggtitle(paste0(toupper(sol@model@name)," model with : ",max(sol@cl)," clusters."))+
    ggplot2::theme_minimal()
}

#' @title Make a matrix of plots with a given data and gmm fitted parameters
#' 
#' @description 
#' Make a matrix of plots with a given data and gmm fitted parameters with ellipses.
#' @param sol a \code{\link{gmm_fit-class}} or \code{\link{diaggmm_fit-class}}
#' @param X the data used for the fit a data.frame or matrix.
#' @return a \code{\link{ggplot2}} graphic 
#' @export
gmmpairs = function(sol,X){
  if(!(methods::is(sol,"gmm_fit") | methods::is(sol,"diaggmm_fit") )){
    stop("Input sol must be a gmm_fit or diaggmm_fit object.",call. = FALSE)
  }
  if(!(methods::is(X,"matrix") | methods::is(X,"data.frame") )){
    stop("Input X must be a matrix or data.frame.",call. = FALSE)
  }
  if(nrow(X)!=length(sol@cl) | ncol(X)!=length(sol@obs_stats$regs[[1]]$m) ){
    stop("Dimension mismatch between the fitted model and the data.",call. = FALSE)
  }
  vnames = names(X)
  plts.df = list()
  ii=1
  nb_ech = 500
  params = coef(sol)
  for(i in 1: ncol(X)){
    for (j in 1:ncol(X)){
      if(i==j){
        #plt = ggally_text(paste("Plot #", i, sep = ""))
        limsi=c(range(X[,vnames[i]])[1]-diff(range(X[,vnames[i]]))*0.1,range(X[,i])[2]+diff(range(X[,vnames[i]]))*0.1)
        x = seq(limsi[1],limsi[2],length.out = nb_ech)
        pdfs    = lapply(1:sol@K,function(k){params$pi[k]*stats::dnorm(x,mean=params$muk[[k]][j],sd=sqrt(params$Sigmak[[k]][j,j]))})
        pdfs.df = data.frame(pdf=do.call(c,pdfs),cl=rep(1:sol@K,each=nb_ech),x=rep(x,sol@K))
        pdfmix  = data.frame(pdf=colSums(do.call(rbind,pdfs)),x=x)
        
        plt = ggplot2::ggplot()+
          ggplot2::geom_line(data=pdfmix,ggplot2::aes_(x=~x,y=~pdf))+
          ggplot2::geom_line(data=pdfs.df,ggplot2::aes_(x=~x,y=~pdf,color=~factor(cl),group=~cl))+
          ggplot2::geom_segment(data=X,ggplot2::aes_(x=as.name(vnames[i]),xend=as.name(vnames[i]),y=0,yend=max(pdfs.df$pdf)*0.05,col=factor(sol@cl)))+
          ggplot2::scale_x_continuous(limits=limsi)
      }else{
        elipses.df = do.call(rbind,lapply(1:sol@K,function(k){ 
          el.df = ellips(params$muk[[k]][c(i,j)],params$Sigmak[[k]][c(i,j),c(i,j)],nb_ech = nb_ech)
          el.df$cl=k
          el.df
          }))
        
        limsi=c(range(X[,vnames[i]])[1]-diff(range(X[,vnames[i]]))*0.1,range(X[,i])[2]+diff(range(X[,vnames[i]]))*0.1)
        limsj=c(range(X[,vnames[j]])[1]-diff(range(X[,vnames[j]]))*0.1,range(X[,j])[2]+diff(range(X[,vnames[j]]))*0.1)
        #elipses.df = elipses.df[elipses.df$x >= limsi[1] & elipses.df$x <= limsi[2] & elipses.df$y >= limsj[1] & elipses.df$y <= limsj[2], ]
        plt = ggplot2::ggplot(X,ggplot2::aes_(as.name(vnames[i]),as.name(vnames[j]),color=factor(sol@cl)))+
          ggplot2::geom_point()+
          ggplot2::geom_path(data=elipses.df,ggplot2::aes_(x=~x,y=~y,color=~factor(cl),group=~cl),na.rm=TRUE)+
          ggplot2::scale_x_continuous(limits=limsi)+
          ggplot2::scale_y_continuous(limits=limsj)
      }
      plts.df[[ii]]= plt
      ii=ii+1
    }
  }
  pm <- GGally::ggmatrix(plts.df,3, 3,xAxisLabels = vnames,yAxisLabels = vnames,byrow = FALSE,legend = c(1,i))+ggplot2::theme_bw()
  pm
  
}



ellips <- function(mu = c(0,0), Sigma=diag(rep(1,2)), nb_ech = 100, l=stats::qchisq(.95, df=2)){
  t <- seq(0, 2*pi, len=nb_ech)
  a <- sqrt(l*eigen(Sigma)$values[1])
  b <- sqrt(l*eigen(Sigma)$values[2])
  x <- a*cos(t)
  y <- b*sin(t)
  X <- cbind(x, y)
  R <- eigen(Sigma)$vectors
  Xf = X%*%t(R)
  data.frame(x=Xf[,1]+mu[1],y=Xf[,2]+mu[2])
}





