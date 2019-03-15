#' @importFrom graphics plot
#' @include models_classes.R fit_classes.R
#' @title Plot a clustering results
#' @description Main methods to explore clusterings results visualy.  
#' @name plot
NULL

#' @rdname plot
#' @param x \code{\link{sbm_fit-class}} object to be ploted
#' @param type 'blocks' for a block matrix view, 'nodelink' for a node link graph, 'tree' for a dendogram, 'path' to show the evolution of -log(alpha) 
#' with respect to the number of clusters.
#' @return A ggplot2 graphics which summarize the results.
#' @export
setMethod(f = "plot", 
          signature = signature("sbm_fit","missing"),
          definition = function(x,type="blocks"){
            switch(type,blocks=graph_blocks(x),nodelink=nodelink(x))
          });


#' @rdname plot
#' @param x \code{\link{dcsbm_fit-class}} object to be ploted
#' @export
setMethod(f = "plot", 
          signature = signature("dcsbm_fit","missing"),
          definition = function(x,type="blocks"){
            switch(type,blocks=graph_blocks(x),nodelink=nodelink(x))
          });

#' @rdname plot
#' @param x \code{\link{mm_fit-class}} object to be ploted
#' @export
setMethod(f = "plot", 
          signature = signature("mm_fit","missing"),
          definition = function(x,type='blocks'){
            mat_blocks(x)      
          });


#' @rdname plot
#' @param x \code{\link{dcsbm_path-class}} object to be ploted
#' @export
setMethod(f = "plot", 
          signature = signature("dcsbm_path","missing"),
          definition = function(x,type='blocks'){
            switch(type,tree = {
              dendo(x)
            },
            path ={
              lapath(x)
            },
            blocks ={
              callNextMethod()
            },
            nodelink={
              callNextMethod()
            })   
          });


#' @rdname plot
#' @param x \code{\link{sbm_path-class}} object to ploted
#' @export
setMethod(f = "plot", 
          signature = signature("sbm_path","missing"),
          definition = function(x,type='blocks'){
            switch(type,tree = {
              dendo(x)
            },
            path ={
              lapath(x)
            },
            blocks ={
              callNextMethod()
            },
            nodelink={
              callNextMethod()
            })   
          });

#' @rdname plot
#' @param x \code{\link{mm_path-class}} object to be plot
#' @export
setMethod(f = "plot", 
          signature = signature("mm_path","missing"),
          definition = function(x,type='blocks'){
            switch(type,tree = {
              dendo(x)
            },
            path ={
              lapath(x)
            },
            blocks ={
              callNextMethod()
            },
            nodelink={
              callNextMethod()
            })   
          });



nodelink = function(sol){
  ij = Matrix::which(sol@obs_stats$x_counts>0,arr.ind = TRUE)
  ld = sol@obs_stats$x_counts
  #/(sol@obs_stats$counts%*%t(sol@obs_stats$counts))
  ij = ij[ij[,1]!=ij[,2],]
  gglink = data.frame(from=ij[,1],to=ij[,2],p=ld[ij])
  ggnode = data.frame(i=1:length(sol@obs_stats$counts),pi=diag(sol@obs_stats$x_counts))
  gl = ggplot2::guide_legend()
  ggplot2::ggplot()+ggplot2::geom_curve(data=gglink,ggplot2::aes(x=from,xend=to,y=ifelse(from<to,-0.3,0.3),yend=ifelse(from<to,-0.3,0.3),size=p,alpha=p),arrow=grid::arrow(length = unit(2,"mm")),curvature = 0.7)+
    ggplot2::scale_x_continuous("",c())+
    ggplot2::scale_y_continuous("",c(),limits = c(-5,5))+
    #scale_size("Ld sizes:",limits=c(0,max(gglink$p)))+
    ggplot2::scale_alpha("Link density:",limits=c(0,max(gglink$p)),guide="none")+
    ggplot2::scale_size_area("Clusters size:",limits=c(0,max(ggnode$pi)),max_size = 4,guide="none")+
    ggplot2::geom_point(data=ggnode,aes(x=i,y=0,size=pi))+
    ggplot2::ggtitle(paste0(toupper(sol@model@name)," model with : ",max(sol@cl)," clusters."))+
    #coord_equal()+
    ggplot2::theme_minimal()
}

graph_blocks = function(x){
  K = length(x@obs_stats$counts)
  gg=data.frame(kc=rep(cumsum(x@obs_stats$counts),each=K),
                lc=rep(cumsum(x@obs_stats$counts),K),
                sizek = rep(x@obs_stats$counts,each=K),
                sizel = rep(x@obs_stats$counts,K), 
                count=as.vector(x@obs_stats$x_counts))
  ggplot2::ggplot(gg[gg$count>0,])+ggplot2::geom_tile(ggplot2::aes(x=kc-sizek/2,y=lc-sizel/2,width=sizek,height=sizel,fill=log(count/(sizek*sizel)),alpha=count/(sizek*sizel)))+
    ggplot2::scale_fill_distiller("Link density",type="seq",direction = 1)+
    ggplot2::scale_alpha("Link density")+
    ggplot2::ggtitle(paste0(toupper(x@model@name)," model with : ",max(x@cl)," clusters."))+
    ggplot2::scale_x_continuous("",breaks=cumsum(x@obs_stats$counts),labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),minor_breaks = NULL)+
    ggplot2::scale_y_continuous("",breaks=cumsum(x@obs_stats$counts),labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),minor_breaks = NULL)+
    ggplot2::coord_fixed()+ggplot2::theme_bw()
  
}

mat_blocks = function(x){
  K = length(x@obs_stats$counts)
  D = dim(x@obs_stats$x_counts)[2]
  gg=data.frame(kc=rep(cumsum(x@obs_stats$counts),D),
                lc=rep(1:D,each=K),
                sizek = rep(x@obs_stats$counts,D),
                sizel = rep(1,K*D), 
                count=as.vector(x@obs_stats$x_counts))
  ggplot2::ggplot(gg)+ggplot2::geom_tile(ggplot2::aes(y=kc-sizek/2,x=lc-sizel/2,height=sizek,width=sizel,fill=count/sizek,alpha=count/sizek))+
    ggplot2::scale_fill_distiller("E[X]",type="seq",direction = 1)+
    ggplot2::scale_alpha("E[X]")+
    ggplot2::ggtitle(paste0("MM Model with : ",max(x@cl)," clusters."))+
    ggplot2::scale_x_continuous("Features",breaks=1:D,labels=rep("",D),minor_breaks = NULL)+
    ggplot2::scale_y_continuous("Clusters",breaks=cumsum(x@obs_stats$counts),labels = paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),minor_breaks = NULL)+
    ggplot2::theme_bw()      
}
   
dendo = function(x){
  ggtree = x@ggtree
  ggplot2::ggplot()+ggplot2::geom_segment(data=ggtree[ggtree$node %in% ggtree$tree,],ggplot2::aes(x=xmin,y=H,xend=xmax,yend=H))+
    ggplot2::geom_segment(data=ggtree[-1,],ggplot2::aes(x=x,y=H,xend=x,yend=Hend))+
    ggplot2::scale_x_continuous("",breaks=c())+
    ggplot2::ylab(expression(paste("-log(",alpha,")")))+
    ggplot2::ggtitle(paste0(toupper(x@model@name)," ",length(x@obs_stats$counts)," clusters, dendogram"))+
    ggplot2::theme_bw()
}


lapath = function(x){
  gg = data.frame(k=sapply(x@path,function(p){length(p$counts)}),logalpha=sapply(x@path,function(p){p$logalpha}))
  gg = rbind(gg,data.frame(k=length(x@obs_stats$counts),logalpha=x@logalpha)) 
  ggplot2::ggplot(data=gg)+ggplot2::geom_line(ggplot2::aes(x=k,y=-logalpha))+
  ggplot2::geom_point(ggplot2::aes(x=k,y=-logalpha))+
  ggplot2::ylab(expression(paste("-log(",alpha,")")))+
  ggplot2::ggtitle(paste0(toupper(x@model@name)," ",length(x@obs_stats$counts)," clusters"))+
  ggplot2::theme_bw()
}


pprint =function(x,M,l){
  K = length(x@obs_stats$counts)
  na = colnames(M)
  D=Matrix::rowSums(M)
  for (k in 1:K){
    print("###########################")
    print(paste0("#########",k, "#########"))
    print("###########################")
    ii=which(x@cl==k)
    topk=order(D[ii],decreasing = TRUE)[1:l]
    print(na[ii[topk]])
  }  
}

spy = function(x){
  ij=Matrix::which(x!=0,arr.ind = TRUE)
  gg=data.frame(i=ij[,1],j=ij[,2],v=x[ij])
  ggplot2::ggplot(gg)+ggplot2::geom_point(aes(y=-i,x=j,size=v))+
    ggplot2::scale_x_continuous("",c())+
    ggplot2::scale_y_continuous("",c())+
    ggplot2::scale_size_area(max_size=1,guide='none')
}

#' @rdname print
#' @param x \code{\link{icl_path-class}} object to print
#' @export
setMethod(f = "print", 
          signature = signature("icl_path"),
          definition = function(x){
            print(paste0("ICL clustering with a ",toupper(x@model@name)," model, ",length(x@obs_stats$counts), " clusters and an icl of ", round(x@icl),"."))
          })

#' @rdname show
#' @param object \code{\link{icl_path-class}} object to print
#' @export
setMethod(f = "show", 
          signature = signature("icl_path"),
          definition = function(object){
            print(object)
          })
