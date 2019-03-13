#' @importFrom graphics plot
#' @include models_classes.R fit_classes.R
#' @importFrom ggraph guide_edge_colourbar
#' @export 
ggraph::guide_edge_colourbar
NULL


#' Plot a clustering results
#' @param sol \code{\link{sbm_fit}} obejct to be ploted
#' @return a ggplot2 graphics which summarize the results 
#' @export
setMethod(f = "plot", 
          signature = signature("sbm_fit","missing"),
          definition = function(x,y,...,type="blocks"){
            switch(type,
                   blocks={
                     K = length(x@obs_stats$counts)
                     gg=data.frame(kc=rep(cumsum(x@obs_stats$counts),each=K),
                                   lc=rep(cumsum(x@obs_stats$counts),K),
                                   sizek = rep(x@obs_stats$counts,each=K),
                                   sizel = rep(x@obs_stats$counts,K), 
                                   count=as.vector(x@obs_stats$x_counts))
                     ggplot2::ggplot(gg[gg$count>0,])+ggplot2::geom_tile(ggplot2::aes(x=kc-sizek/2,y=lc-sizel/2,width=sizek,height=sizel,fill=count/(sizek*sizel),alpha=count/(sizek*sizel)))+
                       ggplot2::scale_fill_distiller("Link density",type="seq",direction = 1)+
                       ggplot2::scale_alpha("Link density")+
                       ggplot2::ggtitle(paste0("SBM Model with : ",max(x@cl)," clusters."))+
                       ggplot2::scale_x_continuous("",breaks=cumsum(x@obs_stats$counts),labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),minor_breaks = NULL)+
                       ggplot2::scale_y_continuous("",breaks=cumsum(x@obs_stats$counts),labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),minor_breaks = NULL)+
                       ggplot2::coord_fixed()+ggplot2::theme_bw()
                   },nodelink=graph(x,'linear'))
          });


#' Plot a clustering results
#' @param sol \code{\link{sbm_fit}} obejct to be ploted
#' @return a ggplot2 graphics which summarize the results 
#' @export
setMethod(f = "plot", 
          signature = signature("dcsbm_fit","missing"),
          definition = function(x,y,...,type="blocks"){
            switch(type,
                   blocks={
            K = length(x@obs_stats$counts)
            gg=data.frame(kc=rep(cumsum(x@obs_stats$counts),each=K),
                          lc=rep(cumsum(x@obs_stats$counts),K),
                          sizek = rep(x@obs_stats$counts,each=K),
                          sizel = rep(x@obs_stats$counts,K), 
                          count=as.vector(x@obs_stats$x_counts))
            ggplot2::ggplot(gg[gg$count>0,])+ggplot2::geom_tile(ggplot2::aes(x=kc-sizek/2,y=lc-sizel/2,width=sizek,height=sizel,fill=count/(sizek*sizel),alpha=count/(sizek*sizel)))+
              ggplot2::scale_fill_distiller("Link density",type="seq",direction = 1)+
              ggplot2::scale_alpha("Link density")+
              ggplot2::ggtitle(paste0("DC-SBM Model with : ",max(x@cl)," clusters."))+
              ggplot2::scale_x_continuous("",breaks=cumsum(x@obs_stats$counts),labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),minor_breaks = NULL)+
              ggplot2::scale_y_continuous("",breaks=cumsum(x@obs_stats$counts),labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),minor_breaks = NULL)+
              ggplot2::coord_fixed()+ggplot2::theme_bw()
                   },nodelink=graph(x,'linear'))
          });

#' Plot a clustering results
#' @param sol \code{\link{mm_fit}} obejct to be ploted
#' @return a ggplot2 graphics which summarize the results
#' @export
setMethod(f = "plot", 
          signature = signature("mm_fit","missing"),
          definition = function(x,y,...,type='blocks'){
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
          });


#' Plot a clustering results
#' @param sol \code{\link{mm_fit}} obejct to be ploted
#' @return a ggplot2 graphics which summarize the results
#' @export
setMethod(f = "plot", 
          signature = signature("dcsbm_path","missing"),
          definition = function(x,y,...,type='tree'){
            switch(type,tree = {
              ggtree = x@ggtree
              ggplot2::ggplot()+ggplot2::geom_segment(data=ggtree[ggtree$node %in% ggtree$tree,],ggplot2::aes(x=xmin,y=H,xend=xmax,yend=H))+
                ggplot2::geom_segment(data=ggtree[-1,],ggplot2::aes(x=x,y=H,xend=x,yend=Hend))+
                ggplot2::scale_x_continuous("",breaks=c())+
                ggplot2::ylab(expression(paste("-log(",alpha,")")))+
                ggplot2::ggtitle(paste0(toupper(x@model@name)," ",length(x@obs_stats$counts)," clusters, dendogram"))+
                ggplot2::theme_bw()
            },
              path ={
                gg = data.frame(k=sapply(x@path,function(p){length(p$counts)}),logalpha=sapply(x@path,function(p){p$logalpha}))
                gg = rbind(gg,data.frame(k=length(x@obs_stats$counts),logalpha=x@logalpha)) 
                ggplot2::ggplot(data=gg)+ggplot2::geom_line(ggplot2::aes(x=k,y=-logalpha))+
                ggplot2::geom_point(ggplot2::aes(x=k,y=-logalpha))+
                ggplot2::ylab(expression(paste("-log(",alpha,")")))+
                ggplot2::ggtitle(paste0(toupper(x@model@name)," ",length(x@obs_stats$counts)," clusters"))+
                ggplot2::theme_bw()
              },
              blocks ={
                callNextMethod()
              },
              nodelink={
                callNextMethod()
              })   
          });


#' Plot a clustering results
#' @param sol \code{\link{mm_fit}} obejct to be ploted
#' @return a ggplot2 graphics which summarize the results
#' @export
setMethod(f = "plot", 
          signature = signature("sbm_path","missing"),
          definition = function(x,y,...,type='tree'){
            switch(type,tree = {
              ggtree = x@ggtree
              ggplot2::ggplot()+ggplot2::geom_segment(data=ggtree[ggtree$node %in% ggtree$tree,],ggplot2::aes(x=xmin,y=H,xend=xmax,yend=H))+
                ggplot2::geom_segment(data=ggtree[-1,],ggplot2::aes(x=x,y=H,xend=x,yend=Hend))+
                ggplot2::scale_x_continuous("",breaks=c())+
                ggplot2::ylab(expression(paste("-log(",alpha,")")))+
                ggplot2::ggtitle(paste0(toupper(x@model@name)," ",length(x@obs_stats$counts)," clusters, dendogram"))+
                ggplot2::theme_bw()
            },
            path ={
              gg = data.frame(k=sapply(x@path,function(p){length(p$counts)}),logalpha=sapply(x@path,function(p){p$logalpha}))
              gg = rbind(gg,data.frame(k=length(x@obs_stats$counts),logalpha=x@logalpha)) 
              ggplot2::ggplot(data=gg)+ggplot2::geom_line(ggplot2::aes(x=k,y=-logalpha))+
                ggplot2::geom_point(ggplot2::aes(x=k,y=-logalpha))+
                ggplot2::ylab(expression(paste("-log(",alpha,")")))+
                ggplot2::ggtitle(paste0(toupper(x@model@name)," ",length(x@obs_stats$counts)," clusters"))+
                ggplot2::theme_bw()
            },
            blocks ={
              callNextMethod()
            },
            nodelink={
              callNextMethod()
            })   
          });

#' Plot a clustering results
#' @param sol \code{\link{mm_fit}} obejct to be ploted
#' @return a ggplot2 graphics which summarize the results
#' @export
setMethod(f = "plot", 
          signature = signature("mm_path","missing"),
          definition = function(x,y,...,type='tree'){
            switch(type,tree = {
              ggtree = x@ggtree
              ggplot2::ggplot()+ggplot2::geom_segment(data=ggtree[ggtree$node %in% ggtree$tree,],ggplot2::aes(x=xmin,y=H,xend=xmax,yend=H))+
                ggplot2::geom_segment(data=ggtree[-1,],ggplot2::aes(x=x,y=H,xend=x,yend=Hend))+
                ggplot2::scale_x_continuous("",breaks=c())+
                ggplot2::ylab(expression(paste("-log(",alpha,")")))+
                ggplot2::ggtitle(paste0(toupper(x@model@name)," ",length(x@obs_stats$counts)," clusters, dendogram"))+
                ggplot2::theme_bw()
            },
            path ={
              gg = data.frame(k=sapply(x@path,function(p){length(p$counts)}),logalpha=sapply(x@path,function(p){p$logalpha}))
              gg = rbind(gg,data.frame(k=length(x@obs_stats$counts),logalpha=x@logalpha)) 
              ggplot2::ggplot(data=gg)+ggplot2::geom_line(ggplot2::aes(x=k,y=-logalpha))+
                ggplot2::geom_point(ggplot2::aes(x=k,y=-logalpha))+
                ggplot2::ylab(expression(paste("-log(",alpha,")")))+
                ggplot2::ggtitle(paste0(toupper(x@model@name)," ",length(x@obs_stats$counts)," clusters"))+
                ggplot2::theme_bw()
            },
            blocks ={
              callNextMethod()
            },
            nodelink={
              callNextMethod()
            })   
          });


graph = function(sol,layout='linear'){
  ig = igraph::graph_from_adjacency_matrix(sol@obs_stats$x_counts/(sol@obs_stats$counts%*%t(sol@obs_stats$counts)),weighted = TRUE)
  igraph::V(ig)$weight = sol@obs_stats$counts
  igraph::V(ig)$ilinks = diag(sol@obs_stats$x_counts)/(sol@obs_stats$counts^2)
  ggraph::ggraph(ig, layout = layout)+
    ggraph::geom_node_point(ggplot2::aes(size=weight,color=ilinks))+
    ggraph::geom_edge_arc(ggplot2::aes(edge_width=weight,edge_alpha=weight,color=weight),curvature = 0.3,arrow = grid::arrow(length = grid::unit(3, 'mm')),start_cap=ggraph::circle(3,'mm'))+
    ggraph::scale_edge_color_distiller("Links density",palette="YlGn",direction = 1)+
    ggplot2::scale_color_distiller(palette="YlGn",guide = "none",direction = 1)+
    ggraph::scale_edge_alpha("Links density",range=c(0,1))+
    ggraph::scale_edge_width("Links density",range=c(0,5))+
    ggplot2::scale_size_area("Clusters size",max_size = 10)+
    ggplot2::scale_y_continuous("",c())+
    ggplot2::scale_x_continuous("",c())+
    ggplot2::ggtitle(paste0(toupper(sol@model@name)," ",length(sol@obs_stats$counts)," clusters"))+
    ggplot2::theme_bw()
}

