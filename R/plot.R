#' @importFrom graphics plot
#' @include models_classes.R fit_classes.R
NULL


#' Plot a clustering results
#' @param sol \code{\link{sbm_fit}} obejct to be ploted
#' @return a ggplot2 graphics which summarize the results 
#' @export
setMethod(f = "plot", 
          signature = signature("sbm_fit","missing"),
          definition = function(x){
            K = length(x@obs_stats$counts)
            gg=data.frame(kc=rep(cumsum(x@obs_stats$counts),each=K),
                          lc=rep(cumsum(x@obs_stats$counts),K),
                          sizek = rep(x@obs_stats$counts,each=K),
                          sizel = rep(x@obs_stats$counts,K), 
                          count=as.vector(x@obs_stats$x_counts))
          ggplot2::ggplot(gg[gg$count>0,])+ggplot2::geom_tile(ggplot2::aes(x=kc-sizek/2,y=lc-sizel/2,width=sizek,height=sizel,fill=log10(count/(sizek*sizel))))+
              ggplot2::scale_fill_distiller("Link density",type="seq",direction = 1)+
              #ggplot2::scale_alpha("Link density")+
              ggplot2::ggtitle(paste0("SBM Model with : ",max(x@cl)," clusters."))+
              ggplot2::scale_x_continuous("",breaks=cumsum(x@obs_stats$counts),labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),minor_breaks = NULL)+
              ggplot2::scale_y_continuous("",breaks=cumsum(x@obs_stats$counts),labels = ifelse(x@obs_stats$counts/sum(x@obs_stats$counts)>0.05,paste0(round(100*x@obs_stats$counts/sum(x@obs_stats$counts)),"%"),""),minor_breaks = NULL)+
              ggplot2::coord_fixed()+ggplot2::theme_bw()
            
          });


#' Plot a clustering results
#' @param sol \code{\link{sbm_fit}} obejct to be ploted
#' @return a ggplot2 graphics which summarize the results 
#' @export
setMethod(f = "plot", 
          signature = signature("dcsbm_fit","missing"),
          definition = function(x){
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
            
          });

#' Plot a clustering results
#' @param sol \code{\link{mm_fit}} obejct to be ploted
#' @return a ggplot2 graphics which summarize the results
#' @export
setMethod(f = "plot", 
          signature = signature("mm_fit"),
          definition = function(x){
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
          definition = function(x,y,type='tree'){
            print(type)   
          });
