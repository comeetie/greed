#' @include models_classes.R fit_classes.R
NULL

#' Plot a clustering results
#' @param sol \code{\link{icl_fit}} obejct to be ploted
#' @return a ggplot2 graphics which summarize the results 
#' @export
setGeneric("clustplot", function(sol) standardGeneric("clustplot")) 

#' Plot a clustering results
#' @param sol \code{\link{sbm_fit}} obejct to be ploted
#' @return a ggplot2 graphics which summarize the results 
#' @export
setMethod(f = "clustplot", 
          signature = signature("sbm_fit"),
          definition = function(sol){
            K = length(sol@counts)
            gg=data.frame(kc=rep(cumsum(sol@counts),each=K),
                          lc=rep(cumsum(sol@counts),K),
                          sizek = rep(sol@counts,each=K),
                          sizel = rep(sol@counts,K), 
                          count=as.vector(sol@x_counts))
            ggplot2::ggplot(gg)+ggplot2::geom_tile(ggplot2::aes(x=kc-sizek/2,y=lc-sizel/2,width=sizek,height=sizel,fill=count/(sizek*sizel),alpha=count/(sizek*sizel)))+
              ggplot2::scale_fill_distiller("Link density",type="seq",direction = 1)+
              ggplot2::scale_alpha("Link density")+
              ggplot2::ggtitle(paste0("SBM Model with : ",max(sol@cl)," clusters."))+
              ggplot2::scale_x_continuous("",breaks=cumsum(sol@counts),labels = paste0(round(100*sol@counts/sum(sol@counts)),"%"))+
              ggplot2::scale_y_continuous("",breaks=cumsum(sol@counts),labels = paste0(round(100*sol@counts/sum(sol@counts)),"%"))+
              ggplot2::coord_fixed()+ggplot2::theme_bw()
            
          });

#' Plot a clustering results
#' @param sol \code{\link{mm_fit}} obejct to be ploted
#' @return a ggplot2 graphics which summarize the results
#' @export
setMethod(f = "clustplot", 
          signature = signature("mm_fit"),
          definition = function(sol){
            K = length(sol@counts)
            D = dim(sol@x_counts)[2]
            gg=data.frame(kc=rep(cumsum(sol@counts),D),
                          lc=rep(1:D,each=K),
                          sizek = rep(sol@counts,D),
                          sizel = rep(1,K*D), 
                          count=as.vector(sol@x_counts))
            ggplot2::ggplot(gg)+ggplot2::geom_tile(ggplot2::aes(y=kc-sizek/2,x=lc-sizel/2,height=sizek,width=sizel,fill=count/sizek,alpha=count/sizek))+
              ggplot2::scale_fill_distiller("E[X]",type="seq",direction = 1)+
              ggplot2::scale_alpha("E[X]")+
              ggplot2::ggtitle(paste0("MM Model with : ",max(sol@cl)," clusters."))+
              ggplot2::scale_x_continuous("Features",breaks=1:D,labels=rep("",D))+
              ggplot2::scale_y_continuous("Individuals",breaks=cumsum(sol@counts),labels = paste0(round(100*sol@counts/sum(sol@counts)),"%"))+
              ggplot2::theme_bw()         
          });
