setGeneric("plot", function(sol) standardGeneric("plot")) 

setMethod(f = "plot", 
          signature = signature("sbm_fit"),
          definition = function(sol){
            gg=data.frame(kc=rep(cumsum(sol@counts),each=10),
                          lc=rep(cumsum(sol@counts),10),
                          sizek = rep(sol@counts,each=10),
                          sizel = rep(sol@counts,10), 
                          count=as.vector(sol@x_counts))
            ggplot(gg)+geom_tile(aes(x=kc-sizek/2,y=lc-sizel/2,width=sizek,height=sizel,fill=count/(sizek*sizel),alpha=count/(sizek*sizel)))+
              scale_fill_distiller("Link density",type="seq",direction = 1)+
              scale_alpha("Link density")+
              ggtitle(paste0("SBM Model with : ",max(sol@cl)," clusters."))+
              scale_x_continuous("",breaks=cumsum(sol@counts),labels = paste0(round(100*sol@counts/sum(sol@counts)),"%"))+
              scale_y_continuous("",breaks=cumsum(sol@counts),labels = paste0(round(100*sol@counts/sum(sol@counts)),"%"))+
              coord_fixed()+theme_bw()
            
          });


setMethod(f = "plot", 
          signature = signature("mm_fit"),
          definition = function(sol){
            K = length(sol@counts)
            D = dim(sol@x_counts)[2]
            gg=data.frame(kc=rep(cumsum(sol@counts),D),
                          lc=rep(1:D,each=K),
                          sizek = rep(sol@counts,D),
                          sizel = rep(1,K*D), 
                          count=as.vector(sol@x_counts))
            ggplot(gg)+geom_tile(aes(y=kc-sizek/2,x=lc-sizel/2,height=sizek,width=sizel,fill=count/sizek,alpha=count/sizek))+
              scale_fill_distiller("E[X]",type="seq",direction = 1)+
              scale_alpha("E[X]")+
              ggtitle(paste0("MM Model with : ",max(sol@cl)," clusters."))+
              scale_x_continuous("Features",breaks=1:D,labels=rep("",D))+
              scale_y_continuous("Individuals",breaks=cumsum(sol@counts),labels = paste0(round(100*sol@counts/sum(sol@counts)),"%"))+
              theme_bw()         
          });
