#' @include models_classes.R fit_classes.R dcsbm.R
NULL



#' @title Degree Corrected Stochastic Block Model for bipartite graph class
#' 
#' @description An S4 class to represent a degree corrected stochastic block model for co_clustering of bipartite graph, extends \code{\link{icl_model-class}} class.
#' Such model can be used to cluster graph vertex, and model a bipartite graph adjacency matrix \eqn{X} with the following generative model :  
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i^r  \sim \mathcal{M}(1,\pi^r)}
#' \deqn{ Z_j^c  \sim \mathcal{M}(1,\pi^c)}
#' \deqn{ \theta_{kl} \sim Exponential(p)}
#' \deqn{ \gamma_i^r\sim \mathcal{U}(S_k)}
#' \deqn{ \gamma_i^c\sim \mathcal{U}(S_l)}
#' \deqn{ X_{ij}|Z_{ik}^cZ_{jl}^r=1 \sim \mathcal{P}(\gamma_i^r\theta_{kl}\gamma_j^c)}
#' The individuals parameters \eqn{\gamma_i^r,\gamma_j^c} allow to take into account the node degree heterogeneity. 
#' These parameters have uniform priors over simplex \eqn{S_k}. This class mainly store the prior parameters value \eqn{\alpha} of this generative model in the following slots (the prior parameter \eqn{p} is estimated from the data as the global average probability of connection between two nodes):
#' @slot alpha Dirichlet parameters for the prior over clusters proportions (default to 1)
#' @slot p Exponential prior parameter (default to Nan, in this case p will be estimated from data as the average intensities of X) 
#' @examples 
#' new("co_dcsbm")
#' new("co_dcsbm", p = 0.1)
#' @export
setClass("co_dcsbm",
         representation = list(p="numeric"),
         contains = "icl_model",
         prototype(name="co_dcsbm",alpha=1,p=NaN))



#' @title Degree corrected stochastic block model for bipartite graph fit results class
#' 
#' @description An S4 class to represent a fit of a degree corrected stochastic block model for co_clustering, extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{co_dcsbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot Krow number of extracted row clusters
#' @slot Kcol number of extracted column clusters
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item din: numeric vector of size K which store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K which store the sums of out-degrees for each clusters 
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters 
#' \item co_x_counts: matrix of size Krow*Kcol with the number of links between each pair of row and column cluster 
#' }
#' @slot clrow a numeric vector with row cluster indexes
#' @slot clcol a numeric vector with column cluster indexes
#' @slot Nrow number of rows
#' @slot Ncol number of columns
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("co_dcsbm_fit",slots = list(model="co_dcsbm",clrow="numeric",clcol="numeric",Krow="numeric",Kcol="numeric",Nrow="numeric",Ncol="numeric"),contains="icl_fit")




#' @title Degree corrected stochastic block model for bipartite graph hierarchical fit results class
#' 
#' 
#' @description An S4 class to represent a fit of a degree corrected stochastic block model for co_clustering, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{co_dcsbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot Krow number of extracted row clusters
#' @slot Kcol number of extracted column clusters
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item din: numeric vector of size K which store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K which store the sums of out-degrees for each clusters 
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters 
#' \item co_x_counts: matrix of size Krow*Kcol with the number of links between each pair of row and column cluster 
#' }
#' @slot clrow a numeric vector with row cluster indexes
#' @slot clcol a numeric vector with column cluster indexes
#' @slot Nrow number of rows
#' @slot Ncol number of columns
#' @slot path a list of size K-1 with each part of the path described by:
#' \itemize{
#' \item icl1: icl value reach with this solution for alpha=1 
#' \item logalpha: log(alpha) value were this solution is better than its parent
#' \item K: number of clusters
#' \item cl: vector of cluster indexes
#' \item k,l: index of the cluster that were merged at this step
#' \item merge_mat: lower triangular matrix of delta icl values 
#' \item obs_stats: a list with the elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item din: numeric vector of size K which store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K which store the sums of out-degrees for each clusters 
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters 
#' \item co_x_counts: matrix of size Krow*Kcol with the number of links between each pair of row and column cluster 
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father  
#' @slot ggtreerow data.frame with complete merge tree of row clusters for easy plotting with \code{ggplot2}
#' @slot ggtreecol data.frame with complete merge tree of column clusters for easy plotting with \code{ggplot2}
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("co_dcsbm_path",slots = list(ggtreerow="data.frame",ggtreecol="data.frame"),contains=c("icl_path","co_dcsbm_fit"))




#' @title method to cut a path solution to a desired number of cluster 
#' 
#' @description this method take a \code{\link{co_dcsbm_path-class}} and an integer K and return the solution from the path with K clusters 
#' @param x A an \code{\link{co_dcsbm_path-class}} solution 
#' @param K Desired number of cluster
#' @return an icl_path object with the desired number of cluster
#' @export
setMethod(f = "cut", 
          signature = signature("co_dcsbm_path"), 
          definition = function(x, K){
            if(K<x@K){
              i = which(sapply(x@path,function(p){p$K})==K)
              x@K = K
              x@logalpha=x@path[[i]]$logalpha
              x@icl = x@path[[i]]$icl
              x@cl = as.vector(x@path[[i]]$cl)
              for(st in names(x@path[[i]]$obs_stats)){
                x@obs_stats[st] = x@path[[i]]$obs_stats[st]
              }
              
              x@path=x@path[(i+1):length(x@path)]
              x=postprocess(x)
            }else{
              warning(paste0("This clustering has ",x@K," clusters and you requested ",K ," clusters. Please provide a value for K smaller than ",x@K,"."),call. = FALSE)
            }
            x
            
          })


#' @title plot a \code{\link{co_dcsbm_fit-class}}
#' 
#' 
#' @param x a \code{\link{co_dcsbm_fit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between row and column clusters
#' \item \code{'nodelink'}: plot a nodelink diagram of the bipartite graph summarizing connections between row and column clusters
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("co_dcsbm_fit","missing"),
          definition = function(x,type="blocks"){
            switch(type,blocks=co_blocks(x),nodelink=co_nodelink(x))
          })

#' @title plot a \code{\link{co_dcsbm_path-class}}
#' 
#' @param x a \code{\link{co_dcsbm_path-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between row and column clusters
#' \item \code{'nodelink'}: plot a nodelink diagram of the bipartite graph summarizing connections between row and column clusters
#' \item \code{'front'}: plot the extracted front ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the associated dendrograms one for the row clusters and one for the column clusters
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("co_dcsbm_path","missing"),
          definition = function(x,type='blocks'){
            switch(type,tree = {
              co_dendo(x)
            },
            path ={
              lapath(x)
            },
            front = {
              plot_front(x)
            },
            blocks ={
              methods::callNextMethod()
            },
            nodelink={
              methods::callNextMethod()
            })   
          })


#' @title Extract parameters from an \code{\link{co_dcsbm_fit-class}} object
#' 
#' @param object a \code{\link{co_dcsbm_fit-class}}
#' @return a list with the model parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pirows'}: row cluster proportions 
#' \item \code{'picols'}: row cluster proportions 
#' \item \code{'thetakl'}: between clusters connection probabilities (matrix of size Krow x Kcol),
#' \item \code{'gammarows'}: rows degree correction parameters (size Nrows),
#' \item \code{'gammacols'}: cols degree correction parameters (size Ncols),
#' }
#' @export 
setMethod(f = "coef", 
          signature = signature(object = "co_dcsbm_fit"),
          definition = function(object){
            sol=object
            pirows=(sol@obs_stats$rows_counts+sol@model@alpha-1)/sum(sol@obs_stats$rows_counts+sol@model@alpha-1)
            picols=(sol@obs_stats$cols_counts+sol@model@alpha-1)/sum(sol@obs_stats$cols_counts+sol@model@alpha-1)
            gammarows = sol@obs_stats_cst$drow/sol@obs_stats$dr[sol@clrow]
            gammacols = sol@obs_stats_cst$dcol/sol@obs_stats$dc[sol@clcol]
            thetakl=(sol@obs_stats$co_x_counts)/(t(t(sol@obs_stats$rows_counts))%*%sol@obs_stats$cols_counts+1/sol@model@p)
            list(pirows=pirows,picols=picols,thetakl=thetakl,gammarows=gammarows,gammacols=gammacols)
          })

setMethod(f = "preprocess", 
          signature = signature("co_dcsbm"), 
          definition = function(model,data){
            X=as.sparse(data)
            if(class(X)!="dgCMatrix"){
              stop(paste0("Unsupported data type :", class(X) ,"for co_dcsbm model."),call. = FALSE)
            }
            if(nrow(X)==ncol(X)){
              stop("Square matrix as input, try dcsbm model instead.",call. = FALSE)
            }
            ij=which(X>0,arr.ind = TRUE)
            if(!all(X[ij]==round(X[ij]))){
              stop("Only integer matrix allowed as input, non integer values found.",call. = FALSE)
            }
            if(length(model@alpha)>1){
              stop("Model prior misspecification, alpha must be of length 1.",call. = FALSE)
            }
            if(is.na(model@alpha)){
              stop("Model prior misspecification, alpha is NA.",call. = FALSE)
            }
            if(model@alpha<=0){
              stop("Model prior misspecification, alpha must be positive.",call. = FALSE)
            }
            if(length(model@p)>1){
              stop("Model prior misspecification, p must be of length 1.",call. = FALSE)
            }
            if(!is.nan(model@p) && model@p<=0){
              stop("Model prior misspecification, p must be positive.",call. = FALSE)
            }
            
            di=dim(X)
            N=sum(di)
            list(X=X,N=N,Nrows=di[1],Ncols=di[2])
          })

setMethod(f = "postprocess", 
          signature = signature("co_dcsbm_path"), 
          definition = function(path,data=NULL){
            sol = path
            if(!is.null(data)){
              sol@Nrow = data$Nrows
              sol@Ncol = data$Ncols
            }
            clusters_type = apply(table(sol@cl,c(rep(1,sol@Nrow),rep(2,sol@Ncol))),1,which.max)
            clust_rows = which(clusters_type==1)
            clust_cols = which(clusters_type==2)
            icol = (sol@Nrow+1):length(sol@cl)
            irow = 1:sol@Nrow

            sol@clrow = as.numeric(factor(sol@cl[irow],levels=clust_rows))
            sol@Krow = max(sol@clrow,na.rm=TRUE)
            sol@clcol = as.numeric(factor(sol@cl[icol],levels=clust_cols))
            sol@Kcol = max(sol@clcol,na.rm=TRUE)
            sol@obs_stats$co_x_counts=sol@obs_stats$x_counts[clust_rows,clust_cols]
            sol@obs_stats$dr=sol@obs_stats$dr[clust_rows]
            sol@obs_stats$dc=sol@obs_stats$dc[clust_cols]
            sol@obs_stats=sol@obs_stats[names(sol@obs_stats)!="x_counts"]
            sol@obs_stats$rows_counts = sol@obs_stats$counts[clust_rows]
            sol@obs_stats$cols_counts = sol@obs_stats$counts[clust_cols]            
            if(!is.null(data)){
              
              tree=sol@ggtree
              root=tree$node[tree$tree==0]
              tree=tree[tree$node!=root,]
              tree$tree[tree$tree==root]=0
              
              xcol=tree[tree$node %in% clust_cols,]$x
              coltree=tree[tree$x>=min(xcol)&tree$x<=max(xcol),]
              xrow=tree[tree$node %in% clust_rows,]$x
              rowtree=tree[tree$x>=min(xrow)&tree$x<=max(xrow),]
              
              
              #cat('-- post-processing --')
              # tree= sol@ggtree[order(sol@ggtree$H,sol@ggtree$node),]
              # coltree = tree[tree$node %in% clust_cols,]
              # coltree$x = seq(1,-1,length.out = length(clust_cols))
              # leafs   =  coltree
              # fathers = unique(leafs$tree)
              # while(length(fathers)>1){
              #   leafs = tree[tree$node %in% fathers,]
              #   coltree = rbind(coltree,leafs)
              #   fathers = setdiff(unique(leafs$tree),coltree$node)
              #   fathers=fathers[fathers!=0]
              # }
              # coltree=coltree[order(coltree$H),]
              # noleaves = (which(coltree$H>0)[1]):nrow(coltree)
              # for (nl in noleaves){
              #   no =coltree$node[nl]
              #   xch = coltree$x[coltree$tree==no]
              #   coltree[nl,"x"]=mean(xch)
              #   coltree[nl,"xmin"]=min(xch)
              #   coltree[nl,"xmax"]=max(xch)
              # }
              # 
              # 
              # rowtree = tree[tree$node %in% clust_rows,]
              # rowtree$x = seq(1,-1,length.out = length(clust_rows))
              # leafs   =  rowtree
              # fathers = unique(leafs$tree)
              # while(length(fathers)>1){
              #   leafs = tree[tree$node %in% fathers,]
              #   rowtree = rbind(rowtree,leafs)
              #   fathers = setdiff(unique(leafs$tree),rowtree$node)
              #   fathers=fathers[fathers!=0]
              # }
              # rowtree=rowtree[order(rowtree$H),]
              # noleaves = (which(rowtree$H>0)[1]):nrow(rowtree)
              # for (nl in noleaves){
              #   no =rowtree$node[nl]
              #   xch = rowtree$x[rowtree$tree==no]
              #   rowtree[nl,"x"]=mean(xch)
              #   rowtree[nl,"xmin"]=min(xch)
              #   rowtree[nl,"xmax"]=max(xch)
              # }
              
              sol@ggtreecol = coltree
              sol@ggtreerow = rowtree 
            }

            sol
          })




 setMethod(f = "sample_cl", 
           signature = signature("co_dcsbm","list","numeric"), 
           definition = function(model,data,K){
             c(sample(1:floor(K/2),data$Nrows,replace = TRUE),sample((floor(K/2)+1):K,data$Ncols,replace = TRUE))
           })

reorder_codcsbm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$dr = obs_stats$dr[or]
  obs_stats$dc = obs_stats$dc[or]
  obs_stats$x_counts = obs_stats$x_counts[or,or]
  obs_stats
}
setMethod(f = "reorder", 
          signature = signature("co_dcsbm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_codcsbm(obs_stats,order)
          })


setMethod(f = "seed", 
          signature = signature("co_dcsbm","list","numeric"), 
          definition = function(model,data, K){
            kmrow = stats::kmeans(data$X,floor(K/2));
            kmcol = stats::kmeans(t(data$X),floor(K/2));
            c(kmrow$cluster,kmcol$cluster+max(kmrow$cluster))
            #spectral(data$X,K)
          })

