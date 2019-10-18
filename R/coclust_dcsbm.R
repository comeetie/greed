#' @include models_classes.R fit_classes.R dcsbm.R
NULL



#' @title Co-clustering with a degree correted stochastick block model class
#' 
#' @description An S4 class to represent a degree corrected stochastick block model for co_clustering, extends \code{\link{icl_model-class}} class.
#' @slot alpha dirichlet parameters for the prior over clusters proportions (default to 1)
#' @examples 
#' new("co_dcsbm")
#' @export
setClass("co_dcsbm",
         representation = list(),
         contains = "icl_model",
         prototype(name="co_dcsbm",alpha=1))



#' @title Co-clustering with a degree correted stochastick block model fit results class
#' 
#' @description An S4 class to represent a fit of a degree corrected stochastick block model for co_clustering, extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{co_dcsbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot Krow number of extracted row clusters
#' @slot Kcol number of extracted column clusters
#' @slot cl a numeric vector with row and clolumns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item din: numeric vector of size K wich store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K wich store the sums of out-degrees for each clusters 
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters 
#' \item co_x_counts: matrix of size Krow*Kcol with the number of links between each pair of row and column cluster 
#' }
#' @slot clrow a numeric vector with row cluster indexes
#' @slot clcol a numeric vector with column cluster indexes
#' @slot Nrow number of rows
#' @slot Ncol number of columns
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history infromation (details depends on the training procedure)
#' @export 
setClass("co_dcsbm_fit",slots = list(model="co_dcsbm",clrow="numeric",clcol="numeric",Krow="numeric",Kcol="numeric",Nrow="numeric",Ncol="numeric"),contains="icl_fit")




#' @title Co-clustering with a degree correted stochastick block model path extraction results class
#' 
#' 
#' @description An S4 class to represent a fit of a degree corrected stochastick block model for co_clustering, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{co_dcsbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot Krow number of extracted row clusters
#' @slot Kcol number of extracted column clusters
#' @slot cl a numeric vector with row and clolumns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item din: numeric vector of size K wich store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K wich store the sums of out-degrees for each clusters 
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
#' \item obs_stats: a list with the same elements
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy ploting with gggplot
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father  
#' @slot ggtreerow data.frame with complete merge tree of row clusters for easy ploting with gggplot
#' @slot ggtreecol data.frame with complete merge tree od column clusters for easy ploting with gggplot
#' @slot train_hist  data.frame with training history infromation (details depends on the training procedure)
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
            i = which(sapply(x@path,function(p){p$K})==K)
            x@K = K
            x@logalpha=x@path[[i]]$logalpha
            x@icl = x@path[[i]]$icl
            x@cl = as.vector(x@path[[i]]$cl)
            for(st in names(x@path[[i]]$obs_stats)){
              x@obs_stats[st] = x@path[[i]]$obs_stats[st]
            }
            
            x@path=x@path[(i+1):length(x@path)]
            postprocess(x)
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
#' \item \code{'path'}: plot the veolution of ICL with repsect to K
#' \item \code{'tree'}: plot the associated dendograms one for the row clustrers and one for the column clusters
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

setMethod(f = "preprocess", 
          signature = signature("co_dcsbm"), 
          definition = function(model,data,K){
            X=as.sparse(data)
            if(class(X)!="dgCMatrix"){
              stop(paste0("Unsupported data type :", class(X) ,"for co_dcsbm model."))
            }
            if(nrow(X)==ncol(X)){
              stop("Square matrix as input, try dcsbm model instead.")
            }
            ij=which(X>0,arr.ind = TRUE)
            if(!all(X[ij]==round(X[ij]))){
              stop("Only integer matrix allowed as input, non integer values found.")
            }

            di=dim(X)
            N=sum(di)
            moves=matrix(0,K,K);
            moves[1:floor(K/2),1:floor(K/2)]=1
            moves[(floor(K/2)+1):K,(floor(K/2)+1):K]=1
            X =  sparseMatrix(i=c(ij[,1],ij[,2]+di[1]),j=c(ij[,2]+di[1],ij[,1]),x = c(X[ij],X[ij]))
            list(X=X,N=nrow(X),Nrows=di[1],Ncols=di[2],moves=as.sparse(moves))
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
            # if(!(max(clust_cols)<min(clust_rows) | max(clust_rows)<min(clust_cols))){
            #   message("Co clustering failed bi partite structure not found")
            #   return(sol)
            # }
            icol = (sol@Nrow+1):length(sol@cl)
            irow = 1:sol@Nrow
            # row_problems = !sol@cl[irow] %in% which(clusters_type==1)
            # col_problems = !sol@cl[icol] %in% which(clusters_type==2)
            # if((sum(row_problems)>0 | sum(col_problems)>0) & !is.null(data)){
            # probas=greed:::post_probs(sol@model,data,sol@cl)
            # if(sum(row_problems)>0){
            #   row_probas = probas[irow,]
            #   sol@cl[which(row_problems)]=clust_rows[apply(matrix(row_probas[which(row_problems),clusters_type==1],nrow = sum(row_problems)),1,which.max)]
            # }
            # if(sum(col_problems)>0){
            #   col_probas = probas[icol,]
            #   sol@cl[which(col_problems)+data$Nrows]=clust_cols[apply(matrix(col_probas[which(col_problems),clusters_type==2],nrow = sum(col_problems)),1,which.max)]
            # }
            # }
            sol@clrow = as.numeric(factor(sol@cl[irow],levels=clust_rows))
            sol@Krow = max(sol@clrow,na.rm=TRUE)
            sol@clcol = as.numeric(factor(sol@cl[icol],levels=clust_cols))
            sol@Kcol = max(sol@clcol,na.rm=TRUE)
            sol@obs_stats$co_x_counts=sol@obs_stats$x_counts[clust_rows,clust_cols]
            
            if(!is.null(data)){
              tree= sol@ggtree[order(sol@ggtree$H,sol@ggtree$node),]
              coltree = tree[tree$node %in% clust_cols,]
              coltree$x = seq(1,-1,length.out = length(clust_cols))
              leafs   =  coltree
              fathers = unique(leafs$tree)
              while(length(fathers)>0){
                leafs = tree[tree$node %in% fathers,]
                coltree = rbind(coltree,leafs)
                fathers = setdiff(unique(leafs$tree),coltree$node)
                fathers=fathers[fathers!=0]
              }
              coltree=coltree[order(coltree$H),]
              noleaves = (which(coltree$H>0)[1]):nrow(coltree)
              for (nl in noleaves){
                no =coltree$node[nl]
                xch = coltree$x[coltree$tree==no]
                coltree[nl,"x"]=mean(xch)
                coltree[nl,"xmin"]=min(xch)
                coltree[nl,"xmax"]=max(xch)
              }
              
              
              rowtree = tree[tree$node %in% clust_rows,]
              rowtree$x = seq(1,-1,length.out = length(clust_rows))
              leafs   =  rowtree
              fathers = unique(leafs$tree)
              while(length(fathers)>0){
                leafs = tree[tree$node %in% fathers,]
                rowtree = rbind(rowtree,leafs)
                fathers = setdiff(unique(leafs$tree),rowtree$node)
                fathers=fathers[fathers!=0]
              }
              rowtree=rowtree[order(rowtree$H),]
              noleaves = (which(rowtree$H>0)[1]):nrow(rowtree)
              for (nl in noleaves){
                no =rowtree$node[nl]
                xch = rowtree$x[rowtree$tree==no]
                rowtree[nl,"x"]=mean(xch)
                rowtree[nl,"xmin"]=min(xch)
                rowtree[nl,"xmax"]=max(xch)
              }
              
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


setMethod(f = "reorder", 
          signature = signature("co_dcsbm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_dcsbm(obs_stats,order)
          })


setMethod(f = "seed", 
          signature = signature("co_dcsbm","list","numeric"), 
          definition = function(model,data, K){
            kmrow = stats::kmeans(data$X[1:data$Nrows,(data$Nrows+1):data$N],floor(K/2));
            kmcol = stats::kmeans(t(data$X[1:data$Nrows,(data$Nrows+1):data$N]),floor(K/2));
            c(kmrow$cluster,kmcol$cluster+max(kmrow$cluster))
            #spectral(data$X,K)
          })
