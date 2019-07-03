#' @include models_classes.R fit_classes.R dcsbm.R
NULL


#' @rdname models-classes
#' @title co_dcsbm
#' 
#' An S4 class to represent a stochastick block model that extends \code{icl_model} class.
#' \itemize{
#' \item slots : \code{name,alpha,a0,b0}
#' }
#' @slot a0 a numeric vector of length 1 which define the parameters of the beta prior over the edges (default to 1)
#' @slot b0 a numeric vector of length 1 which define the parameters of the beta prior over the non-edges (default to 1)
#' @examples 
#' new("co_dcsbm")
#' @export
setClass("co_dcsbm",
         representation = list(),
         contains = "icl_model",
         prototype(name="co_dcsbm",alpha=1))


#' @rdname fits-classes
#' @title co_dcsbm_fit
#' 
#' An S4 class to represent a fit of a stochastick block model that extend \code{icl_fit}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats,model}
#' }
#' @slot model an \code{\link{icl_model}} to store the model fitted
#' @export 
setClass("co_dcsbm_fit",slots = list(model="co_dcsbm",clrow="numeric",clcol="numeric",Krow="numeric",Kcol="numeric",Nrow="numeric",Ncol="numeric"),contains="icl_fit")



#' @rdname fits-classes
#' @title co_dcsbm_path
#' 
#' An S4 class to represent a hierachical path of solutions for a DC-SBM model that extend \code{dcsbm_fit-class} and \code{icl_path-class}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats, model, path, tree, ggtree, logalpha}
#' }
#' @export
setClass("co_dcsbm_path",slots = list(ggtreerow="data.frame",ggtreecol="data.frame"),contains=c("icl_path","co_dcsbm_fit"))




setMethod(f = "reorder", 
          signature = signature("co_dcsbm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_dcsbm(obs_stats,order)
          })


setMethod(f = "seed", 
          signature = signature("co_dcsbm","list","numeric"), 
          definition = function(model,data, K){
            kmrow = kmeans(data$X[1:data$Nrows,(data$Nrows+1):data$N],floor(K/2));
            kmcol = kmeans(t(data$X[1:data$Nrows,(data$Nrows+1):data$N]),floor(K/2));
            c(kmrow$cluster,kmcol$cluster+max(kmrow$cluster))
            #spectral(data$X,K)
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
            if(!(max(clust_cols)<min(clust_rows) | max(clust_rows)<min(clust_cols))){
              message("Co clustering failed bi partite structure not found")
              return(sol)
            }
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
            sol@clrow = sol@cl[irow]-min(sol@cl[irow])+1
            sol@Krow = max(sol@clrow)
            sol@clcol = sol@cl[icol]-min(sol@cl[icol])+1
            sol@Kcol = max(sol@clcol)
            sol@obs_stats$co_x_counts=sol@obs_stats$x_counts[clust_rows,clust_cols]
            
            tree= sol@ggtree
            tree=tree[tree$tree!=0,]
            leafs = tree[tree$H==0,]
            leafs = leafs[order(leafs$x),]
            # cols first
            if(max(clust_cols)<min(clust_rows)){
              xsplit = leafs$x[length(clust_cols)]
              coltree = tree[tree$x<=xsplit,]
              rowtree = tree[tree$x>xsplit,]
            # rows first
            }else{
              xsplit = leafs$x[length(clust_rows)]
              rowtree = tree[tree$x<=xsplit,]
              coltree = tree[tree$x>xsplit,]
            }
            sol@ggtreecol = coltree #coltree[2:nrow(coltree),]
            sol@ggtreerow = rowtree #rowtree[2:nrow(rowtree),]
            sol
          })

 setMethod(f = "sample_cl", 
           signature = signature("co_dcsbm","list","numeric"), 
           definition = function(model,data,K){
             c(sample(1:floor(K/2),data$Nrows,replace = TRUE),sample((floor(K/2)+1):K,data$Ncols,replace = TRUE))
           })


#' @rdname plot
#' @export
setMethod(f = "plot", 
          signature = signature("co_dcsbm_fit","missing"),
          definition = function(x,type="blocks"){
            switch(type,blocks=co_blocks(x),nodelink=co_nodelink(x))
          })

#' @rdname plot
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
