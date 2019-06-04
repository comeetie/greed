#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
#' @importFrom future %<-%
#' @name %<-%
NULL

#' @include models_classes.R fit_classes.R cleanpath.R
#' @import Matrix
NULL

hybrid = function(model,alg,data,K, verbose=FALSE){
            
            fi = function(ncl){ fit_greed(model,data,ncl,type="merge",verbose=verbose) }
            train.hist = data.frame(generation=c(),icl=c(),K=c())

            # multi-start in //
            #future::plan(future::multiprocess)
            
            solutions = listenv::listenv()
            # first generation of solutions
            pop_size = alg@pop_size
            for (i in 1:pop_size){
              solutions[[i]] %<-% fit_greed(model,data,sample(1:K,data$N,replace = TRUE),verbose = verbose)
            }
            solutions = as.list(solutions)
            icls  = sapply(solutions,function(s){s@icl})
            # check for errors 
            solutions=solutions[!is.nan(icls)]
            icls=icls[!is.nan(icls)]
            old_best = -Inf
            best_icl = max(icls)
            nbgen = 1
            # while maximum number of generation // all solutions are equals // no improvements

            while((max(icls)-min(icls))>1 & (best_icl > old_best) & nbgen < alg@nb_max_gen){
              
              
              train.hist=rbind(train.hist,data.frame(generation=nbgen,icl=icls,K=sapply(solutions,function(s){max(s@cl)})))
              # selection keep the top half solutions
              icl_order = order(icls,decreasing = TRUE)
              selected  = icl_order[1:(pop_size/2)]
              # cross_over between the kept solution
              new_solutions = listenv::listenv()
              selected_couples = matrix(selected[sample(1:length(selected),length(selected)*2,replace = TRUE)],ncol=2)
              for (i in 1:nrow(selected_couples)){
                new_solutions[[i]] %<-% full_cross_over(solutions[[selected_couples[i,1]]],solutions[[selected_couples[i,2]]],fi)
              }
              new_solutions = as.list(new_solutions)
              solutions = c(solutions[selected],new_solutions)
              icls = sapply(solutions,function(s){s@icl})
              old_best=best_icl
              best_icl = max(icls)
              nbgen = nbgen + 1;
            }
            

            train.hist=rbind(train.hist,data.frame(generation=nbgen,icl=icls,K=sapply(solutions,function(s){max(s@cl)})))
            #parallel::stopCluster(cl)
            # best solution
            res = solutions[[order(icls,decreasing = TRUE)[1]]]
            # compute merge path
            path = fit_greed_path(data,res)
            # clean the resuts (compute, merge tree,...)
            path = cleanpath(path)
            # store train history
            path@train_hist = train.hist
            # stop future plan
            #oplan <- future::plan()
            #on.exit(future::plan(oplan), add = TRUE)
            path
          }


full_cross_over = function(sol1,sol2,fi){
  # cartesian product on the z of the two solution
  ncl = unclass(factor(paste(sol1@cl,sol2@cl)))
  # greedy merge
  fi(ncl)
}




