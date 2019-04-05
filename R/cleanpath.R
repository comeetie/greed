#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
#' @importFrom future %<-%
#' @name %<-%
NULL

#' @include models_classes.R fit_classes.R
#' @import Matrix
NULL

# clean the merge path 
cleanpath = function(pathsol){
  K=pathsol@K
  pathsol@logalpha = 0
  path=pathsol@path
  tree =c(0)
  xtree=c(0)
  cn  = 1
  lab = c(1)
  xpos = c(0)
  H = rep(0,2*K-1)
  x = 0
  w = 0.5
  K=1
  for (lev in seq(length(path),1)){
    pl = length(path)-lev
    ord = order(xpos)
    path[[lev]]$obs_stats = reorder(pathsol@model,path[[lev]]$obs_stats,ord)
    path[[lev]]$cl  = order(ord)[path[[lev]]$cl] 
    k=path[[lev]]$k
    l=path[[lev]]$l
    tree=c(tree,lab[l],lab[l])
    
    H[lab[l]]=-path[[lev]]$logalpha
    if(tree[lab[l]]!=0 && H[lab[l]]>H[tree[lab[l]]]){
      H[lab[l]]=H[tree[lab[l]]]
    }
    
    lab[l]=cn+1
    fpos = xpos[l]
    xtree=c(xtree,fpos-w^pl,fpos+w^pl)
    # choisir + ou - en fonctionde la taille ?
    xpos[l]=fpos-w^pl
    if(k>K){
      xpos = c(xpos,fpos+w^pl)
      lab=c(lab,cn+2)  
    }else{
      xpos = c(xpos[1:(k-1)],fpos+w^pl,xpos[k:length(lab)])
      lab=c(lab[1:(k-1)],cn+2,lab[k:length(lab)])
    }
    K  = K+1
    cn = cn + 2 
  }
  ggtree=data.frame(H=H,tree=tree,x=xtree,node=1:length(tree),xmin=0,xmax=0)
  # recalculer les x
  leafs = which(ggtree$H==0)
  or = order(ggtree[leafs,"x"])
  ggtree$x[leafs[or]]=seq(-1,1,length.out = length(leafs))
  others = ggtree$node[ggtree$H!=0]
  for(n in others[seq(length(others),1)]){
    sons=which(ggtree$tree==n)
    ggtree$x[n]=mean(ggtree$x[sons])
    ggtree$xmin[n] = min(ggtree$x[sons])
    ggtree$xmax[n] = max(ggtree$x[sons])
  }
  
  
  
  ggtree$Hend = c(-1,ggtree$H[ggtree$tree])
  ggtree$xend = c(-1,ggtree$x[ggtree$tree])
  
  or = order(xpos)
  pathsol@obs_stats = reorder(pathsol@model,pathsol@obs_stats,or)
  pathsol@cl=order(or)[pathsol@cl] 
  pathsol@path = path
  pathsol@tree = tree
  pathsol@ggtree = ggtree 
  pathsol
} 
