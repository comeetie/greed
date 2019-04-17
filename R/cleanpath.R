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
  
  
  # check for possible better solution than init with alpha=1 along the path
  if(length(path)>0){
    icli = sapply(path,function(p){p$icl1})
    if(max(icli)>pathsol@icl){
      im = which.max(icli)
      K = path[[im]]$K
      pathsol@K = K
      pathsol@obs_stats = path[[im]]$obs_stats
      pathsol@icl = path[[im]]$icl1
      pathsol@cl = as.vector(path[[im]]$cl)
      if(im>1){
        path=path[(im-1):1]  
      }else{
        path=list()
      }
      pathsol@path=path
    }
    
    
    
    
    if(length(path)>0){
      Hfront = extract_front_height(pathsol)
      # initialisation
      # vector with the tree information
      tree =c(0)
      # x position of the tree nodes the root is at 0 
      xtree=c(0)
      # current node
      cn  = 1
      
      # vector of length K which cluster name in tree notation
      lab = c(1)
      # vector of length K with cluster x positions
      xpos = c(0)
      # height of the nodes
      H  = rep(0,2*K-1)
      Kc = rep(K,2*K-1)
      x = 0
      w = 0.5
      K=1
      # go back from the tree root
      # is there a tree here ?
      
      best_merge = length(path)
      # for each merge
      for (lev in seq(length(path),1)){
        
        # merge number starting from root
        pl = length(path)-lev
        # get the order form the x position of the cluster at this level
        ord = order(xpos)
        
        # reorder the cluster accordingly
        path[[lev]]$obs_stats = reorder(pathsol@model,path[[lev]]$obs_stats,ord)
        path[[lev]]$cl  = order(ord)[path[[lev]]$cl]
        
        # store in the tree the merge of k and l
        k=path[[lev]]$k
        l=path[[lev]]$l
        tree=c(tree,lab[l],lab[l])
        
        
        #H[lab[l]]=-path[[lev]]$logalpha
        #if(tree[lab[l]]!=0 && H[lab[l]]>H[tree[lab[l]]]){
        #  H[lab[l]]=H[tree[lab[l]]]
        #}
        Kc[lab[l]]=K
        H[lab[l]]=Hfront[K]
        
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
      

      ggtree=data.frame(H=H,tree=tree,x=xtree,node=1:length(tree),xmin=0,xmax=0,K=Kc)
      # recompute the x bottom to top for constant spacing of leafs
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
    }else{
      pathsol@tree=c(0)
      pathsol@ggtree = data.frame(H=0,tree=0,x=0,node=1,xmin=0,max=0)
    }
  }else{
    pathsol@tree=c(0)
    pathsol@ggtree = data.frame(H=0,tree=0,x=0,node=1,xmin=0,max=0)
  }
  pathsol
} 


extract_front_height=function(sol){
  icl=c(sol@icl,sapply(sol@path,function(v){v$icl1}))
  icl = icl[length(icl):1]
  K = 1:length(icl)
  
  H=rep(0,length(icl))
  cdi = Inf
  bestline = 1
  Front = c(1)
  
  for (l in 2:length(icl)){
    di = (icl[l]-icl[bestline])
    din = di/(l-bestline)
    if (di > 0){
      while(din > cdi & length(Front)>1){
        Front=Front[-length(Front)]
        H[bestline]=-1
        
        bestline = Front[length(Front)]
        
        di = (icl[l]-icl[bestline])
        din = di/(l-bestline)
        cdi = H[bestline]
        
      }
      H[Front[length(Front)]]=din
      cdi = din
      bestline = l
      Front=c(Front,l)
      
    }else{
      H[l]= -1
    }
  }
  
  for(l in 2:length(icl)){
    if(H[l]==-1){
      H[l]=H[l-1]
    }
  }
  H
}



