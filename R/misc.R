

dendo =function(resH){
  ggtree=data.frame(node=1:length(resH$tree),father=resH$tree,H=c(rep(0,resH$stats0$K),-resH$logalpha),x=0,xmin=0,xmax=0,k=c(rep(resH$stats0$K,resH$stats0$K),(resH$stats0$K-1):1))
  fathers = c(which(resH$tree==0))
  w = 1/2
  for(l in 1:resH$stats0$K){
    for(f in fathers){
      sons=which(ggtree$father==f)
      ggtree$x[sons[1]]=ggtree$x[f]+w
      ggtree$x[sons[2]]=ggtree$x[f]-w
    }
    fathers=ggtree$node[ggtree$father %in% fathers]
    w =w*0.5
  }
  leafs = which(ggtree$H==0)
  or = order(ggtree[leafs,"x"])
  ggtree$x[leafs[or]]=seq(-1,1,length.out = length(leafs))
  for(n in setdiff(1:nrow(ggtree),leafs)){
    sons=which(ggtree$father==n)
    ggtree$x[n]=mean(ggtree$x[sons])
    ggtree$xmin[n] = min(ggtree$x[sons])
    ggtree$xmax[n] = max(ggtree$x[sons])
  }
  ggleaf=ggtree[ggtree$father!=0,]
  
  #roots = ggtree[ggtree$k>=k,]
  #roots = roots[!(roots$father %in% roots$node),]
  #print(dim(roots))
  ggplot()+geom_segment(aes(x=ggleaf$x,y=ggleaf$H,xend=ggleaf$x,yend=ggtree$H[ggleaf$father]))+
    geom_segment(data=ggtree[ggtree$H>0,],aes(x=xmin,y=H,xend=xmax,yend=H))+
    #geom_point(data=roots,aes(x=x,y=H),color="red")+
    scale_x_continuous("",labels = c(),breaks = c())+
    scale_y_continuous("-log(alpha)")+theme_bw()
}

clust_permutation = function(resH,K){
  ggtree=data.frame(node=1:length(resH$tree),father=resH$tree,H=c(rep(0,resH$stats0$K),-resH$logalpha),k=c(rep(resH$stats0$K,resH$stats0$K),(resH$stats0$K-1):1),x=0,xmin=0,xmax=0)
  fathers = c(which(resH$tree==0))
  w = 1/2
  for(l in 1:resH$stats0$K){
    for(f in fathers){
      sons=which(ggtree$father==f)
      ggtree$x[sons[1]]=ggtree$x[f]+w
      ggtree$x[sons[2]]=ggtree$x[f]-w
    }
    fathers=ggtree$node[ggtree$father %in% fathers]
    w =w*0.5
  }
  roots = ggtree[ggtree$k>=K,]
  roots = roots[!(roots$father %in% roots$node),]
  order(roots$x)
}

plotpath = function(hierch){
  gg=data.frame(icl=sapply(hierch@path,function(p){p$icl}),logalpha=sapply(hierch@path,function(p){p$logalpha}),K=sapply(hierch@path,function(p){p$K}))
  gg$logalpha[1:6]=gg$logalpha[2:7]
  ggplot(gg)+geom_abline(aes(slope=K-1,intercept=icl))
}