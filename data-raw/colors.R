library(Matrix)
library(greed)
N=400
K=6
prop=rep(1/K,K)
lambda  = 0.1
lambda_o = 0.01
Ks=3
mu = bdiag(lapply(1:(K/Ks), function(k){
matrix(lambda_o,Ks,Ks)+diag(rep(lambda,Ks))}))+0.001
sbm = rsbm(N,prop,mu)
sol=greed(sbm$x,model=Sbm())

col_pos = cmdscale(-sol@path[[1]]$merge_mat,2)

r1= range(col_pos[,1])
c1 = ((col_pos[,1]-r1[1])/(r1[2]-r1[1])*0.9+0.05)*200-100

r2= range(col_pos[,2])
c2 = ((col_pos[,2]-r2[1])/(r2[2]-r2[1])*0.9+0.05)*200-100

r3= range(col_pos[,3])
b = ((col_pos[,3]-r3[1])/(r3[2]-r3[1])*0.8+0.1)*200-100


plot(colorspace::LAB(rep(70,6),c1,c2))

plot(sol@ggtree$x,sol@ggtree$H)

H=1828 



tree.df= sol_greed@ggtree
cnode = tree.df$node[tree.df$tree==0]
normfact = 150/tree.df$H[tree.df$node==cnode]
posi = c(0,0)
angles = c(0,pi)

expand_colors  = function(tree.df,cnode,angles,normfact,posi){
  children = tree.df$node[tree.df$tree==cnode]
  print(posi)
  if(length(children)>0){
    children_pos = purrr::map2(children,angles,function(child,angle){
      r = (tree.df$H[tree.df$node==cnode]-tree.df$H[tree.df$node==child])*normfact
      print(r)
      print(angle)
      c(child,posi[1]+r*cos(angle),posi[2]+r*sin(angle))
    })
    newpos = do.call(rbind,children_pos)
    newposc1 = expand_colors(tree.df,children[1],angles+pi/2,normfact,newpos[1,2:3])
    newposc2 = expand_colors(tree.df,children[2],angles+pi/2,normfact,newpos[2,2:3])
    pos = rbind(newpos,newposc1,newposc2)
  }else{
    pos=c()
  }
  pos
}

cols = expand_colors(tree.df,cnode,angles,normfact,c(0,0))

plot(colorspace::LAB(rep(70,6),cols[,2],cols[,3]))

children = sol@ggtree$node[sol@ggtree$tree==cnode]

x=cut(sol_greed,12)
x=sol_dcsbm
K=x@K
tdif = c(0,-x@path[[1]]$merge_mat[cbind(2:K,1:(K-1))],-x@path[[1]]$merge_mat[K,1])
tdifnorm = tdif/sum(tdif)
tdifclip = tdifnorm
pmax= min(12,K)
tdifclip[tdifnorm>1/pmax]=1/pmax
tdifclip[tdifnorm<1/pmax]=tdif[tdifnorm<1/pmax]/sum(tdif[tdifnorm<1/pmax])*(1-sum(tdifnorm>1/pmax)/pmax)
hue = cumsum(tdifclip/sum(tdifclip))*360
plot(colorspace::polarLAB(rep(70,K),rep(100,K),hue[1:K]))



