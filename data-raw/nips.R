library(readr)
library(Matrix)
ijn = read_delim("./data-raw/docword.nips.txt",delim=' ',skip=3,col_names = c("i","j","n"))
X = Matrix::sparseMatrix(ijn$i,ijn$j,x=ijn$n)
words = readLines("./data-raw/vocab.nips.txt")
th=50
Xs = X[,colSums(X)>th]
voc = words[colSums(X)>th]
colnames(Xs)=voc

Nips=Xs

library(greed)
library(future)
plan(multisession)
sol=greed(Nips,model=DcLbm(),alg=Hybrid())


plot(sol)
plot(cut(sol,60))

sol_simple = cut(sol,30)
co=coef(sol_simple)
voc.df=data.frame(voc=voc,cl=sol@clcol,cls=sol_simple@clcol,w=as.numeric(co$gammacols)) %>% arrange(desc(cl),desc(w))
wc=voc.df %>% group_by(cls) %>% top_n(3,wt=w) %>% ungroup()
ggplot(wc)+geom_text(aes(label=voc,x=cl,y=-log(w),size=w,color=factor(cl)))+scale_color_discrete(guide="none")+scale_size(guide="none",range=c(3,6))

sel_cl=14
pcl=colSums(Nips[sol_simple@clrow==sel_cl,]/rowSums(Nips[sol_simple@clrow==sel_cl,]))*(1/sum(sol_simple@clrow==sel_cl))
pref=colSums(Nips/rowSums(Nips))*(1/nrow(Nips))
sort(log(pcl/pref),decreasing = TRUE) %>% head(20)

