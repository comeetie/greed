library(readr)
library(dplyr)
library(Matrix)
library(stringr)
library(sf)
library(cartography)
library(CARTElette)
library(COGugaison)
library(ggplot2)
Mig=read_delim("./data-raw/FD_MIGCOM_2015.txt",del=';')

# Echelle Communale
MigCom=Mig %>% enlever_PLM("COMMUNE",agregation = FALSE) %>% 
  enlever_PLM("DCRAN",agregation = FALSE) %>% 
  filter(IRAN>1) %>% #supprimer les personnes n'ayant pas déménagées
  group_by(COMMUNE,DCRAN) %>% summarise(vol = sum(IPONDI)) %>% rename(to=COMMUNE,from=DCRAN)

COM2017_sf = loadMap(COG=2017,nivsupra="COM")


ComLab = data.frame(idinsee=unique(c(MigCom$from,MigCom$to)),stringsAsFactors = FALSE) %>% 
  filter(idinsee %in% COM2017_sf$INSEE_COM) %>%
  mutate(id=1:n())
Nc = nrow(ComLab)
MigCom$from[!MigCom$from %in% ComLab$idinsee]="99999"
ComLab = rbind(ComLab,data.frame(idinsee="99999",id=Nc+1))

MigComNum = MigCom %>% left_join(ComLab,by=c("from"="idinsee"))%>% rename(ifrom=id) %>%
  left_join(ComLab,by=c("to"="idinsee"))%>% rename(ito=id)
Nc = max(c(MigComNum$ifrom,MigComNum$ito))
Xmigr.com = sparseMatrix(MigComNum$ifrom,MigComNum$ito,x = MigComNum$vol,dims = c(Nc,Nc))
colnames(Xmigr.com)=ComLab$idinsee
devtools::use_data(Xmigr.com,overwrite = TRUE)



# Echelle Departementale clutser regionaux ?
ComDep = COM2017_sf %>% select(INSEE_COM,INSEE_DEP) %>% st_drop_geometry()
ComDep = rbind(ComDep,data.frame(INSEE_COM="99999",INSEE_DEP="99",stringsAsFactors = FALSE))
MigDep = MigCom %>% left_join(ComDep,by=c("from"="INSEE_COM")) %>% rename(from_dep=INSEE_DEP) %>%
  left_join(ComDep,by=c("to"="INSEE_COM")) %>% rename(to_dep=INSEE_DEP) %>%
  group_by(from_dep,to_dep) %>%
  summarise(vol=round(sum(vol))) %>%
  rename(from=from_dep,to=to_dep)

DepLab = data.frame(idinsee=unique(c(MigDep$from,MigDep$to)),stringsAsFactors = FALSE) %>% mutate(id=1:n())
MigDepNum = MigDep %>% left_join(DepLab,by=c("from"="idinsee"))%>% rename(ifrom=id) %>%
  left_join(DepLab,by=c("to"="idinsee"))%>% rename(ito=id)
Nd = max(c(MigDepNum$ifrom,MigDepNum$ito))
MigDepNum = MigDepNum
Xmigr.dep = sparseMatrix(MigDepNum$ifrom,MigDepNum$ito,x = MigDepNum$vol,dims = c(Nd,Nd))
colnames(Xmigr.dep)=DepLab$idinsee

devtools::use_data(Xmigr.dep,overwrite = TRUE)

data("Xmigr.dep")
DepLab=data.frame(idinsee=colnames(Xmigr.dep),id=1:nrow(Xmigr.dep),stringsAsFactors = FALSE)
library(greed)
sol=fit(Xmigr.dep,20,new("dcsbm"))
plot(sol,type='tree')
plot(sol,type='path')
Kf=13
sol.df = data.frame(code_insee=DepLab$idinsee,cl=cut(sol,Kf)@cl,dintern=diag(Xmigr.dep),stringsAsFactors = FALSE)
plot(cut(sol,Kf),type='nodelink')
plot(cut(sol,Kf))
DEP2017_sf = loadMap(COG=2017,nivsupra="DEP")
deps = DEP2017_sf %>% left_join(sol.df,by=c("DEP"="code_insee"))

cols =c("#777777",
        "#a6cee3",
        "#1f78b4",
        "#b2df8a",
        "#33a02c",
        "#fb9a99",
        "#e31a1c",
        "#fdbf6f",
        "#ff7f00",
        "#cab2d6",
        "#6a3d9a",
        "#ffff99",
        "#b15928",
        "#8dd3c7")

typoLayer(x = deps, var = "cl", legend.values.order = 1:max(deps$cl),col = cols[2:(Kf+1)],legend.title.txt = "Clusters :")

#Flowmap

xy = DEP2017_sf %>% st_geometry()%>%st_centroid() %>% st_coordinates()
dep_points = data.frame(x=xy[,1],y=xy[,2],DEP=DEP2017_sf$DEP,name=DEP2017_sf$nom,stringsAsFactors = FALSE)


data("Xmigr.dep")
Xmigr.dep=Xmigr.dep+t(Xmigr.dep)
Xmigr.dep=tril(Xmigr.dep)
ij=which(Xmigr.dep>0,arr.ind = TRUE)
depname=colnames(Xmigr.dep)

cl = cut(sol,Kf)@cl
gg=data.frame(i=depname[ij[,1]],j=depname[ij[,2]],f=Xmigr.dep[ij],col=factor(ifelse(cl[ij[,1]]==cl[ij[,2]],cl[ij[,1]],0),levels=0:max(cl)),stringsAsFactors = FALSE)

ggf = gg %>% left_join(dep_points,by=c("i"="DEP")) %>% rename(from=i,fromx=x,fromy=y) %>%
  left_join(dep_points,by=c("j"="DEP")) %>% rename(to=j,tox=x,toy=y)

names(cols)=0:max(cl)
ggplot(ggf%>%arrange(f,col)%>%filter(f>500))+
  geom_segment(aes(x=fromx,y=fromy,xend=tox,yend=toy,size=f,alpha=log(f),color=col))+
  geom_point(data=dep_points %>% left_join(sol.df,by=c("DEP"="code_insee")),aes(x=x,y=y,col=factor(cl,levels=1:max(cl))),size=1.5)+
  coord_equal()+
  theme_minimal()+
  scale_size_area("Flux",limits = c(0,35000))+
  scale_alpha_continuous(guide="none",range = c(0,1),limits=c(log(500),log(35000)))+
  scale_x_continuous("",breaks = c())+
  scale_y_continuous("",breaks = c())+
  scale_color_manual(values = cols,guide='none')# echelle bassin de vie aggregationdes communes

graph_blocks_degnorm(cut(sol,Kf))
graph_blocks_balance(cut(sol,Kf))

data("Xmigr.com")
BV2017_sf = loadMap(COG=2017,nivsupra="BV2012")
ComLab=data.frame(idinsee=colnames(Xmigr.com),id=1:nrow(Xmigr.com),stringsAsFactors = FALSE)
ComDescr=ComLab %>% left_join(table_supracom_2017%>% select(CODGEO,BV2012,DEP,REG,LIBGEO),by=c("idinsee"="CODGEO")) %>% left_join(sol.df,by=c("DEP"="code_insee"))

# echelle regionale interne
cldep = 8
ComsF = ComDescr %>% filter(cl==cldep,!is.na(cl)) %>% mutate(BV2012=factor(BV2012)) %>% mutate(bvnum=unclass(BV2012))
isel = ComDescr$cl==cldep & !is.na(ComDescr$cl)
X=init(new("dcsbm"),Xmigr.com[isel,isel],ComsF$bvnum)
X=Matrix(X@obs_stats$x_counts,sparse=TRUE)
colnames(X)=levels(ComsF$BV2012)
X=round(X)
diag(X)=0
BVLab = data.frame(BV2012=colnames(X),id=1:nrow(X),din=colSums(X),dout=rowSums(X),stringsAsFactors = FALSE)



library(future)
plan(multisession)
sol=fit(X,20,new("dcsbm"),new("hybrid",pop_size=40))
sol@icl
plot(cut(sol,30),type='path')
plot(sol,type='tree')
sol.c=cut(sol,10)
plot(sol.c)

cl=sol.c@cl
BVLab$cl = cl
bvs = BV2017_sf %>% left_join(BVLab) %>% filter(!is.na(cl))


col20 =c("#777777",
         "#393b79",
         "#5254a3",
         "#6b6ecf",
         "#9c9ede",
         "#637939",
         "#8ca252",
         "#b5cf6b",
         "#cedb9c",
         "#8c6d31",
         "#bd9e39",
         "#e7ba52",
         "#e7cb94",
         "#843c39",
         "#ad494a",
         "#d6616b",
         "#e7969c",
         "#7b4173",
         "#a55194",
         "#ce6dbd",
         "#de9ed6")
cols= col20[1:(max(cl)+1)]
names(cols)=0:max(cl)

typoLayer(x = bvs, var = "cl",legend.values.order = sort(unique(bvs$cl)), col=cols[sort(unique(bvs$cl))])


xy = BV2017_sf %>% st_geometry()%>%st_centroid() %>% st_coordinates()
bv_points = data.frame(x=xy[,1],y=xy[,2],BV=BV2017_sf$BV2012,name=BV2017_sf$nom,stringsAsFactors = FALSE)



ij=which(X>0,arr.ind = TRUE)
depname=colnames(X)

gg=data.frame(i=depname[ij[,1]],j=depname[ij[,2]],f=X[ij],col=factor(ifelse(cl[ij[,1]]==cl[ij[,2]],cl[ij[,1]],0),levels=0:max(cl)),stringsAsFactors = FALSE)

ggf = gg %>% left_join(bv_points,by=c("i"="BV")) %>% rename(from=i,fromx=x,fromy=y) %>%
  left_join(bv_points,by=c("j"="BV")) %>% rename(to=j,tox=x,toy=y)

vm=max(X-diag(diag(X)))
fth = vm*0.01
ggplot(ggf%>%arrange(col,f)%>%filter(f>fth))+
  geom_segment(aes(x=fromx,y=fromy,xend=tox,yend=toy,size=f,alpha=log(f),color=col))+
  geom_point(data=bv_points %>% left_join(BVLab,by=c("BV"="BV2012"))%>%filter(!is.na(cl)),aes(x=x,y=y,col=factor(cl,levels=1:max(cl,na.rm=TRUE))),size=2)+
  coord_equal()+
  theme_minimal()+
  scale_size_area("Flux :",limits = c(0,1.1*vm),max_size=6)+
  scale_alpha_continuous(guide="none",range = c(0,0.8),limits=c(log(fth),log(vm)))+
  scale_x_continuous("",breaks = c())+
  scale_y_continuous("",breaks = c())+
  scale_color_manual(values = cols,guide='none')# echelle bassin de vie aggregationdes communes

graph_blocks_degnorm(sol.c)
graph_blocks_balance(sol.c)


bvs %>% filter(cl==10) %>% st_drop_geometry()

# analyse des résidus et comparaison avec un modèle gravitaire
cmat = sol.c@obs_stats$x_counts
cmatn= cmat/(sol.c@obs_stats$din%*%t(sol.c@obs_stats$dout))
Ep=Matrix(cmatn[sol.c@cl,sol.c@cl]*(colSums(X)%*%t(rowSums(X))))
diag(Ep)=0
ggpredflow = data.frame(prediction=as.vector(Ep),reality=as.vector(X)) 
mv=max(Ep)
sum(abs(Ep-X))
ggplot(ggpredflow)+geom_point(aes(x=reality,y=prediction))+geom_abline(color='red')+scale_x_log10(limits=c(1,1.1*mv))+scale_y_log10(limits=c(1,1.1*mv))


centers=st_centroid(st_geometry(bvs))
dists = st_distance(centers,centers)
metropoles = colnames(X)[order(rowSums(X)+colSums(X),decreasing = TRUE)[1:2]]
metro_fact = rep("metropoles",nrow(X))
metro_fact[!colnames(X) %in% metropoles]="autres"
grav_data = data.frame(flow=as.vector(X),d=as.vector(dists),
                       dout=rep(rowSums(X),nrow(X)),din=rep(colSums(X),each=nrow(X)),
                       cin=rep(metro_fact,each=nrow(X)),cout=rep(metro_fact,nrow(X)))
grav_nodiag = grav_data %>% filter(d>0 & din>0 & dout>0)

lmfit=glm(flow ~ log(d)+log(dout)+log(din),grav_nodiag,family = poisson())
fp =predict(lmfit,grav_nodiag,type = "response")
ggpredflow = data.frame(prediction=fp,reality=grav_nodiag$flow,cin=grav_nodiag$cin,cout=grav_nodiag$cout)
mv= max(grav_nodiag$flow)
ggplot(ggpredflow)+geom_point(aes(x=reality,y=prediction,col=cout))+scale_x_log10(limits=c(1,1.1*mv))+scale_y_log10(limits=c(1,1.1*mv))+geom_abline(color='red')
sum(abs(fp-grav_nodiag$flow))
summary(lmfit)


# echelle bassin de vie descitpiton concentrique + mm
# liste des poles

Xc = c()
Lpoles = ComDescr %>% mutate(vol=(colSums(Xmigr.com)+rowSums(Xmigr.com))) %>% filter(BV2012==idinsee) %>% filter(vol>50000)
for (i in 1:nrow(Lpoles)){
  idpole = Lpoles$idinsee[i]
  namepole = Lpoles$LIBGEO[i]
  reg = Lpoles$REG[i]
  print(namepole)
  levs = c(idpole,paste0("BV",idpole),"REG","FF","99999")
  ComDescrBVl = ComDescr %>% 
    mutate(clbv = factor(case_when(idinsee==idpole ~ idpole, BV2012==idpole ~ paste0("BV",BV2012),REG==reg ~ "REG",idinsee=="99999"~idinsee,TRUE ~ "FF"),levels=levs)) %>%
    mutate(clnum=unclass(clbv))
  
  Xm=init(new("dcsbm"),Xmigr.com,ComDescrBVl$clnum)@obs_stats$x_counts
  x=round(c(Xm[2:5,1],Xm[1,2:4]))
  Xc=rbind(Xc,x)
}
Xc=Matrix(Xc,sparse = TRUE)
sol=fit(Xc,20)
plot(sol,type='tree')
plot(sol,type='path')

sol.c=cut(sol,10)

Lpoles$cl = sol.c@cl
bvs = BV2017_sf %>% left_join(Lpoles) %>% filter(!is.na(cl))
typoLayer(x = bvs, var = "cl")


cmat= sol.c@obs_stats$x_counts
cmat = cmat/rowSums(cmat)
K = nrow(cmat)
labs = c("Bassin de vie -> Pôle","Region -> Pôle","France -> Pôle","Etranger -> Pôle","Pôle -> Bassin de vie","Pôle -> Région","Pôle -> France")
levs = c("Bassin de vie -> Pôle","Pôle -> Bassin de vie","Region -> Pôle","Pôle -> Région","France -> Pôle","Pôle -> France","Etranger -> Pôle")

gg= data.frame(p=as.vector(cmat),type=factor(rep(labs,each=K),levels=levs),K=rep(1:K,7))

ggplot(gg)+geom_bar(aes(x=type,y=p,group=K),stat = "identity")+facet_wrap(~K)+coord_flip()
