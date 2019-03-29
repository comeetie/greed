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

devtools::use_data(Xmigr.com)



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
MigDepNum = MigDepNum %>% filter(ito!=ifrom)
Xmigr.dep = sparseMatrix(MigDepNum$ifrom,MigDepNum$ito,x = MigDepNum$vol,dims = c(Nd,Nd))
colnames(Xmigr.dep)=DepLab$idinsee
devtools::use_data(Xmigr.dep)

data("Xmigr.dep")
DepLab=data.frame(idinsee=colnames(Xmigr.dep),id=1:nrow(Xmigr.dep),stringsAsFactors = FALSE)
library(greed)
sol=fit(Xmigr.dep,20,new("dcsbm"),new("hybrid",pop_size=50))



plot(sol,type='tree')
plot(sol,type='path')
Kf=12
sol.df = data.frame(code_insee=DepLab$idinsee,cl=cut(sol,Kf)@cl,stringsAsFactors = FALSE)
plot(cut(sol,Kf),type='nodelink')

DEP2017_sf = loadMap(COG=2017,nivsupra="DEP")
deps = DEP2017_sf %>% left_join(sol.df,by=c("DEP"="code_insee"))
typoLayer(x = deps, var = "cl", legend.values.order = 1:max(deps$cl))



# echelle bassin de vie aggregationdes communes
data("Xmigr.com")
BV2017_sf = loadMap(COG=2017,nivsupra="BV2012")
ComLab=data.frame(idinsee=colnames(Xmigr.com),id=1:nrow(Xmigr.com),stringsAsFactors = FALSE)
ComDescr=ComLab %>% left_join(table_supracom_2017%>% select(CODGEO,BV2012,DEP,REG,LIBGEO),by=c("idinsee"="CODGEO")) %>% left_join(sol.df,by=c("DEP"="code_insee"))

# echelle regionale interne
cldep = 10
ComsF = ComDescr %>% filter(cl==cldep,!is.na(cl)) %>% mutate(BV2012=factor(BV2012)) %>% mutate(bvnum=unclass(BV2012))
isel = ComDescr$cl==cldep & !is.na(ComDescr$cl)
X=init(new("dcsbm"),Xmigr.com[isel,isel],ComsF$bvnum)
X=Matrix(X@obs_stats$x_counts,sparse=TRUE)
colnames(X)=levels(ComsF$BV2012)
diag(X)=0
X=round(X)

BVLab = data.frame(BV2012=colnames(X),id=1:nrow(X),din=colSums(X),dout=rowSums(X),stringsAsFactors = FALSE)
sol=fit(X,20,new("dcsbm"),new("hybrid",pop_size=10))
sol@icl
plot(sol,type='path')
plot(sol,type='tree')
sol.c=cut(sol,20)
plot(sol.c)

BVLab$cl = sol.c@cl
bvs = BV2017_sf %>% left_join(BVLab) %>% filter(!is.na(cl))
typoLayer(x = bvs, var = "cl",legend.values.order = 1:max(bvs$cl))
bvs %>% filter(cl==3) %>% st_drop_geometry()


# analyse des résidus et comparaison avec un modèle gravitaire
cmat = sol.c@obs_stats$x_counts
cmatn= cmat/(sol.c@obs_stats$din%*%t(sol.c@obs_stats$dout))
Ep=Matrix(cmatn[sol.c@cl,sol.c@cl]*(colSums(X)%*%t(rowSums(X))))
diag(Ep)=0
ggpredflow = data.frame(prediction=Ep[X!=0],reality=X[X!=0])
sum(abs(Ep-X))
ggplot(ggpredflow)+geom_point(aes(x=reality,y=prediction))+scale_x_log10(limits=c(1,1000))+scale_y_log10(limits=c(1,1000))+geom_abline(color='red')


centers=st_centroid(st_geometry(bvs))
dists = st_distance(centers,centers)
metropoles = colnames(X)[order(rowSums(X),decreasing = TRUE)[1:3]]
metro_fact = colnames(X)
metro_fact[!metro_fact %in% metropoles]="autres"
grav_data = data.frame(flow=as.vector(X),d=as.vector(dists),dout=rep(rowSums(X),nrow(X)),din=rep(colSums(X),each=nrow(X)),
                       metr_fact_in=factor(rep(metro_fact,each=nrow(X))),metr_fact_out=factor(rep(metro_fact,nrow(X))))
lmfit=glm(flow ~ .,grav_data,family = poisson())
fp =predict(lmfit,grav_data,type = "response")
ggpredflow = data.frame(prediction=fp[grav_data$flow!=0],reality=grav_data$flow[grav_data$flow!=0])
ggplot(ggpredflow)+geom_point(aes(x=reality,y=prediction))+scale_x_log10(limits=c(1,1000))+scale_y_log10(limits=c(1,1000))+geom_abline(color='red')
sum(abs(fp-grav_data$flow))
summary(lmfit)
# avec les clusters 
grav_data = data.frame(flow=as.vector(X),d=as.vector(dists),dout=rep(rowSums(X),nrow(X)),din=rep(colSums(X),each=nrow(X)),
                       clust_in=factor(rep(sol.c@cl,each=nrow(X))),clust_out=factor(rep(sol.c@cl,nrow(X))))
lmfit=glm(flow ~ log(d)+log(din)+log(dout),grav_data%>% filter(d!=0),family = poisson())
fp =predict(lmfit,grav_data%>% filter(d!=0),type = "response")
ggpredflow = data.frame(prediction=fp[grav_data$flow!=0],reality=grav_data$flow[grav_data$flow!=0])
ggplot(ggpredflow)+geom_point(aes(x=reality,y=prediction))+scale_x_log10(limits=c(1,1000))+scale_y_log10(limits=c(1,1000))+geom_abline(color='red')
sum(abs(fp- (grav_data$flow%>% filter(d!=0))))
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
