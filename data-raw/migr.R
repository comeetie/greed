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



# Echelle Departementale
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

library(greed)
sol=fit(Xmigr.dep,25,new("dcsbm"))

Kf=12
sol.df = data.frame(code_insee=DepLab$idinsee,cl=factor(cut(sol,Kf)@cl,levels=1:Kf),stringsAsFactors = FALSE)
plot(cut(sol,Kf))
plot(sol,type='tree')
plot(cut(sol,Kf),type='nodelink')

DEP2017_sf = loadMap(COG=2017,nivsupra="DEP")
deps = DEP2017_sf %>% left_join(sol.df,by=c("DEP"="code_insee"))
typoLayer(x = deps, var = "cl")


library(future)
plan(multisession)
sol=fit(Xmigr.dep,25,new("dcsbm"))
