library(readr)
library(dplyr)
library(Matrix)
library(stringr)
library(sf)
library(cartography)
Mig=read_delim("./data-raw/FD_MIGCOM_2015.txt",del=';')

# Echelle Communale

library(CARTElette)
library(COGugaison)

# gestion des arrondissement et aggregation
MigCom=Mig %>% 
  enlever_PLM("COMMUNE",agregation=FALSE) %>%
  enlever_PLM("DCRAN",agregation=FALSE) %>%
  group_by(COMMUNE,DCRAN) %>% 
  summarise(vol = sum(IPONDI)) %>% 
  rename(to=COMMUNE,from=DCRAN)


COM2017_sf <- loadMap(COG=2017,nivsupra="COM")
ComLab = data.frame(idinsee=unique(c(MigCom$from,MigCom$to)),stringsAsFactors = FALSE) %>% 
  filter(idinsee %in% COM2017_sf$INSEE_COM) %>% 
  mutate(id=1:n())

Nc = nrow(ComLab)
# gestion des communauté d'outre mer et de l'étranger
Etr = data.frame(idinsee=unique(c(MigCom$from,MigCom$to)),stringsAsFactors = FALSE) %>% filter(!idinsee %in% COM2017_sf$INSEE_COM) %>% mutate(id=Nc+1)
ComLab = ComLab %>% bind_rows(Etr)

MigComNum = MigCom %>% left_join(ComLab,by=c("from"="idinsee"))%>% rename(ifrom=id) %>%
  left_join(ComLab,by=c("to"="idinsee")) %>% rename(ito=id)

Xmigr.com = sparseMatrix(MigComNum$ifrom,MigComNum$ito,x = MigComNum$vol,dims = c(Nc+1,Nc+1))
colnames(Xmigr.com)=c(ComLab$idinsee[1:Nc],"99999")
sol.com=fit(Xmigr.com,50,new("dcsbm"))



# Echelle Departementale
DEP2017_sf <- loadMap(COG=2017,nivsupra="DEP")

MigDep=Mig %>% mutate(to = str_sub(COMMUNE,1,2),from=str_sub(DCRAN,1,2)) %>% group_by(from,to) %>% summarise(vol = sum(IPONDI))
DepLab = data.frame(idinsee=unique(c(MigDep$from,MigDep$to)),stringsAsFactors = FALSE) %>% mutate(id=1:n())
MigDepNum = MigDep %>% left_join(DepLab,by=c("from"="idinsee"))%>% rename(ifrom=id) %>%
  left_join(DepLab,by=c("to"="idinsee"))%>% rename(ito=id)
Nd = max(c(MigDepNum$ifrom,MigDepNum$ito))
Xmigr.dep = sparseMatrix(MigDepNum$ifrom,MigDepNum$ito,x = MigDepNum$vol,dims = c(Nd,Nd))
colnames(Xmigr.dep)=DepLab$idinsee


sol=fit(Xmigr.dep,25,new("dcsbm"))
Kf=12
sol.df = data.frame(code_insee=DepLab$idinsee,cl=factor(cut(sol,Kf)@cl,levels=1:Kf),stringsAsFactors = FALSE)
plot(cut(sol,Kf))
plot(cut(sol,Kf),type='nodelink')
deps=DEP2017_sf %>% left_join(sol.df,by=c("DEP"="code_insee")) 
typoLayer(x = deps, var = "cl")

is = sample(33000,200)
Xs = Xmigr.com[is,is]
na=colnames(Xs)


