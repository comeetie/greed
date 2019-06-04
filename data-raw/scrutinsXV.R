library(rjson)
library(dplyr)
library(tidyr)
library(Matrix)
library(future)
plan(multisession)
scrutins_json = fromJSON(file="./data-raw/Scrutins_XV.json")
scrutins = scrutins_json$scrutin$scrutin
votes = NULL
for (is in 1:length(scrutins)){
  groupes = scrutins[[is]]$ventilationVotes$organe$groupes$groupe
  for (g in 1:length(groupes)){
    if(as.numeric(groupes[[g]]$vote$decompteVoix$pour)>1){
      pours   =  tibble(scrutin = is, groupe= groupes[[g]]$organeRef, acteur=sapply(groupes[[g]]$vote$decompteNominatif$pours$votant,function(v){v$acteurRef}),mandat=sapply(groupes[[g]]$vote$decompteNominatif$pours$votant,function(v){v$mandatRef}), vote='pour')
      votes = votes %>% bind_rows(pours) 
    }
    if(as.numeric(groupes[[g]]$vote$decompteVoix$pour)==1){
      pours   =  tibble(scrutin = is, groupe= groupes[[g]]$organeRef, acteur=groupes[[g]]$vote$decompteNominatif$pours$votant$acteurRef, mandat=groupes[[g]]$vote$decompteNominatif$pours$votant$mandatRef, vote='pour')
      votes = votes %>% bind_rows(pours) 
    }
    if(as.numeric(groupes[[g]]$vote$decompteVoix$contre)>1){
    contres = tibble(scrutin = is, groupe= groupes[[g]]$organeRef,acteur=sapply(groupes[[g]]$vote$decompteNominatif$contres$votant,function(v){v$acteurRef}), mandat=sapply(groupes[[g]]$vote$decompteNominatif$contres$votant,function(v){v$mandatRef}), vote='contre')
    votes = votes %>% bind_rows(contres)
    }
    if(as.numeric(groupes[[g]]$vote$decompteVoix$contre)==1){
      pours   =  tibble(scrutin = is, groupe= groupes[[g]]$organeRef,  acteur=groupes[[g]]$vote$decompteNominatif$contres$votant$acteurRef,mandat=groupes[[g]]$vote$decompteNominatif$contres$votant$mandatRef, vote='contre')
      votes = votes %>% bind_rows(pours) 
    }
  }
}





sumscr = votes %>% group_by(scrutin) %>% summarise(nbv = n(),p=mean(if_else(vote=="pour",1,0)))
vclean = votes %>% select(-groupe,-mandat) 

Xv=spread(vclean,scrutin,vote)
iss = sumscr %>% filter(nbv>100)
ij = which(Xv[,iss$scrutin+1]=="pour",arr.ind = TRUE)
Xxv=sparseMatrix(ij[,1],ij[,2],x = rep(1,nrow(ij)))


meta=fromJSON(file="/home/come/Bureau/AMO10_deputes_actifs_mandats_actifs_organes_XV.json")
meta=fromJSON(file="./data-raw/AMO30_tous_acteurs_tous_mandats_tous_organes_historique.json")


deputes.df=do.call(rbind,lapply(meta$export$acteurs$acteur,
                                function(dep){
                                  tibble(uid=dep$uid[[2]],name=paste(dep$etatCivil$ident[1:3],collapse = " "),
                                         mandat=sapply(dep$mandats$mandat,function(m){m$uid}),
                                         organe=sapply(dep$mandats$mandat,function(m){m$organes$organeRef[[1]]}))}))

orgs.df = do.call(rbind,lapply(meta$export$organes$organe,
                               function(org){
                                 tibble(uid=org$uid,name=org$libelle,
                                        name_abr = org$libelleAbrev,
                                        type=org$codeType)}))


parties.df = deputes.df %>% 
  left_join(orgs.df,by=c("organe"="uid")) %>% 
  filter(type =="PARPOL") %>% 
  arrange(uid,desc(mandat)) %>% 
  filter(!duplicated(uid)) %>%
  mutate(partie=name.y,partie_abr=name_abr) %>%
  select(uid,partie,partie_abr)

groupes.df =deputes.df%>% left_join(orgs.df,by=c("organe"="uid")) %>% 
  filter(type =="GP") %>% 
  arrange(uid,desc(mandat)) %>% 
  filter(!duplicated(uid)) %>%
  mutate(groupe=name.y,groupe_abr=name_abr) %>%
  select(uid,groupe,groupe_abr)

assemblee.df =deputes.df%>% left_join(orgs.df,by=c("organe"="uid")) %>% 
  filter(type =="ASSEMBLEE") %>% 
  arrange(uid,desc(mandat)) %>% 
  filter(!duplicated(uid)) %>%
  mutate(name=name.x,mandat_an=mandat) %>%
  select(uid,name,mandat_an)

dep_meta=assemblee.df %>% left_join(parties.df) %>% left_join(groupes.df)


Xg=sparseMatrix(c(1),c(1),x = c(0),dims = c(593+570,593+570))
Xg[1:593,594:(570+593)]=Xxv
Xg[594:(570+593),1:593]=t(Xxv)


sol=greed(Xxv,15)
plot(sol)
pprobs = post_probs(new("mm"),list(X=Xxv),sol@cl[1:593])

clusters=tibble(uid = unlist(Xv[1:nrow(Xxv),1]),cl=sol@cl[1:593],pr=apply(pprobs,1,max))
clusters.dep.df=clusters %>% left_join(dep_meta,by="uid")

table(clusters.dep.df$cl,clusters.dep.df$groupe)
