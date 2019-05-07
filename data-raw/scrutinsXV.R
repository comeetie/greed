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
      pours   =  tibble(scrutin = is, groupe= groupes[[g]]$organeRef, mandat=sapply(groupes[[g]]$vote$decompteNominatif$pours$votant,function(v){v$mandatRef}), vote='pour')
      votes = votes %>% bind_rows(pours) 
    }
    if(as.numeric(groupes[[g]]$vote$decompteVoix$pour)==1){
      pours   =  tibble(scrutin = is, groupe= groupes[[g]]$organeRef, mandat=groupes[[g]]$vote$decompteNominatif$pours$votant$mandatRef, vote='pour')
      votes = votes %>% bind_rows(pours) 
    }
    if(as.numeric(groupes[[g]]$vote$decompteVoix$contre)>1){
    contres = tibble(scrutin = is, groupe= groupes[[g]]$organeRef, mandat=sapply(groupes[[g]]$vote$decompteNominatif$contres$votant,function(v){v$mandatRef}), vote='contre')
    votes = votes %>% bind_rows(contres)
    }
    if(as.numeric(groupes[[g]]$vote$decompteVoix$contre)==1){
      pours   =  tibble(scrutin = is, groupe= groupes[[g]]$organeRef, mandat=groupes[[g]]$vote$decompteNominatif$contres$votant$mandatRef, vote='contre')
      votes = votes %>% bind_rows(pours) 
    }
  }
}





sumscr = votes %>% group_by(scrutin) %>% summarise(nbv = n(),p=mean(if_else(vote=="pour",1,0)))
vclean = votes %>% select(-groupe) %>% filter(!is.na(mandat))

Xv=spread(vclean,scrutin,vote)
iss = sumscr %>% filter(nbv>100)
ij = which(Xv[,iss$scrutin+1]=="pour",arr.ind = TRUE)
Xxv=sparseMatrix(ij[,1],ij[,2],x = rep(1,nrow(ij)))


meta=fromJSON(file="/home/come/Bureau/AMO10_deputes_actifs_mandats_actifs_organes_XV.json")


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


sol=greed(Xxv)
plot(sol)
pprobs = post_probs(new("mm"),Xxv,sol@cl)

clusters=data.frame(mandat = Xv[,1],cl=sol@cl, p=apply(pprobs,1,max))
mandats=votes %>% filter(!is.na(mandat)) %>% group_by(mandat,groupe) %>% summarise(nbv=n())

groupes_clust = clusters %>% left_join(mandats) %>% filter(p>0.95)

table(groupes_clust$cl,groupes_clust$groupe)
