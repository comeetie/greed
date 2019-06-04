library(rjson)
library(dplyr)
library(tidyr)
library(Matrix)
scrutins = fromJSON(file="./data-raw/Scrutins_XV.json")
scrutins = scrutins$scrutin$scrutin
votes = NULL
for (is in 1:length(scrutins)){
  groupes = scrutins[[is]]$ventilationVotes$organe$groupes$groupe
  for (g in 1:length(groupes)){
    pours   =  tibble(scrutin = is, title=scrutins[[is]]$titre, groupe= groupes[[g]]$organeRef, mandat=sapply(groupes[[g]]$vote$decompteNominatif$pours$votant,function(v){v[2]}), vote='pour')
    contres = tibble(scrutin = is, title=scrutins[[is]]$titre, groupe= groupes[[g]]$organeRef, mandat=sapply(groupes[[g]]$vote$decompteNominatif$contres$votant,function(v){v[2]}), vote='contres')
    votes = votes %>% bind_rows(pours) %>% bind_rows(contres)
  }
}





sumscr = votes %>% group_by(scrutin) %>% summarise(nbv = n(),p=mean(if_else(vote=="pour",1,0)))

vclean = votes %>% select(-groupe) %>% filter(!is.na(mandat))
Xv=spread(vclean,scrutin,vote)


iss = sumscr %>% filter(nbv>100)
ij = which(Xv[,iss$scrutin+1]=="pour",arr.ind = TRUE)

Xxv=sparseMatrix(ij[,1],ij[,2],x = rep(1,nrow(ij)))
rownames(Xxv)=unlist(Xv[,1])

meta=fromJSON('./data-raw/AMO10_deputes_actifs_mandats_actifs_organes_XV.json',encoding = "iso-8859-1")
deputes = lapply(meta$export$acteurs$acteur,function(dep){tibble(uid=dep$uid[2],name = paste(dep$etatCivil$ident[1:3],collapse = " "),mandat=sapply(dep$mandats$mandat,function(mand){mand$uid}))})
deputes.df = do.call(rbind,deputes)

organisations = lapply(meta$export$organes$organe,function(org){tibble(uid=org$uid,name = org$libelle,name_abrv=org$libelleAbrege,type=org$codeType)})
organisations.df = do.call(rbind,organisations)

unique(organisations.df$type)
sol=greed(Xxv)
plot(sol)
pprobs = post_probs(new("mm"),list(X=Xxv),sol@cl)

clusters=data.frame(mandat = Xv[,1],cl=sol@cl, p=apply(pprobs,1,max))
mandats=votes %>% filter(!is.na(mandat)) %>% group_by(mandat,groupe) %>% summarise(nbv=n())
# ! un députés peut avoir plusieurs groupes
groupes_clust = clusters %>% inner_join(mandats,by=c("mandat"="mandat")) %>% 
  left_join(organisations.df,by=c("groupe"="uid")) %>% 
  left_join(deputes.df,by=c("mandat"="mandat")) %>% filter(!duplicated(paste(mandat)))

table(groupes_clust$cl,as.character(groupes_clust$name_abrv),exclude = TRUE)
plot(sol,type='front')
plot(sol,type='tree')
