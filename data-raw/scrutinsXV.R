library(RJSONIO)
library(dplyr)
library(tidyr)
library(Matrix)
scrutins = fromJSON("./data-raw/Scrutins_XV.json")
scrutins = scrutins$scrutin$scrutin
votes = NULL
for (is in 1:length(scrutins)){
  groupes = scrutins[[is]]$ventilationVotes$organe$groupes$groupe
  for (g in 1:length(groupes)){
    pours   =  tibble(scrutin = is, groupe= groupes[[g]]$organeRef, mandat=sapply(groupes[[g]]$vote$decompteNominatif$pours$votant,function(v){v[2]}), vote='pour')
    contres = tibble(scrutin = is, groupe= groupes[[g]]$organeRef, mandat=sapply(groupes[[g]]$vote$decompteNominatif$contres$votant,function(v){v[2]}), vote='contres')
    votes = votes %>% bind_rows(pours) %>% bind_rows(contres)
  }
}





sumscr = votes %>% group_by(scrutin) %>% summarise(nbv = n(),p=mean(if_else(vote=="pour",1,0)))

vclean = votes %>% select(-groupe) %>% filter(!is.na(mandat))
Xv=spread(vclean,scrutin,vote)



iss = sumscr %>% filter(nbv>100)

ij = which(Xv[,iss$scrutin+1]=="pour",arr.ind = TRUE)

Xxv=sparseMatrix(ij[,1],ij[,2],x = rep(1,nrow(ij)))

sol=greed(Xxv)
plot(sol)
pprobs = post_probs(new("mm"),Xxv,sol@cl)

clusters=data.frame(mandat = Xv[,1],cl=sol@cl, p=apply(pprobs,1,max))
mandats=votes %>% filter(!is.na(mandat)) %>% group_by(mandat,groupe) %>% summarise(nbv=n())

groupes_clust = clusters %>% left_join(mandats) %>% filter(p>0.95)

table(groupes_clust$cl,groupes_clust$groupe)
