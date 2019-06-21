library(rjson)
library(dplyr)
library(tidyr)
library(Matrix)
library(future)
plan(multisession)
scrutins_json = fromJSON(file="./data-raw/Scrutins_XV.json")
scrutins = scrutins_json$scrutin$scrutin
votes = NULL
scr_info = NULL
for (is in 1:length(scrutins)){

  scr_info = scr_info %>% bind_rows(tibble(suid = scrutins[[is]]$uid,numero = scrutins[[is]]$numero, type = scrutins[[is]]$typeVote$libelleTypeVote, 
         seance = scrutins[[is]]$seanceRef, date = scrutins[[is]]$dateScrutin, titre = scrutins[[is]]$titre,
         demandeur = ifelse(!is.null(scrutins[[is]]$demandeur$texte),scrutins[[is]]$demandeur$texte,NA),
         objet = scrutins[[is]]$objet$libelle))
  
  groupes = scrutins[[is]]$ventilationVotes$organe$groupes$groupe
  for (g in 1:length(groupes)){
    #pours   =  tibble(scrutin = is, title=scrutins[[is]]$titre, groupe= groupes[[g]]$organeRef, mandat=sapply(groupes[[g]]$vote$decompteNominatif$pours$votant,function(v){v[2]}), vote='pour')
    #contres = tibble(scrutin = is, title=scrutins[[is]]$titre, groupe= groupes[[g]]$organeRef, mandat=sapply(groupes[[g]]$vote$decompteNominatif$contres$votant,function(v){v[2]}), vote='contres')
    #votes = votes %>% bind_rows(pours) %>% bind_rows(contres)
    if(as.numeric(groupes[[g]]$vote$decompteVoix$pour)>1){
      pours   =  tibble(scrutin = scrutins[[is]]$uid, groupe= groupes[[g]]$organeRef, acteur=sapply(groupes[[g]]$vote$decompteNominatif$pours$votant,function(v){v$acteurRef}),mandat=sapply(groupes[[g]]$vote$decompteNominatif$pours$votant,function(v){v$mandatRef}), vote='pour')
      votes = votes %>% bind_rows(pours) 
    }
    if(as.numeric(groupes[[g]]$vote$decompteVoix$pour)==1){
      pours   =  tibble(scrutin = scrutins[[is]]$uid, groupe= groupes[[g]]$organeRef, acteur=groupes[[g]]$vote$decompteNominatif$pours$votant$acteurRef, mandat=groupes[[g]]$vote$decompteNominatif$pours$votant$mandatRef, vote='pour')
      votes = votes %>% bind_rows(pours) 
    }
    if(as.numeric(groupes[[g]]$vote$decompteVoix$contre)>1){
    contres = tibble(scrutin = scrutins[[is]]$uid, groupe= groupes[[g]]$organeRef,acteur=sapply(groupes[[g]]$vote$decompteNominatif$contres$votant,function(v){v$acteurRef}), mandat=sapply(groupes[[g]]$vote$decompteNominatif$contres$votant,function(v){v$mandatRef}), vote='contre')
    votes = votes %>% bind_rows(contres)
    }
    if(as.numeric(groupes[[g]]$vote$decompteVoix$contre)==1){
      pours   =  tibble(scrutin = scrutins[[is]]$uid, groupe= groupes[[g]]$organeRef,  acteur=groupes[[g]]$vote$decompteNominatif$contres$votant$acteurRef,mandat=groupes[[g]]$vote$decompteNominatif$contres$votant$mandatRef, vote='contre')
      votes = votes %>% bind_rows(pours) 
    }
  }
}





sumscr = votes %>% group_by(scrutin) %>% summarise(nbv = n(),p=mean(if_else(vote=="pour",1,0)))
vclean = votes %>% select(-groupe,-mandat) 




meta=fromJSON(file="./data-raw/AMO30_tous_acteurs_tous_mandats_tous_organes_historique.json")

ij = which(Xv[,2:ncol(Xv)]=="pour",arr.ind = TRUE)

Xxv=sparseMatrix(ij[,1],ij[,2],x = rep(1,nrow(ij)))
rownames(Xxv)=unlist(Xv[,1])
colnames(Xxv)=sumscr$scrutin



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

dep_meta = dep_meta %>% semi_join(tibble(uid=rownames(Xxv)))

Xvlegislature = list(X=Xxv,row_meta=dep_meta,col_meta=scr_info)
