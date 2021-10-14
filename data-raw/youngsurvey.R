library(readr)
library(dplyr)
library(stringr)
library(careless)
library(ggrepel)
library(greed)
library(ggplot2)
#https://www.kaggle.com/miroslavsabo/young-people-survey
#survey=read_csv("data-raw/responses.csv")
data("Young_people_survey")

nc  = 76




selected = Young_people_survey %>% 
  select(1:nc) %>%
  mutate(string = longstring(.)) %>%
  mutate(sel = if_else(string <= 10,TRUE,FALSE) ) %>% pull(sel)

X = Young_people_survey %>% 
  select(1:nc) %>%
  filter(selected)  %>% 
  mutate_all(\(x){factor(x,levels=1:5)}) %>% 
  droplevels()

library(future)
plan(multisession)
sol=greed(X)
plot(sol,type='tree')
# plot(cut(sol,7))
params = coef(cut(sol,7))
means_scores = lapply(params$Thetak,\(x){apply(x,1,\(r){sum(r*as.numeric(names(r)))})})
means_scores_long = do.call(rbind,purrr::map2(means_scores,names(means_scores),\(x,y){tibble(cluster=names(x),mean=x,var=y)}))  %>% mutate(var = gsub("\\."," ",var))
means_scores_glob = interests %>% summarise_all(\(x){mean(x,na.rm=TRUE)}) %>% tidyr::pivot_longer(1:nc,names_to = "var",)%>% mutate(var = gsub("[,-/]"," ",var))
gg = means_scores_long %>% left_join(means_scores_glob) %>% mutate(dm=mean-value) 

cluster_profile = function(clid){
  df = gg %>% filter(cluster==paste0("cluster",clid)) %>% filter(abs(dm)>0.3) 
  ggplotGrob(ggplot(df)+
    geom_text_repel(aes(x=0,y=dm,size=abs(dm),color=dm,label=var),
                    direction="x",
                    min.segment.length = 800,    
                    force_pull    = 0,
                    max.time      = 0.5,
                    max.iter      = 1e5,
                    max.overlaps  = Inf,
                    segment.color = NA,
                    point.padding = NA)+
    scale_color_gradient2(guide="none")+
    scale_size_area(guide="none",max_size = 5)+
    theme_minimal()+
    scale_x_discrete("",breaks=c(),expand = c(0.2,1))+
    scale_y_continuous("cluster score difference with average",limits = c(-1.5,1.2))+
    ggtitle(paste("Cluster",clid),subtitle = paste(round(params$pi[clid]*100),"% of the respondents")))
}

plts_list = lapply(1:7,\(cli){cluster_profile(cli)})


grid::grid.draw(gridExtra::arrangeGrob(grobs=plts_list, nrow = 3,ncol =3))






