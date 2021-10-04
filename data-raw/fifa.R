library(readr)
library(dplyr)
library(stringr)
# data from https://www.kaggle.com/stefanoleone992/fifa-20-complete-player-dataset?select=players_20.csv
fifa = read_csv("./data-raw/players_20.csv")


positions = fifa$player_positions %>% str_remove_all(" ") %>% str_split(",") %>% unlist() %>% unique()
fifa=fifa %>% mutate(player_positions_vec = str_split(str_remove_all(player_positions," "),","))
create_pos_feature <- function(data, position) {
  position_var = enquo(position)
  data %>% mutate(!!position :=  factor(sapply(player_positions_vec,\(pvec){position %in% pvec}))) 
}


for( p in positions){
  fifa = create_pos_feature(fifa,p)
}



X = fifa%>% select(short_name,nationality,preferred_foot,age,height_cm,weight_kg,value_eur,pace,shooting,passing,dribbling,defending,physic,!!!syms(positions)) %>% 
  filter(!is.na(pace)) %>% mutate(preferred_foot=factor(preferred_foot)) 

library(greed)

data(Fifa)
Ns=5000
sol=greed(Fifa[sample(nrow(X),Ns),-c(1:2)],model=new("mmm"))
library(xml2)

team_pos=xml2::read_xml("/home/come/Bureau/team_position.svg")
team_pos_xy = tibble(team_position=xml_find_all(team_pos, "//svg:circle") %>% xml_attr("id"),
y=xml_find_all(team_pos, "//svg:circle") %>% xml_attr("cy") %>% as.numeric(),
x=xml_find_all(team_pos, "//svg:circle") %>% xml_attr("cx")%>% as.numeric()) %>% mutate(x=x-5.7452493,y=100-(y-18.658329))


pos_clust = data.frame(do.call(rbind,lapply(1:sol@K,\(k){
  sapply(params$Thetak[2:length(params$Thetak)],\(x){1-x[k,1]})
})))
pos_clust$cluster=factor(1:13)

pos_clust_long = pos_clust %>% pivot_longer(cols = -16,names_to = "position",values_to = "p") %>% 
  mutate(position=tolower(position)) %>% 
  left_join(team_pos_xy,by=c("position"="team_position"))


library(ggplot2)

ggplot(pos_clust_long %>% filter(cluster==13))+geom_point(aes(x=x,y=y,size=p))+coord_fixed(ratio=0.73)+scale_size_area("p",limits=c(0,1))

pos_clust_mean = pos_clust_long %>% group_by(cluster) %>% summarise(x=weighted.mean(x,p),y=weighted.mean(y,p))

img <- png::readPNG("/home/come/Bureau/foot_field.png")
ggplot(pos_clust_mean)+background_image(img)+geom_text(aes(x=x,y=y,label=cluster))+
  coord_fixed(ratio=1)+
  scale_x_continuous(limits=c(0,203.2),expand = c(0,1))+
  scale_y_continuous(limits=c(0,101.6),expand = c(0,1))+theme_void()


ggplot(team_pos_xy)+background_image(img)+geom_text(aes(x=x,y=y,label=team_position))+
  coord_fixed(ratio=1)+
  scale_x_continuous(limits=c(0,203.2),expand = c(0,1))+
  scale_y_continuous(limits=c(0,101.6),expand = c(0,1))


clust_means = data.frame(do.call(rbind,params$muk))
clust_means$cluster=factor(1:13)

clust_means_long = clust_means[,5:11] %>% pivot_longer(cols = -7,names_to = "feat",values_to = "score") 

ggplot(clust_means_long) + geom_line(aes(x=feat,y=score,group=cluster,color=cluster))+coord_polar()


ggplot(clust_means_long) + geom_point(aes(x=factor(feat,levels=c("physic","defending","pace","passing","dribling","shooting")),y=cluster,color=score-50,size=abs(score-50)),shape=43)+scale_color_distiller(palette = "RdBu",limits=c(-35,35))
