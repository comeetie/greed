library(readr)
library(dplyr)
# data from https://www.kaggle.com/stefanoleone992/fifa-20-complete-player-dataset?select=players_20.csv
fifa = read_csv("data-raw/players_20.csv")

X = fifa %>% dplyr::select(short_name,nationality,preferred_foot,team_position,age,height_cm,weight_kg,value_eur,pace,shooting,passing,dribbling,defending,physic) %>% 
  filter(!is.na(pace),!is.na(team_position)) %>% mutate(preferred_foot=factor(preferred_foot),team_position=factor(team_position)) 

library(greed)

sol=greed(X[,-c(1:2)],model=new("mmm"))


team_pos=xml2::read_xml("/home/come/Bureau/team_position.svg")
team_pos_xy = tibble(team_position=xml_find_all(team_pos, "//svg:circle") %>% xml_attr("id"),
y=xml_find_all(team_pos, "//svg:circle") %>% xml_attr("cy") %>% as.numeric(),
x=xml_find_all(team_pos, "//svg:circle") %>% xml_attr("cx")%>% as.numeric()) %>% mutate(y=100-y)

library(ggplot2)

ggplot(team_pos_xy)+geom_text(aes(x=x,y=y,label=team_position))+coord_fixed(ratio=0.73)

res_cl =tibble(team_position=X$team_position,cl=sol@cl) %>% count(team_position,cl) %>% mutate(team_position=tolower(team_position))%>% left_join(team_pos_xy)


ggplot(res_cl %>% filter(cl==1) )+geom_point(aes(x=x,y=y,size=n))+geom_text(aes(x=x,y=y,label=team_position))+coord_fixed(ratio=0.73)
