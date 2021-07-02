library(readr)
library(dplyr)
# data from https://www.kaggle.com/stefanoleone992/fifa-20-complete-player-dataset?select=players_20.csv
fifa = read_csv("data-raw/players_20.csv")

X = fifa %>% dplyr::select(short_name,nationality,preferred_foot,team_position,age,height_cm,weight_kg,value_eur,pace,shooting,passing,dribbling,defending,physic) %>% 
  filter(!is.na(pace),!is.na(team_position)) %>% mutate(preferred_foot=factor(preferred_foot),team_position=factor(team_position)) 

library(greed)
sol=greed(X[,-c(1,2)],model=new("mmm"))
