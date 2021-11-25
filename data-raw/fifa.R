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

data("Fifa")
poscolnames = tolower(colnames(Fifa)[14:28])
poscolnames
data("Fifa_positions")

pos = Fifa_positions$positions[match(poscolnames, Fifa_positions$positions$team_position),]

Xp=matrix(rep(pos$x,each=nrow(Fifa)),nrow(Fifa))
Yp=matrix(rep(pos$y,each=nrow(Fifa)),nrow(Fifa))
Fposbin=as.matrix(Fifa[,toupper(poscolnames)] %>% mutate_all(\(f){if_else(f==TRUE,1,0)}))

mean_pos_x =rowSums(Fposbin*Xp)/rowSums(Fposbin)
mean_pos_y =rowSums(Fposbin*Yp)/rowSums(Fposbin)
Fifa$pos_y=mean_pos_y
Fifa$pos_x=mean_pos_x
