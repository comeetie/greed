library(readr)
library(tidyr)
library(jsonlite)
library(dplyr)
spec=cols(
  X1 = col_datetime(format = ""),
  X2 = col_character(),
  X3 = col_character(),
  X4 = col_character(),
  X5 = col_character(),
  X6 = col_character(),
  X7 = col_double(),
  X8 = col_double(),
  X9 = col_character(),
  X10 = col_character(),
  X11 = col_character(),
  X12 = col_double(),
  X13 = col_character()
)




meta = fromJSON("https://api.hearthstonejson.com/v1/latest/enUS/cards.json")
collect = meta[,c("id","collectible")]




library(lubridate)
deb = as.POSIXct("2019-04-09")
cur = Sys.time()
base="http://files.hearthscry.com/collectobot/"
pmonth_urls=paste0(year(deb),sprintf("-%02d.zip",month(deb):(month(cur)-1)))
cmonth_urls = paste0(year(cur),sprintf("-%02d",month(cur)),sprintf("-%02d.zip",1:(day(cur)-1)))
urls = c(pmonth_urls,cmonth_urls)


nbg = 1
for(u in urls){
  download.file(paste0(base,u),paste0("./data-raw/hs/",u))
  data=fromJSON(unzip(paste0("./data-raw/hs/",u)))
  drank = data$games[data$games$mode=="ranked",]
  for (i in 1:nrow(drank)){
    game = drank[i,"card_history"][[1]]
    if(nrow(game)>5){
      gg=data.frame(date=drank$added[i],hero = drank$hero[i], deck =drank$hero_deck[i], opp=drank$opponent[i],opp_deck=drank$opponent_deck[i], uid = drank$user_hash[i],rank = drank$rank[i],gid=nbg,player=game$player,card=game$card$name,result=drank$result[i],turn=game$turn,cardid=game$card$id)
      write.table(gg,"./data-raw/matchs2019S2.csv",append=TRUE,col.names=FALSE,row.names = FALSE,sep=',')
      nbg=nbg+1
    }
  }
}

cdat = read_csv("./data-raw/matchs2019S2.csv",col_names = FALSE,col_types = spec)
names(cdat)=c("date","hero","deck","opp","opp_deck","uid","rank","gid","player","card","result","turn","cardid")


Hunter=cdat %>% filter(date>=as.POSIXct("2019-04-09")) %>% left_join(collect,by=c("cardid"="id")) %>% 
  filter((hero=="Hunter" & player=="me")|(opp=="Hunter" & player=="opponent"),!is.na(collectible)) %>% 
  mutate(gdid = paste(gid,player,sep="_"),res=if_else(player=="me",result=="win",result=="loss")) %>% 
  select(gdid,card,cardid,res) %>% group_by(gdid,card,cardid,res) %>% summarise(n=n()) %>% ungroup()

Hmat=Hunter %>% select(gdid,cardid,n) %>% spread(cardid,n,fill=0)
Hmatnum= as.matrix(Hmat[,-1])

ij=which(Hmatnum>0,arr.ind = TRUE)
library(Matrix)
Hspmat = sparseMatrix(ij[,1],ij[,2],x=Hmatnum[ij])
fit=greed(Hspmat[,colSums(Hspmat)>25])


clusts = tibble(gdid=Hmat$gdid,cluster=fit@cl)

resultats=Hunter%>% group_by(gdid,res) %>% summarise(n=n()) %>% ungroup()

clusts %>% left_join(resultats) %>% group_by(cluster) %>% summarize(p=mean(res))

cards = Hunter %>% group_by(card,cardid) %>% summarise(n=sum(n))
lp=fit@obs_stats$x_counts[9,]
lpn = lp/sum(lp)*30
s=0.5

cn = colnames(Hmatnum)[colSums(Hspmat)>25]
dcards=tibble(cardid=cn[lpn>s],p=lpn[lpn>s],nbc=lp[lpn>s]) %>% left_join(cards) 
dcards



