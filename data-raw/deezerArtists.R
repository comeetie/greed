library(dplyr)
library(purrr)
library(stringr)
library(Matrix)

con = file("./data-raw/playlists2.txt", "r")
ArtistsLinks = list()
i=0
while ( TRUE ) {
  i=i+1
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  if(i %% 100 ==0){
    print(i)
  }
  line_tab = str_split(line,',')[[1]]
  artists = unique(line_tab[seq(2,length(line_tab),by=2)])
  ArtistsLinks[[i]]=data.frame(from=rep(artists,length(artists)),to=rep(artists,each=length(artists)),stringsAsFactors = FALSE)
}

close(con)

ArtistsLinks.df = do.call(rbind,ArtistsLinks)
ArtistsLinksVal = ArtistsLinks.df %>% group_by(from,to) %>% summarize(v=n()) %>% filter(from!=to)
ArtistsScores = ArtistsLinksVal %>% group_by(from) %>% summarise(v=sum(v)) %>% rename(did=from)


ArtistsInfo = read.csv("./data-raw/Artists.csv") %>% mutate(did=paste0("artist_",id))
ArtistsFiltered = ArtistsScores %>% filter(v>500) %>% mutate(num=1:n()) %>% left_join(ArtistsInfo %>% select(did,name))

ArtistsLinksValNum = ArtistsLinksVal %>% left_join(ArtistsFiltered %>% select(did,num),by=c("from"="did")) %>% rename(from.num=num) %>% 
  left_join(ArtistsFiltered %>% select(did,num),by=c("to"="did")) %>% rename(to.num=num) %>% 
  filter(!is.na(from.num),!is.na(to.num))


X = sparseMatrix(ArtistsLinksValNum$from.num,ArtistsLinksValNum$to.num,x = ArtistsLinksValNum$v)
colnames(X) = ArtistsFiltered$name

XdeezerArtist = X
devtools::use_data(XdeezerArtist)
