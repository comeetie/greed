library(dplyr)
library(Matrix)
data("Marvel")
Heroes = Marvel %>% group_by(heroes) %>% summarise(v=n()) %>% arrange(desc(v))
Heroes_f = Heroes %>% filter(v>100)
Marvel_f =Marvel %>% filter(heroes %in% Heroes_f$heroes)
XMarvel.df=Marvel_f %>% rename(heroes1=heroes) %>% 
  left_join(Marvel_f %>% rename(heroes2=heroes)) %>% 
  group_by(heroes1,heroes2) %>% 
  summarise(v=n()) %>% filter(heroes1!=heroes2)



HLab = data.frame(name=unique(XMarvel.df$heroes1),stringsAsFactors = FALSE) %>% mutate(id=1:n())
XmarvelNum = XMarvel.df %>% left_join(HLab,by=c("heroes1"="name")) %>% rename(ifrom=id) %>%
  left_join(HLab,by=c("heroes2"="name")) %>% rename(ito=id)
XMarvel=sparseMatrix(XmarvelNum$ifrom,XmarvelNum$ito,x = XmarvelNum$v)
colnames(XMarvel)=HLab$name

devtools::use_data(XMarvel)

library(greed)

sol=fit(XMarvel,25,new("dcsbm"))
