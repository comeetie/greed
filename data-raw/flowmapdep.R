DEP2017_sf = loadMap(COG=2017,nivsupra="DEP")
xy = DEP2017_sf %>% st_geometry()%>%st_centroid() %>% st_coordinates()
dep_points = data.frame(x=xy[,1],y=xy[,2],DEP=DEP2017_sf$DEP,name=DEP2017_sf$nom,stringsAsFactors = FALSE)
ggplot(dep_points)+geom_text(aes(x=x,y=y,label=name))+coord_equal()

data("Xmigr.dep")
Xmigr.dep=Xmigr.dep+t(Xmigr.dep)
Xmigr.dep=tril(Xmigr.dep)
ij=which(Xmigr.dep>0,arr.ind = TRUE)
depname=colnames(Xmigr.dep)
gg=data.frame(i=depname[ij[,1]],j=depname[ij[,2]],f=Xmigr.dep[ij],stringsAsFactors = FALSE)

ggf = gg %>% left_join(dep_points,by=c("i"="DEP")) %>% rename(from=i,fromx=x,fromy=y) %>%
  left_join(dep_points,by=c("j"="DEP")) %>% rename(to=j,tox=x,toy=y)


ggplot(ggf%>%arrange(f)%>%filter(f>500))+
  geom_segment(aes(x=fromx,y=fromy,xend=tox,yend=toy,size=f,alpha=log(f)))+
  geom_point(data=dep_points,aes(x=x,y=y),size=1.5)+
  coord_equal()+
  theme_minimal()+
  scale_size_continuous("Flux",range=c(0,6),limits = c(0,35000))+
  scale_alpha_continuous(guide="none",range = c(0,1),limits=c(log(500),log(35000)))+
  scale_x_continuous("",breaks = c())+
  scale_y_continuous("",breaks = c())
  
