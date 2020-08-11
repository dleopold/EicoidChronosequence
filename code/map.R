# Map of chronosequence sites

library(tidyverse)
library(magrittr)
library(ggmap)
library(ggrepel)
library(cowplot)

#Needs an active Google Maps API key
#register_google("##########################")
SHP <- rgdal::readOGR("data/Coastline-shp/Coastline.shp") %>% 
  fortify(SHP) 

points <- read.csv("data/Sites.csv")

map <- make_bbox(points$Long,points$Lat,f=0.2) %>%
  get_map(zoom=7, source = "google", maptype = "satellite", crop=T, force=T) 
map_attributes <- attributes(map)
map_transparent <- matrix(adjustcolor(map, alpha.f = 0.74), 
                                    nrow = nrow(map))
map_transparent <- matrix(colorspace::lighten(map,0.3), 
                          nrow = nrow(map))
attributes(map_transparent) <- map_attributes

set.seed(13)
ggmap(map_transparent) +
  #geom_point(data=points, aes(Long, Lat))+
  geom_polygon(data=SHP, aes(long, lat, group=group),fill=NA, color="black") +
  geom_label_repel(data=points, aes(Long, Lat, label = lab), size = 4.75,
                   box.padding = 0.8, fill="#ffffff99") +
  lims(y=c(18.75,22.4), x=c(-160.2,-154.2)) +
  labs(x="Longitude", y="Latitude") +
  geom_richtext(data=tibble(label=c("1 - Thurston<span style='font-size:11pt'> (3.0&times;10^2 yr)</span><br>2 - Olaa<span style='font-size:11pt'> (2.1&times;10^3 yr)</span><br>3 - Laupahoehoe<span style='font-size:11pt'> (2.0&times;10^4 yr)</span><br>4 - Kohala<span style='font-size:11pt'> (1.5&times;10^5 yr)</span><br>5 - Kokee<span style='font-size:11pt'> (4.1&times;10^6 yr)</span>")), 
                inherit.aes = F, hjust=0, vjust=0, 
                size=4.5, fill="#ffffff99",
                aes(x=-160.2,y=18.85, label=label))
ggsave("output/figs/map.pdf",width=14,height=10,units = "cm")  



