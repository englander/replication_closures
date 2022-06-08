rm(list=ls())

#Peru time
Sys.setenv(TZ='America/Lima')

options(scipen=999)

library(ggplot2); library(rworldmap); library(geosphere)
library(sf); library(dplyr); library(cowplot); library(purrr)
library(viridis); library(lubridate)

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

myThemeStuff <- theme(panel.background = element_rect(fill = NA),
                      panel.border = element_rect(fill = NA, color = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.ticks = element_line(color = "gray5",size=.35),
                      axis.text = element_text(color = "black", size = 5.5, family="sans"),
                      axis.title = element_text(color = "black", size = 6.5, family = "sans"),
                      #axis.title.y.right = element_text(angle = 90,hjust=0),
                      axis.title.y = element_text(hjust = .5),
                      legend.key = element_blank(),
                      plot.title = element_text(hjust = 0.5, size = 8), 
                      legend.text=element_text(size=6.5, family = "sans"),
                      legend.title = element_text(size=6.5, family = "sans"),
                      plot.margin = unit(c(0.01,0.01,0.01,0.01),"in"),
                      plot.tag = element_text(family = "sans", size = 9)
)


#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Make sf
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

besf <- st_sf(besf, fullbe)

#Shape of Peru land
bbox <- c(-84, -16, -74, -3.4) 
names(bbox) <- names(st_bbox(besf))

#First make a plot of South America to add as an inset for first season 2017 plot
southamerica <- getMap(resolution='low') %>% st_as_sf()

southamerica <- southamerica %>%
  filter(Stern=="South America")

southamerica <- st_union(southamerica)

#Crop out islands of South America
southamerica <- st_crop(southamerica, xmin = -82, xmax = -34, ymin=-56, ymax=13.5)

#Just peru
peru <- getMap(resolution='high') %>% st_as_sf()

peru <- filter(peru, ADMIN=="Peru")

#Crop Peru to bounding box
peru <- st_crop(peru,bbox)
southamerica <- st_crop(southamerica, bbox)

#Plot juvenile catch in thousands
besf <- mutate(besf, numtjuv = numjuv/10^3)

#Make an inset of zoomed in area
xmin <- -79; ymin <- -9; width <- .1
xmax <- xmin + width; ymax <- ymin + width 

croparea <- c(xmin, ymin, xmax, ymax)
names(croparea) <- names(st_bbox(besf))

#Crop besf to this area
croppedbesf <- st_crop(besf, croparea)

insetplot <- ggplot() + myThemeStuff + 
  geom_sf(data=croppedbesf, aes(fill=numtjuv,col=numtjuv),
          alpha=.5,size=.5) + 
  scale_x_continuous(name=NULL,expand = c(0,0),
                     breaks = c(xmin, (xmin+xmax)/2, xmax)) +
  scale_y_continuous(name=NULL,expand = c(0,0),
                     breaks = c(ymin, (ymin+ymax)/2, ymax)) +
  scale_color_viridis("Juvenile catch (thousands of juveniles)", trans='log', 
                      na.value = '#440154ff',
                      limits = c(min(besf$numtjuv[besf$numtjuv>0], na.rm=T), max(besf$numtjuv)),
                      breaks = c(1,20,400,8000),labels = c("1","20","400","8,000"),
                      guide = guide_colorbar(direction = "horizontal",
                                             position="bottom",
                                             title.position="top")) + 
  scale_fill_viridis("Juvenile catch (thousands of juveniles)", trans='log', 
                     breaks = c(1,20,400,8000),labels = c("1","20","400","8,000"),
                     na.value = '#440154ff',
                     limits = c(min(besf$numtjuv[besf$numtjuv>0], na.rm=T), max(besf$numtjuv)),
                     guide = guide_colorbar(direction = "horizontal",
                                            position="bottom",
                                            title.position="top")) + 
  guides(fill=FALSE,color=FALSE)  +
  theme(panel.background = element_rect(fill='white'),
        panel.border = element_rect(fill = NA, color = NA),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,0,.0),"in"))

redbox <- data.frame(x = c(croparea["xmin"],croparea["xmax"],croparea["xmax"],croparea["xmin"],croparea["xmin"]), 
                     y=c(croparea["ymax"],croparea["ymax"],croparea["ymin"],croparea["ymin"],croparea["ymax"]))


#Make plot of BE data next to Peru coast
beplot <- ggplot() + myThemeStuff + 
  geom_sf(data=southamerica,fill='grey93',col=NA) +
  geom_sf(data=peru,fill='grey85',col=NA) + 
  geom_sf(data=besf, aes(fill=numtjuv,col=numtjuv),
          alpha=.05,size=.15) + 
  scale_x_continuous(name=NULL,expand = c(0,0), 
                     breaks = seq(from=-83,to=-75,by=2)) + 
  scale_y_continuous(name=NULL,expand = c(0,0)) + 
  scale_color_viridis("Juvenile catch (thousands of juveniles)", trans='log', 
                      na.value = '#440154ff',
                      limits = c(min(besf$numtjuv[besf$numtjuv>0], na.rm=T), max(besf$numtjuv)),
                      breaks = c(1,20,400,8000),labels = c("1","20","400","8,000"),
                      guide = guide_colorbar(direction = "horizontal",
                                             position="bottom",
                                             title.position="top")) + 
  scale_fill_viridis("Juvenile catch (thousands of juveniles)", trans='log', 
                     breaks = c(1,20,400,8000),labels = c("1","20","400","8,000"),
                     na.value = '#440154ff',
                     limits = c(min(besf$numtjuv[besf$numtjuv>0], na.rm=T), max(besf$numtjuv)),
                     guide = guide_colorbar(direction = "horizontal",
                                            position="bottom",
                                            title.position="top")) + 
  theme(legend.position = "bottom", 
        legend.key.width=unit(1.5,"cm"),
        legend.title = element_text(vjust = 1, hjust = .5, size = 7.5)) + 
  geom_polygon(data=redbox,aes(x=x,y=y),fill=NA,col='darkorange1',size=.5) + 
  geom_polygon(data=data.frame(y=c(croparea["ymax"],-7.65),x=c(croparea["xmax"],-78.3)), aes(x=x,y=y),col='darkorange1') + 
  geom_polygon(data=data.frame(y=c(-7.65,-6.65),x=c(-78.3,-78.3)), aes(x=x,y=y),col='darkorange1') +
  geom_polygon(data=data.frame(y=c(-7.65,-7.65),x=c(-78.3,-77.3)), aes(x=x,y=y),col='darkorange1')


#Add inset map
beinset <- ggdraw() +
  draw_plot(beplot) +
  draw_plot(insetplot, x = .6, y = .65, width = 0.385, height = 0.385)

ggsave(beinset, file=paste0("Output/Figures/figure4.png"),
       w=4,h=6, units = "in", dpi = 900)

#Number of sets
nrow(besf)

#Number of unique vessels
unique(besf$Matricula) %>% length()

##Sets per day
nrow(besf) / unique(date(besf$FechaInicioCala)) %>% length()

#Average set catches 
mean(besf$numjuv,na.rm=T) #570103.1
mean(besf$betons) #50.22304
#Average number of adults
mean(besf$numadults, na.rm=T) #2508788

#Calculate distance between each set and the coast
disttocoast <- st_distance(besf, peru)

distnum <- as.numeric(disttocoast)
quantile(distnum, probs=seq(from=0,to=1,by=.01))
#95% of sets are within 79759.164 m of coast

sessionInfo()
