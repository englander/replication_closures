rm(list=ls())
setwd("C:/Users/englander/Documents/replication_closures/")
library(dplyr); library(readxl); library(ggplot2)
library(sf); library(gridExtra); library(lwgeom)
library(grid); library(geosphere); library(viridis)
library(sp); library(cowplot); library(rworldmap)

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

options(scipen=999)

myThemeStuff <- theme(panel.background = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_line(color = 'black'),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.ticks = element_line(color = "gray5",size=.35),
                      axis.text = element_text(color = "black", size = 6.5, family="sans"),
                      axis.title = element_text(color = "black", size = 6.5, family = "sans"),
                      #axis.title.y.right = element_text(angle = 90,hjust=0),
                      axis.title.y = element_text(hjust = .5),
                      legend.key = element_blank(),
                      plot.title = element_text(hjust = 0.5, size = 9), 
                      legend.text=element_text(size=6.5, family = "sans"),
                      legend.title = element_text(size=6.5, family = "sans"),
                      plot.margin = unit(c(0,0,0,0),"in"),
                      plot.tag = element_text(family = "sans", size = 9),
                      plot.tag.position = c(.02, .99)
)

#Peru time
Sys.setenv(TZ='America/Lima')

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Make sf
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

fullbe <- st_sf(geometry = besf, fullbe)

fullbe <- rename(fullbe, calatime = FechaInicioCala)

rm(besf)

#Load potential closures created in 4. make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

myseason <- "s1_2019"

#Plot first five days of season
firstfive <- seq(from=min(rddf$start[rddf$season==myseason & rddf$bin=="active_in"],na.rm=T),
                 to=min(rddf$start[rddf$season==myseason & rddf$bin=="active_in"],na.rm=T)+4*24*3600,
                 by=24*3600)


#Grab active_in rectangles that begin on these five days
fiverect <- filter(rddf, bin=="active_in" & season==myseason)

#Bounding box for plot
bbox <- c(-76.7, -14.7, -75.5,-13.7) 
names(bbox) <- names(st_bbox(fiverect))

#Just peru
peru <- getMap(resolution='high') %>% st_as_sf()

peru <- filter(peru, ADMIN=="Peru")

#Crop
cropperu <- st_crop(peru, bbox)

#Crop potential closures and sets to this bbox as well
cropbe <- st_crop(fullbe, bbox)

croprect <- st_crop(fiverect, bbox)

#1. Plot raw BE data and resulting clusters

#Filter to BE observations between 24 and 9 hours before closure begins
#(sets that generated potential closure)
today <- filter(cropbe, (firstfive[1]-24*3600)<=calatime & 
                 calatime<=(firstfive[1]-9*3600)) 
  
#Make clusters
{
xy <- SpatialPointsDataFrame(
  matrix(c(today$lon,today$lat), ncol=2), data.frame(bepjhat = today$bepjhat),
  proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

# use the distm function to generate a geodesic distance matrix in meters
mdist <- distm(xy)

# cluster all points using a friend-of-friends approach
hc <- hclust(as.dist(mdist), method="single")

# define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
xy$clust <- cutree(hc, h=5*1.852*1000) #Nautical miles

#Convert xy to sf so can visualize
xy <- st_as_sf(xy)
xy$clust <- as.factor(xy$clust)

#Drop clusters made up of 3 or fewer points
numobs <- group_by(xy, clust) %>% 
  summarise(n = n())

numobs <- numobs$clust[numobs$n<=3]

xy <- filter(xy, clust %not in% numobs) 

#Also drop any cluster if it only has one unique point
distpts <- st_coordinates(xy) %>% as.data.frame() %>% 
  mutate(clust = xy$clust) %>% 
  distinct() %>% 
  group_by(clust) %>% 
  summarise(n = n())

if(1 %in% distpts$n){
  drop1 <- distpts$clust[distpts$n==1]
  xy <- filter(xy, clust %not in% drop1)
}

#Create convex polygon for each cluster
cpy <- group_by(xy, clust) %>% 
  summarize(geometry = st_union(geometry)) %>%
  st_convex_hull()

}

#Make plot area a data frame
bbsf <- data.frame(x = c(bbox["xmin"],bbox["xmax"],bbox["xmax"],bbox["xmin"]),
                   y=c(bbox["ymax"],bbox["ymax"],bbox["ymin"],bbox["ymin"]))

#Inset map
(insetmap <- ggplot() + 
  geom_sf(data=peru,fill='grey85',col='grey85') + 
  geom_polygon(data=bbsf,aes(x=x,y=y),fill=NA,col='black')+ 
  theme(panel.border = element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,0,.0),"in"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white')))




#Plot rectangles and declared closures
#Filter to active rectangles
myrect <- filter(croprect, start <= firstfive[1] & firstfive[1] <= end) 

#Create df of text for labeling potential closures
textdf <- data.frame(x=c(-76.57,-76.58),y=c(-13.75,-14.5),
                     label=c("Potential\nClosure 1","Potential Closure 2")) 

textdf$label <- as.character(textdf$label)

#Part a plot
(firstplot_present <- ggplot() + 
    #Plot closures that begin at midnight or at 6 am of next day
    geom_sf(data=today,alpha=.3,size=2,col='black',fill='black') +
    myThemeStuff+ 
    theme(panel.border = element_rect(color = 'black',fill=NA),
          legend.position = c(.87,.66), plot.margin = unit(c(0,0,0,0),"in"),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank()) + 
    geom_sf(data=cpy, col='red',fill=NA) +
    ggtitle("a. Sets between midnight and 3 pm on April 28, 2019\n") + 
    geom_sf(data=cropperu) + 
    scale_x_continuous(name=NULL,expand = c(0,0), limits = c(bbox["xmin"]-.05,bbox['xmax'])) + 
    scale_y_continuous(name=NULL,expand = c(0,0), limits = c(bbox["ymin"], bbox["ymax"])))


textdf$x[textdf$label=="Potential Closure 2"] <- -76.62

rectplot_present <-  ggplot() + 
    geom_sf(data=myrect, fill="dodgerblue2",col="dodgerblue2",alpha=.4,size=.15) + 
    myThemeStuff+ 
    ggtitle("b. Potential closures begin midnight April 29\nand end at 11:59 PM on May 1, 2019") + 
    geom_sf(data=cropperu) + 
    scale_x_continuous(name=NULL,expand = c(0,0), limits = c(bbox["xmin"]-.05,bbox['xmax'])) + 
    scale_y_continuous(name=NULL,expand = c(0,0),limits = c(bbox["ymin"], bbox["ymax"])) + 
    theme(plot.margin = unit(c(0,0,0,0),"in"),
          panel.border = element_rect(color = 'black',fill=NA),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank()) + 
    geom_text(data=textdf, aes(x=x,y=y,label=label),hjust=0,vjust=1,size=2.5)

#Add inset to second plot
rectplot_present <- ggdraw() +
  draw_plot(rectplot_present) + 
  draw_plot(insetmap, x = 0.57, y = 0.35, width = 0.5, height = 0.5) + 
  theme(plot.margin = unit(c(0,0,0,0),"in"))

stackplot_present <- plot_grid(firstplot_present, rectplot_present, ncol=2,nrow=1,rel_widths = c(1,1),rel_heights = c(1,1))


ggsave(stackplot_present, file="Output/Figures/figure5.png",
       w=6,h=3, units = "in", dpi=1200)

sessionInfo()
