rm(list=ls())

#Peru time
Sys.setenv(TZ='America/Lima')

options(scipen=999)

library(ggplot2); library(rworldmap); library(geosphere)
library(sf); library(dplyr); library(cowplot); library(purrr)

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

myThemeStuff <- theme(panel.background = element_blank(),
                      panel.border = element_blank(),
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
                      plot.margin = unit(c(0,0,0.2,0),"in"),
                      plot.tag = element_text(family = "sans", size = 9)
)

#Load closed areas created in make_closures_df.R
load("Output/Data/closed.Rdata")

#Only want to use closures declared by PRODUCE from BE data
closed <- filter(closed, days <= 5)

#Filter to closures during study period
closed <- filter(closed, !is.na(season) & 
                   (season=="s1_2017" | season=="s2_2017" | 
                      season=="s1_2018" | season=="s2_2018" | season=="s1_2019" | 
                      season=="s2_2019")
)

#Crop to North-Central zone
bbox <- c(-84, -15.99999, -73, -3.5) 
names(bbox) <- names(st_bbox(closed))

closed <- st_crop(closed, bbox)

nrow(closed) #410

#Calculate area
areas <- st_area(closed)
areas <- as.numeric(areas) / 10^6

min(areas) #169.6286
max(areas) # 12512.12
mean(areas) #1328.252
median(areas) #1211.445

filter(closed, days==3) %>% nrow() / nrow(closed) #0.4829268
filter(closed, days==4) %>% nrow() / nrow(closed) #0.1536585
filter(closed, days==5) %>% nrow() / nrow(closed) #0.3609756

#What % of closures begin at 6 am? (during study period)
grep("\\.25",closed$start_raw) %>% length() / nrow(closed) #0.09268293

#How many closures are active per day on average during a fishing season? (during fishing seasons, excluding exploratory days)
datevec <- c(
  seq(from=as.POSIXct("2017-04-26 00:00:00 -05"),to=as.POSIXct("2017-07-31 00:00:00 -05"),by=3600*24),
  seq(from=as.POSIXct("2017-11-27 00:00:00 -05"),to=as.POSIXct("2018-01-25 00:00:00 -05"),by=3600*24),
  seq(from=as.POSIXct("2018-04-12 00:00:00 -05"),to=as.POSIXct("2018-08-08 00:00:00 -05"), by=3600*24),
  seq(from=as.POSIXct("2018-11-15 00:00:00 -05"),to=as.POSIXct("2019-04-03 00:00:00 -05"), by=3600*24),
  seq(from=as.POSIXct("2019-05-04 00:00:00 -05"), to=as.POSIXct("2019-07-30 00:00:00 -05"), by=3600*24),
  seq(from=as.POSIXct("2019-11-16 00:00:00 -05"),to=as.POSIXct("2020-01-14 00:00:00 -05"),by=3600*24)
)

nrow(closed) / length(datevec) #.73 closures on average per active day

##What is the distance to the nearest next closure within-season?
#Given communicado name, find next communicado within season. 
#For each closure in given communicado, find minimum distance to closures in next communicado
comDist <- function(comm){
  
  #Closures in given communicado
  curcl <- filter(closed, `Comunicado Name` == comm)
  
  #Next communicado
  nextcl <- filter(closed, season==unique(curcl$season)) %>% 
    mutate(nextcom = lead(`Comunicado Name`, 1)) %>%
    filter(`Comunicado Name`==comm)
  
  nextcl <- nextcl[nrow(nextcl),"nextcom"]
  
  nextcl <- filter(closed, `Comunicado Name`==nextcl$nextcom)
  
  #Minimum distance to next closure
  curcl$mindisttonext <- as.numeric(NA)

  if(nrow(nextcl)>0){
    for(i in 1:nrow(curcl)){
      
      mindist <- st_distance(curcl[i,],nextcl)
      
      mindist <- as.numeric(mindist)
      
      curcl$mindisttonext[i] <- min(mindist)
    }
  }
    
  out <- data.frame(closeid = curcl$closeid, mindisttonext = curcl$mindisttonext)
  
  return(out)
}

mindistdf <- map_df(unique(closed$`Comunicado Name`), function(x){
  comDist(x)
})

#What % of closures are right next to a closure announced in immediately preceding announcement?
#(within 100 m since sometimes distance is a few meters but in practice they are right next to each other)
#26%
filter(mindistdf, mindisttonext < 100) %>% nrow() / filter(mindistdf, !is.na(mindisttonext)) %>% nrow()

#Mean is 178 km and median is 92 km
summary(mindistdf$mindisttonext)


#What % of closure announcements contain multiple closures?
clperann <- as.data.frame(closed) %>% dplyr::select(-geometry) %>% 
  group_by(`Comunicado Name`) %>% 
  summarise(count = n())

#43%
filter(clperann, count > 1) %>% nrow() / nrow(clperann) #0.4285714

#Average closure announcement creates 1.58 closures
mean(clperann$count)

#Shape of Peru land
bbox <- c(-84, -16, -74, -3.4) 
names(bbox) <- names(st_bbox(closed))

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

#Make plot area a data frame
bbsf <- data.frame(x = c(-84,-74,-74,-84), y=c(-3.4,-3.4,-16,-16))

saplot <- ggplot() + geom_sf(data=southamerica,fill='grey93',col=NA) + 
  geom_sf(data=peru,fill='grey85',col='grey85') + 
    geom_polygon(data=bbsf,aes(x=x,y=y),fill=NA,col='black')+ 
  myThemeStuff + 
  theme(panel.border = element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,0,.0),"in"),
        panel.background = element_rect(fill='white'))

#Crop Peru to bounding box
peru <- st_crop(peru,bbox)
southamerica <- st_crop(southamerica, bbox)

#Determine correct proportions for closed area plots
myheight <- distGeo(c(-79,-3.4),c(-79,-16))
mywidth <- distGeo(c(-84,-9.7),c(-74,-9.7))


seasonPlot <- function(myseason,includeinset=NULL){
  
  #Title
  if(length(grep("s1",myseason))>0){
    seasontext <- paste0(
      "First season ",
      gsub("*s1_","",myseason)
    )
  } else{
    seasontext <- paste0(
      "Second season ",
      gsub("*s2_","",myseason)
    )
  }
  
  plot <- ggplot() + myThemeStuff + 
    geom_sf(data=southamerica,fill='grey93',col=NA) +
    geom_sf(data=peru,fill='grey85',col=NA) + 
    geom_sf(data=filter(closed, season==myseason),
            fill='red',alpha=.3,size=.15,col='red') + 
    scale_x_continuous(name=NULL,expand = c(0,0), limits = c(-82,-74),
                       breaks = seq(from=-81,to=-75,by=2)) + 
    scale_y_continuous(name=NULL,expand = c(0,0)) + 
    ggtitle(seasontext)
  
  return(plot)
}


plotlist <- lapply(unique(closed$season), function(x){
  seasonPlot(x)
})

##Add inset onto first season 2017
f17 <- ggdraw() +
  draw_plot(plotlist[[1]]) + 
    draw_plot(saplot, x = 0.545, y = 0.59, width = 0.34, height = 0.34)


#3 x 2 plot
tbt <- plot_grid(f17, plotlist[[2]], plotlist[[3]], plotlist[[4]],
                 plotlist[[5]], plotlist[[6]], ncol=3, nrow=2)

ggsave(tbt, file=paste0("Output/Figures/figure1.png"),
       w=7,h=(7/3)*(myheight/mywidth)*2, units = "in", dpi=1200)



##Make a similar plot for potential closures
#Load potential closures created in 4. make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

rddf <- filter(rddf, bin=="active_in")

rddf$season <- as.character(rddf$season)

seasonPlot <- function(myseason){
  
  #Title
  if(length(grep("s1",myseason))>0){
    seasontext <- paste0(
      "First season ",
      gsub("*s1_","",myseason)
    )
  } else{
    seasontext <- paste0(
      "Second season ",
      gsub("*s2_","",myseason)
    )
  }
  
  plot <- ggplot() + myThemeStuff + 
    geom_sf(data=southamerica,fill='grey93',col=NA) +
    geom_sf(data=peru,fill='grey85',col=NA) + 
    geom_sf(data=filter(rddf, season==myseason),
            fill='dodgerblue2',alpha=.4,size=.15,col='dodgerblue2') + 
    scale_x_continuous(name=NULL,expand = c(0,0), limits = c(-82,-74),
                       breaks = seq(from=-81,to=-75,by=2)) + 
    scale_y_continuous(name=NULL,expand = c(0,0)) + 
    ggtitle(seasontext)
  
  return(plot)
}


plotlist_potcl <- lapply(unique(closed$season), function(x){
  seasonPlot(x)
})


##Add inset onto first season 2017
p17 <- ggdraw() +
  draw_plot(plotlist_potcl[[1]]) + 
  draw_plot(saplot, x = 0.535, y = 0.58, width = 0.35, height = 0.35)


#3 x 2 plot
tbt <- plot_grid(p17, plotlist_potcl[[2]], plotlist_potcl[[3]], plotlist_potcl[[4]],
                 plotlist_potcl[[5]], plotlist_potcl[[6]], ncol=3, nrow=2)

ggsave(tbt, file="Output/Figures/figure4.png",
       w=7,h=(7/3)*(myheight/mywidth)*2, units = "in", dpi=1200)


sessionInfo()

