rm(list=ls())

#Peru time
Sys.setenv(TZ='America/Lima')

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

library(ggplot2); library(rworldmap); library(cowplot)
library(sf); library(dplyr); library(geosphere)

myThemeStuff <- theme(panel.background = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_line(color = 'black'),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.ticks = element_line(color = "gray5",size=.35),
                      axis.text = element_text(color = "black", size = 5.5, family="sans"),
                      axis.title = element_text(color = "black", size = 6.5, family = "sans"),
                      #axis.title.y.right = element_text(angle = 90,hjust=0),
                      axis.title.y = element_text(hjust = .5),
                      legend.key = element_blank(),
                      plot.title = element_text(hjust = 0.5, size = 7.5), 
                      legend.text=element_text(size=6.5, family = "sans"),
                      legend.title = element_text(size=6.5, family = "sans"),
                      plot.margin = unit(c(0,0,0,0),"in"),
                      plot.tag = element_text(family = "sans", size = 9)
)

#Median-sized potential closure
{
  #Load potential closures created in 4. make_rddf.R
  load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

  rddf <- filter(rddf, bin=="active_in")
  
  #To calculate areas, project
  medarea <- st_transform(rddf, st_crs("+proj=laea +lon_0=-76.5"))
  
  #Area
  medarea <- mutate(medarea, area_km2 = st_area(medarea)/10^6) %>% 
    mutate(area_km2 = as.numeric(area_km2)) %>%
    mutate(areadif = abs(area_km2 - median(area_km2)))
  
  #min, mean, median, and max area of potential closures (for comparison to actual closures)
  min(medarea$area_km2) 
  mean(medarea$area_km2) #956.7765
  median(medarea$area_km2) 
  max(medarea$area_km2) 
  
  #How many potential closures are active per day during the fishing season? (during fishing seasons, excluding exploratory days)
  datevec <- c(
    seq(from=as.POSIXct("2017-04-26 00:00:00 -05"),to=as.POSIXct("2017-07-31 00:00:00 -05"),by=3600*24),
    seq(from=as.POSIXct("2017-11-27 00:00:00 -05"),to=as.POSIXct("2018-01-25 00:00:00 -05"),by=3600*24),
    seq(from=as.POSIXct("2018-04-12 00:00:00 -05"),to=as.POSIXct("2018-08-08 00:00:00 -05"), by=3600*24),
    seq(from=as.POSIXct("2018-11-15 00:00:00 -05"),to=as.POSIXct("2019-04-03 00:00:00 -05"), by=3600*24),
    seq(from=as.POSIXct("2019-05-04 00:00:00 -05"), to=as.POSIXct("2019-07-30 00:00:00 -05"), by=3600*24),
    seq(from=as.POSIXct("2019-11-16 00:00:00 -05"),to=as.POSIXct("2020-01-14 00:00:00 -05"),by=3600*24)
  )
  
  nrow(medarea) / length(datevec) #1.718085
  
  #What is the average number of potential closures per two-week-of-sample by two-degree grid cell?
  nrow(medarea) / unique(medarea$twoweek_cellid_2p) %>% length() #4.384615
  
  
  medclosed <- as.data.frame(medarea) %>% dplyr::select(areadif) %>% 
    summarise(min(areadif)) %>% as.matrix() %>% as.numeric()
  
  #There is a tie, so randomly pick one
  set.seed(20200422)
  
  medclosed <- filter(medarea, areadif==medclosed) %>% 
    sample_n(1)
  
  ##Create annuli
  #Function of buffer distance
  buF <- function(myb){
    #dist given in m,  
    buffered <- st_buffer(medclosed, dist = myb*1000) %>%
      #but I want bdist variable to be in km
      mutate(bdist = myb)
    
    return(buffered)
  }
  
  #Apply over buffer distances
  buffers <- lapply(seq(from=10,to=50,by=10), function(x){
    buF(x)
  })
  
  buffers <- do.call("rbind", buffers)
  
  #Bind with medclosed
  medclosed <- rbind(mutate(medclosed, bdist = 0), buffers)
  
  #Iteratively difference to create annuli
  annuli <- dplyr::select(medclosed, bdist)
  
  for(i in 5:1){
    annuli[i+1,] <- st_difference(annuli[i+1,],annuli[i,]) %>% dplyr::select(bdist)
  }
  
  potcl_annuli <- annuli
  
  potcl_annuli <- st_transform(potcl_annuli, st_crs(rddf))
  
  rm(buffers, medarea, medclosed, i, rddf, buF, annuli)
}

#Side note: 89% of actual closures intersect potential closures
{
#Load potential closures created in 4. make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

rddf <- filter(rddf, bin=="active_in")

#Load closed areas. Created in make_closures_df.R
load("Output/Data/closed.Rdata")

#Filter to closures during study period
closed <- filter(closed, !is.na(season) & 
                   (season=="s1_2017" | season=="s2_2017" | 
                      season=="s1_2018" | season=="s2_2018" | season=="s1_2019" | 
                      season=="s2_2019")
)

bbox <- c(-84, -15.9999, -74, -3.5) 
names(bbox) <- names(st_bbox(closed))

closed <- st_crop(closed, bbox)

#Only plot closures declared by PRODUCE
closed <- filter(closed, days <= 5)

#Given row of closed, return 1 if it intersects a potential closure
closed$intersectspotcl <- 0

for(i in 1:nrow(closed)){
  
  row <- closed[i,]
  
  #Filter potential closures to same time period
  potcl <- filter(rddf, start <= row$end & end >= row$start)
  
  inter <- st_intersects(row, potcl)
  
  if(length(inter[[1]])>0){
    closed$intersectspotcl[i] <- 1
  }
  
}

mean(closed$intersectspotcl) #0.8902439

rm(i, row, potcl, inter)

#Also calculate area of overlap of each closure with potential closure(s)
#as a fraction (relative to area of actual closure)
closed$intersect_area_frac <- 0

#Since want areas, project both
rddf <- st_transform(rddf, st_crs("+proj=laea +lon_0=-76.5"))
closed <- st_transform(closed, st_crs("+proj=laea +lon_0=-76.5"))

for(i in 1:nrow(closed)){
  
  row <- closed[i,]
  
  #Filter potential closures to same time period
  potcl <- filter(rddf, start <= row$end & end >= row$start)
  
  inter <- st_intersection(row, potcl)
  
  if(nrow(inter) > 0){
    closed$intersect_area_frac[i] <- as.numeric(sum(st_area(inter)) / st_area(row))
    
  }
  
}

#intersect_area_frac exceeds 1 when overlapped by multiple potential closures
#In these cases, replace it with 1 because fraction cannot exceed 1
closed$intersect_area_frac[closed$intersect_area_frac > 1] <- 1

summary(closed$intersect_area_frac)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.2030  0.5001  0.4890  0.7500  1.0000 

rm(i, row, potcl, inter)

}


#Make a single plot of potential closure and rings
medpotcl <- ggplot() + geom_sf(data=filter(potcl_annuli, bdist==0), fill='dodgerblue2') + 
  geom_sf(data=filter(potcl_annuli, bdist!=0),fill=NA) + 
  myThemeStuff

#Plot without axes or title
noax <- medpotcl + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(margin=margin(0,0,0,0)))

##Try making a version of this plot with larger font and single row
unfilled <- ggplot() + geom_sf(data=filter(potcl_annuli, bdist==0), fill=NA) + 
  geom_sf(data=filter(potcl_annuli, bdist!=0),fill=NA) + 
  myThemeStuff + 
  theme(plot.title = element_text(hjust = 0.5, size = 11, margin=margin(0,0,0,0)),
        axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())

p1 <- unfilled + ggtitle("Announcement")
p2 <- noax + ggtitle("Closure period") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11))
p3 <- unfilled + ggtitle("1 day after") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11))
p4 <- unfilled + ggtitle("2 days after") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11))
p5 <- unfilled + ggtitle("3 days after") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11))
p6 <- unfilled + ggtitle("4 days after") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11))

tbt <- plot_grid(p1, p2, p3, p4, p5, p6, ncol=6, nrow=1)

ggsave(tbt, file=paste0("Output/Figures/figure5.pdf"),
       w=7,h=2.8, units = "in", dpi=1200)

sessionInfo()
