rm(list=ls())

library(dplyr); library(readxl); library(ggplot2)
library(rworldmap); library(sf); library(lwgeom)
library(rgdal); library(geosphere); library(sp)
library(purrr); library(lubridate); library(glmnet)
library(lfe); library(Formula); library(smoothr)
library(parallel); library(xtable); library(latex2exp)
library(cowplot); library(tidyr)

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

options(scipen=999)
options(lfe.threads=24)

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

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
                      plot.margin = unit(c(0,0.04,0,0),"in"),
                      plot.tag = element_text(family = "sans", size = 9)
)

#Peru time
Sys.setenv(TZ='America/Lima')

#Load full BE from PRODUCE where I have imputed size distribution of non-SNP observations
#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

fullbe <- rename(fullbe, calatime = FechaInicioCala)

#Load closed areas created in make_closures_df.R
load("Output/Data/closed.Rdata")

#Only want to use closures declared by PRODUCE from BE data
closed <- filter(closed, days <= 5)

#Given closeday and distance threshold, create clusters
makeClusts <- function(closeday, distthresh){
  
  today <- filter(fullbe, (closeday-24*3600)<=calatime & 
                    calatime<=(closeday-9*3600))
  
  #Cannot have a cluster unless more than 1 point
  if(nrow(today)>1){
    xy <- SpatialPointsDataFrame(
      matrix(c(today$lon,today$lat), ncol=2), 
      data.frame(bepjhat = today$bepjhat, bepj = today$bepj, betons = today$betons, calatime=today$calatime, 
                 numindivids=today$numindivids,numjuv=today$numjuv,numadults=today$numadults,
                 prop3hat = today$prop3hat, prop3p5hat = today$prop3p5hat, prop4hat = today$prop4hat, prop4p5hat = today$prop4p5hat,
                 prop5hat = today$prop5hat, prop5p5hat = today$prop5p5hat, prop6hat = today$prop6hat, prop6p5hat = today$prop6p5hat,
                 prop7hat = today$prop7hat, prop7p5hat = today$prop7p5hat, prop8hat = today$prop8hat, prop8p5hat = today$prop8p5hat,
                 prop9hat = today$prop9hat, prop9p5hat = today$prop9p5hat, prop10hat = today$prop10hat, prop10p5hat = today$prop10p5hat,
                 prop11hat = today$prop11hat, prop11p5hat = today$prop11p5hat, prop12hat = today$prop12hat, prop12p5hat = today$prop12p5hat,
                 prop13hat = today$prop13hat, prop13p5hat = today$prop13p5hat, prop14hat = today$prop14hat, prop14p5hat = today$prop14p5hat,
                 prop15hat = today$prop15hat, prop15p5hat = today$prop15p5hat, prop16hat = today$prop16hat, prop16p5hat = today$prop16p5hat,
                 prop17hat = today$prop17hat, prop17p5hat = today$prop17p5hat, prop18hat = today$prop18hat, prop18p5hat = today$prop18p5hat),
      proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
    
    # use the distm function to generate a geodesic distance matrix in meters
    mdist <- distm(xy)
    
    hc <- hclust(as.dist(mdist), method="single")
    
    # define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
    xy$clust <- cutree(hc, h=distthresh*1.852*1000) #Nautical miles
    
    #Convert xy to sf so can visualize
    xy <- st_as_sf(xy)
    
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
    
    
    if(nrow(xy)>0){
      #Create convex polygon for each cluster
      cpy <- group_by(xy, clust) %>% 
        summarize(geometry = st_union(geometry)) %>%
        st_convex_hull()
      
      out <- list(xy, cpy)
      names(out) <- c("xy","cpy")
    } else{
      out <- "NA"
    }
  }
  else{
    out <- "NA"
  }
  return(out)
}

#Round a point to nearest 5-minute interval
#Want to round right and top coordinates up
roundUp <- function(num){
  num <- num*12
  num <- ceiling(num)
  num <- num / 12
}

#Want to round left and bottom coordinates down
roundDown <- function(num){
  num <- num*12
  num <- floor(num)
  num <- num / 12
}

#Small function to make a rectangular polygon from 4 points
Rect <- function(top, bottom, right, left){
  
  lt <- c(left, top)
  rt <- c(right, top)
  rb <- c(right, bottom)
  lb <- c(left, bottom)
  
  #Make a polygon
  poly <- rbind(lt, rt, rb, lb, lt) %>%
    list() %>%
    st_polygon() %>%
    st_sfc(crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% st_sf()
  
  names(poly)[length(poly)] <- "geometry"
  st_geometry(poly) <- "geometry"
  
  poly <- st_transform(poly, st_crs("+proj=laea +lon_0=-76.5"))
  
  return(poly)
}


#Round bounding box up to nearest 5-minute interval and draw rectangles
#Iteratively union intersecting rectangles
drawRects <- function(closeday, distthresh){
  
  #Make clusters
  clusts <- makeClusts(closeday, distthresh)
  
  #Draw rectangles, rounding each coordinate to nearest 5-minute interval
  if(length(clusts)>1){
    myrects <- lapply(1:nrow(clusts[["cpy"]]), function(x){
      
      Rect(
        roundUp(st_bbox(clusts[["cpy"]][x,])["ymax"]),
        roundDown(st_bbox(clusts[["cpy"]][x,])["ymin"]),
        roundUp(st_bbox(clusts[["cpy"]][x,])["xmax"]),
        roundDown(st_bbox(clusts[["cpy"]][x,])["xmin"])
      ) %>% 
        st_transform(st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
      
    })
    
    myrects <- do.call("rbind", myrects)
    
    #If any rectangles have 0 area, drop them and corresponding convex hull
    if(0 %in% st_area(myrects)=="TRUE"){
      #Which rectangles have 0 area
      zeroarea <- which(as.numeric(st_area(myrects))==0)
      myrects <- myrects[-zeroarea,]
      clusts[["cpy"]] <- clusts[["cpy"]][-zeroarea,]
    }
    
    #Union intersecting rectangles
    unrect <- st_union(myrects)
    
    unrect <- st_cast(unrect, "POLYGON") %>% 
      st_sf()
    
    #Draw new rectangles
    newrects <- lapply(1:nrow(unrect), function(x){
      
      Rect(
        st_bbox(unrect[x,])["ymax"],
        st_bbox(unrect[x,])["ymin"],
        st_bbox(unrect[x,])["xmax"],
        st_bbox(unrect[x,])["xmin"]
      ) %>% 
        st_transform(st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
      
    })
    
    newrects <- do.call("rbind", newrects)
    
    #Re-make cluster variable
    newclust <- st_intersects(clusts[["cpy"]], newrects)
    
    out <- clusts[["xy"]] %>% 
      mutate(newclust = as.integer(NA))
    
    for(i in 1:length(newclust)){
      
      oldclust <- clusts[["cpy"]]$clust[i]
      
      #In one case, clusts[["cpy"]][i,] is a linestring and does not intersect its newrects counterpart
      #(clust 8 on closeday=="2019-12-11 -05"). In that case, give clusts[["cpy"]] a tiny buffer so it 
      #intersects with the right newrects
      if(length(newclust[[i]])==0){
        newclust[[i]] <- st_intersects(st_buffer(clusts[["cpy"]][i,],.0001),newrects)[[1]]
      }
      
      out$newclust[out$clust==oldclust] <- newclust[[i]]
    }
    
    #Output underlying BE values and drawn closures
    out <- list(out, 
                mutate(newrects, newclust = 1:nrow(newrects))
    )
    names(out) <- c("be","rects")
  } else{
    out <- "NA"
  }
  return(out)
  
}

#Compute % of individuals in each size bin for a cluster
numindSize <- function(clustdat){
  
  out <- as.data.frame(clustdat) %>% dplyr::select(-geometry) %>%
    #Filter out NAs. These are observations that have location, tons, and % juv,
    #but they occurred in a twoweek by two degree grid cell that did not have SNP observations
    #that I could use to impute the resulting size distribution
    filter(!is.na(prop12hat)) %>%
    summarise(prop3hat = sum(prop3hat*numindivids)/sum(numindivids),
              prop3p5hat = sum(prop3p5hat*numindivids)/sum(numindivids),
              prop4hat = sum(prop4hat*numindivids)/sum(numindivids),
              prop4p5hat = sum(prop4p5hat*numindivids)/sum(numindivids),
              prop5hat = sum(prop5hat*numindivids)/sum(numindivids),
              prop5p5hat = sum(prop5p5hat*numindivids)/sum(numindivids),
              prop6hat = sum(prop6hat*numindivids)/sum(numindivids),
              prop6p5hat = sum(prop6p5hat*numindivids)/sum(numindivids),
              prop7hat = sum(prop7hat*numindivids)/sum(numindivids),
              prop7p5hat = sum(prop7p5hat*numindivids)/sum(numindivids),
              prop8hat = sum(prop8hat*numindivids)/sum(numindivids),
              prop8p5hat = sum(prop8p5hat*numindivids)/sum(numindivids),
              prop9hat = sum(prop9hat*numindivids)/sum(numindivids),
              prop9p5hat = sum(prop9p5hat*numindivids)/sum(numindivids),
              prop10hat = sum(prop10hat*numindivids)/sum(numindivids),
              prop10p5hat = sum(prop10p5hat*numindivids)/sum(numindivids),
              prop11hat = sum(prop11hat*numindivids)/sum(numindivids),
              prop11p5hat = sum(prop11p5hat*numindivids)/sum(numindivids),
              prop12hat = sum(prop12hat*numindivids)/sum(numindivids),
              prop12p5hat = sum(prop12p5hat*numindivids)/sum(numindivids),
              prop13hat = sum(prop13hat*numindivids)/sum(numindivids),
              prop13p5hat = sum(prop13p5hat*numindivids)/sum(numindivids),
              prop14hat = sum(prop14hat*numindivids)/sum(numindivids),
              prop14p5hat = sum(prop14p5hat*numindivids)/sum(numindivids),
              prop15hat = sum(prop15hat*numindivids)/sum(numindivids),
              prop15p5hat = sum(prop15p5hat*numindivids)/sum(numindivids),
              prop16hat = sum(prop16hat*numindivids)/sum(numindivids),
              prop16p5hat = sum(prop16p5hat*numindivids)/sum(numindivids),
              prop17hat = sum(prop17hat*numindivids)/sum(numindivids),
              prop17p5hat = sum(prop17p5hat*numindivids)/sum(numindivids),
              prop18hat = sum(prop18hat*numindivids)/sum(numindivids),
              prop18p5hat = sum(prop18p5hat*numindivids)/sum(numindivids))
  
  return(out)
}


#Compute summary statistics and quantiles for each cluster
sumandQuants <- function(closeday, distthresh){
  
  #Make clusters and draw rectangles
  mymake <- drawRects(closeday, distthresh)
  
  if(length(mymake)>1){
    be <- mymake[["be"]]
    
    #Compute % of individuals in each size bin
    persize <- map_df(unique(be$newclust), function(x){
      mysize <- numindSize(
        filter(be, newclust==x)) %>% 
        mutate(newclust = x)
    })
    
    #total tons caught, number of sets in cluster, number of individuals caught in cluster,
    #and mean percentage juveniles in cluster, using (bepjhat or bepj) X (weighted by tons or unweighted)
    out <- left_join(persize,
                     as.data.frame(be) %>% 
                       dplyr::select(newclust, betons, numindivids, numjuv, numadults, bepjhat, bepj) %>% 
                       mutate(weightpj = bepj*numindivids, weightpjhat = bepjhat*numindivids) %>%
                       group_by(newclust) %>%
                       summarise(tons = sum(betons,na.rm=T),
                                 nobs = n(),
                                 numindivids=sum(numindivids,na.rm=T),
                                 numjuv=sum(numjuv,na.rm=T),
                                 numadults=sum(numadults,na.rm=T),
                                 meanpjhat = mean(bepjhat, na.rm=T),
                                 meanpj = mean(bepj, na.rm=T),
                                 meanpjhat_weighted = sum(weightpjhat,na.rm=T)/sum(numindivids,na.rm=T),
                                 meanpj_weighted = sum(weightpj,na.rm=T)/sum(numindivids,na.rm=T)),
                     by='newclust'
    )
    
    #Add closeday to out
    out <- mutate(out, closeday = closeday)
    
    #Join this df onto rectangles
    out <- left_join(mymake[["rects"]], out, by = 'newclust')
  } else{
    out <- "NA"
  }
  return(out)
}

#Apply over all days in season; not just days that had a closure
#Start 24 hours after first day in season with fishing
datevec <- c(
  seq(from=as.POSIXct("2017-04-23 00:00:00 -05"),to=as.POSIXct("2017-07-31 00:00:00 -05"),by=3600*24),
  seq(from=as.POSIXct("2017-11-24 00:00:00 -05"),to=as.POSIXct("2018-01-25 00:00:00 -05"),by=3600*24),
  seq(from=as.POSIXct("2018-04-07 00:00:00 -05"),to=as.POSIXct("2018-08-08 00:00:00 -05"), by=3600*24),
  seq(from=as.POSIXct("2018-11-15 00:00:00 -05"),to=as.POSIXct("2019-04-03 00:00:00 -05"), by=3600*24),
  seq(from=as.POSIXct("2019-04-29 00:00:00 -05"), to=as.POSIXct("2019-07-30 00:00:00 -05"), by=3600*24),
  #Last day of BE data I have is January 7
  seq(from=as.POSIXct("2019-11-07 00:00:00 -05"),to=as.POSIXct("2020-01-06 00:00:00 -05"),by=3600*24)
)



#Parallel apply over datevec
(myCores <- detectCores())

cl <- makeCluster(myCores - 10)

clusterExport(cl, "fullbe")
clusterExport(cl, "datevec")
clusterExport(cl, "drawRects")
clusterExport(cl, "makeClusts")
clusterExport(cl, "Rect")
clusterExport(cl, "roundDown")
clusterExport(cl, "roundUp")
clusterExport(cl, "numindSize")
clusterExport(cl, "sumandQuants")
clusterExport(cl, "%not in%")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(rgdal))
clusterEvalQ(cl, library(geosphere))
clusterEvalQ(cl, library(sp))
clusterEvalQ(cl, library(purrr))


rdflist <- parLapply(cl = cl,
                     datevec,
                     function(x){
                       
                       sumandQuants(x, 5)
                       
                     })

stopCluster(cl)
rm(cl, myCores)

#Some days have no clusters so need to be dropped before row binding
dropel <- sapply(1:length(rdflist), function(x){
  ifelse("data.frame" %in% class(rdflist[[x]]),1,0)
})
dropel <- which(dropel==0)

rdf <- rdflist[-dropel]
rdf <- do.call("rbind", rdf)

#Areas of rectangles I've created
areas <- st_area(rdf) 

rdf <- mutate(rdf, area_km2 = as.numeric(areas)/10^6)

#Create season indicator
rdf <- mutate(rdf, season = "s1_2017") %>%
  rename(start = closeday)

rdf$season[rdf$start>=as.POSIXct("2017-11-24 00:00:00 -05") & 
             rdf$start<=as.POSIXct("2018-01-25 00:00:00 -05")] <- "s2_2017"

rdf$season[rdf$start>=as.POSIXct("2018-04-07 00:00:00 -05") & 
             rdf$start<=as.POSIXct("2018-08-08 00:00:00 -05")] <- "s1_2018"

rdf$season[rdf$start>=as.POSIXct("2018-11-15 00:00:00 -05") & 
             rdf$start<=as.POSIXct("2019-04-03 00:00:00 -05")] <- "s2_2018"

rdf$season[rdf$start>=as.POSIXct("2019-04-29 00:00:00 -05") & 
             rdf$start<=as.POSIXct("2019-07-30 00:00:00 -05")] <- "s1_2019"

rdf$season[rdf$start>=as.POSIXct("2019-11-07 00:00:00 -05") & 
             rdf$start<=as.POSIXct("2020-01-06 00:00:00 -05")] <- "s2_2019"

#Filter to closures during season
closed <- filter(closed, !is.na(season) & 
                   (season=="s1_2017" | season=="s2_2017" | 
                      season=="s1_2018" | season=="s2_2018" | season=="s1_2019" | 
                      season=="s2_2019")
)

#Calculate minimum area each season
areas <- st_area(closed) 

closed <- mutate(closed, area_km2 = as.numeric(areas)/10^6)

minareas <- as.data.frame(closed) %>% 
  dplyr::select(season, area_km2) %>% 
  group_by(season) %>% 
  summarise(minarea = min(area_km2))

#Drop rectangles if less than minimum closed area that season
rdf <- left_join(rdf, minareas, by = 'season') %>%
  filter(area_km2 >= minarea)

rm(areas, dropel)

closed$days <- as.character(closed$days) %>% as.integer()

#Set end for rectangles equal to 5 days
rdf <- mutate(rdf, end = start + 3600*24*5-1)

#Moving forward in time, modify area if intersecting with already-created rectangle
#First season 2017
for(i in unique(rdf$start[rdf$season=="s1_2017"])[-1]){
  
  #Rectangles that have already been created
  already <- filter(rdf, start < i & i < end)
  
  if(nrow(already)>0){
    
    #This day's rectangles
    thisday <- filter(rdf, start==i)
    
    #Intersections with this day's rectangles
    inter <- st_intersects(thisday, already)
    
    for(k in 1:length(inter)){
      if(length(inter[[k]])>0){
        
        #Difference with first overlapping already existing rectangle
        diff <- st_difference(st_set_precision(st_buffer(st_make_valid(thisday[k,]), dist = 0), 1000000),
                              st_set_precision(st_buffer(st_make_valid(already[inter[[k]][1],]), dist = 0), 1000000)) #Don't want any hanging lines
        
        #Iteratively difference if more than one already existing rectangle
        if(length(inter[[k]])>1){
          for(j in 2:length(inter[[k]])){
            diff <- st_difference(st_set_precision(st_buffer(st_make_valid(diff), dist = 0), 1000000), 
                                  st_set_precision(st_buffer(st_make_valid(already[inter[[k]][j],]), dist = 0), 1000000))
          }
        }
        
        #Drop polygon if it is smaller than minimum closed area created this season
        diff <- drop_crumbs(diff, units::set_units(minareas$minarea[minareas$season=="s1_2017"],km^2))
        
        #Which row does this today's rectangle correspond to in rdf?
        whichrow <- which(rdf$start==i & rdf$newclust==thisday$newclust[k])
        
        #If thisday's rectangle is already completed overlapped, drop it from rdf
        if(nrow(diff)==0){
          rdf <- rdf[-whichrow,]
        }
        
        if(nrow(diff)>0){
          #Is this new polygon convex? Yes if area of convex hull within 1% of polygon
          ch <- st_convex_hull(diff) %>% st_area() #Area of convex hull
          polyarea <- st_area(diff)
          
          #Drop if not
          if(
            ((as.numeric(ch)/10^6 - as.numeric(polyarea)/10^6) / (as.numeric(polyarea)/10^6)) > .01){
            rdf <- rdf[-whichrow,]
          } else{
            #Otherwise replace geometry in rdf (if spatial overlap is partial)
            rdf$geometry[whichrow] <- diff$geometry
          }
        }
      }
    }
    
  }
}

rm(diff, inter, thisday, already, i, k, whichrow, j, ch, polyarea)

#Second season 2017
for(i in unique(rdf$start[rdf$season=="s2_2017"])[-1]){
  
  #Rectangles that have already been created
  already <- filter(rdf, start < i & i < end)
  
  if(nrow(already)>0){
    
    #This day's rectangles
    thisday <- filter(rdf, start==i)
    
    #Intersections with this day's rectangles
    inter <- st_intersects(thisday, already)
    
    for(k in 1:length(inter)){
      if(length(inter[[k]])>0){
        
        #Difference with first overlapping already existing rectangle
        diff <- st_difference(st_set_precision(st_buffer(st_make_valid(thisday[k,]), dist = 0), 1000000),
                              st_set_precision(st_buffer(st_make_valid(already[inter[[k]][1],]), dist = 0), 1000000)) #Don't want any hanging lines
        
        #Iteratively difference if more than one already existing rectangle
        if(length(inter[[k]])>1){
          for(j in 2:length(inter[[k]])){
            diff <- st_difference(st_set_precision(st_buffer(st_make_valid(diff), dist = 0), 1000000), 
                                  st_set_precision(st_buffer(st_make_valid(already[inter[[k]][j],]), dist = 0), 1000000))
          }
        }
        
        #Drop polygon if it is smaller than minimum closed area created this season
        diff <- drop_crumbs(diff, units::set_units(minareas$minarea[minareas$season=="s2_2017"],km^2))
        
        #Which row does this today's rectangle correspond to in rdf?
        whichrow <- which(rdf$start==i & rdf$newclust==thisday$newclust[k])
        
        #If thisday's rectangle is already completed overlapped, drop it from rdf
        if(nrow(diff)==0){
          rdf <- rdf[-whichrow,]
        }
        
        if(nrow(diff)>0){
          #Is this new polygon convex? Yes if area of convex hull within 1% of polygon
          ch <- st_convex_hull(diff) %>% st_area() #Area of convex hull
          polyarea <- st_area(diff)
          
          #Drop if not
          if(
            ((as.numeric(ch)/10^6 - as.numeric(polyarea)/10^6) / (as.numeric(polyarea)/10^6)) > .01){
            rdf <- rdf[-whichrow,]
          } else{
            #Otherwise replace geometry in rdf (if spatial overlap is partial)
            rdf$geometry[whichrow] <- diff$geometry
          }
        }
      }
    }
    
  }
}

rm(diff, inter, thisday, already, i, k, whichrow, j, ch, polyarea)


#First season 2018
for(i in unique(rdf$start[rdf$season=="s1_2018"])[-1]){
  
  #Rectangles that have already been created
  already <- filter(rdf, start < i & i < end)
  
  if(nrow(already)>0){
    
    #This day's rectangles
    thisday <- filter(rdf, start==i)
    
    #Intersections with this day's rectangles
    inter <- st_intersects(thisday, already)
    
    for(k in 1:length(inter)){
      if(length(inter[[k]])>0){
        
        #Difference with first overlapping already existing rectangle
        diff <- st_difference(st_set_precision(st_buffer(st_make_valid(thisday[k,]), dist = 0), 1000000),
                              st_set_precision(st_buffer(st_make_valid(already[inter[[k]][1],]), dist = 0), 1000000)) #Don't want any hanging lines
        
        #Iteratively difference if more than one already existing rectangle
        if(length(inter[[k]])>1){
          for(j in 2:length(inter[[k]])){
            diff <- st_difference(st_set_precision(st_buffer(st_make_valid(diff), dist = 0), 1000000), 
                                  st_set_precision(st_buffer(st_make_valid(already[inter[[k]][j],]), dist = 0), 1000000))
          }
        }
        
        #Drop polygon if it is smaller than minimum closed area created this season
        diff <- drop_crumbs(diff, units::set_units(minareas$minarea[minareas$season=="s1_2018"],km^2))
        
        #Which row does this today's rectangle correspond to in rdf?
        whichrow <- which(rdf$start==i & rdf$newclust==thisday$newclust[k])
        
        #If thisday's rectangle is already completed overlapped, drop it from rdf
        if(nrow(diff)==0){
          rdf <- rdf[-whichrow,]
        }
        
        if(nrow(diff)>0){
          #Is this new polygon convex? Yes if area of convex hull within 1% of polygon
          ch <- st_convex_hull(diff) %>% st_area() #Area of convex hull
          polyarea <- st_area(diff)
          
          #Drop if not
          if(
            ((as.numeric(ch)/10^6 - as.numeric(polyarea)/10^6) / (as.numeric(polyarea)/10^6)) > .01){
            rdf <- rdf[-whichrow,]
          } else{
            #Otherwise replace geometry in rdf (if spatial overlap is partial)
            rdf$geometry[whichrow] <- diff$geometry
          }
        }
      }
    }
    
  }
}

rm(diff, inter, thisday, already, i, k, whichrow, j, ch, polyarea)

#Second season 2018
for(i in unique(rdf$start[rdf$season=="s2_2018"])[-1]){
  
  #Rectangles that have already been created
  already <- filter(rdf, start < i & i < end)
  
  if(nrow(already)>0){
    
    #This day's rectangles
    thisday <- filter(rdf, start==i)
    
    #Intersections with this day's rectangles
    inter <- st_intersects(thisday, already)
    
    for(k in 1:length(inter)){
      if(length(inter[[k]])>0){
        
        #Difference with first overlapping already existing rectangle
        diff <- st_difference(st_set_precision(st_buffer(st_make_valid(thisday[k,]), dist = 0), 1000000),
                              st_set_precision(st_buffer(st_make_valid(already[inter[[k]][1],]), dist = 0), 1000000)) #Don't want any hanging lines
        
        #Iteratively difference if more than one already existing rectangle
        if(length(inter[[k]])>1){
          for(j in 2:length(inter[[k]])){
            diff <- st_difference(st_set_precision(st_buffer(st_make_valid(diff), dist = 0), 1000000), 
                                  st_set_precision(st_buffer(st_make_valid(already[inter[[k]][j],]), dist = 0), 1000000))
          }
        }
        
        #Drop polygon if it is smaller than minimum closed area created this season
        diff <- drop_crumbs(diff, units::set_units(minareas$minarea[minareas$season=="s2_2018"],km^2))
        
        #Which row does this today's rectangle correspond to in rdf?
        whichrow <- which(rdf$start==i & rdf$newclust==thisday$newclust[k])
        
        #If thisday's rectangle is already completed overlapped, drop it from rdf
        if(nrow(diff)==0){
          rdf <- rdf[-whichrow,]
        }
        
        if(nrow(diff)>0){
          #Is this new polygon convex? Yes if area of convex hull within 1% of polygon
          ch <- st_convex_hull(diff) %>% st_area() #Area of convex hull
          polyarea <- st_area(diff)
          
          #Drop if not
          if(
            ((as.numeric(ch)/10^6 - as.numeric(polyarea)/10^6) / (as.numeric(polyarea)/10^6)) > .01){
            rdf <- rdf[-whichrow,]
          } else{
            #Otherwise replace geometry in rdf (if spatial overlap is partial)
            rdf$geometry[whichrow] <- diff$geometry
          }
        }
      }
    }
    
  }
}

rm(diff, inter, thisday, already, i, k, whichrow, j, ch, polyarea)

#First season 2019
for(i in unique(rdf$start[rdf$season=="s1_2019"])[-1]){
  
  #Rectangles that have already been created
  already <- filter(rdf, start < i & i < end)
  
  if(nrow(already)>0){
    
    #This day's rectangles
    thisday <- filter(rdf, start==i)
    
    #Intersections with this day's rectangles
    inter <- st_intersects(thisday, already)
    
    for(k in 1:length(inter)){
      if(length(inter[[k]])>0){
        
        #Difference with first overlapping already existing rectangle
        diff <- st_difference(st_set_precision(st_buffer(st_make_valid(thisday[k,]), dist = 0), 1000000),
                              st_set_precision(st_buffer(st_make_valid(already[inter[[k]][1],]), dist = 0), 1000000)) #Don't want any hanging lines
        
        #Iteratively difference if more than one already existing rectangle
        if(length(inter[[k]])>1){
          for(j in 2:length(inter[[k]])){
            diff <- st_difference(st_set_precision(st_buffer(st_make_valid(diff), dist = 0), 1000000), 
                                  st_set_precision(st_buffer(st_make_valid(already[inter[[k]][j],]), dist = 0), 1000000))
          }
        }
        
        #Drop polygon if it is smaller than minimum closed area created this season
        diff <- drop_crumbs(diff, units::set_units(minareas$minarea[minareas$season=="s1_2019"],km^2))
        
        #Which row does this today's rectangle correspond to in rdf?
        whichrow <- which(rdf$start==i & rdf$newclust==thisday$newclust[k])
        
        #If thisday's rectangle is already completed overlapped, drop it from rdf
        if(nrow(diff)==0){
          rdf <- rdf[-whichrow,]
        }
        
        if(nrow(diff)>0){
          #Is this new polygon convex? Yes if area of convex hull within 1% of polygon
          ch <- st_convex_hull(diff) %>% st_area() #Area of convex hull
          polyarea <- st_area(diff)
          
          #Drop if not
          if(
            ((as.numeric(ch)/10^6 - as.numeric(polyarea)/10^6) / (as.numeric(polyarea)/10^6)) > .01){
            rdf <- rdf[-whichrow,]
          } else{
            #Otherwise replace geometry in rdf (if spatial overlap is partial)
            rdf$geometry[whichrow] <- diff$geometry
          }
        }
      }
    }
    
  }
}

rm(diff, inter, thisday, already, i, k, whichrow, j, ch, polyarea)

#Second season 2019
for(i in unique(rdf$start[rdf$season=="s2_2019"])[-1]){
  
  #Rectangles that have already been created
  already <- filter(rdf, start < i & i < end)
  
  if(nrow(already)>0){
    
    #This day's rectangles
    thisday <- filter(rdf, start==i)
    
    #Intersections with this day's rectangles
    inter <- st_intersects(thisday, already)
    
    for(k in 1:length(inter)){
      if(length(inter[[k]])>0){
        
        #Difference with first overlapping already existing rectangle
        diff <- st_difference(st_set_precision(st_buffer(st_make_valid(thisday[k,]), dist = 0), 1000000),
                              st_set_precision(st_buffer(st_make_valid(already[inter[[k]][1],]), dist = 0), 1000000)) #Don't want any hanging lines
        
        #Iteratively difference if more than one already existing rectangle
        if(length(inter[[k]])>1){
          for(j in 2:length(inter[[k]])){
            diff <- st_difference(st_set_precision(st_buffer(st_make_valid(diff), dist = 0), 1000000), 
                                  st_set_precision(st_buffer(st_make_valid(already[inter[[k]][j],]), dist = 0), 1000000))
          }
        }
        
        #Drop polygon if it is smaller than minimum closed area created this season
        diff <- drop_crumbs(diff, units::set_units(minareas$minarea[minareas$season=="s2_2019"],km^2))
        
        #Which row does this today's rectangle correspond to in rdf?
        whichrow <- which(rdf$start==i & rdf$newclust==thisday$newclust[k])
        
        #If thisday's rectangle is already completed overlapped, drop it from rdf
        if(nrow(diff)==0){
          rdf <- rdf[-whichrow,]
        }
        
        if(nrow(diff)>0){
          #Is this new polygon convex? Yes if area of convex hull within 1% of polygon
          ch <- st_convex_hull(diff) %>% st_area() #Area of convex hull
          polyarea <- st_area(diff)
          
          #Drop if not
          if(
            ((as.numeric(ch)/10^6 - as.numeric(polyarea)/10^6) / (as.numeric(polyarea)/10^6)) > .01){
            rdf <- rdf[-whichrow,]
          } else{
            #Otherwise replace geometry in rdf (if spatial overlap is partial)
            rdf$geometry[whichrow] <- diff$geometry
          }
        }
      }
    }
    
  }
}

rm(diff, inter, thisday, already, i, k, whichrow, j, ch, polyarea)

#Calculate new areas
areas <- st_area(rdf) %>% as.numeric()/10^6
rdf <- mutate(rdf, area_km2 = areas) %>% 
  mutate(logarea = log(area_km2))
rm(areas)

#Drop rectangles if less than minimum closed area that season
rdf <- filter(rdf, area_km2 >= minarea)

#Create bins
rdf <- rdf %>% mutate(bin = "active_in",bdist=0,tvar=0)

#Create buffer bins
BufFun <- function(rowind, bmin, bmax){
  
  row <- rdf[rowind, ]
  
  #Project 
  row <- st_transform(row, st_crs("+proj=laea +lon_0=-76.5"))
  
  #If invalid make it valid
  if(st_is_valid(row)==FALSE){
    row <- st_make_valid(row)
  }
  
  #Buffer (in km)
  buf <- st_buffer(row, bmax*1000)
  
  #Subtract inside
  if(bmin==0){
    inbuf <- row
  } else{
    inbuf <- st_buffer(row, bmin*1000)
  }
  
  out <- st_difference(buf, inbuf)
  
  #Drop all duplicate columns
  out <- out[,-grep("\\.1",names(out))]
  
  #Replace bdist column with buffer used
  out$bdist <- bmax
  
  #Also correct bin
  out$bin <- paste0(gsub("_.*","",out$bin),"_",bmax)
  
  #Transform back into lon lat
  out <- st_transform(out, st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  return(out)
}

#Function to make buffer for each rectangle
applyBufFun <- function(bmin, bmax){
  out <- lapply(1:nrow(rdf), function(x){
    BufFun(x, bmin, bmax)
  })
  
  out <- do.call("rbind",out)
  
  return(out)
}

#Apply over all buffers
(myCores <- detectCores())

cl <- makeCluster(myCores - 10)

clusterExport(cl, "rdf")
clusterExport(cl, "BufFun")
clusterExport(cl, "applyBufFun")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lwgeom))

#Try creating buffers
rbufs <- parLapply(cl=cl,
                   list(
                     c(0,10),c(10,20),c(20,30),c(30,40),c(40,50)                   
                   ), function(x){
                     applyBufFun(bmin=x[[1]],bmax=x[[2]])
                   })

stopCluster(cl)
rm(cl, myCores)

rbufs <- do.call("rbind",rbufs)

rdf <- rbind(rdf, rbufs)

rm(rbufs)

#Duplicate rectangles for leads and lags
rdf <- rbind(
  rdf,
  rdf %>% mutate(end=start-1) %>% 
    mutate(start=start-9*3600,tvar=-1),
  rdf %>% mutate(start = end+1,end=end+24*3600*1,tvar=1),
  rdf %>% mutate(start = end+1+24*3600*1,end=end+24*3600*2,tvar=2),
  rdf %>% mutate(start = end+1+24*3600*2,end=end+24*3600*3,tvar=3),
  rdf %>% mutate(start = end+1+24*3600*3,end=end+24*3600*4,tvar=4)
)

#Remake bin variable
rdf$bin[rdf$tvar==-1] <- gsub("active","lead9hours",rdf$bin[rdf$tvar==-1])
rdf$bin[rdf$tvar==1] <- gsub("active","lag1",rdf$bin[rdf$tvar==1])
rdf$bin[rdf$tvar==2] <- gsub("active","lag2",rdf$bin[rdf$tvar==2])
rdf$bin[rdf$tvar==3] <- gsub("active","lag3",rdf$bin[rdf$tvar==3])
rdf$bin[rdf$tvar==4] <- gsub("active","lag4",rdf$bin[rdf$tvar==4])

#Have to make corresponding closed buffers and time lags
closed <- mutate(closed, bin = "active_in",bdist=0,tvar=0)

#Make closed buffers
BufFun_closed <- function(rowind, bmin, bmax){
  
  row <- closed[rowind, ]
  
  #Project 
  row <- st_transform(row, st_crs("+proj=laea +lon_0=-76.5"))
  
  #If invalid make it valid
  if(st_is_valid(row)==FALSE){
    row <- st_make_valid(row)
  }
  
  #Buffer (in km)
  buf <- st_buffer(row, bmax*1000)
  
  #Subtract inside
  if(bmin==0){
    inbuf <- row
  } else{
    inbuf <- st_buffer(row, bmin*1000)
  }
  
  out <- st_difference(buf, inbuf)
  
  #Drop all duplicate columns
  out <- out[,-grep("\\.1",names(out))]
  
  #Replace bdist column with buffer used
  out$bdist <- bmax
  
  #Also correct bin
  out$bin <- paste0(gsub("_.*","",out$bin),"_",bmax)
  
  #Transform back into lon lat
  out <- st_transform(out, st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  return(out)
}

#Function to make buffer for each rectangle
applyBufFun_closed <- function(bmin, bmax){
  out <- lapply(1:nrow(closed), function(x){
    BufFun_closed(x, bmin, bmax)
  })
  
  out <- do.call("rbind",out)
  
  return(out)
}

#Apply over all buffers
(myCores <- detectCores())

cl <- makeCluster(myCores - 10)

clusterExport(cl, "closed")
clusterExport(cl, "BufFun_closed")
clusterExport(cl, "applyBufFun_closed")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lwgeom))

cbufs <- parLapply(cl=cl,
                   list(
                     c(0,10),c(10,20),c(20,30),c(30,40),c(40,50)                   
                   ), function(x){
                     applyBufFun_closed(bmin=x[[1]],bmax=x[[2]])
                   })

stopCluster(cl)
rm(cl, myCores)

cbufs <- do.call("rbind",cbufs)

names(cbufs)[names(cbufs)!="geometry"] <- names(closed)[names(closed)!="geometry"]

closed <- rbind(closed, cbufs)

rm(cbufs)

#Indicator for whether closure begins at 6 am
closed$six <- 0
closed$six[grep("\\.25",closed$start_raw)] <- 1

#Duplicate rectangles for leads and lags
closed <- rbind(
  closed,
  #Preperiod is 9 hours before if begins at midnight; 12 hours before if begins at 6 am
  rbind(
    filter(closed, six==0) %>% mutate(end=start-1) %>% 
      mutate(start=start-9*3600,tvar=-1),
    filter(closed, six==1) %>% mutate(end=start-1) %>% 
      mutate(start=start-12*3600,tvar=-1)
  ),
  closed %>% mutate(start = end+1,end=end+24*3600*1,tvar=1),
  closed %>% mutate(start = end+1+24*3600*1,end=end+24*3600*2,tvar=2),
  closed %>% mutate(start = end+1+24*3600*2,end=end+24*3600*3,tvar=3),
  closed %>% mutate(start = end+1+24*3600*3,end=end+24*3600*4,tvar=4))

#Remake bin variable
closed$bin[closed$tvar==-1] <- gsub("active","lead9hours",closed$bin[closed$tvar==-1])
closed$bin[closed$tvar==1] <- gsub("active","lag1",closed$bin[closed$tvar==1])
closed$bin[closed$tvar==2] <- gsub("active","lag2",closed$bin[closed$tvar==2])
closed$bin[closed$tvar==3] <- gsub("active","lag3",closed$bin[closed$tvar==3])
closed$bin[closed$tvar==4] <- gsub("active","lag4",closed$bin[closed$tvar==4])


#If any rectangles are invalid, make them valid
if(st_is_valid(rdf) %>% sum < nrow(rdf)){
  rdf <- lapply(1:nrow(rdf), function(x){
    if(st_is_valid(rdf[x,])==FALSE){
      st_make_valid(rdf[x,])
    } else{
      rdf[x,]
    }
  })
  
  rdf <- do.call("rbind", rdf)
  
}

#Calculate rectangle overlap with actual closures
#Function of start and bin
treatVar <- function(mystart, mybin){
  
  #Rectangles
  myrect <- filter(rdf, start==mystart & bin==mybin) %>% 
    mutate(treatfrac = 0) #Space-time fraction overlapping with actual closed area
  
  #When these rectangles end
  myend <- unique(myrect$end)
  
  #Filter to active closures
  #Only care about closures declared by PRODUCE
  actclosed <- filter(closed, days<=5 & bin==mybin & 
                        start <= myend & 
                        mystart <= end)
  
  if(nrow(myrect)>0 & nrow(actclosed)>0){
    
    inter <- st_intersects(myrect, actclosed)
    
    for(i in 1:length(inter)){
      if(length(inter[[i]])>0){
        for(j in inter[[i]]){
          #Intersection area of rect i with each actual closed area j
          myintersection <- st_intersection(myrect[i,],actclosed[j,]) %>%
            st_area()
          #Fraction of rectangle area
          areafrac <- as.numeric(myintersection) / 
            st_area(myrect[i,]) %>% as.numeric()
          if(areafrac>1){
            print("areafrac > 1")
          }
          #Now calculate time overlap as well. Fraction of overlapping hours
          recthours <- seq(from=myrect$start[i],to=myrect$end[i]+1,by=3600)
          acthours <- seq(from=actclosed$start[j],to=actclosed$end[j]+1,by=3600)
          timefrac <- sum(recthours %in% acthours) / length(recthours)
          
          #Add to treatfrac
          myrect$treatfrac[i] <- myrect$treatfrac[i] + areafrac*timefrac
        }
      }
    }
  }
  return(myrect)
}

#Apply over all start values for a bin
treatBin <- function(mybin){
  mylist <- lapply(unique(rdf$start[rdf$bin==mybin]), function(x){
    
    treatVar(x, mybin)
  })
  
  out <- do.call("rbind",mylist)
  
  return(out)
}

#Apply over all bins
(myCores <- detectCores())

cl <- makeCluster(myCores - 10)

clusterExport(cl, "rdf")
clusterExport(cl, "closed")
clusterExport(cl, "treatVar")
clusterExport(cl, "treatBin")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lubridate))

rdf <- parLapply(cl = cl,
                 unique(rdf$bin),
                 function(x){
                   
                   treatBin(x)
                   
                 })

stopCluster(cl)
rm(cl, myCores)

rdf <- do.call("rbind",rdf)

rddf <- mutate(rdf, startdate = as.Date(start) %>% as.factor(),
               season = as.factor(season))

#Create unique identifier for each rectangle
rid <- as.data.frame(rddf) %>% 
  dplyr::select(newclust, tons, nobs, area_km2, season) %>%
  distinct()

rid <- mutate(rid, rid = 1:nrow(rid))

#Join onto rddf
rddf <- left_join(rddf, rid)

rm(rid)

rddf <- dplyr::select(rddf, -minarea, -logarea)

#Calculate centroid for each unique closure
centr <- filter(rddf, bin=="active_in") %>% 
  st_centroid()

#Calculate centroid distance to coast
peru <- getMap(resolution = "high")
peru <- st_as_sf(peru) %>% filter(ADMIN.1 == "Peru")

dist2coast <- st_distance(centr, peru)

kmtocoast <- as.numeric(dist2coast) / 10^3

centr <- mutate(centr, kmtocoast = kmtocoast)

#Create 2 degree grid, created in 1. match_be_landings.R
load("Output/Data/grid2p.Rdata")

#Which cell is each centroid inside of
insidecell <- st_intersects(centr, grid2p)

insidecell <- as.numeric(as.character(insidecell))

#Add cellid column onto centr
centr <- mutate(centr, cellid = insidecell)

#Join cellid onto rddf
rddf <- left_join(rddf, 
                  as.data.frame(centr) %>% dplyr::select(rid, cellid, kmtocoast),
                  by='rid')

rm(insidecell, centr, grid2p, kmtocoast, dist2coast, peru)

#Create two-week of sample variable. 
twoweek <- as.data.frame(rddf) %>% dplyr::select(rid, bin, start)

twoweek <- mutate(twoweek, week = week(start), year = year(start))

twowk <- distinct(twoweek, week, year) %>% arrange(year, week)

twowk$twowk <- as.integer(NA); twowk$twowk[1] <- 1
counter <- 1
for(i in 2:nrow(twowk)){
  if((counter%%2==1 & 
      ((twowk$year[i]==twowk$year[i-1] & 
        twowk$week[i]==twowk$week[i-1]+1) | 
       twowk$year[i]==twowk$year[i-1]+1 & 
       twowk$week[i]==1 & twowk$week[i-1]>=52
      ))){
    twowk$twowk[i] <- twowk$twowk[i-1]
  } else{
    twowk$twowk[i] <- twowk$twowk[i-1]+1
  }
  counter <- counter+1
  #If gap, reset counter
  if(twowk$week[i]!=twowk$week[i-1]+1 & 
     twowk$week[i]!=1 & twowk$week[i-1]!=53){
    counter <- 1
  }
}

#Join
twoweek <- left_join(twoweek, twowk, by = c("week","year"))

#Join onto rddf
rddf <- left_join(rddf, dplyr::select(twoweek, rid, bin, twowk), by = c('rid','bin'))

rddf <- rename(rddf, cellid_2p = cellid)

#Now create variable to cluster my standard errors on
rddf <- mutate(rddf,twoweek_cellid_2p = paste0(twowk, "_", cellid_2p))

#Clean up
rm(twoweek, twowk, counter, i)

##Calculate catch outcomes inside each element of rddf
#Make bedat an sf object
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

besf <- st_sf(besf, fullbe)

outcomesFun <- function(rdrow){
  
  row <- rddf[rdrow,]
  
  #Filter BE observations to within same time period
  mybesf <- filter(besf, row$start<=calatime & calatime<row$end)
  
  if(nrow(mybesf)>0){
    
    #BE observations that are spatially inside element of rddf
    inter <- st_intersects(mybesf, row)
    
    inter <- as.numeric(as.character(inter))  
    
    if(length(which(!is.na(inter))>0)){
      tons <- sum(mybesf$betons[which(!is.na(inter))],na.rm=T)
      sdtons <- sd(mybesf$betons[which(!is.na(inter))],na.rm=T)
      numindivids <- sum(mybesf$numindivids[which(!is.na(inter))],na.rm=T)
      numjuv <- sum(mybesf$numjuv[which(!is.na(inter))],na.rm=T)
      tonsjuv <- sum(mybesf$tonsjuv[which(!is.na(inter))],na.rm=T)
      tonsadult <- sum(mybesf$tonsadult[which(!is.na(inter))],na.rm=T)
      nobs <- nrow(mybesf[which(!is.na(inter)),])
      uniquevessels <- length(unique(mybesf$Matricula[which(!is.na(inter))]))
      avgweightg <- (tons / numindivids) * 10^6 #convert to g
    } else{
      tons <- 0; numindivids <- 0; numjuv <- 0; nobs <- 0; sdtons = as.numeric(NA)
      uniquevessels <- 0; tonsjuv <- 0; tonsadult <- 0; avgweightg <- as.numeric(NA)
    }
  } else{
    tons <- 0; numindivids <- 0; numjuv <- 0; nobs <- 0; sdtons = as.numeric(NA)
    uniquevessels <- 0; tonsjuv <- 0; tonsadult <- 0; avgweightg <- as.numeric(NA)
  }
  
  #Output df with bin and rid so I can join onto rddf
  out <- data.frame(bin=row$bin, rid=row$rid,tons=tons,
                    #Round number of individuals to integer
                    numindivids=round(numindivids), numjuv=round(numjuv),
                    tonsjuv = tonsjuv, tonsadult = tonsadult,
                    nobs = nobs, sdtons = sdtons, uniquevessels = uniquevessels,
                    avgweightg = avgweightg)
  
  return(out)
}

#Apply over all rows
(myCores <- detectCores())

cl <- makeCluster(myCores - 10)

clusterExport(cl, "rddf")
clusterExport(cl, "besf")
clusterExport(cl, "outcomesFun")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lubridate))

#Apply over rows of rddf
myoutcomes <- parLapply(cl = cl,
                        1:nrow(rddf),
                        function(x){
                          
                          outcomesFun(x)
                          
                        })

stopCluster(cl)
rm(cl, myCores)

myoutcomes <- do.call("rbind",myoutcomes)

#Need to first rename other tons variables for clarity
rddf <- rename(rddf, clusttons = tons, 
               clustnumindivids = numindivids, clustnumjuv = numjuv, 
               clustnumadults = numadults, clustnobs = nobs, 
               clustarea_km2 = area_km2)

#Join onto rddf
rddf <- left_join(rddf, myoutcomes, by = c("bin",'rid'))

#Calculate number of caught adults in bin
rddf <- mutate(rddf, numadults = numindivids - numjuv)


##Now move to estimating effect of policy
rddf <- as.data.frame(rddf) %>% dplyr::select(-geometry)

rddf <- arrange(rddf, tvar, bdist)

#Millions of juveniles
rddf <- mutate(rddf, nummjuv = numjuv/10^6) %>%
  #Tons caught per set
  mutate(tonsperset = tons/nobs) %>%
  #Cluster tons/set and tons/area
  mutate(clusttonsperset = clusttons/clustnobs, clusttonsperarea = clusttons/clustarea_km2) %>%
  #clust millions of juveniles, individuals, and adults
  mutate(clustmjuv = clustnumjuv/10^6, clustmindivids = clustnumindivids/10^6,
         clustmadults = clustnumadults/10^6)

rddf <- mutate(rddf, asinhnummjuv = asinh(nummjuv), 
               asinhtons = asinh(tons))

rddf$bin <- as.factor(rddf$bin)
rddf$bin <- relevel(rddf$bin, ref="active_in")

rddf$twoweek_cellid_2p <- as.factor(rddf$twoweek_cellid_2p)
rddf$twowk <- as.factor(rddf$twowk)
rddf$cellid_2p <- as.factor(rddf$cellid_2p)

#Drop potential closures that have NA for size distribution
rddf <- filter(rddf, !is.na(prop12hat))

#Given variable, interact it with bin indicators, giving
interVars <- function(var){

  mydf <- rddf
  names(mydf)[names(mydf)==var] <- "myvar"

  #Want to manually interact var with bin indicator so I can look at each bin's coefficient
  #relative to 0 (rather than relative to omitted category)
  bininds <- model.matrix(~bin,data=rddf) %>% as.data.frame()

  #Drop intercept column and manually create active_in indicator column
  bininds <- dplyr::select(bininds, -`(Intercept)`)

  bininds <- bind_cols(
    dplyr::select(mydf, bin) %>% mutate(binactive_in = if_else(bin=="active_in",1,0)) %>%
      dplyr::select(-bin),
    bininds
  )

  #Create instrument columns
  outdf <- lapply(names(bininds), function(x){

    #Bind bin indicator column with instr
    mycols <- dplyr::select(bininds, x) %>%
      bind_cols(dplyr::select(mydf, myvar))

    #Rename bin indicator column so can refer to it directly
    names(mycols)[1] <- "mybin"

    #Create instrument for bin
    mycols <- mutate(mycols, inter = mybin*myvar) %>%
      #And only keep this column
      dplyr::select(inter)

    #Rename instrument column
    names(mycols) <- paste0(substr(x,4,nchar(x)),"_",var)

    mycols
  })

  outdf <- do.call("cbind",outdf)

  return(outdf)
}

rddf <- bind_cols(
  rddf,
  interVars("treatfrac")
)


#Realized juveniles caught
juvcatch <- felm(
  as.Formula(paste0(
  "asinhnummjuv", "~ ", 
  paste0(grep("_treatfrac",names(rddf),value=T),collapse="+"),
  " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
  paste0(grep("prop",names(rddf),value=T),collapse="+"),
  "| bin + twowk:cellid_2p + startdate",
  " | 0 | twoweek_cellid_2p")),
  data =rddf)

juvcatch$N #29664

jvtab <- summary(juvcatch)[["coefficients"]]

jvtab <- mutate(as.data.frame(jvtab), bin = rownames(jvtab))

jvtab <- jvtab[grep("treatfrac",jvtab$bin),]

jvtab$bin <- gsub("_treatfrac","",jvtab$bin)

#Calculate juv1 in potential closures, 
#then calculate juv0. Then scale up both by ratio of nummjuv in potential closures 
#to total numjuv.
toteffect_juv <- group_by(rddf, bin, tvar, bdist) %>% 
  summarise(juv1 = sum(nummjuv)) %>%
  left_join(jvtab, by = 'bin') %>% 
  mutate(juv0 = juv1/(exp(Estimate))) %>% 
  mutate(chmjuv = juv1 - juv0) %>%
  arrange(tvar, bdist) %>% ungroup()

#Percent effect, not accounting for reallocation
sum(toteffect_juv$chmjuv) / sum(toteffect_juv$juv0) 

#Scaled total effect see above comment. 
#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Total effect, not accounting for reallocation
#I guess potential closures nummjuv larger because of double-counting.
(sum(fullbe$numjuv,na.rm=T) / sum(rddf$numjuv, na.rm=T)) * sum(toteffect_juv$chmjuv)

toteffect_juv <- dplyr::select(toteffect_juv, -`t value`, -`Pr(>|t|)`)

toteffect_juv <- rename(toteffect_juv, juvcoef = Estimate, juvse = `Cluster s.e.`)

#Calculate delta method standard errors for change in millions of juveniles caught
library(msm)

toteffect_juv <- mutate(toteffect_juv, chmjuvse = as.numeric(NA), 
                        chmjuv_scaled = chmjuv * (sum(fullbe$numjuv,na.rm=T) / sum(rddf$numjuv, na.rm=T)),
                        chmjuvse_scaled = as.numeric(NA), 
                        juv0se = as.numeric(NA))

scaleconstant <- (sum(fullbe$numjuv,na.rm=T) / sum(rddf$numjuv, na.rm=T))

for(mybin in toteffect_juv$bin){
  
  #Filter toteffect_juv to mybin
  mydf <- filter(toteffect_juv, bin==mybin)
  
  #mypj <- mydf$meanpj/100 %>% as.character() %>% as.numeric()
  myjuv1 <- mydf$juv1
  myjuvcoef <- mydf$juvcoef
  myjuvvcov <- mydf$juvse^2

  #Delta se for juv0. Random variable being transformed is myjuvcoef
  juv0_delta <- deltamethod(~ (myjuv1 / (exp(x1))), myjuvcoef, myjuvvcov, ses=T)
  
  #Plug this value into toteffect_juv
  toteffect_juv$juv0se[toteffect_juv$bin==mybin] <- juv0_delta
  
  #Delta se for chmjuvse. Random variable being transformed is myjuvcoef
  chmjuv_delta <- deltamethod( ~ myjuv1 - (myjuv1/(exp(x1))),
                            myjuvcoef,
                            myjuvvcov,
                            ses=T
  )
  
  #Plug this value into toteffect_juv
  toteffect_juv$chmjuvse[toteffect_juv$bin==mybin] <- chmjuv_delta
  
  #Delta se for scaled chmjuvse. Random variable being transformed is myjuvcoef
  chmjuv_scaled_delta <- deltamethod( ~ (myjuv1 - (myjuv1/(exp(x1))))*scaleconstant,
                               myjuvcoef,
                               myjuvvcov,
                               ses=T
  )
  
  #Plug this value into toteffect_juv
  toteffect_juv$chmjuvse_scaled[toteffect_juv$bin==mybin] <- chmjuv_scaled_delta
  
}


##Calculate total effect, not accounting for reallocation of tons caught. (in millions)
changejuv <- (sum(fullbe$numjuv,na.rm=T) / sum(rddf$numjuv, na.rm=T)) * sum(toteffect_juv$chmjuv)

#Closures cannot increases tons caught because of TAC, so account for this reallocation
#by estimating how policy affects tons caught
tonscaught <- felm(
  as.Formula(paste0(
    "asinhtons", "~ ", 
    "treatfrac ",
    " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(rddf),value=T),collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data =rddf)

#Increase tons by exp(coef) - 1
summary(tonscaught)[["coefficients"]]["treatfrac",]

tonscoef <- summary(tonscaught)[["coefficients"]]["treatfrac","Estimate"]

#Don't reallocate tons from 2017 second season or 2019 second season because they were both 
#shut down well before TAC was reached. 

#Tons caught in state of world I observe in seasons where TAC hit
tons1 <- sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"])

#Change in tons because of policy
ctons <- tons1 - tons1/exp(tonscoef)

#Average pj outside of treatment window
avgpjoutside <- filter(fullbe, lead_0==0 & lead_10==0 & lead_20==0 & lead_30==0 & lead_40==0 & lead_50==0 & 
                         active_0==0 & active_10==0 & active_20==0 & active_30==0 & active_40==0 & active_50==0 & 
                         lag1_0==0 & lag1_10==0 & lag1_20==0 & lag1_30==0 & lag1_40==0 & lag1_50==0 & 
                         lag2_0==0 & lag2_10==0 & lag2_20==0 & lag2_30==0 & lag2_40==0 & lag2_50==0 & 
                         lag3_0==0 & lag3_10==0 & lag3_20==0 & lag3_30==0 & lag3_40==0 & lag3_50==0 & 
                         lag4_0==0 & lag4_10==0 & lag4_20==0 & lag4_30==0 & lag4_40==0 & lag4_50==0 & 
                         !is.na(numindivids) & !is.na(bepjhat) & Temporada!="2017-II" & Temporada!="2019-II") %>%
  #Weight by tons
  mutate(pjweighted = bepjhat*numindivids) %>%  
  summarise(perjuv = sum(pjweighted)/sum(numindivids)) %>% as.numeric() / 100 #0.09045436


#Avg weight of individual caught outside of treatment window
avgweightoutside <- filter(fullbe, lead_0==0 & lead_10==0 & lead_20==0 & lead_30==0 & lead_40==0 & lead_50==0 & 
                             active_0==0 & active_10==0 & active_20==0 & active_30==0 & active_40==0 & active_50==0 & 
                             lag1_0==0 & lag1_10==0 & lag1_20==0 & lag1_30==0 & lag1_40==0 & lag1_50==0 & 
                             lag2_0==0 & lag2_10==0 & lag2_20==0 & lag2_30==0 & lag2_40==0 & lag2_50==0 & 
                             lag3_0==0 & lag3_10==0 & lag3_20==0 & lag3_30==0 & lag3_40==0 & lag3_50==0 & 
                             lag4_0==0 & lag4_10==0 & lag4_20==0 & lag4_30==0 & lag4_40==0 & lag4_50==0 & 
                             !is.na(numindivids) & !is.na(avgweightg) & Temporada!="2017-II" & Temporada!="2019-II") %>%
  #Weight by tons
  mutate(weightweighted = avgweightg*numindivids) %>% 
  summarise(avgweightg = sum(weightweighted)/sum(numindivids)) %>% as.numeric()


#Decrease in individuals caught outside of treatment window in millions
#(converting tons to g cancels out conversion to millions)
chindividsoutside <- -ctons/avgweightoutside

chjuvsoutside <- chindividsoutside*avgpjoutside

#Now can calculate change in juvenile catch due to policy, accounting for reallocation
(chmjuvsstart <- changejuv + chjuvsoutside) #55969.1

#How many juveniles are caught during my sample period in total?
#F(1)*pj*individuals/VMS fishing obs
juv1 <- sum(fullbe$numjuv, na.rm=T) / 10^6 

#chmjuvsstart = juv(1) - juv(0)
juv0 <- juv1 - chmjuvsstart 

#Then increase in juvenile catch as a percentage is 
chmjuvsstart / juv0 # 0.6602906

#Calculate standard error on total change in juvenile catch and in total percentage change
mycoefs <- toteffect_juv$chmjuv_scaled
mybigvcov <- diag(toteffect_juv$chmjuvse_scaled^2)

#This includes reallocation, so this is what I want: 
(changebillionsse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                          x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                          x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                          x31 + x32 + x33 + x34 + x35 + x36) + chjuvsoutside) / 
                      1000, mycoefs, mybigvcov, ses=T))
#3.776392

#Now get SE on total percentage change
(totperse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                          x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                          x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                          x31 + x32 + x33 + x34 + x35 + x36) + chjuvsoutside) / 
                      juv0, mycoefs, mybigvcov, ses=T))
#0.04455166


finaldf <- toteffect_juv

finaldf$tvar[finaldf$tvar==-1] <- "lead"
finaldf$tvar[finaldf$tvar==0] <- "active"
finaldf$tvar[finaldf$tvar%in%c(1,2,3,4)] <- paste0("lag",finaldf$tvar[finaldf$tvar%in%c(1,2,3,4)])

#Plot treatment coefficients for consistent comparison with other appendix figures
finaldf <- rename(finaldf, juvcoef_5day = juvcoef, 
                  juvse_5day = juvse)

finaldf <- mutate(finaldf, 
                  juvcoef_5day_ub = juvcoef_5day + juvse_5day*qnorm(.975),
                  juvcoef_5day_lb = juvcoef_5day - juvse_5day*qnorm(.975))


#Function of one event day and dependent variable
singlePlot <- function(myvar, mytvar, ylab){
  
  #Title
  if(mytvar=="lead"){
    tit <- "After announcement, before closure"
  } else if(mytvar=="active"){
    tit <- "Closure period"
  } else if(mytvar=="lag1"){
    tit <- "1 day after"
  } else{
    tit <- paste0(gsub("lag","",mytvar), " days after")
  }
  
  #Rename desired variable
  usedf <- finaldf
  names(usedf)[names(usedf)==myvar] <- "plotvar"
  
  #Rename lb and ub so can refer to directly
  names(usedf)[names(usedf)==paste0(myvar,"_lb")] <- "lb"
  names(usedf)[names(usedf)==paste0(myvar,"_ub")] <- "ub"
  
  #Want consistent y range across given variable
  myymin <- min(usedf$lb)
  myymax <- max(usedf$ub)
  
  plot <- ggplot(data=filter(usedf, tvar==mytvar),aes(x=bdist)) + 
    geom_hline(aes(yintercept=0)) + 
    geom_point(aes(y=plotvar)) + 
    scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                       labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) +
    scale_y_continuous("",limits = c(myymin,myymax),
                       breaks = scales::pretty_breaks(n=10)) + 
    myThemeStuff + 
    ggtitle(tit) + 
    geom_errorbar(data=filter(usedf, tvar==mytvar), aes(x=bdist,ymin=lb,ymax=ub),width=0)
  
  #Add ylabel for leftmost plots
  if(mytvar=="lead"|mytvar=="lag2"){
    plot <- plot + scale_y_continuous(ylab,limits = c(myymin,myymax),
                       breaks = scales::pretty_breaks(n=10))
  }

  #There also used to be code here to add arrows for CI that got off, but again no longer necessary
  
  return(plot)
}




#Make a 2x3 plot for paper with labeled panels
paperFig <- function(myvar, ylab){
  leadplot <- singlePlot(myvar, "lead", ylab) + 
    labs(tag = "a") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  activeplot <- singlePlot(myvar, "active", ylab)+ 
    labs(tag = "b") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  lag1plot <- singlePlot(myvar, "lag1", ylab)+ 
    labs(tag = "c") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  lag2plot <- singlePlot(myvar, "lag2", ylab)+ 
    labs(tag = "d") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  lag3plot <- singlePlot(myvar, "lag3", ylab)+ 
    labs(tag = "e") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  lag4plot <- singlePlot(myvar, "lag4", ylab)+ 
    labs(tag = "f") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  
  tbt <- plot_grid(leadplot, activeplot, lag1plot, 
                   lag2plot, lag3plot, lag4plot, nrow=2, ncol=3, 
                   rel_widths = c(1.01,1,1))
  
  ggsave(tbt, file=paste0("Output/Figures/figureA8.png"),
         w=7,h=(7/1.69)*2, units = "in", dpi=1200)
}


paperFig("juvcoef_5day", TeX("$\\beta_{st}$ coefficient and 95% confidence interval (Equation 1)"))

sessionInfo()
