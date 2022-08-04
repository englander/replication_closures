rm(list=ls())

library(dplyr); library(readxl); library(ggplot2)
library(sf); library(gridExtra); library(lwgeom)
library(grid); library(geosphere); library(viridis)
library(sp); library(cowplot); library(rworldmap)
library(tidyr); library(rgdal)
library(purrr); library(lubridate); library(glmnet)
library(lfe); library(Formula); library(smoothr)
library(furrr); library(collapse)
library(scales)

sf::sf_use_s2(FALSE) 

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
                      plot.tag.position = c(.1, .98)
)

#Peru time
Sys.setenv(TZ='America/Lima')


#I have the length distribution for potential closures, but I need to recover 
#the percentage juvenile draws of the sets that generated them.
#This requires modifying 4. make_rddf.R, so that there is an additional column 
#with the percentage juvenile draws that I can unpack 
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
    
    #This is the new part from 4. make_rddf
    pjdrawsdf <- data.frame(newclust = unique(be$newclust), 
                          pjdraws = as.character(NA))
    
    for(i in unique(be$newclust)){
      pjdrawsdf$pjdraws[pjdrawsdf$newclust == i] <- 
        paste0(be$bepj[be$newclust == i], collapse = ";")
    }
    
    #Join onto out
    out <- left_join(out, pjdrawsdf, by = 'newclust')
    
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
plan(multisession, workers = 4)

rdflist <- future_map(datevec,
                     function(x){
                       
                       sumandQuants(x, 5)
                       
                     })

#Some days have no clusters so need to be dropped before row binding
dropel <- sapply(1:length(rdflist), function(x){
  ifelse("data.frame" %in% class(rdflist[[x]]),1,0)
})
dropel <- which(dropel==0)

rdf <- rdflist[-dropel]
rdf <- bind_rows(rdf)

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

#Set end for rectangles equal to 3 days since that is modal closure period length
rdf <- mutate(rdf, end = start + 3600*24*3-1)

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
#Function of start
treatVar <- function(mystart){
  
  #Rectangles
  myrect <- filter(rdf, start==mystart) %>% 
    mutate(treatfrac = 0) #Space-time fraction overlapping with actual closed area
  
  #When these rectangles end
  myend <- unique(myrect$end)
  
  #Filter to active closures
  #Only care about closures declared by PRODUCE
  actclosed <- filter(closed, days<=5 & 
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

rdf <- lapply(unique(rdf$start), function(x){
  treatVar(x)
})

rddf <- bind_rows(rdf)

#Temporarily save
#save(rddf, file = "G:/My Drive/Downloads/figure_lengthdist_example.Rdata")
load("G:/My Drive/Downloads/figure_lengthdist_example.Rdata")

#Don't need geometry anymore
rddf <- as.data.frame(rddf) %>% dplyr::select(-geometry) %>% 
  as_tibble()

#Potential closures with treatfrac equal 1
treatment <- filter(rddf, treatfrac > 0.9999) #R numerical sensitivity
  
#Potential closures with treatfrac equal 0
control <- filter(rddf, treatfrac == 0)

#Drop control potcl that are missing length distribution
control <- fsubset(control, !is.na(prop12hat))

#Create unique ids
treatment <- mutate(treatment, treatid = 1:nrow(treatment))

control <- mutate(control, controlid = 1:nrow(control))

#Given id of treated potcl and id of control potcl, calculate distance between length distribution, 
#mean reported pj, and number of sets
dissimFun <- function(mytreatid, mycontrolid){
  
  mytreat <- fsubset(treatment, treatid == mytreatid)
  mycontrol <- fsubset(control, controlid == mycontrolid)
  
  #Sum of absolute differences in length distribution for most common length intervals
  lengthdif <- sapply(names(mytreat)[grep("prop8hat", names(mytreat)):grep("prop15hat", names(mytreat))], 
                      function(x){
                        abs(
                          dplyr::select(mytreat, all_of(x)) %>% as.matrix() %>% as.numeric() - 
                            dplyr::select(mycontrol, all_of(x)) %>% as.matrix() %>% as.numeric()
                        )
                      }) %>% 
    sum()
  
  #Difference in (unweighted) mean reported pj
  pjdif <- abs(mytreat$meanpj - mycontrol$meanpj)
  
  #Difference in number of sets
  setsdif <- abs(mytreat$nobs - mycontrol$nobs)
  
  #Return df with id of treated potcl, control potcl, and three elements calculated above
  out <- data.frame(treatid = mytreatid, controlid = mycontrolid, lengthdif = lengthdif, 
                    pjdif = pjdif, setsdif = setsdif)
  
  return(out)
}

#Apply over control potcl
applyControl <- function(mytreatid){
  
  out <- map_df(unique(control$controlid), function(x){
    dissimFun(mytreatid, x)
  }) 
  
  return(out)
}

#Parallel apply over treated potcl
plan(multisession, workers = 6)

dsdf <- future_map(unique(treatment$treatid), function(x){
  try(applyControl(x))
})

#Temporarily save
#save(dsdf, file = 'G:/My Drive/Downloads/figure_lengthdist_dsdf.Rdata')
load('G:/My Drive/Downloads/figure_lengthdist_dsdf.Rdata')

dsdf <- bind_rows(dsdf)

#Given treatid and controlid, make a 2x2 plot of 
#length distribution and reported percentage juvenile
plotFun <- function(mytreatid, mycontrolid){
  
  #Select treated and control potential closure
  mytreat <- filter(treatment, treatid == mytreatid)
  
  mycontrol <- filter(control, controlid == mycontrolid)
  
  #Reshape length distribution columns
  treatlength <- mytreat[,grep("prop",names(mytreat))]
  names(treatlength) <- gsub("prop","",names(treatlength))
  names(treatlength) <- gsub("hat","",names(treatlength))
  names(treatlength) <- gsub("p","\\.",names(treatlength))
  treatlength <- pivot_longer(treatlength, everything()) %>% 
    rename(lengthbin = name, prop = value) %>% 
    mutate_all(as.numeric)
  
  controllength <- mycontrol[,grep("prop",names(mycontrol))]
  names(controllength) <- gsub("prop","",names(controllength))
  names(controllength) <- gsub("hat","",names(controllength))
  names(controllength) <- gsub("p","\\.",names(controllength))
  controllength <- pivot_longer(controllength, everything()) %>% 
    rename(lengthbin = name, prop = value) %>% 
    mutate_all(as.numeric)
  
  #12 cm bin really means individuals between 12 and 12.5 cm, 
  #And geom_col takes x coordinate as right interval
  #so shift bins so they reflect correct interval in plot
  treatlength <- mutate(treatlength, lengthbin = lengthbin + .5)
  controllength <- mutate(controllength, lengthbin = lengthbin + .5)
  
  #Label juveniles as less than 12 cm
  labtext <- data.frame(x=7,y=.1,text="Juveniles are < 12 cm")
  
  #Part a
  treatlengthdist <- ggplot() + 
    geom_col(data=treatlength, aes(x=lengthbin,y=prop)) +
    scale_x_continuous("Half-cm length interval (cm)",breaks=c(seq(from=3.25,to=17.25,by=2),18.75),
                       labels = c(seq(from=3,to=17,by=2),18.5)) + 
    scale_y_continuous("Mean proportion of individuals in length interval (corrected)",
                       limits = c(0, max(treatlength$prop, controllength$prop)),
                       breaks = breaks_extended(n = 5)) +
    myThemeStuff + 
    geom_vline(aes(xintercept=12.25),col='red') + 
    geom_text(data=labtext,aes(x=x,y=y,label=text),col='red',size=2.5) + 
    theme(axis.text = element_text(color = "black", size = 8, family="sans"),
          axis.title = element_text(color = "black", size = 8, family = "sans"),
          plot.margin = unit(c(.05, .1, .1, .01), 'in'), 
          plot.tag.position = c(.08, .985)) + 
    ggtitle("Length distribution of treated potential closure") + 
    labs(tag = "a")
  
  
  #Part b is length distribution of control
  controllengthdist <- ggplot() + 
    geom_col(data=controllength, aes(x=lengthbin,y=prop)) +
    scale_x_continuous("Half-cm length interval (cm)",breaks=c(seq(from=3.25,to=17.25,by=2),18.75),
                       labels = c(seq(from=3,to=17,by=2),18.5)) + 
    scale_y_continuous("Mean proportion of individuals in length interval (corrected)",
                       limits = c(0, max(treatlength$prop, controllength$prop)),
                       breaks = breaks_extended(n = 5)) +
    myThemeStuff + 
    geom_vline(aes(xintercept=12.25),col='red') + 
    theme(axis.text = element_text(color = "black", size = 8, family="sans"),
          axis.title = element_text(color = "black", size = 8, family = "sans"),
          plot.margin = unit(c(.05, .01, .1, .1), 'in'),
          plot.tag.position = c(.08, .985)) + 
    ggtitle("Length distribution of control potential closure") + 
    labs(tag = "b")
  
  #Now need to unpack reported percentage juvenile values
  treatpjdf <- data.frame(bepj = strsplit(mytreat$pjdraws, ";") %>% 
                            unlist()) %>% 
    mutate_all(as.character) %>% 
    mutate_all(as.numeric)
  
  controlpjdf <- data.frame(bepj = strsplit(mycontrol$pjdraws, ";") %>% 
                              unlist()) %>% 
    mutate_all(as.character) %>% 
    mutate_all(as.numeric)
  
  #Maximum count value across both dfs
  maxcountval <- max(
    count(treatpjdf, bepj) %>% 
      summarise(max(n)) %>% as.matrix() %>% as.numeric(),
    count(controlpjdf, bepj) %>% 
      summarise(max(n)) %>% as.matrix() %>% as.numeric()
  )
  
  #part c is reported percentage juvenile values for treated potential closure
  treatedpjdraws <- ggplot() + 
    geom_histogram(data = treatpjdf,
                   aes(x=bepj), bins = 400) + 
    myThemeStuff + 
    scale_x_continuous("Percentage juvenile (reported to regulator)",
                       breaks = seq(from = 0, to = 100, by = 20),
                       labels=paste0(seq(from=0,to=100,by=20),"%"),
                       limits = c(-0.25,100)) +
    scale_y_continuous("Count", limits = c(0, maxcountval),
                       breaks = breaks_extended(n = 5)) + 
    theme(axis.text = element_text(color = "black", size = 8, family="sans"),
          axis.title = element_text(color = "black", size = 8, family = "sans"),
          plot.margin = unit(c(.1, .1, .05, .01), 'in'),
          plot.tag.position = c(.06, .98)) + 
    ggtitle("% juvenile values for treated potential closure") + 
    labs(tag = "c")
  
  controlpjdraws <- ggplot() + 
    geom_histogram(data = controlpjdf,
                   aes(x=bepj), bins = 400, closed = 'left') + 
    myThemeStuff + 
    scale_x_continuous("Percentage juvenile (reported to regulator)",
                       breaks = seq(from = 0, to = 100, by = 20),
                       labels=paste0(seq(from=0,to=100,by=20),"%"),
                       limits = c(-0.25, 100)) + 
    scale_y_continuous("Count", limits = c(0, maxcountval),
                       breaks = breaks_extended(n = 5)) +
    theme(axis.text = element_text(color = "black", size = 8, family="sans"),
          axis.title = element_text(color = "black", size = 8, family = "sans"),
          plot.margin = unit(c(.1, .01, .05, .1), 'in'),
          plot.tag.position = c(.06, .98)) + 
    ggtitle("% juvenile values for control potential closure") + 
    labs(tag = "d")
  
  tbt <- plot_grid(treatlengthdist, controllengthdist, 
                   treatedpjdraws, controlpjdraws, 
                   nrow = 2, ncol=2, 
                   rel_heights = c(1.25, 1))
  
  ggsave(tbt, file=paste0("G:/My Drive/Peru/Output/PaperFigures/plot_treat",
                          mytreatid, "_control", mycontrolid, 
                          ".png"),
         w=7,h=(7/1.69)*(4/3), units = "in", dpi=900)
  
}

arrange(dsdf, lengthdif) %>% head()


plotFun(12, 164) #control potcl has many more zero draws. good in that same dist, within 5 days of each other, and similar area
plotFun(19, 222) #maybe perfect: visually identical length distributions, can see treated has slightly more high pj draws, similar number of sets. within two weeks of each other, treated is about 1/2 size
plotFun(22, 470) #not good. pj draws don't look different enough
plotFun(12, 150) #not good. difference in sets too great.
plotFun(12, 208) #not good. length dist look a little different
plotFun(15, 181) #ok. length dist look a little different. but can see higher pj draws for treated


filter(treatment, treatid == 19) %>% as.data.frame()
filter(control, controlid == 222) %>% as.data.frame()

#Summary statistics for 19 and 222
pj19 <- strsplit(treatment$pjdraws[19], ";") %>% 
  unlist() %>% as.numeric()

pj222 <- strsplit(control$pjdraws[222], ";") %>% 
  unlist() %>% as.numeric()

which(pj19 > 25) %>% length() #28

which(pj222 > 25) %>% length() #14


#For ref 3 comment 4, what % of potential closures have at least one set 
#that reports 10% or more juveniles?
g <- sapply(1:nrow(rddf), function(x){
  
  thedraws <- strsplit(rddf$pjdraws[x], ";") %>% 
    unlist() %>% as.numeric()
  
  ifelse(max(thedraws) >= 10, 1, 0)
})

which(g == 1) %>% length()
670/973
