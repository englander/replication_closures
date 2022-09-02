#In preparation for synthetic controls, calculate control variables 
#among sets that generate actual closures
rm(list=ls())

library(dplyr); library(readxl); library(ggplot2)
library(rworldmap); library(sf); library(lwgeom)
library(rgdal); library(geosphere); library(sp)
library(purrr); library(lubridate); library(glmnet)
library(lfe); library(Formula); library(smoothr)
library(furrr); library(tidyr)

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

options(scipen=999)

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

#Peru time
Sys.setenv(TZ='America/Lima')

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

fullbe <- rename(fullbe, calatime = FechaInicioCala)

#Make sf
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

fullbe <- st_sf(geometry = besf, fullbe)

rm(besf)

#Load closed areas created in make_closures_df.R
load("Output/Data/closed.Rdata")

#Only want to use closures declared by PRODUCE from BE data
closed <- filter(closed, days <= 5)

#North-Central closures during season
closed <- filter(closed, !is.na(season) & 
                   (season=="s1_2017" | season=="s2_2017" | 
                      season=="s1_2018" | season=="s2_2018" | season=="s1_2019" | 
                      season=="s2_2019")
)


bbox <- c(-84, -15.99999, -73, -3.5) 
names(bbox) <- names(st_bbox(closed))

closed <- st_crop(closed, bbox)

#Calculate area of each closure
clustarea_km2 <- st_area(closed)
clustarea_km2 <- as.numeric(clustarea_km2) / 10^6

closed <- mutate(closed, clustarea_km2 = clustarea_km2)

#Calculate centroid distance to coast
peru <- getMap(resolution = "high")
peru <- st_as_sf(peru) %>% filter(ADMIN.1 == "Peru")

centr <- st_centroid(closed)

kmtocoast <- st_distance(centr, peru)
kmtocoast <- as.numeric(kmtocoast) / 10^3

closed <- mutate(closed, kmtocoast = kmtocoast)

rm(clustarea_km2, centr, kmtocoast, bbox, peru)

#Indicator for whether closure begins at 6 am
closed$six <- 0
closed$six[grep("\\.25",closed$start_raw)] <- 1

#For instances when there are no sets with length distribution inside closures 9-24 hours before,
#make a df of all zeros for use in below function
lengthzero <- grep("hat",names(fullbe),value=T)
lengthzero <- grep("prop",lengthzero,value=T)
lzdf <- matrix(0,ncol=length(lengthzero),nrow=1) %>% as.data.frame()
names(lzdf) <- lengthzero
lengthzero <- lzdf; rm(lzdf)


#Calculate same control variables for actual closures as for potential closures
#using sets that occur inside actual closures 9-24 hours before they begin
closedControls <- function(rowind){
  
  myrow <- closed[rowind,]
  
  if(myrow$six==0){
  #Sets 9-24 hours before 
  mysets <- filter(fullbe, calatime >= (myrow$start-24*3600) & 
                     calatime <= (myrow$start - 9*3600))
  } else{
    #Then do 12-27 hours before
    mysets <- filter(fullbe, calatime >= (myrow$start-27*3600) & 
                       calatime <= (myrow$start - 12*3600))
  }
  #Sets inside
  inter <- st_intersects(myrow,mysets)
  
  mysets <- mysets[inter[[1]],]
  
  if(nrow(mysets)>0){
    
    #Length distribution
    if(filter(mysets, !is.na(prop12hat))%>%nrow()>0){
      
      out <- as.data.frame(mysets) %>% dplyr::select(-geometry) %>%
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

    } else{
      out <- lengthzero
    }
    
    #Number of sets, tons caught by them, tons per set, and tonsperarea
    out <- mutate(out, clustnobs = nrow(mysets), clusttons = sum(mysets$betons)) %>% 
      mutate(clusttonsperset = clusttons / clustnobs, 
             clusttonsperarea = clusttons / myrow$clustarea_km2)
    
  } else{
    out <- mutate(lengthzero, clustnobs = 0, clusttons = 0, 
                  clusttonsperset = 0, clusttonsperarea = 0)
  }
  
  #Start date of closure and closureid
  out <- mutate(out, startdate = date(myrow$start), 
                closeid = myrow$closeid)
  
  return(out)
}

#Parallel apply over closures
plan(multisession, workers = 14)

closedlist <- future_map(1:nrow(closed),
                     function(x){
                       
                       closedControls(x)
                       
                     })

closedlist <- bind_rows(closedlist)

#Join onto closed
closed <- left_join(closed, closedlist, by = 'closeid')

rm(closedlist, lengthzero)

#What is the percentage difference in tons per set 
#comparing above-median to below-median closures?
perdif <- as.data.frame(closed) %>% dplyr::select(-geometry) %>% 
  mutate(medsize = median(clustarea_km2)) %>% 
  mutate(abovemed = if_else(clustarea_km2 > medsize, 1, 0)) %>% 
  group_by(abovemed) %>% 
  summarise(clusttonsperset = mean(clusttonsperset)) %>% 
  t()
  
#16% higher
(perdif[,2] - perdif[,1]) / perdif[,1]

rm(perdif)

##Create closure treatment bins
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
  
  out <- bind_rows(out)
  
  return(out)
}

#Apply over all buffers
plan(multisession, workers = 14)

cbufs <- future_map(list(
                     c(0,10),c(10,20),c(20,30),c(30,40),c(40,50)                   
                   ), function(x){
                     applyBufFun_closed(bmin=x[[1]],bmax=x[[2]])
                   })

cbufs <- bind_rows(cbufs)

names(cbufs)[names(cbufs)!="geometry"] <- names(closed)[names(closed)!="geometry"]

closed <- bind_rows(closed, cbufs)

rm(cbufs)

#Duplicate closures for leads and lags
closed <- bind_rows(
  closed,
  #Preperiod is 9 hours before if begins at midnight; 12 hours before if begins at 6 am
  bind_rows(
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


##Calculate realized juvenile catch inside closure-treatment bins
outcomesFun <- function(rowind){
  
  row <- closed[rowind,]
  
  #Filter BE observations to within same time period
  mybesf <- filter(fullbe, row$start<=calatime & calatime<row$end)
  
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
    } else{
      tons <- 0; numindivids <- 0; numjuv <- 0; nobs <- 0; sdtons = as.numeric(NA)
      uniquevessels <- 0; tonsjuv <- 0; tonsadult <- 0
    }
  } else{
    tons <- 0; numindivids <- 0; numjuv <- 0; nobs <- 0; sdtons = as.numeric(NA)
    uniquevessels <- 0; tonsjuv <- 0; tonsadult <- 0
  }
  
  #Output df with bin and closeid so I can join onto closed
  out <- data.frame(bin=row$bin, closeid=row$closeid,tons=tons,
                    #Round number of individuals to integer
                    numindivids=round(numindivids), numjuv=round(numjuv),
                    tonsjuv = tonsjuv, tonsadult = tonsadult,
                    nobs = nobs, sdtons = sdtons, uniquevessels = uniquevessels)
  
  return(out)
}

#Apply over all bins
plan(multisession, workers = 14)

#Apply over rows of closed
myoutcomes <- future_map(1:nrow(closed),
                        function(x){
                          
                          outcomesFun(x)
                          
                        })

myoutcomes <- bind_rows(myoutcomes)

#Join on closeid and bin
closed <- left_join(closed, myoutcomes, by = c("bin","closeid"))

#Create 2 degree grid, created in 1. match_be_landings.R
load("Output/Data/grid2p.Rdata")

centr <- filter(closed, bin=="active_in") %>% 
  st_centroid()

#Which cell is each centroid inside of
insidecell <- st_intersects(centr, grid2p)

insidecell <- as.numeric(as.character(insidecell))

#Add cellid column onto centr
centr <- mutate(centr, cellid_2p = insidecell)

#Join cellid onto closed
closed <- left_join(closed, 
                  as.data.frame(centr) %>% dplyr::select(closeid, cellid_2p),
                  by='closeid')

#Treatfrac for all rows is 1
closed <- mutate(closed, treatfrac=1)


##Re-define two-week-of-sample for rddf and define it for closed so they are consistently defined
#Created in 4. make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

#Create two-week of sample variable. 
twowk <- bind_rows(
  as.data.frame(rddf) %>% dplyr::select(start) %>%
    mutate(week = week(start), year = year(start)),
  as.data.frame(closed) %>% dplyr::select(start) %>%
    mutate(week = week(start), year = year(start)))%>% 
  distinct(week,year)

twowk <- arrange(twowk, year, week)

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

#Join onto rddf and closed
rddf <- mutate(rddf, week = week(start), year = year(start)) %>% 
  dplyr::select(-twowk, -twoweek_cellid_2p) %>% 
  left_join(twowk, by = c("week","year")) %>% 
  dplyr::select(-week,-year)

closed <- mutate(closed, week = week(start), year = year(start)) %>% 
  left_join(twowk, by = c("week","year")) %>% 
  dplyr::select(-week,-year)

#For clustering
rddf <- mutate(rddf,twoweek_cellid_2p = paste0(twowk, "_", cellid_2p))
closed <- mutate(closed,twoweek_cellid_2p = paste0(twowk, "_", cellid_2p))

rm(twowk, counter, i, insidecell, centr, grid2p)


##Row bind closed and rddf
rddf <- dplyr::select(rddf, -newclust, -clustnumindivids, -clustnumjuv, -clustnumadults,
                      -meanpjhat, -meanpj, -meanpjhat_weighted, -meanpj_weighted, -avgweightg) %>% 
  mutate(startdate = date(start))

closed <- dplyr::select(closed, -`Comunicado Name`, -start_raw, -end_raw, 
                        -`top (S)`, -`bottom (S)`, -`right (W)`, -`left (W)`, 
                        -closepercom, -days, -six)

closed <- rename(closed, id = closeid)
rddf <- rename(rddf, id = rid)

rddf <- mutate(rddf, clusttonsperset = clusttons / clustnobs,
                 clusttonsperarea = clusttons / clustarea_km2)

closed <- mutate(closed, numadults = numindivids - numjuv)

#Row bind
regdf <- bind_rows(
  mutate(closed, closuretype = "actual"),
  mutate(rddf, closuretype='potential', 
         id = as.character(id))
)

save(regdf, file = "Output/TempData/actualclosure_regressioncontrol.Rdata")