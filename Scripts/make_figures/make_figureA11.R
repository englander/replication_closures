rm(list=ls())

#Peru time
Sys.setenv(TZ='America/Lima')

library(dplyr); library(ggplot2); library(sf)
library(lubridate); library(Formula)
library(purrr); library(parallel); library(lfe)
library(cowplot); library(latex2exp)

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

options(lfe.threads=12)

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

#Load actual closures created in make_closures_df.R and create treatment bins
{
load("Output/Data/closed.Rdata")

#Only want to use closures declared by PRODUCE from BE data
closed <- filter(closed, days <= 5)

#Filter to closures during season
closed <- filter(closed, !is.na(season) & 
                   (season=="s1_2017" | season=="s2_2017" | 
                      season=="s1_2018" | season=="s2_2018" | season=="s1_2019" | 
                      season=="s2_2019")
)

#Indicator for whether closure begins at 6 am
closed$six <- 0
closed$six[grep("\\.25",closed$start_raw)] <- 1

closed <- dplyr::select(closed, start, end, closeid, season, six)

#Create treatment bins
#Have to make corresponding closed buffers and time lags
closed <- mutate(closed, bin = "active_0",bdist=0,tvar="active")

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

cl <- makeCluster(12)

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


#Duplicate rectangles for leads and lags
closed <- rbind(
  closed,
  #Preperiod is 9 hours before if begins at midnight; 12 hours before if begins at 6 am
  rbind(
    filter(closed, six==0) %>% mutate(end=start-1) %>% 
      mutate(start=start-9*3600,tvar="lead"),
    filter(closed, six==1) %>% mutate(end=start-1) %>% 
      mutate(start=start-12*3600,tvar="lead")
  ),
  closed %>% mutate(start = end+1,end=end+24*3600*1,tvar="lag1"), 
  closed %>% mutate(start = end+1+24*3600*1,end=end+24*3600*2,tvar="lag2"),
  closed %>% mutate(start = end+1+24*3600*2,end=end+24*3600*3,tvar="lag3"),
  closed %>% mutate(start = end+1+24*3600*3,end=end+24*3600*4,tvar="lag4"))

#Remake bin variable
closed$bin[closed$tvar=="lead"] <- gsub("active","lead",closed$bin[closed$tvar=="lead"])
closed$bin[closed$tvar=="lag1"] <- gsub("active","lag1",closed$bin[closed$tvar=="lag1"])
closed$bin[closed$tvar=="lag2"] <- gsub("active","lag2",closed$bin[closed$tvar=="lag2"])
closed$bin[closed$tvar=="lag3"] <- gsub("active","lag3",closed$bin[closed$tvar=="lag3"])
closed$bin[closed$tvar=="lag4"] <- gsub("active","lag4",closed$bin[closed$tvar=="lag4"])
}


#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

fullbe <- rename(fullbe, calatime = FechaInicioCala)

#Make sf; only need a few variables
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))
besf <- st_sfc(besf) %>% st_cast('POINT')
st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
besf <- st_sf(geometry = besf, dplyr::select(fullbe, calatime, numjuv, betons, Temporada))

besf <- mutate(besf, season = paste0(ifelse(nchar(Temporada)==6,"s1_","s2_"),substr(Temporada,1,4))) %>% 
  dplyr::select(-Temporada)

rm(fullbe)

#Create .05 degree grid
#Load peru EEZ. Downloaded from https://www.marineregions.org/downloads.php on July 2, 2018.
#Version 2 - 2012 (4.96 MB) [Known issues] (created from EEZ version 7)
eez <- st_read("C:/Users/englander/Box Sync/VMS/Data/Intersect_IHO_EEZ_v2_2012/eez.shp") %>%
  filter(EEZ == "Peruvian Exclusive Economic Zone (disputed - Peruvian point of view)")

grid <- st_make_grid(st_bbox(eez), cellsize = .05, what = 'polygons') %>%
  st_sf()

#Keep only grid cells that intersect eez
inter <- st_intersects(eez, grid) %>% unlist()

grid <- grid[inter,]

#Also crop to those in North-Central zone
grid <- st_crop(grid, c(ymin = -16, ymax = -3, xmin=-85, xmax=-70))

rm(inter)

#Add centroid of grid as variable
centr <- st_centroid(grid)

centrcoords <- st_coordinates(centr)

#Add coordinates of centroid as columns for joining between two objects
grid <- mutate(grid, lon_centr = centrcoords[,1], 
               lat_centr = centrcoords[,2])

centr <- mutate(centr, lon_centr = centrcoords[,1], 
                lat_centr = centrcoords[,2])

rm(centrcoords)

#Create two-degree grid cell variable
#Create 2 degree grid, created in 1. match_be_landings.R
load("Output/Data/grid2p.Rdata")

twog <- st_intersects(centr, grid2p)
twog <- as.numeric(as.character(twog))

grid <- mutate(grid, cellid_2p = twog)
centr <- mutate(centr, cellid_2p = twog)

grid <- mutate(grid, gridcellrow = 1:nrow(grid))

rm(twog, grid2p, eez)

#Create single-row containing bin indicator variables
fillinds <- matrix(0,ncol=36,nrow=1) %>% as.data.frame()
names(fillinds) <- unique(closed$bin)

#Small function. Given row that I know is inside treatment window, calculate 
#all treatment bins the row is inside
smallTreat <- function(centrrow, threeclosed_nosf, inter, centr_nosf, mythree){
  
  giveninter <- inter[[centrrow]]
  
  #Rows of threeclosed row is inside
  insideclosed <- threeclosed_nosf[giveninter,] %>% 
    #Can do this because spatial rings are rings, not buffers
    distinct(tvar, bdist)
  
  #Mash tvar_bdist together 
  charvec <- mutate(insideclosed, yeszone = paste0(tvar,"_",bdist)) %>% 
    dplyr::select(yeszone) %>% as.matrix() %>% as.character()
  
  #Zones that should be 1 because observation is in these zones distances
  yeszones <- which(names(fillinds) %in% charvec)
  
  out <- fillinds
  
  out[,yeszones] <- 1
  
  out <- bind_cols(centr_nosf[centrrow,] %>% mutate(time=mythree), out)
  
  return(out)
}

#Calculate treatment for each cell-time and each of the 36 treatment bins
#Function of three-hour unit 
treatThree <- function(mythree){
  
  #Filter myclosed to this time period
  threeclosed <- filter(myclosed, start <= (mythree + 3600*3) & 
                          mythree <= end)
  
  #Cells in treatment window
  inter <- st_intersects(centr, threeclosed)
  
  #Non-missing
  treatwindow <- sapply(1:nrow(centr), function(x){
    ifelse(length(inter[[x]])>0,1,0)
  })
  
  #Centroids object without geometry column
  centr_nosf <- as.data.frame(centr) %>% dplyr::select(-geometry)
  
  #Three closed without geometry column
  threeclosed_nosf <- as.data.frame(threeclosed) %>% dplyr::select(-geometry)
  
  #Apply smallTreat over rows in treatwindow
  inwindow <- map_df(which(treatwindow==1), function(x){
    smallTreat(x, threeclosed_nosf, inter, centr_nosf, mythree)
  })
  
  #Other rows have 0 for all treatment bins
  notwindow <- map_df(which(treatwindow==0), function(x){
    bind_cols(centr_nosf[x,] %>% mutate(time = mythree), fillinds)
  })
  
  out <- bind_rows(inwindow, notwindow)
  
  return(out)
  
}


#Parallel apply one season at a time

useseason <- "s1_2017"
{
#Filter to season
mybe <- filter(besf, season==useseason)
myclosed <- filter(closed, season==useseason)

#First and last day in season (either first fishing or first closure)
firstday <- as.Date(min(mybe$calatime, myclosed$start))
lastday <- as.Date(max(mybe$calatime, myclosed$end)) 

#Three hour time blocks
timevec <- seq(from=as.POSIXct(paste0(firstday, "00:00:00 -05")), 
               to = as.POSIXct(paste0(lastday + 1, "00:00:00 -05")),
               by = 3600*3)

#Calculate outcome variable in each cell-time by cutting by timevec to get time, 
#then intersecting with grid to get cell
mybe$threehour <- cut(mybe$calatime, timevec, right = FALSE)

inter <- st_intersects(mybe, grid)

inter <- as.numeric(as.character(inter))

mybe <- mutate(mybe, gridcellrow = inter)

#Get coordinates of centroid for .05 degree grid cell
mybe <- left_join(mybe, 
                  as.data.frame(grid) %>% dplyr::select(-geometry),
                  by = 'gridcellrow')

#Outcome variable is summed juvenile catch, tons, and number of sets
myoutcome <- as.data.frame(mybe) %>% dplyr::select(-geometry) %>% 
  group_by(threehour, lon_centr, lat_centr) %>% 
  summarise(nummjuv = sum(numjuv)/10^6, tons = sum(betons), nobs = n()) %>% 
  group_by()

myoutcome$threehour <- as.character(myoutcome$threehour) %>% as.POSIXct()

rm(inter)

#Calculate treatment fraction for each grid-time
#Apply over all buffers
(myCores <- detectCores())

cl <- makeCluster(12)

clusterExport(cl, "timevec")
clusterExport(cl, "myclosed")
clusterExport(cl, "centr")
clusterExport(cl, "smallTreat")
clusterExport(cl, "treatThree")
clusterExport(cl, "fillinds")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(purrr))
clusterEvalQ(cl, library(lubridate))


treatlist <- parLapply(cl=cl,
                      timevec, function(x){
                     treatThree(x)
                   })

stopCluster(cl)
rm(cl, myCores)

paneldf <- bind_rows(treatlist)

#Join outcomes onto paneldf
paneldf <- left_join(
  rename(paneldf, threehour = time), myoutcome, 
  by = c("lon_centr","lat_centr","threehour"))

#If missing outcome, it is 0
paneldf$nummjuv[is.na(paneldf$nummjuv)] <- 0
paneldf$tons[is.na(paneldf$tons)] <- 0
paneldf$nobs[is.na(paneldf$nobs)] <- 0


#Calculate two-week-of-sample and join onto df
twoweek <- data.frame(threehour = timevec)

twoweek <- mutate(twoweek, week = week(threehour), year = year(threehour))

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

#Join onto paneldf
paneldf <- left_join(paneldf, dplyr::select(twoweek, threehour, twowk), by = c("threehour"))

#Add season
paneldf <- mutate(paneldf, season = useseason)

#Make sure two week of sample values are different for each season
paneldf <- mutate(paneldf, twowk = paste0(twowk,"_",season))

#Create two-week-of-sample by two-degree grid cell variable
paneldf <- mutate(paneldf, twoweek_cellid_2p = paste0(twowk, "_", cellid_2p))

save(paneldf, file = paste0("Output/TempData/rasterjuv_",useseason,".Rdata"))

#Clean up
rm(useseason, paneldf, twoweek, firstday, lastday, 
   timevec, myclosed, mybe, myoutcome, treatlist, twowk, i, counter)
}

useseason <- "s2_2017"
{
  #Filter to season
  mybe <- filter(besf, season==useseason)
  myclosed <- filter(closed, season==useseason)
  
  #First and last day in season (either first fishing or first closure)
  firstday <- as.Date(min(mybe$calatime, myclosed$start))
  lastday <- as.Date(max(mybe$calatime, myclosed$end)) 
  
  #Three hour time blocks
  timevec <- seq(from=as.POSIXct(paste0(firstday, "00:00:00 -05")), 
                 to = as.POSIXct(paste0(lastday + 1, "00:00:00 -05")),
                 by = 3600*3)
  
  #Calculate outcome variable in each cell-time by cutting by timevec to get time, 
  #then intersecting with grid to get cell
  mybe$threehour <- cut(mybe$calatime, timevec, right = FALSE)
  
  inter <- st_intersects(mybe, grid)
  
  inter <- as.numeric(as.character(inter))
  
  mybe <- mutate(mybe, gridcellrow = inter)
  
  #Get coordinates of centroid for .05 degree grid cell
  mybe <- left_join(mybe, 
                    as.data.frame(grid) %>% dplyr::select(-geometry),
                    by = 'gridcellrow')
  
  #Outcome variable is summed juvenile catch, tons, and number of sets
  myoutcome <- as.data.frame(mybe) %>% dplyr::select(-geometry) %>% 
    group_by(threehour, lon_centr, lat_centr) %>% 
    summarise(nummjuv = sum(numjuv)/10^6, tons = sum(betons), nobs = n()) %>% 
    group_by()
  
  myoutcome$threehour <- as.character(myoutcome$threehour) %>% as.POSIXct()
  
  rm(inter)
  
  #Calculate treatment fraction for each grid-time
  #Apply over all buffers
  (myCores <- detectCores())
  
  cl <- makeCluster(12)
  
  clusterExport(cl, "timevec")
  clusterExport(cl, "myclosed")
  clusterExport(cl, "centr")
  clusterExport(cl, "smallTreat")
  clusterExport(cl, "treatThree")
  clusterExport(cl, "fillinds")
  clusterEvalQ(cl, library(dplyr))
  clusterEvalQ(cl, library(sf))
  clusterEvalQ(cl, library(purrr))
  clusterEvalQ(cl, library(lubridate))
  
  
  treatlist <- parLapply(cl=cl,
                         timevec, function(x){
                           treatThree(x)
                         })
  
  stopCluster(cl)
  rm(cl, myCores)
  
  paneldf <- bind_rows(treatlist)
  
  #Join outcomes onto paneldf
  paneldf <- left_join(
    rename(paneldf, threehour = time), myoutcome, 
    by = c("lon_centr","lat_centr","threehour"))
  
  #If missing outcome, it is 0
  paneldf$nummjuv[is.na(paneldf$nummjuv)] <- 0
  paneldf$tons[is.na(paneldf$tons)] <- 0
  paneldf$nobs[is.na(paneldf$nobs)] <- 0
  
  
  #Calculate two-week-of-sample and join onto df
  twoweek <- data.frame(threehour = timevec)
  
  twoweek <- mutate(twoweek, week = week(threehour), year = year(threehour))
  
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
  
  #Join onto paneldf
  paneldf <- left_join(paneldf, dplyr::select(twoweek, threehour, twowk), by = c("threehour"))
  
  #Add season
  paneldf <- mutate(paneldf, season = useseason)
  
  #Make sure two week of sample values are different for each season
  paneldf <- mutate(paneldf, twowk = paste0(twowk,"_",season))
  
  #Create two-week-of-sample by two-degree grid cell variable
  paneldf <- mutate(paneldf, twoweek_cellid_2p = paste0(twowk, "_", cellid_2p))
  
  save(paneldf, file = paste0("Output/TempData/rasterjuv_",useseason,".Rdata"))
  
  #Clean up
  rm(useseason, paneldf, twoweek, firstday, lastday, 
     timevec, myclosed, mybe, myoutcome, treatlist, twowk, i, counter)
}

useseason <- "s1_2018"
{
  #Filter to season
  mybe <- filter(besf, season==useseason)
  myclosed <- filter(closed, season==useseason)
  
  #First and last day in season (either first fishing or first closure)
  firstday <- as.Date(min(mybe$calatime, myclosed$start))
  lastday <- as.Date(max(mybe$calatime, myclosed$end)) 
  
  #Three hour time blocks
  timevec <- seq(from=as.POSIXct(paste0(firstday, "00:00:00 -05")), 
                 to = as.POSIXct(paste0(lastday + 1, "00:00:00 -05")),
                 by = 3600*3)
  
  #Calculate outcome variable in each cell-time by cutting by timevec to get time, 
  #then intersecting with grid to get cell
  mybe$threehour <- cut(mybe$calatime, timevec, right = FALSE)
  
  inter <- st_intersects(mybe, grid)
  
  inter <- as.numeric(as.character(inter))
  
  mybe <- mutate(mybe, gridcellrow = inter)
  
  #Get coordinates of centroid for .05 degree grid cell
  mybe <- left_join(mybe, 
                    as.data.frame(grid) %>% dplyr::select(-geometry),
                    by = 'gridcellrow')
  
  #Outcome variable is summed juvenile catch, tons, and number of sets
  myoutcome <- as.data.frame(mybe) %>% dplyr::select(-geometry) %>% 
    group_by(threehour, lon_centr, lat_centr) %>% 
    summarise(nummjuv = sum(numjuv)/10^6, tons = sum(betons), nobs = n()) %>% 
    group_by()
  
  myoutcome$threehour <- as.character(myoutcome$threehour) %>% as.POSIXct()
  
  rm(inter)
  
  #Calculate treatment fraction for each grid-time
  #Apply over all buffers
  (myCores <- detectCores())
  
  cl <- makeCluster(12)
  
  clusterExport(cl, "timevec")
  clusterExport(cl, "myclosed")
  clusterExport(cl, "centr")
  clusterExport(cl, "smallTreat")
  clusterExport(cl, "treatThree")
  clusterExport(cl, "fillinds")
  clusterEvalQ(cl, library(dplyr))
  clusterEvalQ(cl, library(sf))
  clusterEvalQ(cl, library(purrr))
  clusterEvalQ(cl, library(lubridate))
  
  
  treatlist <- parLapply(cl=cl,
                         timevec, function(x){
                           treatThree(x)
                         })
  
  stopCluster(cl)
  rm(cl, myCores)
  
  paneldf <- bind_rows(treatlist)
  
  #Join outcomes onto paneldf
  paneldf <- left_join(
    rename(paneldf, threehour = time), myoutcome, 
    by = c("lon_centr","lat_centr","threehour"))
  
  #If missing outcome, it is 0
  paneldf$nummjuv[is.na(paneldf$nummjuv)] <- 0
  paneldf$tons[is.na(paneldf$tons)] <- 0
  paneldf$nobs[is.na(paneldf$nobs)] <- 0
  
  
  #Calculate two-week-of-sample and join onto df
  twoweek <- data.frame(threehour = timevec)
  
  twoweek <- mutate(twoweek, week = week(threehour), year = year(threehour))
  
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
  
  #Join onto paneldf
  paneldf <- left_join(paneldf, dplyr::select(twoweek, threehour, twowk), by = c("threehour"))
  
  #Add season
  paneldf <- mutate(paneldf, season = useseason)
  
  #Make sure two week of sample values are different for each season
  paneldf <- mutate(paneldf, twowk = paste0(twowk,"_",season))
  
  #Create two-week-of-sample by two-degree grid cell variable
  paneldf <- mutate(paneldf, twoweek_cellid_2p = paste0(twowk, "_", cellid_2p))
  
  save(paneldf, file = paste0("Output/TempData/rasterjuv_",useseason,".Rdata"))
  
  #Clean up
  rm(useseason, paneldf, twoweek, firstday, lastday, 
     timevec, myclosed, mybe, myoutcome, treatlist, twowk, i, counter)
}

useseason <- "s2_2018"
{
  #Filter to season
  mybe <- filter(besf, season==useseason)
  myclosed <- filter(closed, season==useseason)
  
  #First and last day in season (either first fishing or first closure)
  firstday <- as.Date(min(mybe$calatime, myclosed$start))
  lastday <- as.Date(max(mybe$calatime, myclosed$end)) 
  
  #Three hour time blocks
  timevec <- seq(from=as.POSIXct(paste0(firstday, "00:00:00 -05")), 
                 to = as.POSIXct(paste0(lastday + 1, "00:00:00 -05")),
                 by = 3600*3)
  
  #Calculate outcome variable in each cell-time by cutting by timevec to get time, 
  #then intersecting with grid to get cell
  mybe$threehour <- cut(mybe$calatime, timevec, right = FALSE)
  
  inter <- st_intersects(mybe, grid)
  
  inter <- as.numeric(as.character(inter))
  
  mybe <- mutate(mybe, gridcellrow = inter)
  
  #Get coordinates of centroid for .05 degree grid cell
  mybe <- left_join(mybe, 
                    as.data.frame(grid) %>% dplyr::select(-geometry),
                    by = 'gridcellrow')
  
  #Outcome variable is summed juvenile catch, tons, and number of sets
  myoutcome <- as.data.frame(mybe) %>% dplyr::select(-geometry) %>% 
    group_by(threehour, lon_centr, lat_centr) %>% 
    summarise(nummjuv = sum(numjuv)/10^6, tons = sum(betons), nobs = n()) %>% 
    group_by()
  
  myoutcome$threehour <- as.character(myoutcome$threehour) %>% as.POSIXct()
  
  rm(inter)
  
  #Calculate treatment fraction for each grid-time
  #Apply over all buffers
  (myCores <- detectCores())
  
  cl <- makeCluster(12)
  
  clusterExport(cl, "timevec")
  clusterExport(cl, "myclosed")
  clusterExport(cl, "centr")
  clusterExport(cl, "smallTreat")
  clusterExport(cl, "treatThree")
  clusterExport(cl, "fillinds")
  clusterEvalQ(cl, library(dplyr))
  clusterEvalQ(cl, library(sf))
  clusterEvalQ(cl, library(purrr))
  clusterEvalQ(cl, library(lubridate))
  
  
  treatlist <- parLapply(cl=cl,
                         timevec, function(x){
                           treatThree(x)
                         })
  
  stopCluster(cl)
  rm(cl, myCores)
  
  paneldf <- bind_rows(treatlist)
  
  #Join outcomes onto paneldf
  paneldf <- left_join(
    rename(paneldf, threehour = time), myoutcome, 
    by = c("lon_centr","lat_centr","threehour"))
  
  #If missing outcome, it is 0
  paneldf$nummjuv[is.na(paneldf$nummjuv)] <- 0
  paneldf$tons[is.na(paneldf$tons)] <- 0
  paneldf$nobs[is.na(paneldf$nobs)] <- 0
  
  
  #Calculate two-week-of-sample and join onto df
  twoweek <- data.frame(threehour = timevec)
  
  twoweek <- mutate(twoweek, week = week(threehour), year = year(threehour))
  
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
  
  #Join onto paneldf
  paneldf <- left_join(paneldf, dplyr::select(twoweek, threehour, twowk), by = c("threehour"))
  
  #Add season
  paneldf <- mutate(paneldf, season = useseason)
  
  #Make sure two week of sample values are different for each season
  paneldf <- mutate(paneldf, twowk = paste0(twowk,"_",season))
  
  #Create two-week-of-sample by two-degree grid cell variable
  paneldf <- mutate(paneldf, twoweek_cellid_2p = paste0(twowk, "_", cellid_2p))
  
  save(paneldf, file = paste0("Output/TempData/rasterjuv_",useseason,".Rdata"))
  
  #Clean up
  rm(useseason, paneldf, twoweek, firstday, lastday, 
     timevec, myclosed, mybe, myoutcome, treatlist, twowk, i, counter)
}

useseason <- "s1_2019"
{
  #Filter to season
  mybe <- filter(besf, season==useseason)
  myclosed <- filter(closed, season==useseason)
  
  #First and last day in season (either first fishing or first closure)
  firstday <- as.Date(min(mybe$calatime, myclosed$start))
  lastday <- as.Date(max(mybe$calatime, myclosed$end)) 
  
  #Three hour time blocks
  timevec <- seq(from=as.POSIXct(paste0(firstday, "00:00:00 -05")), 
                 to = as.POSIXct(paste0(lastday + 1, "00:00:00 -05")),
                 by = 3600*3)
  
  #Calculate outcome variable in each cell-time by cutting by timevec to get time, 
  #then intersecting with grid to get cell
  mybe$threehour <- cut(mybe$calatime, timevec, right = FALSE)
  
  inter <- st_intersects(mybe, grid)
  
  inter <- as.numeric(as.character(inter))
  
  mybe <- mutate(mybe, gridcellrow = inter)
  
  #Get coordinates of centroid for .05 degree grid cell
  mybe <- left_join(mybe, 
                    as.data.frame(grid) %>% dplyr::select(-geometry),
                    by = 'gridcellrow')
  
  #Outcome variable is summed juvenile catch, tons, and number of sets
  myoutcome <- as.data.frame(mybe) %>% dplyr::select(-geometry) %>% 
    group_by(threehour, lon_centr, lat_centr) %>% 
    summarise(nummjuv = sum(numjuv)/10^6, tons = sum(betons), nobs = n()) %>% 
    group_by()
  
  myoutcome$threehour <- as.character(myoutcome$threehour) %>% as.POSIXct()
  
  rm(inter)
  
  #Calculate treatment fraction for each grid-time
  #Apply over all buffers
  (myCores <- detectCores())
  
  cl <- makeCluster(12)
  
  clusterExport(cl, "timevec")
  clusterExport(cl, "myclosed")
  clusterExport(cl, "centr")
  clusterExport(cl, "smallTreat")
  clusterExport(cl, "treatThree")
  clusterExport(cl, "fillinds")
  clusterEvalQ(cl, library(dplyr))
  clusterEvalQ(cl, library(sf))
  clusterEvalQ(cl, library(purrr))
  clusterEvalQ(cl, library(lubridate))
  
  
  treatlist <- parLapply(cl=cl,
                         timevec, function(x){
                           treatThree(x)
                         })
  
  stopCluster(cl)
  rm(cl, myCores)
  
  paneldf <- bind_rows(treatlist)
  
  #Join outcomes onto paneldf
  paneldf <- left_join(
    rename(paneldf, threehour = time), myoutcome, 
    by = c("lon_centr","lat_centr","threehour"))
  
  #If missing outcome, it is 0
  paneldf$nummjuv[is.na(paneldf$nummjuv)] <- 0
  paneldf$tons[is.na(paneldf$tons)] <- 0
  paneldf$nobs[is.na(paneldf$nobs)] <- 0
  
  
  #Calculate two-week-of-sample and join onto df
  twoweek <- data.frame(threehour = timevec)
  
  twoweek <- mutate(twoweek, week = week(threehour), year = year(threehour))
  
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
  
  #Join onto paneldf
  paneldf <- left_join(paneldf, dplyr::select(twoweek, threehour, twowk), by = c("threehour"))
  
  #Add season
  paneldf <- mutate(paneldf, season = useseason)
  
  #Make sure two week of sample values are different for each season
  paneldf <- mutate(paneldf, twowk = paste0(twowk,"_",season))
  
  #Create two-week-of-sample by two-degree grid cell variable
  paneldf <- mutate(paneldf, twoweek_cellid_2p = paste0(twowk, "_", cellid_2p))
  
  save(paneldf, file = paste0("Output/TempData/rasterjuv_",useseason,".Rdata"))
  
  #Clean up
  rm(useseason, paneldf, twoweek, firstday, lastday, 
     timevec, myclosed, mybe, myoutcome, treatlist, twowk, i, counter)
}

useseason <- "s2_2019"
{
  #Filter to season
  mybe <- filter(besf, season==useseason)
  myclosed <- filter(closed, season==useseason)
  
  #First and last day in season (either first fishing or first closure)
  firstday <- as.Date(min(mybe$calatime, myclosed$start))
  lastday <- as.Date(max(mybe$calatime, myclosed$end)) 
  
  #Three hour time blocks
  timevec <- seq(from=as.POSIXct(paste0(firstday, "00:00:00 -05")), 
                 to = as.POSIXct(paste0(lastday + 1, "00:00:00 -05")),
                 by = 3600*3)
  
  #Calculate outcome variable in each cell-time by cutting by timevec to get time, 
  #then intersecting with grid to get cell
  mybe$threehour <- cut(mybe$calatime, timevec, right = FALSE)
  
  inter <- st_intersects(mybe, grid)
  
  inter <- as.numeric(as.character(inter))
  
  mybe <- mutate(mybe, gridcellrow = inter)
  
  #Get coordinates of centroid for .05 degree grid cell
  mybe <- left_join(mybe, 
                    as.data.frame(grid) %>% dplyr::select(-geometry),
                    by = 'gridcellrow')
  
  #Outcome variable is summed juvenile catch, tons, and number of sets
  myoutcome <- as.data.frame(mybe) %>% dplyr::select(-geometry) %>% 
    group_by(threehour, lon_centr, lat_centr) %>% 
    summarise(nummjuv = sum(numjuv)/10^6, tons = sum(betons), nobs = n()) %>% 
    group_by()
  
  myoutcome$threehour <- as.character(myoutcome$threehour) %>% as.POSIXct()
  
  rm(inter)
  
  #Calculate treatment fraction for each grid-time
  #Apply over all buffers
  (myCores <- detectCores())
  
  cl <- makeCluster(12)
  
  clusterExport(cl, "timevec")
  clusterExport(cl, "myclosed")
  clusterExport(cl, "centr")
  clusterExport(cl, "smallTreat")
  clusterExport(cl, "treatThree")
  clusterExport(cl, "fillinds")
  clusterEvalQ(cl, library(dplyr))
  clusterEvalQ(cl, library(sf))
  clusterEvalQ(cl, library(purrr))
  clusterEvalQ(cl, library(lubridate))
  
  
  treatlist <- parLapply(cl=cl,
                         timevec, function(x){
                           treatThree(x)
                         })
  
  stopCluster(cl)
  rm(cl, myCores)
  
  paneldf <- bind_rows(treatlist)
  
  #Join outcomes onto paneldf
  paneldf <- left_join(
    rename(paneldf, threehour = time), myoutcome, 
    by = c("lon_centr","lat_centr","threehour"))
  
  #If missing outcome, it is 0
  paneldf$nummjuv[is.na(paneldf$nummjuv)] <- 0
  paneldf$tons[is.na(paneldf$tons)] <- 0
  paneldf$nobs[is.na(paneldf$nobs)] <- 0
  
  
  #Calculate two-week-of-sample and join onto df
  twoweek <- data.frame(threehour = timevec)
  
  twoweek <- mutate(twoweek, week = week(threehour), year = year(threehour))
  
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
  
  #Join onto paneldf
  paneldf <- left_join(paneldf, dplyr::select(twoweek, threehour, twowk), by = c("threehour"))
  
  #Add season
  paneldf <- mutate(paneldf, season = useseason)
  
  #Make sure two week of sample values are different for each season
  paneldf <- mutate(paneldf, twowk = paste0(twowk,"_",season))
  
  #Create two-week-of-sample by two-degree grid cell variable
  paneldf <- mutate(paneldf, twoweek_cellid_2p = paste0(twowk, "_", cellid_2p))
  
  save(paneldf, file = paste0("Output/TempData/rasterjuv_",useseason,".Rdata"))
  
  #Clean up
  rm(useseason, paneldf, twoweek, firstday, lastday, 
     timevec, myclosed, mybe, myoutcome, treatlist, twowk, i, counter)
}

#Load and bind data from each season
load("Output/TempData/rasterjuv_s1_2017.Rdata")
regdf <- paneldf

loadvec <- c("s2_2017","s1_2018","s2_2018","s1_2019","s2_2019")

for(i in loadvec){

  load(paste0("Output/TempData/rasterjuv_",i,".Rdata"))

  regdf <- bind_rows(regdf, paneldf)

}

rm(i, paneldf, loadvec)

regdf$lon_centr <- as.factor(regdf$lon_centr)
regdf$lat_centr <- as.factor(regdf$lat_centr)
regdf$threehour <- as.factor(regdf$threehour)
regdf$twoweek_cellid_2p <- as.factor(regdf$twoweek_cellid_2p)

regdf <- mutate(regdf, asinhnummjuv = asinh(nummjuv))

formlfe <- as.Formula(paste0(
  "asinhnummjuv", "~ ", 
  paste0(grep("lead_|active_|lag",names(regdf),value=T),collapse="+"),
  " | lon_centr:lat_centr + threehour + twoweek_cellid_2p | 0 | twoweek_cellid_2p"))

juvcatch <- felm(formlfe, data=regdf)

#Plot treatment coefficients and standard errors
jvtab <- summary(juvcatch)[["coefficients"]]

jvtab <- mutate(as.data.frame(jvtab), bin = rownames(jvtab))

jvtab <- dplyr::select(jvtab, -`t value`, -`Pr(>|t|)`)

jvtab <- rename(jvtab, rasterjuv = Estimate, rasterjuvse = `Cluster s.e.`)

#Confidence intervals
jvtab <- mutate(jvtab, rasterjuv_ub = rasterjuv + rasterjuvse*qnorm(.975), 
       rasterjuv_lb = rasterjuv - rasterjuvse*qnorm(.975))

#Separate tvar and bdist variables
jvtab$bdist <- gsub(".*_","",jvtab$bin)
jvtab$bdist <- as.numeric(jvtab$bdist)
jvtab$tvar <- gsub("_.*","",jvtab$bin)

#Function of one time period and dependent variable
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
  usedf <- jvtab
  names(usedf)[names(usedf)==myvar] <- "plotvar"
  
  #Rename lb and ub so can refer to directly
  names(usedf)[names(usedf)==paste0(myvar,"_lb")] <- "lb"
  names(usedf)[names(usedf)==paste0(myvar,"_ub")] <- "ub"
  
  #Want consistent y range across given variable
  myymin <- min(usedf$lb)
  myymax <- max(usedf$ub)
  
  
  if(mytvar=="lead"|mytvar=="lag2"){
    
    plot <- ggplot(data=filter(usedf, tvar==mytvar),aes(x=bdist)) + 
      geom_hline(aes(yintercept=0)) + 
      geom_point(aes(y=plotvar)) + 
      scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                         labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) +
      scale_y_continuous(ylab,limits = c(myymin,myymax),
                         breaks = scales::pretty_breaks(n=10)) + 
      myThemeStuff + 
      ggtitle(tit) + 
      geom_errorbar(data=filter(usedf, tvar==mytvar), aes(x=bdist,ymin=lb,ymax=ub),width=0)

  } else{
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
  }

  
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
  
  ggsave(tbt, file=paste0("Output/Figures/figureA11.png"),
         w=7,h=(7/1.69)*2, units = "in", dpi=1200)
}

paperFig("rasterjuv",TeX("$\\beta_{st}$ coefficient and 95% confidence interval"))


sessionInfo()