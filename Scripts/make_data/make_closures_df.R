#Make an sf data frame with the id, start and end date and time, and geometry of temporarily closed area
rm(list=ls())

#Peru time
Sys.setenv(TZ='America/Lima')

library(dplyr); library(ggplot2); library(sf)
library(readr); library(rworldmap); library(rworldxtra)
library(lwgeom); library(lubridate)
library(parallel); library(readxl)

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
                      plot.margin = unit(c(0,0,0,0),"in"),
                      plot.tag = element_text(family = "sans", size = 9)
)

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

#Load closures data
closed <- read_xlsx("Data/closures_2014to2019.xlsx") %>% 
  arrange(Start) %>% 
  mutate_all(as.character)

#Want to create start and end dttm columns
#Current Start and End columns are dates and hand-coded from PRODUCE pdfs 
closed <- rename(closed, start_raw = Start, end_raw = End) %>%
  mutate(start_raw = as.character(start_raw), end_raw = as.character(end_raw)) %>%
  mutate(start = as.POSIXct(NA), end = as.POSIXct(NA))

#start_raw and end_raw that contain .25 are closures that begin on 6 am of start_raw date and end on 6 am of end_raw date
closed$start[grep("\\.25",closed$start_raw)] <- paste0(
  substr(closed$start_raw[grep("\\.25",closed$start_raw)],1,8),
  " 06:00:00"
) %>% strptime(format = "%Y%m%d %H:%M:%S", tz = 'America/Lima')

closed$end[grep("\\.25",closed$end_raw)] <- paste0(
  substr(closed$end_raw[grep("\\.25",closed$end_raw)],1,8),
  " 06:00:00"
) %>% strptime(format = "%Y%m%d %H:%M:%S", tz = 'America/Lima')

#Remaining start_raw and end_raw closed periods begin on start_raw date at midnight and end at 11:59 pm of end_raw date
closed$start[is.na(closed$start)] <- paste0(
  substr(closed$start_raw[is.na(closed$start)], 1, 8),
  " 00:00:00"
) %>% strptime(format = "%Y%m%d %H:%M:%S", tz = 'America/Lima')

closed$end[is.na(closed$end)] <- paste0(
  substr(closed$end_raw[is.na(closed$end)], 1, 8),
  " 23:59:59"
) %>% strptime(format = "%Y%m%d %H:%M:%S", tz = 'America/Lima')

#Small function to convert lat or lon from degree minutes to degree decimals
#Also special because lat and lon are always negative in this case
Dec <- function(num){
  
  num <- strsplit(num, " ") %>% unlist()
  
  if(length(num)==2){
    num <- -as.numeric(num[1]) - as.numeric(num[2])/60
  } else{
    num <- -as.numeric(num)
  }
  
  return(num)
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


#Shape of Peru land
peru <- getMap(resolution = "high")
peru <- st_as_sf(peru) %>%
  filter(NAME_SORT == "Peru")

#Load peru EEZ. Downloaded from https://www.marineregions.org/downloads.php on July 2, 2018.
#Version 2 - 2012 (4.96 MB) [Known issues] (created from EEZ version 7)
eez <- st_read("Data/Intersect_IHO_EEZ_v2_2012/eez.shp") %>%
  filter(EEZ == "Peruvian Exclusive Economic Zone (disputed - Peruvian point of view)")

#Project both.
peru <- st_transform(peru, st_crs("+proj=laea +lon_0=-76.5"))

eez <- st_transform(eez, st_crs("+proj=laea +lon_0=-76.5"))


#Another function to make a shape for closed area defined by lon and distance from shore
Buf <- function(rowind){
  
  row <- closed[rowind,]
  
  #Make an outer buffer from Peru land
  outbuf <- gsub(" nm","",row$`left (W)`) %>% as.numeric()
  
  #Convert nautical miles to meters
  outbuf <- st_buffer(peru, dist = outbuf*1852)
  
  #Keep intersection with eez
  outbuf <- st_intersection(outbuf, eez) %>%
    #Don't keep column variables
    dplyr::select()
  
  #Create inner buffer
  inbuf <- gsub(" nm","",row$`right (W)`) %>% as.numeric()
  
  #If closure extends from coast, just keep outbuf
  if(inbuf==0){
    buf <- outbuf %>%
      st_transform("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  } else{
    #Convert nautical miles to meters
    inbuf <- st_buffer(peru, dist = inbuf*1852)
    
    #Keep intersection with eez
    inbuf <- st_intersection(inbuf, eez) %>%
      #Don't keep column variables
      dplyr::select()
    
    #Subtract inbuf from outbuf to get restricted area in lon direction
    buf <- st_difference(outbuf, inbuf) %>%
      st_transform("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  }
  
  if(st_is_valid(buf)==FALSE){
    buf <- st_make_valid(buf)
  }
  
  #Now crop in lat direction. Create a rectangle to crop with.
  croprect <- Rect(top = Dec(row$`top (S)`), bottom = Dec(row$`bottom (S)`),
                             right = -70, left = -85)
  
  poly <- st_crop(
    buf,
    st_transform(croprect, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  )
  
  return(poly)
  
}

#Function that checks whether "nm" in right (W) variable. Execute first function if not, second function if yes
closedShape <- function(rowind){
  
  row <- closed[rowind,]
  
  if(grep(" nm",row$`right (W)`) %>% length==1){
    poly <- Buf(rowind)
  } else{
    poly <- Rect(top=Dec(row$`top (S)`),bottom=Dec(row$`bottom (S)`),
                 right=Dec(row$`right (W)`),left=Dec(row$`left (W)`)) %>%
      #Transform back to lat lon
      st_transform("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  }
  
  #Add data columns
  poly <- bind_cols(poly, row)
  
  return(poly)
}

#Apply over rows
closed <- lapply(1:nrow(closed), function(x){

  closedShape(x)

})

closed <- do.call("rbind", closed)

#Add a unique id for each closure event
closed <- mutate(closed, closeid = "NA")

#Record number of closures per Comunicado Name
closed <- left_join(closed,
                    group_by(closed, `Comunicado Name`) %>%
                      summarise(closepercom = n()) %>% ungroup() %>%
                      as.data.frame() %>% dplyr::select(-geometry))

closed <- lapply(unique(closed$`Comunicado Name`), function(x){

  df <- filter(closed, `Comunicado Name`==x)

  df$closeid <- paste0(df$`Comunicado Name`, "-", 1:unique(df$closepercom))

  df

})

closed <- do.call("rbind",closed)

#Calculate days closure takes
closed <- mutate(closed, days = difftime(end, start,units='days') %>% 
                   round() %>% as.numeric())

#Create season indicators corresponding to BE data
closed$season <- as.character(NA)

#2017 Season 1 begins April 26, 2017 and ends July 31, 2017
#2017 Season 2 begins November 27, 2017 and ends January 26, 2018 (last active day is January 25)
#2018 season 1 begins April 12, 2018 and ends August 9, 2018 (last active day is August 8)
#2018 season 2 begins November 15, 2018 and ends April 4, 2019 (last active day is April 3)
#2019 season 1 begins May 4, 2018 and ends July 31, 2019 (last active day is July 30, 2019)
#First BE observation for 2019 season 2 is Nov 6, 2019 and have coded closures up to Jan 9, 2020
#Start 24 hours after first BE observation in season (manually checking that this observation occurs at midnight)
closed$season[closed$start>=as.POSIXct("2017-04-23 00:00:00 -05") & 
                closed$start <= as.POSIXct("2017-07-31 23:59:59 -05")] <- "s1_2017"
closed$season[closed$start>=as.POSIXct("2017-11-24 00:00:00 -05") & 
                closed$start <= as.POSIXct("2018-01-25 23:59:59 -05")] <- "s2_2017"
closed$season[closed$start>=as.POSIXct("2018-04-08 00:00:00 -05") & 
                closed$start <= as.POSIXct("2018-08-08 23:59:59 -05")] <- "s1_2018"
closed$season[closed$start>=as.POSIXct("2018-11-15 00:00:00 -05") & 
                closed$start <= as.POSIXct("2019-04-03 23:59:59 -05")] <- "s2_2018"
closed$season[closed$start>=as.POSIXct("2019-04-29 00:00:00 -05") & 
                closed$start <= as.POSIXct("2019-07-30 23:59:59 -05")] <- "s1_2019"
closed$season[closed$start>=as.POSIXct("2019-11-07 00:00:00 -05") & 
                closed$start <= as.POSIXct("2020-01-14 23:59:59 -05")] <- "s2_2019"

save(closed, file = "Output/Data/closed.Rdata")

#What percent of closures start at 6 AM?
closed$six <- 0
closed$six[grep("\\.25",closed$start_raw)] <- 1
mean(closed$six) # 0.06867846

sessionInfo()