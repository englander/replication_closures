#This is the third script cleaning the BE data. It follows 2. impute_size_be.R
#In this script, I will calculate average trip-level pj in the BE data, 
#calculate bepjhat--what pj should be to make individuals-weighted pj match landing level pj--
#shift length distribution so that it implies bepjhat and calculate updated numindivids, numjuv etc. 
#I will then repeat procedure for the sets not matched to a landing event, 
#using matched sets in same twoweek-cell group

rm(list=ls())

library(dplyr); library(readxl); library(ggplot2)
library(sf); library(purrr); library(lubridate)
library(parallel); library(tidyr)

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

options(scipen=999)

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

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

#Peru time
Sys.setenv(TZ='America/Lima')

#Created in impute_size_be*.R
load("Output/Data/pbe_imp_uncorrected.Rdata")

#Calculate average be at landid level, weighted by number of individuals, and join onto fullbe
fullbe <- left_join(fullbe, 
                       mutate(fullbe, weightpj = bepj*numindivids) %>% 
                         filter(!is.na(landid) & !is.na(landpj)) %>%
                         group_by(landid) %>% 
                         summarise(avgbepj = sum(weightpj)/sum(numindivids), sumnumindivids = sum(numindivids), maxpj = max(bepj)) %>%
                         ungroup() %>% 
                         #If numindivids = 0 and all bepj = 0, avgbepj should equal 0 too
                         mutate(avgbepj = if_else(sumnumindivids==0 & maxpj==0, 0, avgbepj)) %>% 
                         dplyr::select(landid, avgbepj),
                       by = 'landid')

#Inflate or deflate reported BE % juv by ratio of land pj to average pj 
fullbe <- mutate(fullbe, bepjhat = if_else(landpj==0 & avgbepj==0,0,bepj*(landpj/avgbepj)))

#If avgbepj==0 & landpj > 0, add (landpj - avgbepj) (i.e. set it equal to landpj)
fullbe$bepjhat[fullbe$avgbepj==0 & fullbe$landpj>0 & !is.na(fullbe$landpj) & !is.na(fullbe$avgbepj)] <- 
  fullbe$landpj[fullbe$avgbepj==0 & fullbe$landpj>0 & !is.na(fullbe$landpj) & !is.na(fullbe$avgbepj)]

#If landpj==100, bepjhat should equal 100 as well (as long as tons are not 0)
fullbe$bepjhat[fullbe$landpj==100 & !is.na(fullbe$landpj) & 
                    fullbe$numindivids > 0] <- 100

#For landid trips with some sets = 0 and others > 0, inflation above gives correct avg pj, but can't give 
#some sets pj > 100. In those cases, recalculate bepjhat to  so that  
#landpj equals weighted average bepjhat 
greater100 <- filter(fullbe, bepjhat > 100.000000000001) %>% distinct(landid) #For some reason a value that seems to equal 100 is giving TRUE for being > 100

for(i in greater100$landid){
  
  whichrows <- which(fullbe$landid==i)
  
  maxval <- max(fullbe$bepjhat[whichrows])
  
  #Reset bepjhat for these rows
  fullbe$bepjhat[whichrows] <- fullbe$bepj[whichrows]
  
  #constant to add to bepj to construct new bepjhat
  constant <- unique(fullbe$landpj[whichrows]) - 
    #Weighted average bepjhat
    as.numeric(fullbe$bepjhat[whichrows]%*%fullbe$numindivids[whichrows]) / 
    sum(fullbe$numindivids[whichrows])
  
  #Add this constant to bepjhat
  fullbe$bepjhat[whichrows] <- fullbe$bepjhat[whichrows] + constant
  
  #Check whether can stop
  maxval <- max(fullbe$bepjhat[whichrows])
  
  while(maxval > 100.000000000001){
    
    #Average of rows with bepjhat >= 100
    avgexc <- mean(fullbe$bepjhat[fullbe$landid==i & fullbe$bepjhat >= 100 & 
                                       !is.na(fullbe$bepjhat)])
    
    #Difference with landpj, divided by number of rows
    landpjdiff <- (avgexc - unique(fullbe$landpj[whichrows])) / length(whichrows)
    
    #Add to rows with bepjhat < 100
    fullbe$bepjhat[fullbe$landid==i & fullbe$bepjhat < 100 & 
                        !is.na(fullbe$bepjhat)]  <- 
      fullbe$bepjhat[fullbe$landid==i & fullbe$bepjhat < 100 & 
                          !is.na(fullbe$bepjhat)] + 
      landpjdiff
    
    #Subtract from rows with bepjhat >= 100
    fullbe$bepjhat[fullbe$landid==i & fullbe$bepjhat >= 100 & 
                        !is.na(fullbe$bepjhat)]  <- 
      fullbe$bepjhat[fullbe$landid==i & fullbe$bepjhat >= 100 & 
                          !is.na(fullbe$bepjhat)] - 
      landpjdiff
    
    #Make weighted average bepjhat equal landpj again
    correctfac <- unique(fullbe$landpj[whichrows]) - 
      (sum(fullbe$bepjhat[whichrows]%*%fullbe$numindivids[whichrows]) / 
         sum(fullbe$numindivids[whichrows]))
    
    fullbe$bepjhat[whichrows] <- fullbe$bepjhat[whichrows] + correctfac
    
    #Check whether can exit while loop
    maxval <- max(fullbe$bepjhat[whichrows]) 
  }
}

rm(i, maxval, constant, whichrows, correctfac, greater100, landpjdiff, avgexc)


#Some observations could not be matched to landing or the % juvenile in landings was missing
#For these, want to update bepj using non-missing observations in same twoweek-cell group
#Calculate average landpj and bepj at two-week cell group level and join onto fullbe
fullbe <- left_join(fullbe, 
                       filter(fullbe, !is.na(landid) & !is.na(landpj)) %>% 
                         mutate(lweightpj = landpj*numindivids, bweightpj = avgbepj*numindivids) %>%
                         group_by(twoweek_cellid_2p_landidlevel) %>%
                         summarise(grouplandpj = sum(lweightpj)/sum(numindivids),
                                   groupbepj = sum(bweightpj)/sum(numindivids)) %>% ungroup(),
                       by='twoweek_cellid_2p_landidlevel')

#If missing bepjhat, use bepj * ratio of grouplandpj to groupbepj
fullbe <- mutate(fullbe, bepjhat = if_else(is.na(bepjhat) & !is.na(groupbepj) & !is.na(grouplandpj) & 
                                                   groupbepj!=0 & grouplandpj!=0,
                                                 bepj*(grouplandpj/groupbepj), bepjhat))

#If groupbepj==0 & grouplandpj > 0, add (grouplandpj - groupbepj) (i.e. set it equal to grouplandpj)
fullbe$bepjhat[is.na(fullbe$bepjhat) & fullbe$groupbepj==0 & fullbe$grouplandpj>0 & 
                    !is.na(fullbe$grouplandpj) & !is.na(fullbe$groupbepj)] <- 
  fullbe$grouplandpj[is.na(fullbe$bepjhat) & fullbe$groupbepj==0 & fullbe$grouplandpj>0 & 
                          !is.na(fullbe$grouplandpj) & !is.na(fullbe$groupbepj)]

#If groupbepj>=0 & grouplandpj == 0 and bepjhat missing, set it equal to 0
fullbe$bepjhat[is.na(fullbe$bepjhat) & fullbe$groupbepj>=0 & fullbe$grouplandpj==0 & 
                    !is.na(fullbe$grouplandpj) & !is.na(fullbe$groupbepj)] <- 0

#For remaining 148 observations with missing bepjhat, just use bepj
fullbe$bepjhat[is.na(fullbe$bepjhat)] <- fullbe$bepj[is.na(fullbe$bepjhat)]

#Several hundred observations > 100 (because of multiplying by ratio above). 
#Not worth it to correct with while script like I did above; Just use bepj for these observations
fullbe$bepjhat[fullbe$bepjhat > 100] <- fullbe$bepj[fullbe$bepjhat > 100]


#Use values from IMARPE (2019)
lengthweight <- function(length){
  #length in cm; weight in g
  weight <- .0036*length^3.238
  return(weight)
}


#For calculating weight of individuals in subsequent function
sizedf <- data.frame(lengthname = c("3hat","3p5hat","4hat","4p5hat","5hat","5p5hat","6hat","6p5hat",
                                    "7hat","7p5hat","8hat","8p5hat","9hat","9p5hat","10hat","10p5hat",
                                    "11hat","11p5hat","12hat","12p5hat","13hat","13p5hat","14hat","14p5hat",
                                    "15hat","15p5hat","16hat","16p5hat","17hat","17p5hat","18hat","18p5hat"))
sizedf$lengthnum <- gsub("hat","",sizedf$lengthname)
sizedf$lengthnum <- gsub("p",".",sizedf$lengthnum)
sizedf$lengthnum <- as.numeric(sizedf$lengthnum)

sizedf <- mutate(sizedf, weight = lengthweight(lengthnum))

sizedf <- arrange(sizedf, lengthnum)

sizedf <- mutate(sizedf, charlength = gsub("hat","",lengthname))

#For grabbing columns from fullbe in below function
sizedf <- mutate(sizedf, proplengthname = paste0("prop",lengthname))

#Names of variables that are less than 12 cm
juvvars <- sizedf$proplengthname[sizedf$lengthnum < 12] %>% as.character()

#Names of variables that are 12 cm or above
adultvars <- sizedf$proplengthname[sizedf$lengthnum >= 12] %>% as.character()

##Update length distribution for all rows to match bepjhat
#And re-calculate numjuv etc.
shiftDist <- function(rowind){
  
  myrow <- fullbe[rowind,]
  
  if(!is.na(myrow$prop12hat)){
    
    #New columns to hold shifted size distribution values
    holddf <- dplyr::select(myrow,juvvars,adultvars) %>% 
      gather("length","prop") %>%
      mutate(lengthnum = gsub("prop","",length)) %>% 
      mutate(lengthnum = gsub("hat","",lengthnum)) %>%
      mutate(lengthnum = gsub("p",".",lengthnum) %>% as.numeric())
    
    #If prop sums to less than 1, inflate all positive rows so that prop adds up to 1
    if(sum(holddf$prop) < 1){
      holddf <- mutate(holddf, prop = if_else(prop > 0, prop*(1/sum(holddf$prop)), prop))
    }
    
    holddf <- mutate(holddf, propshift = prop) 
    
    #Difference between bepjhat and bepj
    pjdif <- myrow$bepjhat - myrow$bepj
    
    #Is there positive weight in the 3 cm column? 
    #If so, that's fine, but if not, I will tell while loop to stop when distribution is shifted 
    #so far that there is positive weight in the 3 cm column. This avoids the three instances
    #when bepjhat is so high and distribution so diffuse that it is impossible to shift 
    #distribution enough to match bepjhat, because positive values eventually get shifted out of the data faster
    #than they are shifted in
    pos3 <- ifelse(myrow$prop3hat>0,1,0)
    
    if(pjdif > 0.00001){
      uppj <- myrow$bepj
      
      if(pos3==1){
        #In case of numeric sensitivity (one row has uppj = 100 but returns true for uppj < myrow$bepjhat (=100))
        while(uppj < (myrow$bepjhat-.0000000001)){
          
          #Shift length distribution down one half cm interval and calculate new pj
          holddf <- mutate(holddf, propshift = lead(propshift, 1))
          
          #New pj 
          uppj <- sum(holddf$propshift[holddf$lengthnum < 12], na.rm=T)*100
          
        }
      } else{
        while(uppj < (myrow$bepjhat-.0000000001) & holddf$propshift[holddf$lengthnum==3] == 0){
          
          #Shift length distribution down one half cm interval and calculate new pj
          holddf <- mutate(holddf, propshift = lead(propshift, 1))
          
          #New pj 
          uppj <- sum(holddf$propshift[holddf$lengthnum < 12], na.rm=T)*100
          
        }
      }
      
      #Shift length distribution up one half cm interval (up = smaller)
      #These two shifted distributions bound the "true" pj (bepjhat)
      #Number of shifts 
      numshifts <- is.na(holddf$propshift) %>% sum()
      
      #So shift one less than this 
      holddf <- mutate(holddf, propshift2 = lead(prop, numshifts - 1))
      
      #Use propshift that is closest to bepjhat
      if( abs(sum(holddf$propshift[holddf$lengthnum < 12], na.rm=T) - myrow$bepjhat/100) <= 
          abs(myrow$bepjhat/100 - sum(holddf$propshift2[holddf$lengthnum < 12], na.rm=T))
      ){
        holddf <- rename(holddf, prophat = propshift)
      } else{
        holddf <- rename(holddf, prophat = propshift2)
      }
      
      holddf <- dplyr::select(holddf, length, prophat) %>% 
        mutate(prophat = if_else(is.na(prophat),0,prophat)) #%>% 
      #spread(length, prophat)
      
    }else if(pjdif < -0.00001){
      uppj <- myrow$bepj
      
      while(uppj > myrow$bepjhat){
        
        #Shift length distribution down one half cm interval and calculate new pj
        holddf <- mutate(holddf, propshift = lag(propshift, 1))
        
        #New pj 
        uppj <- sum(holddf$propshift[holddf$lengthnum < 12], na.rm=T)*100
        
      }
      
      #Shift length distribution down one half cm interval (down = bigger)
      #These two shifted distributions bound the "true" pj (bepjhat)
      #Number of shifts 
      numshifts <- is.na(holddf$propshift) %>% sum()
      #So shift one less than this 
      holddf <- mutate(holddf, propshift2 = lag(prop, numshifts - 1))
      
      #Use propshift that is closest to bepjhat
      if( abs(sum(holddf$propshift[holddf$lengthnum < 12], na.rm=T) - myrow$bepjhat/100) <= 
          abs(myrow$bepjhat/100 - sum(holddf$propshift2[holddf$lengthnum < 12], na.rm=T))
      ){
        holddf <- rename(holddf, prophat = propshift)
      } else{
        holddf <- rename(holddf, prophat = propshift2)
      }
      
      holddf <- dplyr::select(holddf, length, prophat) %>% 
        mutate(prophat = if_else(is.na(prophat),0,prophat)) #%>% 
      #spread(length, prophat)
      
    } else{
      holddf <- dplyr::select(holddf, length, prop) %>% 
        rename(prophat = prop) #%>% 
      #spread(length, prophat)
    }
    
    #If length distribution sums to less than 1 (because a value was shifted out of distribution),
    #inflate all positive rows so that prop adds up to 1
    if(sum(holddf$prophat) < 1){
      holddf <- mutate(holddf, prophat = if_else(prophat > 0, prophat*(1/sum(holddf$prophat)), prophat))
    }
    
    #Replace length distribution values
    spreadholddf <- spread(holddf, length, prophat)

    myrow <- myrow[,-grep("prop",names(myrow))] %>% 
      bind_cols(spreadholddf)
    
    #Replace avgweightg, numindivids, numjuv, and numadults for this row
    myrow$avgweightg <- sapply(holddf$length, function(x){
      myrow[,x]*sizedf$weight[sizedf$proplengthname==x]
    }) %>% sum()
    
    myrow <- myrow %>%
      mutate(numindivids = (betons*10^6)/avgweightg) %>%
      #Number of juveniles
      mutate(numjuv = numindivids*bepjhat/100) %>%
      #Number of adults
      mutate(numadults = numindivids-numjuv) 
    
    ##Also calculate juvbiomass and adultbiomass
    #Number of individuals in each length interval
    holddf <- dplyr::select(myrow, juvvars, adultvars) %>% 
      gather(length, prophat) %>% 
      mutate(numindivids_l = prophat*myrow$numindivids) %>% 
      mutate(lengthnum = gsub("hat","",length)) %>% 
      mutate(lengthnum = gsub("prop","",lengthnum)) %>% 
      mutate(lengthnum = gsub("p","\\.",lengthnum)) %>% 
      mutate(lengthnum = as.numeric(lengthnum)) %>% 
      #Weight of individuals in length interval in tons
      mutate(weightindivids_l = (numindivids_l*lengthweight(lengthnum))/10^6)
    
    myrow$tonsjuv <- sum(holddf$weightindivids_l[holddf$lengthnum < 12])
    
    myrow$tonsadult <- sum(holddf$weightindivids_l[holddf$lengthnum >= 12])
    
  } else{
    myrow <- mutate(myrow, tonsjuv = as.numeric(NA), tonsadult = as.numeric(NA))
  }

  return(myrow)
}

#Apply over rows
#Parallel apply over datevec
(myCores <- detectCores())

cl <- makeCluster(22)

clusterExport(cl, "fullbe")
clusterExport(cl, "sizedf")
clusterExport(cl, "adultvars")
clusterExport(cl, "juvvars")
clusterExport(cl, "lengthweight")
clusterExport(cl, "shiftDist")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(purrr))
clusterEvalQ(cl, library(tidyr))


mylist <- parLapply(cl = cl,
                    1:nrow(fullbe),
                    function(x){
                      
                      try(shiftDist(x))
                      
                    })

fullbe <- bind_rows(mylist)

stopCluster(cl)
rm(cl, myCores)

#Make sure all Southern zone observations are dropped 
fullbe <- filter(fullbe, lat >  -16) #(this only drops an additional 6 observations)

#Also calculate which treatment bins of actual closure each set is in.
#Load closed areas and create spatial-temporal spillover zones
{
  #created in make_closures_df.R
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
  
  #Only need start and end date of closure and season and closeid
  closed <- dplyr::select(closed, start, end, season, closeid, six)
  
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
  
  #Function to make buffer for each closed area
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
  closed$bin[closed$tvar==-1] <- gsub("active","lead",closed$bin[closed$tvar==-1])
  closed$bin[closed$tvar==1] <- gsub("active","lag1",closed$bin[closed$tvar==1])
  closed$bin[closed$tvar==2] <- gsub("active","lag2",closed$bin[closed$tvar==2])
  closed$bin[closed$tvar==3] <- gsub("active","lag3",closed$bin[closed$tvar==3])
  closed$bin[closed$tvar==4] <- gsub("active","lag4",closed$bin[closed$tvar==4])
  
  #Make tvar match formatting in bin variable
  closed$tvar <- as.character(closed$tvar)
  closed$tvar[closed$tvar=="-1"] <- "lead"
  closed$tvar[closed$tvar=="0"] <- "active"
  closed$tvar[closed$tvar!="lead" & closed$tvar!="active"] <- paste0("lag",closed$tvar[closed$tvar!="lead" & closed$tvar!="active"])

}

##Determine which treatment zones each set is in
#Create single-row containing zone indicator variables
fillinds <- matrix(0,ncol=36,nrow=1) %>% as.data.frame()
names(fillinds) <- unique(closed$bin)

#Except rename _in to _0 for use in inBuf function. Can change back to _in afterwards
names(fillinds) <- gsub("_in","_0",names(fillinds))

#Re-make fullbe as sf
mp <- st_multipoint(cbind(fullbe$lon, fullbe$lat))
mp <- st_sfc(mp) %>% st_cast("POINT")
st_crs(mp) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
fullbe <- st_sf(geometry = mp, fullbe)

#Given row and temporal lead or lag, see which buffer distances it is in
inBuf <- function(myrowind){
  
  out <- fillinds
  
  row <- fullbe[myrowind,] 
  
  #Filter to active closure zones
  actclosed <- filter(closed, start<=row$FechaInicioCala & row$FechaInicioCala<=end)
  
  inter <- st_intersects(row, actclosed)
  
  #Filter to closure zones that row is spatially inside
  insideclosed <- actclosed[unlist(inter),] %>%
    as.data.frame() %>% dplyr::select(-geometry) %>%
    group_by(closeid, tvar) %>%    #there will only be one value of tvar for a given closeid
    #Take minimum buffer distance for each closure, to get the smallest buffer obs is inside
    summarise(bdist = min(bdist)) %>% ungroup() %>%
    distinct(tvar, bdist)
  
  if(nrow(insideclosed)>0){
    #Mash tvar_bdist together 
    charvec <- mutate(insideclosed, yeszone = paste0(tvar,"_",bdist)) %>% 
      dplyr::select(yeszone) %>% as.matrix() %>% as.character()
    
    #Zones that should be 1 because observation is in these zones distances
    yeszones <- which(names(fillinds) %in% charvec)
    
    out[,yeszones] <- 1
  }
  
  
  #Bind with row
  out <- bind_cols(as.data.frame(row) %>% dplyr::select(-geometry), out)
  
  return(out)
}

rm(mp)

(myCores <- detectCores())

cl <- makeCluster(24)

clusterExport(cl, "fullbe")
clusterExport(cl, "closed")
clusterExport(cl, "inBuf")
clusterExport(cl, "fillinds")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lubridate))


belist <- parLapply(cl = cl,
                    1:nrow(fullbe),
                    function(x){
                      
                      try(inBuf(x))
                      
                    })

fullbe <- bind_rows(belist)

stopCluster(cl)
rm(cl, myCores)

save(fullbe, file = "Output/Data/pbe_imp.Rdata")

#What % of tons landed do SNP vessels represent?
filter(fullbe, !is.na(stons)) %>% 
  summarise(sum(betons)) %>% as.matrix() %>% as.numeric() / sum(fullbe$betons) #0.5645073

#Calculate weighted average percentage juvenile in both datasets
filter(fullbe, !is.na(bepj) & !is.na(numindivids)) %>% 
  mutate(weightpj = bepj*numindivids) %>% 
  summarise(sum(weightpj) / sum(numindivids)) #11.0004

filter(fullbe, !is.na(landpj) & !is.na(numindivids)) %>% 
  mutate(weightpj = landpj*numindivids) %>% 
  summarise(sum(weightpj) / sum(numindivids))  #18.34531

#Percentage difference 
(11.0004 - 18.34531) / 18.34531

#Calculate tons reported by sets matched to landing events
#and compare the two 
filter(fullbe, !is.na(betons) & !is.na(landtons)) %>% 
  summarise(sum(betons)) / filter(fullbe, !is.na(betons) & !is.na(landtons)) %>% 
     distinct(landtons, landid) %>% 
     summarise(sum(landtons))
#1.08391


##Median length anchoveta caught is between 13 and 13.5 cm
lengthvars <- grep("hat",names(fullbe),value=T)
lengthvars <- grep("prop",lengthvars,value=T)

medlength <- dplyr::select(fullbe, numindivids, lengthvars) %>% 
  filter(!is.na(numindivids)) %>% 
  gather(length,prop,-numindivids) %>% 
  mutate(individs_l = numindivids*prop) %>% 
  group_by(length) %>% 
  summarise(numindivids = sum(individs_l)) 

medlength <- mutate(medlength, prop = numindivids/sum(medlength$numindivids)) 

#Get numeric length
medlength <- mutate(medlength, lengthnum = gsub("prop","",length)) %>% 
  mutate(lengthnum = gsub("hat","",lengthnum)) %>% 
  mutate(lengthnum = gsub("p",".",lengthnum)) %>% 
  mutate(lengthnum = as.numeric(lengthnum)) %>% 
  arrange(lengthnum)

filter(medlength, lengthnum < 13) %>% 
  summarise(sum(prop))

filter(medlength, lengthnum > 13) %>% 
  summarise(sum(prop))

filter(medlength, lengthnum == 13)


#What % of non-SNP sets occur in two week of sample by two degree grid cell that has SNP sets with length 
#distribution?
#two week of sample by two degree grid cell that has SNP sets
hassets <- filter(fullbe, !is.na(`12`)) %>% distinct(twoweek_cellid_2p)

nonsnp <- filter(fullbe, is.na(`12`))

filter(nonsnp, twoweek_cellid_2p %in% hassets$twoweek_cellid_2p) %>% 
  nrow() / nrow(nonsnp) #0.964659

#What percent of non-SNP sets don't get a length distribution
filter(nonsnp, is.na(prop12hat)) %>% nrow() / nrow(nonsnp)  #0.0003643399

sessionInfo()
