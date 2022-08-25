#This is the second script cleaning BE data. It follows 1. match_be_landings.R
#Impute length distribution and calculate numjuv for each set, 
#In next script, calculate trip-level pj and compare to landing pj and correct pj, length distribution, etc. if applicable

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

#Load SNP Bitacora Electronica data
sbe <- read_csv("Data/bedat_snp.csv")

#Given length in cm, return weight in g
#Use values from IMARPE (2019)
lengthweight <- function(length){
  #length in cm; weight in g
  weight <- .0036*length^3.238
  return(weight)
}

#Size bins
sizebins <- c(names(sbe)[nchar(names(sbe)) %in% c(1,2)],
              grep("p5",names(sbe),value=T)
)

#Don't want ID column
sizebins <- sizebins[-grep("ID",sizebins)]

#Multiply proportion of individuals in each length bin by weight in grams 
#of individuals of that length
weightpropcols <- map(sizebins, function(x){
  
  mycol <- dplyr::select(sbe, x) 
  
  #Rename so can refer to directly
  names(mycol) <- "mybin"
  
  #If bin contains p, replace with "."
  if(grep("p", x) %>% length() > 0){
    
    binchar <- gsub("p","\\.",x)
    lengthbin <- as.numeric(binchar)
    
  }else{
    lengthbin <- as.numeric(x) 
  }
  
  mycol <- mutate(mycol, weightprop = mybin * lengthweight(lengthbin)) %>% 
    dplyr::select(weightprop)
  
  names(mycol) <- paste0("weightprop_", x)
  
  mycol
})

weightpropcols <- bind_cols(weightpropcols)

#Calculate average weight in grams of each row and add as new column to sbe
sbe <- mutate(sbe, avgweightg = apply(weightpropcols, 1, sum)) %>% 
  #Calculate number of individuals caught by each set
  mutate(numindivids = (tons*10^6)/avgweightg) %>%
    #Number of juveniles
    mutate(numjuv = numindivids*perjuv/100) %>%
    #Number of adults
    mutate(numadults = (numindivids*(100-perjuv))/100)

#Clean up
rm(weightpropcols, sizebins)

#Created in 1. match_be_landings.R
load("Output/Data/matched_be_landings_belevel.Rdata")

#Id column
matchedbe <- mutate(matchedbe, pbeid = 1:nrow(matchedbe))

sbe <- rename(sbe, slat = lat, slon = lon)
sbe <- rename(sbe, spj = perjuv, stons = tons)

sbe <- dplyr::select(sbe, -ID)

sbe <- distinct(sbe)

matchedbe <- mutate(matchedbe, caladate = date(FechaInicioCala))
sbe <- mutate(sbe, caladate = date(calatime))

#Make both sf objects
besf <- st_multipoint(cbind(matchedbe$lon, matchedbe$lat))
besf <- st_sfc(besf) %>% st_cast("POINT")
st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
matchedbe <- st_sf(geometry = besf, matchedbe)

ssf <- st_multipoint(cbind(sbe$slon, sbe$slat))
ssf <- st_sfc(ssf) %>% st_cast("POINT")
st_crs(ssf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
sbe <- st_sf(geometry = ssf, sbe)

rm(ssf, besf)


#Make a match function. Given row of matchedbe, find closest row in sbe to add length distribution values from
matchFun <- function(rowind){
  
  myrow <- matchedbe[rowind,]
  
  #Best match is same vessel-name and same calatime
  firstmatch <- filter(sbe, nave == myrow$Embarcacion & 
                         calatime>=myrow$FechaInicioCala & calatime<=myrow$FechaFinCala)
  
  #If match, get row to add length distribution values from
  if(nrow(firstmatch) > 0){
    
    #If multiple matches, pick one
    if(nrow(firstmatch) > 1){
      
      #Matching tons
      keepmatch <- filter(firstmatch, stons==myrow$betons)
      
      #If still multiple, pick the closer one
      if(nrow(keepmatch) > 1){
        
        dists <- st_distance(keepmatch, myrow) %>% as.numeric()
        minind <- which(dists==min(dists))
        keepmatch <- keepmatch[minind,]
        
        #If there are multiple sets of same minimum distance, pick the one with nearest percentage juvenile
        if(nrow(keepmatch)> 1){
          
          keepmatch <- mutate(keepmatch, pjdif = abs(spj - myrow$bepj))
          keepmatch <- keepmatch[which(keepmatch$pjdif==min(keepmatch$pjdif)),]
          
          #If still multiple rows, randomly pick one
          if(nrow(keepmatch) > 1){
            keepmatch <- sample_n(keepmatch, 1)
          }
          
        }
        
      } else if(nrow(keepmatch)==0){
        #If multiple same vessel-name and same time matches, but none are exact tons matches, 
        #keep the row that has closest tons
        keepmatch <- mutate(firstmatch, tonsdif = abs(stons - myrow$betons))
        keepmatch <- keepmatch[which(keepmatch$tonsdif==min(keepmatch$tonsdif)),]
        
        #If there are ties, pick the closer one
        if(nrow(keepmatch) > 1){
          
          dists <- st_distance(keepmatch, myrow) %>% as.numeric()
          minind <- which(dists==min(dists))
          keepmatch <- keepmatch[minind,]
          
          #If there are multiple sets of same minimum distance, pick the one with nearest percentage juvenile
          if(nrow(keepmatch)> 1){
            
            keepmatch <- mutate(keepmatch, pjdif = abs(spj - myrow$bepj))
            keepmatch <- keepmatch[which(keepmatch$pjdif==min(keepmatch$pjdif)),]
            
            #If still multiple rows, randomly pick one
            if(nrow(keepmatch) > 1){
              keepmatch <- sample_n(keepmatch, 1)
            }
            
          }
        }
      }
    }else if(nrow(firstmatch)==1){
      keepmatch <- firstmatch
    }
    
  } else{
    
    #Next match is same vessel-name and same date
    firstmatch <- filter(sbe, nave == myrow$Embarcacion & caladate==myrow$caladate)
    
    #If set in SNP data by same vessel-name on same day
    if(nrow(firstmatch) > 0){
      
      #If multiple matches, pick one
      if(nrow(firstmatch) > 1){
        
        #Matching tons
        keepmatch <- filter(firstmatch, stons==myrow$betons)
        
        #If still multiple, pick the closer one
        if(nrow(keepmatch) > 1){
          
          dists <- st_distance(keepmatch, myrow) %>% as.numeric()
          minind <- which(dists==min(dists))
          keepmatch <- keepmatch[minind,]
          
          #If there are multiple sets of same minimum distance, pick the one with nearest percentage juvenile
          if(nrow(keepmatch)> 1){
            
            keepmatch <- mutate(keepmatch, pjdif = abs(spj - myrow$bepj))
            keepmatch <- keepmatch[which(keepmatch$pjdif==min(keepmatch$pjdif)),]
            
            #If still multiple rows, randomly pick one
            if(nrow(keepmatch) > 1){
              keepmatch <- sample_n(keepmatch, 1)
            }
            
          }
          
        } else if(nrow(keepmatch)==0){
          #If multiple same vessel-name and same time matches, but none are exact tons matches, 
          #keep the row that is within 10 tons
          keepmatch <- mutate(firstmatch, tonsdif = abs(stons - myrow$betons))
          keepmatch <- keepmatch[which(keepmatch$tonsdif==min(keepmatch$tonsdif, na.rm=T)),]
          
          if(min(keepmatch$tonsdif) > 10){
            
            #If no matching rows within 10 tons, check whether any matching rows within 20 km
            dists <- st_distance(firstmatch, myrow) %>% as.numeric()
            minind <- which(dists==min(dists))
            keepmatch <- firstmatch[minind,]
            
            #If no sets within 20 km, don't match row to any SNP set
            if(min(dists)>20000){
              
              keepmatch <- rep(as.numeric(NA), 34) %>% t() %>% as.data.frame()
              names(keepmatch) <- c("stons","spj",grep("^[0-9]",names(sbe),value=T))
              
            } else if(nrow(keepmatch) > 1){
              #If there are multiple same-distance sets within 20 km and same tonsdif, pick the one with nearest percentage juvenile
              keepmatch <- mutate(keepmatch, pjdif = abs(spj - myrow$bepj))
              keepmatch <- keepmatch[which(keepmatch$pjdif==min(keepmatch$pjdif)),]
              
              #If still multiple rows, randomly pick one
              if(nrow(keepmatch) > 1){
                keepmatch <- sample_n(keepmatch, 1)
              }
              
              
            }
            
          } else if(nrow(keepmatch) > 1 & min(keepmatch$tonsdif, na.rm=T) <= 10){#If tonsdif less than 10 and there are ties, pick the closer one
            dists <- st_distance(keepmatch, myrow) %>% as.numeric()
            minind <- which(dists==min(dists))
            keepmatch <- keepmatch[minind,]
            
            #If there are multiple sets of same minimum distance, pick the one with nearest percentage juvenile
            if(nrow(keepmatch)> 1){
              
              keepmatch <- mutate(keepmatch, pjdif = abs(spj - myrow$bepj))
              keepmatch <- keepmatch[which(keepmatch$pjdif==min(keepmatch$pjdif)),]
              
              #If still multiple rows, randomly pick one
              if(nrow(keepmatch) > 1){
                keepmatch <- sample_n(keepmatch, 1)
              }
              
            }
          }
        }
      }else if(nrow(firstmatch)==1){
        keepmatch <- firstmatch
      }
      
    } else{
      #Add NA rows because no matching set in SNP data for same vessel on same day
      keepmatch <- rep(as.numeric(NA), 34) %>% t() %>% as.data.frame()
      names(keepmatch) <- c("stons","spj",grep("^[0-9]",names(sbe),value=T))
    }
    
  }
  
  
  #Bind columns of keepmatch onto myrow, dropping geometry columns
  keepmatch <- as.data.frame(keepmatch) %>% dplyr::select("stons","spj",grep("^[0-9]",names(sbe),value=T))
  
  out <- as.data.frame(myrow) %>% dplyr::select(-geometry) %>% bind_cols(keepmatch)
  
  return(out)
}

(myCores <- detectCores())

cl <- makeCluster(18)

clusterSetRNGStream(cl, 20200615)

clusterExport(cl, "matchedbe")
clusterExport(cl, "sbe")
clusterExport(cl, "matchFun")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))


matchlist <- parLapply(cl = cl,
                       1:nrow(matchedbe),
                       function(x){

                         try(matchFun(x))

                       })

stopCluster(cl)
rm(cl, myCores)

matched <- bind_rows(matchlist)

#Several rows in matched have NA spj but from looking at size distribution I can see there are no juveniles
#Set spj for these rows equal to 0
matched$spj[is.na(matched$spj) & !is.na(matched$stons)] <- 0

#Several rows have NA spj and NA stons, but I can calculate from their length distribution what spj should be
for(i in which(is.na(matched$spj) & !is.na(matched$`12`))){
  rowspj <- matched[i,] %>%
    dplyr::select(c(`3`,`3p5`,`4`,`4p5`,`5`,`5p5`,`6`,`6p5`,`7`,`7p5`,`8`,`8p5`,`9`,`9p5`,
                    `10`,`10p5`,`11`,`11p5`)) %>%
    as.matrix() %>% as.numeric() %>% sum()

  #Multiply by 100 to go from proportion to percentage (other rows coded as 0 to 100; percentage)
  matched$spj[i] <- rowspj*100

  #And assign stons equal to betons for these rows
  matched$stons[i] <- matched$betons[i]
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

#Row without hat for grabbing columns in matched
sizedf <- mutate(sizedf, charlength = gsub("hat","",lengthname))

#Calculate avgweightg, numindivids, numjuv, and numadults for matched rows
calcIndivids <- function(rowind){
  
  myrow <- matched[rowind,]
  
  if(!is.na(myrow$spj)){
    
    myrow$avgweightg <- sapply(sizedf$charlength, function(x){
        myrow[,x]*sizedf$weight[sizedf$charlength==x]
      }) %>% sum()
    
    myrow <- myrow %>%
      mutate(numindivids = (betons*10^6)/avgweightg) %>%
      #Number of juveniles
      mutate(numjuv = numindivids*spj/100) %>%
      #Number of adults
      mutate(numadults = numindivids-numjuv) 
    
  } else{
    
    myrow <- mutate(myrow, avgweightg = as.numeric(NA), numindivids = as.numeric(NA),
                    numjuv = as.numeric(NA), numadults = as.numeric(NA))
    
  }
  
  return(myrow)
}

#Apply over rows
(myCores <- detectCores())

cl <- makeCluster(18)

clusterExport(cl, "matched")
clusterExport(cl, "sizedf")
clusterExport(cl, "calcIndivids")
clusterEvalQ(cl, library(dplyr))

matchlist <- parLapply(cl = cl,
                       1:nrow(matched),
                       function(x){

                         try(calcIndivids(x))

                       })

stopCluster(cl)
rm(cl, myCores)

matched <- bind_rows(matchlist)

rm(matchlist)

#One row was matched to an spj set with 0 in each length interval column; set these values = NA for this row
zerorow <- which(matched$avgweightg==0)

matched[zerorow,sizedf$charlength] <- as.numeric(NA)
matched[zerorow,c("stons","spj","avgweightg","numindivids","numjuv","numadults")] <- as.numeric(NA)

rm(zerorow)

#There are two rows with 3 slightly greater than 0 that have bepj = 0 and spj > 0
#This causes problems in shiftSNP because length distribution rows keep shifting until 
#updated pj is 0. This results in most of the length distribution being shifted out
#To avoid this problem, set 3 = 0 for these two rows. Actual proportion in this length interval
#for these two rows is both less than .007 and the other length distribution rows will be inflated
#so that length distribution adds up to 1
matched$`3`[matched$`3`>0] <- 0

#Before imputing size distribution for observations without one, update SNP distribution if bepj != spj
shiftSNP <- function(rowind){
  
  myrow <- matched[rowind,]
  
  if(!is.na(myrow$spj)){
  
  #New columns to hold shifted size distribution values
  holddf <- dplyr::select(myrow,"3","3p5","4","4p5","5","5p5","6","6p5","7","7p5","8","8p5",
                          "9","9p5","10","10p5","11","11p5","12","12p5","13","13p5","14","14p5",
                          "15","15p5","16","16p5","17","17p5","18","18p5") %>% 
    gather("length","prop") %>%
    mutate(lengthnum = gsub("p",".",length) %>% as.numeric())
  
  #If prop sums to less than 1, inflate all positive rows so that prop adds up to 1
  if(sum(holddf$prop) < 1){
    holddf <- mutate(holddf, prop = if_else(prop > 0, prop*(1/sum(holddf$prop)), prop))
  }

  holddf <- mutate(holddf, propshift = prop) 
  
  #Difference between bepj and spj
  pjdif <- myrow$bepj - myrow$spj
  
  if(pjdif > 0.00001){
    uppj <- myrow$spj
    
    #First condition is for weird numeric sensitivity where one row has uppj = 100 but returns true for uppj < myrow$bepj (=100)
    #Second condition is for another row which gets shifted so far down that one positive value is lost. 
    #Would rather this not happen at cost of uppj not equaling bepj for this row
    while(uppj < (myrow$bepj-.0000000001) & holddf$propshift[holddf$lengthnum==3] == 0){
      
      #Shift length distribution down one half cm interval and calculate new pj
      holddf <- mutate(holddf, propshift = lead(propshift, 1))
      
      #New pj 
      uppj <- sum(holddf$propshift[holddf$lengthnum < 12], na.rm=T)*100
      
    }
    
    #Shift length distribution up one half cm interval (up = smaller)
    #These two shifted distributions bound the "true" pj (bepj)
    #Number of shifts 
    numshifts <- is.na(holddf$propshift) %>% sum()
    
    #So shift one less than this 
    holddf <- mutate(holddf, propshift2 = lead(prop, numshifts - 1))
    
    #Use propshift that is closest to bepj
    if( abs(sum(holddf$propshift[holddf$lengthnum < 12], na.rm=T) - myrow$bepj/100) <= 
        abs(myrow$bepj/100 - sum(holddf$propshift2[holddf$lengthnum < 12], na.rm=T))
    ){
      holddf <- rename(holddf, prophat = propshift)
    } else{
      holddf <- rename(holddf, prophat = propshift2)
    }
    
    holddf <- mutate(holddf, length = paste0(length,"hat")) %>%
      dplyr::select(length, prophat) %>% 
      mutate(prophat = if_else(is.na(prophat),0,prophat)) #%>% 
    #spread(length, prophat)
    
  }else if(pjdif < -0.00001){
    uppj <- myrow$spj
    
    while(uppj > myrow$bepj){
      
      #Shift length distribution down one half cm interval and calculate new pj
      holddf <- mutate(holddf, propshift = lag(propshift, 1))
      
      #New pj 
      uppj <- sum(holddf$propshift[holddf$lengthnum < 12], na.rm=T)*100
      
    }
    
    #Shift length distribution down one half cm interval (down = bigger)
    #These two shifted distributions bound the "true" pj (bepj)
    #Number of shifts 
    numshifts <- is.na(holddf$propshift) %>% sum()
    #So shift one less than this 
    holddf <- mutate(holddf, propshift2 = lag(prop, numshifts - 1))
    
    #Use propshift that is closest to bepj
    if( abs(sum(holddf$propshift[holddf$lengthnum < 12], na.rm=T) - myrow$bepj/100) <= 
        abs(myrow$bepj/100 - sum(holddf$propshift2[holddf$lengthnum < 12], na.rm=T))
    ){
      holddf <- rename(holddf, prophat = propshift)
    } else{
      holddf <- rename(holddf, prophat = propshift2)
    }
    
    holddf <- mutate(holddf, length = paste0(length,"hat")) %>%
      dplyr::select(length, prophat) %>% 
      mutate(prophat = if_else(is.na(prophat),0,prophat)) #%>% 
    #spread(length, prophat)
    
  } else{
    holddf <- mutate(holddf, length = paste0(length,"hat")) %>%
      dplyr::select(length, prop) %>% 
      rename(prophat = prop) #%>% 
    #spread(length, prophat)
  }
  
  #Bind onto myrow
  myrow <- bind_cols(myrow, 
                     spread(holddf, length, prophat)
  )
  
  #Now calculate avgweightg, numindivids, numjuv, and numadults for this row
  myrow$avgweightg <- sapply(holddf$length, function(x){
      myrow[,x]*sizedf$weight[sizedf$lengthname==x]
    }) %>% sum()
  
  myrow <- myrow %>%
    mutate(numindivids = (betons*10^6)/avgweightg) %>%
    #Number of juveniles
    mutate(numjuv = numindivids*bepj/100) %>%
    #Number of adults
    mutate(numadults = numindivids-numjuv) 
  
  } else{
    
    #Add NA columns for length distribution for rows not matched to SNP set
    nacols <- matrix(as.numeric(NA),ncol=length(sizedf$lengthname), nrow=1) %>% 
      as.data.frame()
    names(nacols) <- sizedf$lengthname
    
    myrow <- bind_cols(myrow, nacols)
  }
  
  return(myrow)
}

#Apply over rows of matched observations
(myCores <- detectCores())

cl <- makeCluster(18)

clusterExport(cl, "matched")
clusterExport(cl, "sizedf")
clusterExport(cl, "shiftSNP")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(tidyr))


mlist <- parLapply(cl = cl,
                   1:nrow(matched),
                   function(x){
                     
                     try(shiftSNP(x))
                     
                   })

stopCluster(cl)
rm(cl, myCores)

matched <- bind_rows(mlist)

rm(mlist)

#Drop sets that occur in Southern Zone
matched <- filter(matched, Zona != "Sur")

###For non-matched observations, impute size distribution from matched SNP observations
#at level of two-week-of-sample by two degree grid cell. 

#Make matched and nonmatched sf
mp <- st_multipoint(cbind(matched$lon, matched$lat))
mp <- st_sfc(mp) %>% st_cast("POINT")
st_crs(mp) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
matched <- st_sf(geometry = mp, matched)

rm(mp)

#Created in 1. match_be_landings.R
load("Output/Data/grid2p.Rdata")

#Intersect 
insidecell <- st_intersects(matched, grid2p)

insidecell <- as.numeric(as.character(insidecell))

#Add cellid column onto matched
matched <- mutate(matched, cellid_2p = insidecell)

rm(insidecell)

#Create two-week of sample variable.
matched <- dplyr::select(matched, -twowk)

twoweek <- mutate(as.data.frame(matched) %>% dplyr::select(-geometry),
                  week = week(caladate), year = year(caladate))

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

#Bind onto matched
matched <- mutate(matched, twowk = twoweek$twowk)

#Now create variable to group on
matched <- mutate(matched, twoweek_cellid_2p = paste0(twowk, "_", cellid_2p))

rm(counter, i, twowk, twoweek)

#Compute empirical size distribution of matched observations 
#(proportion of individuals in each size interval for each two-week by two degree grid cell)

#Add prop to beginnning of name of size bins
names(matched)[names(matched) %in% sizedf$lengthname] <- 
  paste0("prop",
         names(matched)[names(matched) %in% sizedf$lengthname]
  )

#Names of variables that are less than 12 cm
juvvars <- paste0("prop",sizedf$lengthname[sizedf$lengthnum < 12] %>% as.character())

#Names of variables that are 12 cm or above
adultvars <- paste0("prop",sizedf$lengthname[sizedf$lengthnum >= 12] %>% as.character())

#Drop geometry column from matched for now (will add it back in later)
matched <- as.data.frame(matched) %>% dplyr::select(-geometry)

makeEmp <- function(mygroup){
  
  #Filter to matched in this group
  myobs <- filter(matched, !is.na(prop12hat) & twoweek_cellid_2p ==mygroup)
  
  out <- myobs %>%
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
              prop18p5hat = sum(prop18p5hat*numindivids)/sum(numindivids)) %>%
    mutate(twoweek_cellid_2p = mygroup)
  
  #Also calculate average pj
  out <- mutate(out, avgpj = sum(out[,juvvars])*100)
  
  return(out)
}

#Apply over groups
empdf <- map_df(unique(matched$twoweek_cellid_2p), function(x){
  makeEmp(x)
})

#Calculate population-level distribution for use when nonmatched observations has pj > 0 
#in a group that has pj = 0
popdist <- filter(matched, !is.na(prop12hat)) %>%
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
  

#Also calculate average pj
popdist <- mutate(popdist, avgpj = sum(popdist[,juvvars])*100)

#Calculate two-week level distribution when grid cell in that two week does not have any SNP observations
makeEmptwoweek <- function(mytwoweek){
  
  #Filter to matched in this group
  myobs <- filter(matched, !is.na(prop12hat) & twowk == mytwoweek)
  
  out <- myobs %>%
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
              prop18p5hat = sum(prop18p5hat*numindivids)/sum(numindivids)) %>%
    mutate(twowk = mytwoweek)
  
  #Also calculate average pj
  out <- mutate(out, avgpj = sum(out[,juvvars])*100)
  
  return(out)
}

twoweekdist <- map_df(unique(matched$twowk[!is.na(matched$prop12hat)]), function(x){
  makeEmptwoweek(x)
})



#Now given group, filter to non-matched rows and impute size distribution from 
#percentage juvenile using size distribution of matched rows in that group
imputeGroup <- function(mygroup){
  
  #Empirical distribution for this group
  empdist <- filter(empdf, twoweek_cellid_2p==mygroup)
  
  #Non-matched rows
  nonmatched <- filter(matched, is.na(prop12hat) & twoweek_cellid_2p==mygroup)
  
  if(nrow(nonmatched)>0 & !is.na(empdist$prop12hat)){
    
    #Apply imputeFun over rows
    out <- map_df(1:nrow(nonmatched), function(x){
      imputeFun(x, nonmatched, empdist)
    })
  }else if(nrow(nonmatched)>0 & is.na(empdist$prop12hat)){
    #If no SNP in that twoweek-cell, use distribution for twoweek
    if(unique(nonmatched$twowk) %in% twoweekdist$twowk){
      out <- map_df(1:nrow(nonmatched), function(x){
        imputeFun(x, nonmatched, filter(twoweekdist, twowk==unique(nonmatched$twowk)))
      })
    } else{
      #Leave values as NA if no SNP in that twowk
      out <- nonmatched
    }
  } else{
    out <- nonmatched
  }
  
  return(out)
}


imputeFun <- function(nonmatchedrow, nonmatched, empdist){
  
  myrow <- nonmatched[nonmatchedrow,]
  
  #Re-weight empirical distribtion by actual percentage juvenile
  if(empdist$avgpj>0){
    
    myrow[,juvvars] <- empdist[,juvvars] * myrow$bepj/empdist$avgpj
    
    myrow[,adultvars] <- empdist[,adultvars] * 
      (100-myrow$bepj)/(100-empdist$avgpj)
    
  } else if(empdist$avgpj==0 & myrow$bepj==0){
    
    #Don't need to reweight in this case
    myrow[,juvvars] <- empdist[,juvvars] 
    
    myrow[,adultvars] <- empdist[,adultvars] 
    
  } else{
    #If bepj>0 and avgpj=0, use population-level distribution
    myrow[,juvvars] <- popdist[,juvvars] * myrow$bepj/popdist$avgpj
    
    myrow[,adultvars] <- popdist[,adultvars] * 
      (100-myrow$bepj)/(100-popdist$avgpj)
    
  }
  
  
  #Now calculate avgweightg, numindivids, numjuv, and numadults for this row
  myrow$avgweightg <- sapply(c(juvvars, adultvars), function(x){
      myrow[,x]*sizedf$weight[sizedf$lengthname==gsub("prop","",x)]
    }) %>% sum()
  
  myrow <- myrow %>%
    mutate(numindivids = (betons*10^6)/avgweightg) %>%
    #Number of juveniles
    mutate(numjuv = numindivids*bepj/100) %>%
    #Number of adults
    mutate(numadults = numindivids-numjuv)
  
  return(myrow)
}

sizedf$lengthname <- as.character(sizedf$lengthname)

#Apply over groups
#Parallel apply over datevec
(myCores <- detectCores())

cl <- makeCluster(20)

clusterExport(cl, "empdf")
clusterExport(cl, "matched")
clusterExport(cl, "sizedf")
clusterExport(cl, "popdist")
clusterExport(cl, "twoweekdist")
clusterExport(cl, "adultvars")
clusterExport(cl, "juvvars")
clusterExport(cl, "lengthweight")
clusterExport(cl, "imputeFun")
clusterExport(cl, "imputeGroup")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(purrr))
clusterEvalQ(cl, library(tidyr))


nolist <- parLapply(cl = cl,
                    unique(empdf$twoweek_cellid_2p),
                    function(x){
                      
                      imputeGroup(x)
                      
                    })

nodf <- bind_rows(nolist)

stopCluster(cl)
rm(cl, myCores)

#Bind with matched rows
fullbe <- bind_rows(filter(matched, !is.na(prop12hat)),
                    nodf) 

nrow(matched)==nrow(fullbe)

save(fullbe, file = "Output/Data/pbe_imp_uncorrected.Rdata")
#That's enough for one script! Next script is 3. correct_be.R

sessionInfo()
