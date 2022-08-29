#Instead of imputing length distribution of non-SNP sets at the 
#two-week-of-sample by two-degree grid cell, do so at the 
#one week of sample by one degree grid cell. Then 
#re-run main estimation to see how effects of policy change.

#This script is going to be very long. 
#I can start from the output of 1. match_be_landings.R
#But then I need to modify 2. impute_size_be.R, 
#and in order run 3. correct_be.R, 4. make_rddf.R, and 
#make_figure7.R

#modify 2. impute_size_be.R
{
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
  
  cl <- makeCluster(4)
  
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
  
  cl <- makeCluster(4)
  
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
  
  cl <- makeCluster(4)
  
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
  #at level one week by one degree level, rather than two by two
  
  #Make matched and nonmatched sf
  mp <- st_multipoint(cbind(matched$lon, matched$lat))
  mp <- st_sfc(mp) %>% st_cast("POINT")
  st_crs(mp) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
  matched <- st_sf(geometry = mp, matched)
  
  rm(mp)
  
  eez <- st_read("Data/Intersect_IHO_EEZ_v2_2012/eez.shp") %>%
    filter(EEZ == "Peruvian Exclusive Economic Zone (disputed - Peruvian point of view)") %>% 
    st_transform(crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  grid1p <- st_make_grid(
    st_sfc(st_polygon(list(rbind(
      c(st_bbox(eez)["xmin"],st_bbox(eez)["ymax"]), 
      c(st_bbox(eez)["xmax"],st_bbox(eez)["ymax"]), 
      c(st_bbox(eez)["xmax"],st_bbox(eez)["ymin"]), 
      c(st_bbox(eez)["xmin"],st_bbox(eez)["ymin"]),
      c(st_bbox(eez)["xmin"],st_bbox(eez)["ymax"])
    ))),
    crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")),
    cellsize = 1, what = 'polygons') %>%
    st_sf()
  
  #Keep only grid cells that intersect eez
  inter <- st_intersects(eez, grid1p) %>% unlist()
  
  grid1p <- grid1p[inter,]
  
  rm(inter)
  
  #Add a cellid column
  grid1p <- mutate(grid1p, cellid = 1:nrow(grid1p))
  
  #Intersect 
  insidecell <- st_intersects(matched, grid1p)
  
  insidecell <- as.numeric(as.character(insidecell))
  
  #Add cellid column onto matched
  matched <- mutate(matched, cellid_1p = insidecell)
  
  rm(insidecell)
  
  #Also need two degree grid cell for clustering standard errors
  #Created in 1. match_be_landings.R
  load("Output/Data/grid2p.Rdata")
  
  #Intersect 
  insidecell <- st_intersects(matched, grid2p)
  
  insidecell <- as.numeric(as.character(insidecell))
  
  #Add cellid column onto matched
  matched <- mutate(matched, cellid_2p = insidecell)
  
  rm(insidecell)

  #Create two-week of sample variable.
  #Still need this later for clustering standard errors
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
  
  #Now create one-week-of-sample variable
  matched <- mutate(matched, 
                    onewk = paste0(week(caladate), "_", year(caladate)))
  
  #one week by one degree variable
  matched <- mutate(matched, oneweek_cellid_1p = paste0(onewk, "_", cellid_1p))
  
  
  
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
    myobs <- filter(matched, !is.na(prop12hat) & oneweek_cellid_1p == mygroup)
    
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
      mutate(oneweek_cellid_1p = mygroup)
    
    #Also calculate average pj
    out <- mutate(out, avgpj = sum(out[,juvvars])*100)
    
    return(out)
  }
  
  #Apply over groups
  empdf <- map_df(unique(matched$oneweek_cellid_1p), function(x){
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
  
  #Calculate one-week level distribution when grid cell in that one week does not have any SNP observations
  makeEmponeweek <- function(myoneweek){
    
    #Filter to matched in this group
    myobs <- filter(matched, !is.na(prop12hat) & onewk == myoneweek)
    
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
      mutate(onewk = myoneweek)
    
    #Also calculate average pj
    out <- mutate(out, avgpj = sum(out[,juvvars])*100)
    
    return(out)
  }
  
  oneweekdist <- map_df(unique(matched$onewk[!is.na(matched$prop12hat)]), function(x){
    makeEmponeweek(x)
  })
  
  
  
  #Now given group, filter to non-matched rows and impute size distribution from 
  #percentage juvenile using size distribution of matched rows in that group
  imputeGroup <- function(mygroup){
    
    #Empirical distribution for this group
    empdist <- filter(empdf, oneweek_cellid_1p==mygroup)
    
    #Non-matched rows
    nonmatched <- filter(matched, is.na(prop12hat) & oneweek_cellid_1p==mygroup)
    
    if(nrow(nonmatched)>0 & !is.na(empdist$prop12hat)){
      
      #Apply imputeFun over rows
      out <- map_df(1:nrow(nonmatched), function(x){
        imputeFun(x, nonmatched, empdist)
      })
    }else if(nrow(nonmatched)>0 & is.na(empdist$prop12hat)){
      #If no SNP in that oneweek-cell, use distribution for oneweek
      if(unique(nonmatched$onewk) %in% oneweekdist$onewk){
        out <- map_df(1:nrow(nonmatched), function(x){
          imputeFun(x, nonmatched, filter(oneweekdist, onewk==unique(nonmatched$onewk)))
        })
      } else{
        #Leave values as NA if no SNP in that onewk
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
  
  cl <- makeCluster(2)
  
  clusterExport(cl, "empdf")
  clusterExport(cl, "matched")
  clusterExport(cl, "sizedf")
  clusterExport(cl, "popdist")
  clusterExport(cl, "oneweekdist")
  clusterExport(cl, "adultvars")
  clusterExport(cl, "juvvars")
  clusterExport(cl, "lengthweight")
  clusterExport(cl, "imputeFun")
  clusterExport(cl, "imputeGroup")
  clusterEvalQ(cl, library(dplyr))
  clusterEvalQ(cl, library(purrr))
  clusterEvalQ(cl, library(tidyr))
  
  
  nolist <- parLapply(cl = cl,
                      unique(empdf$oneweek_cellid_1p),
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
  
  save(fullbe, file = "Output/TempData/pbe_imp_uncorrected_onebyone.Rdata")

}


#Now run 3. correct_be.R with pbe_imp_uncorrected_onebyone.Rdata created in previous chunk
{
  rm(list=ls())
  
  library(dplyr); library(readxl); library(ggplot2)
  library(sf); library(purrr); library(lubridate)
  library(parallel); library(tidyr); library(furrr)
  
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
  load("Output/TempData/pbe_imp_uncorrected_onebyone.Rdata")
  
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
  
  cl <- makeCluster(4)
  
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

      out <- bind_rows(out)

      return(out)
    }

    #Apply over all buffers
    (myCores <- detectCores())

    cl <- makeCluster(4)

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

    cbufs <- bind_rows(cbufs)

    names(cbufs)[names(cbufs)!="geometry"] <- names(closed)[names(closed)!="geometry"]

    closed <- bind_rows(closed, cbufs)

    rm(cbufs)

    #Duplicate rectangles for leads and lags
    closed <- bind_rows(
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

  cl <- makeCluster(4)

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
  
  save(fullbe, file = "Output/TempData/pbe_imp_onebyone.Rdata")
  
}


#Now run 4. make_rddf.R with pbe_imp_onebyone.Rdata created in previous chunk
{#Make potential closures
  
  rm(list=ls())
  
  library(dplyr); library(readxl); library(ggplot2)
  library(rworldmap); library(sf); library(lwgeom)
  library(rgdal); library(geosphere); library(sp)
  library(purrr); library(lubridate); library(glmnet)
  library(lfe); library(Formula); library(smoothr)
  library(parallel); library(furrr)
  
  options(scipen=999)
  
  `%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
  
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
  
  #Peru time
  Sys.setenv(TZ='America/Lima')
  
  #Load full BE from PRODUCE where I have imputed size distribution of non-SNP observations
  #Created in 3. correct_be.R
  load("Output/TempData/pbe_imp_onebyone.Rdata")
  
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
      
      myrects <- bind_rows( myrects)
      
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
      
      newrects <- bind_rows( newrects)
      
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
  
  cl <- makeCluster(4)
  
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
  rdf <- bind_rows( rdf)
  
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
    
    out <- bind_rows(out)
    
    return(out)
  }
  
  #Apply over all buffers
  (myCores <- detectCores())
  
  cl <- makeCluster(4)
  
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
  
  rbufs <- bind_rows(rbufs)
  
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
    
    out <- bind_rows(out)
    
    return(out)
  }
  
  #Apply over all buffers
  (myCores <- detectCores())
  
  cl <- makeCluster(4)
  
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
  
  cbufs <- bind_rows(cbufs)
  
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
    
    rdf <- bind_rows( rdf)
    
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
    
    out <- bind_rows(mylist)
    
    return(out)
  }
  
  rdf <- lapply(unique(rdf$bin),
                function(x){
                  
                  try(treatBin(x))
                  
                })
  
  rdf <- bind_rows(rdf)
  
  rddf <- mutate(rdf, startdate = as.Date(start) %>% as.factor(),
                 season = as.factor(season))
  
  #Create unique identifier for each rectangle
  rid <- as.data.frame(rddf) %>% 
    dplyr::select(newclust, tons, nobs, area_km2, season) %>%
    distinct()
  
  rid <- mutate(rid, rid = 1:nrow(rid))
  
  #potential closures
  nrow(rid)
  
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
  
  #Number of clusters: 
  unique(rddf$twoweek_cellid_2p) %>% length() #255
  
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
  
plan(multisession, workers = 6)
  
  #Apply over rows of rddf
  myoutcomes <- future_map(
                          1:nrow(rddf),
                          function(x){
                            
                            outcomesFun(x)
                            
                          })
  
  stopCluster(cl)
  rm(cl, myCores)
  
  myoutcomes <- bind_rows(myoutcomes)
  
  #Need to first rename other tons variables for clarity
  rddf <- rename(rddf, clusttons = tons, 
                 clustnumindivids = numindivids, clustnumjuv = numjuv, 
                 clustnumadults = numadults, clustnobs = nobs, 
                 clustarea_km2 = area_km2)
  
  #Join onto rddf
  rddf <- left_join(rddf, myoutcomes, by = c("bin",'rid'))
  
  #Calculate number of caught adults in bin
  rddf <- mutate(rddf, numadults = numindivids - numjuv)
  
  #Save rddf
  save(rddf, file = "Output/TempData/rddf_onebyone.Rdata")
}


#Now finally estimate Equation 1 and calculate change in total juvenile catch from policy
#make_figure7.R
{
  rm(list=ls())
  
  library(dplyr); library(ggplot2)
  library(sf);library(lubridate); library(glmnet)
  library(lfe); library(Formula); library(xtable)
  library(purrr); library(cowplot); library(tidyr)
  library(latex2exp); library(scales)
  
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
  
  #Load rddf created in 4. make_rddf.R
  load("Output/TempData/rddf_onebyone.Rdata")
  
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
  
  #How many clusters are there
  unique(rddf$twoweek_cellid_2p) %>% length() #255
  
  #How many observations
  nrow(rddf) #34,164 (so dropped 21 potential closures)
  
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
  sum(toteffect_juv$chmjuv) / sum(toteffect_juv$juv0) #0.746729
  
  #Scaled total effect see above comment. 
  #Created in 3. correct_be.R
  load("Output/TempData/pbe_imp_onebyone.Rdata")
  
  #Total effect, not accounting for reallocation
  #Rescale (potential closures nummjuv larger because of double-counting.)
  (sum(fullbe$numjuv,na.rm=T) / sum(rddf$numjuv, na.rm=T)) * sum(toteffect_juv$chmjuv) #60162.01
  
  toteffect_juv <- dplyr::select(toteffect_juv, -`t value`, -`Pr(>|t|)`)
  
  toteffect_juv <- rename(toteffect_juv, juvcoef = Estimate, juvse = `Cluster s.e.`)
  
  #Calculate delta method standard errors for change in millions of juveniles caught
  library(msm)
  
  toteffect_juv <- mutate(toteffect_juv, chmjuvse = as.numeric(NA), 
                          chmjuv_scaled = chmjuv * (sum(fullbe$numjuv,na.rm=T) / sum(rddf$numjuv, na.rm=T)),
                          chmjuvse_scaled = as.numeric(NA), 
                          juv0se = as.numeric(NA))
  
  scaleconstant <- (sum(fullbe$numjuv,na.rm=T) / sum(rddf$numjuv, na.rm=T)) #.394
  
  for(mybin in toteffect_juv$bin){
    
    #Filter toteffect_juv to mybin
    mydf <- filter(toteffect_juv, bin==mybin)
    
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
  
  #Decomposing effects:
  #Direct effect
  toteffect_juv$chmjuv_scaled[toteffect_juv$bin=="active_in"]  #-2176.79
  
  #Inside, post announcement
  toteffect_juv$chmjuv_scaled[toteffect_juv$bin=="lead9hours_in"] #1182.788
  
  
  #Contemporaneous spatial spillovers
  toteffect_juv$chmjuv_scaled[toteffect_juv$bin=="active_10" | toteffect_juv$bin=="active_20" | toteffect_juv$bin=="active_30" | 
                                toteffect_juv$bin=="active_40" | toteffect_juv$bin=="active_50"] %>% sum() #44478.46
  
  #Inside, day after
  toteffect_juv$chmjuv_scaled[toteffect_juv$bin=="lag1_in"]  #2262.906
  
  #Inside, two days after
  toteffect_juv$chmjuv_scaled[toteffect_juv$bin=="lag2_in"]  #2258.567
  
  #Total effect, not accounting for reallocation
  sum(toteffect_juv$chmjuv_scaled) #60163.7
  
  #Closures cannot increases tons caught because of TAC, so account for this reallocation
  #by estimating how policy affects tons caught on average across treatment bins in treatment window
  tonscaught <- felm(
    as.Formula(paste0(
      "asinhtons", "~ ", 
      "treatfrac ",
      " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
      paste0(grep("prop",names(rddf),value=T),collapse="+"),
      "| bin + twowk:cellid_2p + startdate",
      " | 0 | twoweek_cellid_2p")),
    data =rddf)
  
  #Increase tons by 35% (exp( 0.29933719 )-1)
  summary(tonscaught)[["coefficients"]]["treatfrac",]
  
  tonscoef <- summary(tonscaught)[["coefficients"]]["treatfrac","Estimate"]
  
  #Don't reallocate tons from 2017 second season or 2019 second season because they were both 
  #shut down well before TAC was reached. 
  
  #Tons caught in state of world I observe in seasons where TAC hit
  tons1 <- sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"])
  
  #Change in tons because of policy
  (ctons <- tons1 - tons1/exp(tonscoef)) #2722546
  
  #Average pj outside of treatment window
  (avgpjoutside <- filter(fullbe, lead_0==0 & lead_10==0 & lead_20==0 & lead_30==0 & lead_40==0 & lead_50==0 & 
                            active_0==0 & active_10==0 & active_20==0 & active_30==0 & active_40==0 & active_50==0 & 
                            lag1_0==0 & lag1_10==0 & lag1_20==0 & lag1_30==0 & lag1_40==0 & lag1_50==0 & 
                            lag2_0==0 & lag2_10==0 & lag2_20==0 & lag2_30==0 & lag2_40==0 & lag2_50==0 & 
                            lag3_0==0 & lag3_10==0 & lag3_20==0 & lag3_30==0 & lag3_40==0 & lag3_50==0 & 
                            lag4_0==0 & lag4_10==0 & lag4_20==0 & lag4_30==0 & lag4_40==0 & lag4_50==0 & 
                            !is.na(numindivids) & !is.na(bepjhat) & Temporada!="2017-II" & Temporada!="2019-II") %>%
      #Weight by number of individuals
      mutate(pjweighted = bepjhat*numindivids) %>%  
      summarise(perjuv = sum(pjweighted)/sum(numindivids)) %>% as.numeric() / 100) #0.09045436
  
  
  #Avg weight of individual caught outside of treatment window
  (avgweightoutside <- filter(fullbe, lead_0==0 & lead_10==0 & lead_20==0 & lead_30==0 & lead_40==0 & lead_50==0 & 
                                active_0==0 & active_10==0 & active_20==0 & active_30==0 & active_40==0 & active_50==0 & 
                                lag1_0==0 & lag1_10==0 & lag1_20==0 & lag1_30==0 & lag1_40==0 & lag1_50==0 & 
                                lag2_0==0 & lag2_10==0 & lag2_20==0 & lag2_30==0 & lag2_40==0 & lag2_50==0 & 
                                lag3_0==0 & lag3_10==0 & lag3_20==0 & lag3_30==0 & lag3_40==0 & lag3_50==0 & 
                                lag4_0==0 & lag4_10==0 & lag4_20==0 & lag4_30==0 & lag4_40==0 & lag4_50==0 & 
                                !is.na(numindivids) & !is.na(avgweightg) & Temporada!="2017-II" & Temporada!="2019-II") %>%
      #Weight by number of individuals
      mutate(weightweighted = avgweightg*numindivids) %>% 
      summarise(avgweightg = sum(weightweighted)/sum(numindivids)) %>% as.numeric())
  
  
  #Decrease in individuals caught outside of treatment window in millions
  #(converting tons to g cancels out conversion to millions)
  chindividsoutside <- -ctons/avgweightoutside
  
  chjuvsoutside <- chindividsoutside*avgpjoutside
  
  #Now can calculate change in juvenile catch due to policy, accounting for reallocation
  (chmjuvsstart <- changejuv + chjuvsoutside) #46369.62, compared to 46829.39
  
  #How many juveniles are caught during my sample period in total?
  #F(1)*pj*individuals/VMS fishing obs
  juv1 <- sum(fullbe$numjuv, na.rm=T) / 10^6 
  
  #chmjuvsstart = juv(1) - juv(0)
  juv0 <- juv1 - chmjuvsstart 
  
  #Then increase in juvenile catch as a percentage is 
  chmjuvsstart / juv0 #0.4906038, compared to 0.4986941
  
  #Calculate standard error on total change in juvenile catch and in total percentage change
  mycoefs <- toteffect_juv$chmjuv_scaled
  mybigvcov <- diag(toteffect_juv$chmjuvse_scaled^2)
  
  #This includes reallocation, so this is what I want: 
  (changebillionsse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                                         x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                                         x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                                         x31 + x32 + x33 + x34 + x35 + x36) + chjuvsoutside) / 
                                     1000, mycoefs, mybigvcov, ses=T))
  #5.371539, compared to 5.147492
  
  #Now get SE on total percentage change
  (totperse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                                 x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                                 x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                                 x31 + x32 + x33 + x34 + x35 + x36) + chjuvsoutside) / 
                             juv0, mycoefs, mybigvcov, ses=T))
  #0.05683242, compared to 0.0548165
  
  
}