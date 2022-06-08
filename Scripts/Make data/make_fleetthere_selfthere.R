#Do vessels respond less to closures if a vessel in their fleet was already there before it is declared?
#Script used to be called "semidirect_infomechanism*.R"

rm(list=ls())

library(dplyr); library(ggplot2); library(lfe)
library(lubridate); library(sf); library(readxl)
library(purrr); library(parallel); 
library(xtable); library(lwgeom); library(tidyr)
library(Formula); library(msm); library(cowplot)

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

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


`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

options(lfe.threads=24)

#Peru time
Sys.setenv(TZ='America/Lima')

#Full electronic logbook data
be <- read_excel("Data/BE_2017to2019_allvessels.xlsx")

be <- rename(be, Embarcacion=Embarcación, Matricula=Matrícula)

be$Zona[be$Zona=="Norte Centro"] <- "Norte-Centro"

be <- mutate(be, lon = LongitudGrados + LongitudMinutos/60 + LongitudSegundos/3600,
             lat = LatitudGrados + LatitudMinutos/60 + LatitudSegundos/3600)

#Should always be negative
be$lon <- -be$lon
be$lat <- -be$lat

#Load ownership information
load("Data/owndf.Rdata")

#Keep columns of interest for this exercise
be <- dplyr::select(be, Matricula, Temporada, Zona, FechaInicioCala, lon, lat)

#Join ownership information onto be
be <- left_join(be, 
                dplyr::select(owndf, Matricula, Temporada, Armador, numowned, casco, eslora),
                by = c("Matricula","Temporada"))

#Make be sf
besf <- dplyr::select(be, lon, lat) %>% as.matrix() %>% 
  st_multipoint() %>% st_sfc(crs = "+proj=longlat +datum=WGS84 +no_defs") %>% 
  st_cast('POINT')

be <- st_sf(be, geometry = besf)

rm(besf)

#Load potential closures created in 4. make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

#Drop outcome variables; Going to redefine them
rddf <- dplyr::select(rddf, -tons, -numindivids, -numjuv, -nobs, -sdtons, -numadults)


#Create lead bin day before closure announcement would occur
#(24-48 hours before potential closure begins, because closure announcement occurs 9-24 hours before closure begins)
rddf <- rbind(
  rddf, 
  filter(rddf, bdist==0 & tvar==0) %>% 
    mutate(end = start-1-24*3600) %>% 
    mutate(start = start-24*3600*2, tvar=-3)
)


#Remake bin variable
rddf$bin[rddf$tvar==-3] <- gsub("active","lead3",rddf$bin[rddf$tvar==-3])
rddf$bin[rddf$tvar==-1] <- gsub("lead9hours","lead1",rddf$bin[rddf$tvar==-1])


#Make tvar match formatting in bin variable
rddf$tvar <- as.character(rddf$tvar)
rddf$tvar[rddf$tvar=="-3"] <- "lead3"
rddf$tvar[rddf$tvar=="-1"] <- "lead1"
rddf$tvar[rddf$tvar=="0"] <- "active"
rddf$tvar[grep("lag",rddf$bin)] <- paste0("lag",rddf$tvar[grep("lag",rddf$bin)])


#Make season in potential closures match format in be
rddf <- mutate(rddf, season = as.character(season)) %>% 
  mutate(Temporada = paste0(substr(season, 4,7),"-",
                            ifelse(substr(season,2,2)=="1","I","II"))) %>%
  dplyr::select(-season)

#Since outcome variable will be juvenile catch in each treatment bin, rather than whether vessel 
#fished in treatment bin of a potential closure,
#I only need to determine whether vessel fished in lead3_in for each potential closure

#Check whether given vessel, in given season, fished inside given potential closure's given treatment bin
fishThere <- function(myvesseason, mypotclseason, myrid){
  
  #Bin for this closure
  mypotcl <- filter(mypotclseason, rid == myrid)
  
  #Same time
  intime <- filter(myvesseason, FechaInicioCala>=mypotcl$start & FechaInicioCala<=mypotcl$end)
  
  if(nrow(intime)>0){
    #Inside?
    inbin <- st_intersects(intime, mypotcl)
    
    inbin <- as.numeric(as.character(inbin))
    
    inbin <- ifelse(1 %in% inbin,1,0)
  } else{
    inbin <- 0
  }
  
  out <- data.frame(bin = "lead3_in", fishedthere = inbin, rid = myrid)
  
  return(out)
}

#Apply fishThere over given season for given Matricula
applySeason <- function(vesbe, myseason){
  
  myvesseason <- filter(vesbe, Temporada==myseason)
  
  mypotclseason <- filter(rddf, Temporada==myseason & bin=="lead3_in")
  
  out <- map_df(unique(mypotclseason$rid), function(x){
    
    fishThere(myvesseason, mypotclseason, x)
    
  }) %>% 
    #Add season
    mutate(Temporada=myseason)
  
  return(out)
}

#Apply applySeason over seasons for given Matricula
applySeasons <- function(mymatricula){
  
  #BE for this matricula
  vesbe <- filter(be, Matricula==mymatricula)
  
  #Don't need to include rows for seasons where this vessel didn't fish at all
  out <- map_df(unique(vesbe$Temporada), function(x){
    applySeason(vesbe, x)
  }) %>% 
    #Add Matricula
    mutate(Matricula=mymatricula)
  
  return(out)
}

#Apply over vessels
(myCores <- detectCores())

cl <- makeCluster(24)

clusterExport(cl, "be")
clusterExport(cl, "rddf")
clusterExport(cl, "applySeason")
clusterExport(cl, "applySeasons")
clusterExport(cl, "fishThere")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(purrr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lubridate))


fishtherelist <- parLapply(cl = cl,
                      unique(be$Matricula),
                      function(x){

                        try(applySeasons(x))
                      })

ftdf <- bind_rows(fishtherelist)

stopCluster(cl)
rm(cl, myCores)

#Want to separate response to own information from response to fleet information
#Already know whether vessel has fished in lead3_in
selfalready <- dplyr::select(ftdf, -bin) %>% 
  rename(self_already3=fishedthere)

#Use selfalready in defining fleet activity. Add ownership column
fleetdf <- left_join(selfalready, dplyr::select(owndf, Matricula, Temporada, Armador), by = c("Temporada","Matricula"))

#Now want fleet_already3 column that I will join onto selfalready by rid, Temporada, and Matricula
#Given rid and same-fleet, same-season BE for given matricula, return whether vessel from fleet
#other than Matricula fished in lead3_in

#Check whether given vessel, in given season, fished inside given potential closure day before closure announcement would occur
fleetThere <- function(myfleetseason, myrid){
  
  #Fleet rows for this closure
  fleetrows <- filter(myfleetseason, rid == myrid)
  
  fleet3 <- max(fleetrows$self_already3)
  
  out <- data.frame(rid = myrid, fleet_already3 = fleet3)
  
  return(out)
}

#Apply fleetThere over closures in season
fleetSeason <- function(bigfleet, fleetbe, myseason){
  
  myfleetseason <- filter(fleetbe, Temporada==myseason)
  
  if(nrow(myfleetseason)>0){
    out <- map_df(unique(myfleetseason$rid), function(x){
      
      fleetThere(myfleetseason, x)
      
    }) %>% 
      #Add season
      mutate(Temporada=myseason)
  } else{
    out <- data.frame(rid = unique(bigfleet$rid[bigfleet$Temporada==myseason]),
                      fleet_already3 = as.numeric(NA), Temporada=myseason)
  }
  return(out)
}

#Apply fleetSeason over Matricula
fleetMatricula <- function(bigfleet, mymatricula){
  
  fleetbe <- filter(bigfleet, Matricula != mymatricula)
  
  if(nrow(fleetbe)>0){
    #Apply over seasons where mymatricula is active
    out <- map_df(unique(bigfleet$Temporada[bigfleet$Matricula==mymatricula]), function(x){
      fleetSeason(bigfleet, fleetbe, x)
    }) %>% 
      #Add Matricula
      mutate(Matricula=mymatricula)
  } else{
    out <- data.frame(rid = unique(bigfleet$rid), fleet_already3 = as.numeric(NA),
                      Matricula=mymatricula)
    
    out <- map_df(unique(bigfleet$Temporada), function(x){
      mutate(out, Temporada=x)
    })
  }
  
  return(out)
}

#Apply fleetMatricula over all vessels in a fleet
fleetArmador <- function(myarm){
  
  bigfleet <- filter(fleetdf, Armador==myarm)
  
  out <- map_df(unique(bigfleet$Matricula), function(x){
    fleetMatricula(bigfleet, x)
  }) %>% 
    mutate(Armador = myarm)
  
  return(out)
}


#Apply over owners
(myCores <- detectCores())

cl <- makeCluster(24)

clusterExport(cl, "fleetdf")
clusterExport(cl, "fleetThere")
clusterExport(cl, "fleetSeason")
clusterExport(cl, "fleetMatricula")
clusterExport(cl, "fleetArmador")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(purrr))


fleetlist <- parLapply(cl = cl,
                       unique(fleetdf$Armador),
                       function(x){
                         
                         try(fleetArmador(x))                       
                       })

fleettheredf <- bind_rows(fleetlist)

stopCluster(cl)
rm(cl, myCores)


fleetdf <- left_join(fleetdf, fleettheredf, by = c("rid","Temporada", "Matricula", "Armador"))


#Now construct outcome variable: realized juvenile catch in each rid-bin for each category

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

fullbe <- rename(fullbe, calatime = FechaInicioCala)

##Calculate tons and individuals caught in each element of rddf
#Make bedat an sf object
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

besf <- st_sf(geometry = besf, fullbe)

#Only keep variables I need
besf <- dplyr::select(besf, Matricula, Temporada, calatime, numindivids, 
                      numjuv, betons)

#Can now drop lead3_in bin from rddf
rddf <- filter(rddf, bin!="lead3_in")

#Make data frame for each type-value combination
tvdf <- data.frame(type = c("self_already3","fleet_already3"))

tvdf <- bind_rows(
  mutate(tvdf, value = 0),
  mutate(tvdf, value = 1)
)

tvdf$type <- as.character(tvdf$type)

#Only need a few columns of rddf in next function
userddf <- dplyr::select(rddf, rid, bin, start, end)

#Given row of rddf, 
#calculate outcomes for each type-value combination (4)
outFun <- function(rowind){
  
  row <- userddf[rowind,]
  
  #Filter BE observations to within same time period
  mybesf <- filter(besf, row$start<=calatime & calatime<row$end)
  
  if(nrow(mybesf)){
    
    #BE observations that are spatially inside element of rddf
    inter <- st_intersects(mybesf, row)
    
    inter <- as.numeric(as.character(inter))  
    
    if(length(which(!is.na(inter))>0)){
      
      #Redefine mybesf to rows that are spatially inside
      mybesf <- mybesf[which(!is.na(inter)),]
      
      #Add rid onto mybesf
      mybesf <- mutate(mybesf, rid = row$rid)
      
      #Join fleetdf onto mybesf
      joineddf <- left_join(as.data.frame(mybesf) %>% dplyr::select(-geometry),
                            fleetdf, 
                            by = c("Matricula","Temporada","rid"))
      
      #Apply function over each type-value combination
      out <- map_df(1:nrow(tvdf), function(x){
        singleVal(joineddf, x)
      }) %>% 
        mutate(rid = row$rid, bin = row$bin)
      
    } else{
      out <- data.frame(rid = row$rid, bin = row$bin)
      
      out <- map_df(1:nrow(tvdf), function(x){
        mutate(out, type = tvdf$type[x], value=tvdf$value[x], tons = 0,
               numindivids = 0, numjuv = 0, nobs = 0)
      })
    }
  } else{
    out <- data.frame(rid = row$rid, bin = row$bin)
    
    out <- map_df(1:nrow(tvdf), function(x){
      mutate(out, type = tvdf$type[x], value=tvdf$value[x], tons = 0,
             numindivids = 0, numjuv = 0, nobs = 0)
    })
  }
  
  return(out)
}


#Given sets inside a potential-closure bin, compute outcome for single type-value
singleVal <- function(joineddf, tvrowind){
  
  #Select type column from joinneddf combination
  filtdf <- dplyr::select(joineddf, numindivids, numjuv, betons, tvdf$type[tvrowind])
  
  #Rename last column so can refer to directly
  names(filtdf)[names(filtdf)==tvdf$type[tvrowind]] <- "typecol"
  
  #Filter to rows with given value
  filtdf <- filter(filtdf, typecol==tvdf$value[tvrowind])
  
  if(nrow(filtdf)>0){
    tons <- sum(filtdf$betons,na.rm=T)
    numindivids <- sum(filtdf$numindivids,na.rm=T)
    numjuv <- sum(filtdf$numjuv,na.rm=T)
    nobs <- nrow(filtdf)
  } else{
    tons <- 0; numindivids <- 0; numjuv <- 0; nobs <- 0
  }
  
  singleout <- data.frame(type = tvdf$type[tvrowind], value = tvdf$value[tvrowind],
                          tons=tons, numindivids = numindivids, numjuv = numjuv, nobs = nobs)
  
  return(singleout) 
}

#Apply over all bins
(myCores <- detectCores())

cl <- makeCluster(myCores - 10)

clusterExport(cl, "userddf")
clusterExport(cl, "besf")
clusterExport(cl, "tvdf")
clusterExport(cl, "fleetdf")
clusterExport(cl, "outFun")
clusterExport(cl, "singleVal")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lubridate))
clusterEvalQ(cl, library(purrr))


#Apply over rows of rddf
myoutcomes <- parLapply(cl = cl,
                        1:nrow(userddf),
                        function(x){

                          outFun(x)

                        })

myoutcomes <- bind_rows(myoutcomes)

stopCluster(cl)
rm(cl, myCores)

#Join rddf onto myoutcomes
heterodf <- left_join(myoutcomes,
                      as.data.frame(rddf) %>% dplyr::select(-geometry),
                      by = c("bin",'rid'))

save(heterodf, file = "Output/Data/fleetthere_selfthere.Rdata")

sessionInfo()
