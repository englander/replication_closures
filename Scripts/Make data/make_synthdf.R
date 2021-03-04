rm(list=ls())
setwd("C:/Users/englander/Documents/replication_closures/")
library(dplyr); library(readxl); library(ggplot2)
library(rworldmap); library(sf); library(lwgeom)
library(rgdal); library(geosphere); library(sp)
library(purrr); library(lubridate); library(glmnet)
library(lfe); library(Formula); library(smoothr)
library(parallel); library(Synth); library(tidyr)
library(xtable); library(cowplot); library(latex2exp)

#Peru time
Sys.setenv(TZ='America/Lima')

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

#Created in make_actualclosure_regressioncontrol.R
load("Output/TempData/actualclosure_regressioncontrol.Rdata")

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

fullbe <- rename(fullbe, calatime = FechaInicioCala)

#Make sf
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

fullbe <- st_sf(geometry = besf, fullbe)

rm(besf)

#Actual closures that begin at six am
sixids <- regdf$id[regdf$closuretype=="actual" & regdf$bin=="active_in"][grep("06:00:00",regdf$start[regdf$closuretype=="actual" & regdf$bin=="active_in"])]

regdf <- mutate(regdf, six = if_else(id %in% sixids, 1, 0))

rm(sixids)

#Need to calculate juvenile catch in preperiod for both actual and potential closures
predf <- rbind(
  filter(regdf, tvar==0) %>% 
    mutate(end = start - 24*3600-1) %>% 
    mutate(start = if_else(six==0,start - 9*3600, start - 12*3600), 
           tvar = -2, bin = paste0("lead2_",bdist)),
  filter(regdf, tvar==0) %>% 
    mutate(end = start - 24*3600-1) %>% 
    mutate(start = start - 48*3600, tvar = -3, bin = paste0("lead3_",bdist)),
  filter(regdf, tvar==0) %>% 
    mutate(end = start - 48*3600-1) %>% 
    mutate(start = start - 72*3600, tvar = -4, bin = paste0("lead4_",bdist)),
  filter(regdf, tvar==0) %>% 
    mutate(end = start - 72*3600-1) %>% 
    mutate(start = start - 96*3600, tvar=-5, bin=paste0("lead5_",bdist)),
  filter(regdf, tvar==0) %>% 
    mutate(end = start - 96*3600-1) %>% 
    mutate(start = start - 120*3600, tvar=-6, bin=paste0("lead6_",bdist)),
  filter(regdf, tvar==0) %>% 
    mutate(end = start - 120*3600-1) %>% 
    mutate(start = start - 144*3600, tvar=-7, bin=paste0("lead7_",bdist)),
  filter(regdf, tvar==0) %>% 
    mutate(end = start - 144*3600-1) %>% 
    mutate(start = start - 168*3600, tvar=-8, bin=paste0("lead8_",bdist)),
  filter(regdf, tvar==0) %>% 
    mutate(end = start - 168*3600-1) %>% 
    mutate(start = start - 192*3600, tvar=-9, bin=paste0("lead9_",bdist))
)

#Make bin format like that in regdf
predf$bin <- gsub("_0","_in",predf$bin)

#Drop outcome variables and recalculate below
#Drop outcome variables and treatment fraction; going to redefine them for pre-bins
predf <- dplyr::select(predf, -tons, -numindivids, -numjuv, -nobs, -sdtons, -numadults) %>% 
  mutate(treatfrac = as.numeric(NA))

st_crs(predf) <- "+proj=longlat +datum=WGS84"
st_crs(fullbe) <- "+proj=longlat +datum=WGS84"

#Calculate preperiod juvenile catch
outcomesFun <- function(rdrow){
  
  row <- predf[rdrow,]
  
  #Filter BE observations to within same time period
  mybesf <- filter(fullbe, row$start<=calatime & calatime<row$end)
  
  if(nrow(mybesf)>0){
    
    #BE observations that are spatially inside element of predf
    inter <- st_intersects(mybesf, row)
    
    inter <- as.numeric(as.character(inter))  
    
    if(length(which(!is.na(inter))>0)){
      tons <- sum(mybesf$betons[which(!is.na(inter))],na.rm=T)
      sdtons <- sd(mybesf$betons[which(!is.na(inter))],na.rm=T)
      numindivids <- sum(mybesf$numindivids[which(!is.na(inter))],na.rm=T)
      numjuv <- sum(mybesf$numjuv[which(!is.na(inter))],na.rm=T)
      nobs <- nrow(mybesf[which(!is.na(inter)),])
    } else{
      tons <- 0; numindivids <- 0; numjuv <- 0; nobs <- 0; sdtons = as.numeric(NA)
    }
  } else{
    tons <- 0; numindivids <- 0; numjuv <- 0; nobs <- 0; sdtons = as.numeric(NA)
  }
  
  #Output df with bin and rid so I can join onto predf
  out <- data.frame(bin=row$bin, id=row$id,tons=tons,
                    #Round number of individuals to integer
                    numindivids=round(numindivids), numjuv=round(numjuv),
                    nobs = nobs, sdtons = sdtons)
  
  return(out)
}


#Apply over all bins
(myCores <- detectCores())

cl <- makeCluster(12)

clusterExport(cl, "predf")
clusterExport(cl, "fullbe")
clusterExport(cl, "outcomesFun")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lubridate))

#Apply over rows of predf
myoutcomes <- parLapply(cl = cl,
                        1:nrow(predf),
                        function(x){
                          
                          outcomesFun(x)
                          
                        })

stopCluster(cl)
rm(cl, myCores)

myoutcomes <- bind_rows(myoutcomes)

#Join onto predf
predf <- left_join(predf, myoutcomes, by = c("bin",'id'))

predf <- mutate(predf, numadults = numindivids - numjuv)

##For each actual closure active_in, calculate weights for potential closure active_in that 
#have treatfrac=0 and non-missing length distribution
#using preperiod juvenile catch and control variables from 9-24 hours before beginning

#Possible control units
possiblecontrols <- filter(regdf, closuretype=="potential" & bin=="active_in" & treatfrac==0) %>% 
  as.data.frame() %>% dplyr::select(-geometry) %>% filter(!is.na(prop12hat)) %>% 
  mutate(id = as.integer(id)) %>%
  arrange(id)

controlvars <- c("clustnobs","clusttons","clusttonsperset","clusttonsperarea","clustarea_km2","kmtocoast",
                 grep("prop",names(regdf),value=T))

#Actual or possible control units
preoutcomes <- filter(predf, bdist==0 & (closuretype=="actual" | id %in% possiblecontrols$id)) %>% 
  as.data.frame() %>% dplyr::select(-geometry) %>% 
  #asinh(nummjuv)
  mutate(asinhnummjuv = asinh(numjuv/10^6))

#Possible control units as matrix where nrow=number of preperiods and ncol=number of possible control units
preoutcontrols <- filter(preoutcomes, closuretype=='potential') %>% 
  dplyr::select(id, tvar, asinhnummjuv) %>% 
  mutate(id = as.integer(id)) %>%
  arrange(id, desc(tvar)) %>% 
  spread(key = id, value = asinhnummjuv) %>% 
  dplyr::select(-tvar) %>% as.matrix()
  
actualactive_in <- filter(regdf, closuretype=='actual' & bin=="active_in") %>% 
  as.data.frame() %>% dplyr::select(-geometry)

#Synthetic control for given actual closure
synthFun <- function(rowind){
  
  synthobj <- synth(
    X1 = dplyr::select(actualactive_in[rowind,], all_of(controlvars)) %>% t(),
    X0 = dplyr::select(possiblecontrols, all_of(controlvars)) %>% t(),
    Z1 = filter(preoutcomes, bdist==0 & id==actualactive_in$id[rowind]) %>% 
      arrange(tvar) %>% dplyr::select(asinhnummjuv) %>% as.matrix(),
    Z0 = preoutcontrols
  )
  
  synthobj <- c(synthobj, actualactive_in$id[rowind])
  
  names(synthobj)[length(synthobj)] <- "id"
  
  return(synthobj)
}

#Apply over actual closures
(myCores <- detectCores())

cl <- makeCluster(10)

clusterExport(cl, "actualactive_in")
clusterExport(cl, "possiblecontrols")
clusterExport(cl, "controlvars")
clusterExport(cl, "preoutcomes")
clusterExport(cl, "preoutcontrols")
clusterExport(cl, "synthFun")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(Synth))


#Apply over rows of predf
bigsynthlist <- parLapply(cl = cl,
                        1:nrow(actualactive_in),
                        function(x){

                          try(synthFun(x))

                        })

#Save
save(bigsynthlist, file = 'Output/TempData/bigsynthlist.Rdata')

stopCluster(cl)
rm(cl, myCores)

#Some failed
whichsuccess <- sapply(1:410, function(x){
  ifelse(class(bigsynthlist[[x]])=="list",1,0)
})

#Number of failures
which(whichsuccess==0) %>% length()

#Add success indicator to actualactive_in
actualactive_in <- mutate(actualactive_in, synth = whichsuccess)

#So proceed with synthetic control on actual closures synth() worked for

#Create treated observations and control observations objects for pre and post period
treatdf <- bind_rows(
  filter(regdf, closuretype=='actual') %>% as.data.frame() %>% dplyr::select(-geometry),
  filter(predf, closuretype=='actual') %>% as.data.frame() %>% dplyr::select(-geometry)
) %>% 
  #Actual closures that synthetic control worked on
  filter(id %in% actualactive_in$id[actualactive_in$synth==1])

controldf <- bind_rows(
  filter(regdf, id %in% possiblecontrols$id)%>% as.data.frame() %>% dplyr::select(-geometry),
  filter(predf, id %in% possiblecontrols$id)%>% as.data.frame() %>% dplyr::select(-geometry)
) %>% 
  mutate(id = as.character(id))

#Second function. Given actual closure, return treated and weighted control asinhnummjuv
#Also calculate tons caught for use when calculating reallocation in tons caught
synthOuts <- function(synthobj){
  
  actualout <- filter(treatdf, id==synthobj$id) %>% 
    mutate(nummjuv = numjuv/10^6, asinhnummjuv = asinh(numjuv/10^6),
           asinhtons = asinh(tons)) %>% 
    dplyr::select(bdist, tvar, nummjuv, asinhnummjuv, asinhtons, id) %>% 
    rename(treatid = id)
  
  #Join weights onto controldf
  controlout <- 
    left_join(
      controldf %>% mutate(id = as.integer(id), asinhtons = asinh(tons)),
      data.frame(id = possiblecontrols$id, synthobj$solution.w),
      by = 'id') %>% 
    mutate(weightnummjuv = (numjuv/10^6)*w.weight,
           weightasinhnummjuv = asinh(numjuv/10^6)*w.weight,
           weightasinhtons = asinhtons*w.weight) %>% 
    group_by(bdist, tvar) %>%
    summarise(nummjuv = sum(weightnummjuv),
              asinhnummjuv = sum(weightasinhnummjuv),
              asinhtons = sum(weightasinhtons)) %>% ungroup() %>% 
    mutate(treatid = synthobj$id)
  
  #Join controlout onto actualout so can substract asinhnummjuv
  out <- left_join(
    rename(actualout, asinhnummjuv_treat = asinhnummjuv, nummjuv_treat = nummjuv,
           asinhtons_treat = asinhtons),
    rename(controlout, asinhnummjuv_control = asinhnummjuv, nummjuv_control = nummjuv,
           asinhtons_control = asinhtons),
    by =c("bdist","tvar","treatid")
  )
  
  
  return(out)
}

#Apply synthOuts over elements of bigsynth (each actual closure id)
synthdf <- map_df(1:length(bigsynthlist), function(x){
  if(whichsuccess[x]==1){
    synthOuts(bigsynthlist[[x]])
  }
})

save(synthdf, file="Output/Data/synthdf.Rdata")


sessionInfo()

