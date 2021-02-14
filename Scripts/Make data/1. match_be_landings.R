#Match PRODUCE Bitacora Electronica to PRODUCE landings data
#This is the first script cleaning BE data. 
#The next one is impute_size_be.R

rm(list=ls())
setwd("C:/Users/englander/Documents/replication_closures/")

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

#Peru time
Sys.setenv(TZ='America/Lima')

library(dplyr); library(readr); library(purrr); library(readxl)
library(lubridate); library(ggplot2); library(lfe); library(car)
library(sf); library(parallel)

#Full 2017 to 2019 landings data from Mariano
land <- read_excel("Data/landings_2017to2019.xlsx")

#Full BE from Mariano
be <- read_excel("Data/BE_2017to2019_allvessels.xlsx")

be <- rename(be, Embarcacion=Embarcación, Matricula=Matrícula)

be$Zona[be$Zona=="Norte Centro"] <- "Norte-Centro"

be <- mutate(be, lon = LongitudGrados + LongitudMinutos/60 + LongitudSegundos/3600,
             lat = LatitudGrados + LatitudMinutos/60 + LatitudSegundos/3600)

#Should always be negative
be$lon <- -be$lon
be$lat <- -be$lat

#Convert both datasets to be in Peru time
land$FechaHoraInicioDescarga <- force_tz(land$FechaHoraInicioDescarga,tz="America/Lima")
land$FechaHoraTerminoDescarga <- force_tz(land$FechaHoraTerminoDescarga,tz="America/Lima")

be$FechaInicioFaena <- force_tz(be$FechaInicioFaena, tz="America/Lima")
be$FechaFinFaena <- force_tz(be$FechaFinFaena, tz="America/Lima")
be$FechaInicioCala <- force_tz(be$FechaInicioCala, tz="America/Lima")
be$FechaFinCala <- force_tz(be$FechaFinCala, tz="America/Lima")

#Match be to landings
land <- arrange(land, Matricula, FechaHoraInicioDescarga)

land <- mutate(land, landid=1:nrow(land)) %>% dplyr::select(-Textbox4)

#Consolidate FechaFinFaena across observations if they 
#have different values but seem to belong to same trip. e.g.:
bd <- distinct(be, Matricula, FechaInicioFaena, FechaFinFaena) %>% 
  arrange(Matricula, FechaInicioFaena, FechaFinFaena)

#If next be trip starts before previous one has ended, set FechaInicioFaena
#equal to the minimum of the two and set FechaFinFaena equal to the maximum of the two
bd <- mutate(bd, FechaInicioFaenahat = FechaInicioFaena, FechaFinFaenahat = FechaFinFaena)

#Join onto be
be <- left_join(be, bd, by = c("Matricula","FechaInicioFaena","FechaFinFaena"))

for(i in 2:nrow(bd)){
  if(bd$Matricula[i]==bd$Matricula[i-1]){
    if(bd$FechaInicioFaena[i] < bd$FechaFinFaena[i-1]){
      bd$FechaInicioFaenahat[i] <- min(bd$FechaInicioFaena[i],bd$FechaInicioFaena[i-1])
      bd$FechaInicioFaenahat[i-1] <- min(bd$FechaInicioFaena[i],bd$FechaInicioFaena[i-1])
      bd$FechaFinFaenahat[i] <- max(bd$FechaFinFaena[i],bd$FechaFinFaena[i-1])
      bd$FechaFinFaenahat[i-1] <- max(bd$FechaFinFaena[i],bd$FechaFinFaena[i-1])
    }
  }
}

rm(bd, i)

#183 sets take place outside of start and end of trip, so update start and end of trip for these rows as well
outsidecalas <- which(be$FechaInicioCala < be$FechaInicioFaenahat | be$FechaFinCala > be$FechaFinFaenahat)

for(i in outsidecalas){
  if(be$FechaInicioCala[i] < be$FechaInicioFaenahat[i]){
    be$FechaInicioFaenahat[i] <- be$FechaInicioCala[i]
  }
  if(be$FechaFinCala[i] > be$FechaFinFaenahat[i]){
    be$FechaFinFaenahat[i] <- be$FechaFinCala[i]
  }
}

rm(i, outsidecalas)

#Drop any BE trips that last more than two weeks
be <- mutate(be, tripdays = difftime(FechaFinFaenahat,FechaInicioFaenahat,units='days') %>% as.numeric())

be <- filter(be, tripdays < 14) %>% dplyr::select(-tripdays)

#Create trip id in be data
betripid <- distinct(be, Temporada, Zona, Matricula, Embarcacion, FechaInicioFaenahat, FechaFinFaenahat)
betripid <- mutate(betripid, betripid = 1:nrow(betripid))

#Join onto be
be <- left_join(be, betripid)

rm(betripid)


#Add date column
land <- land %>% mutate(enddate = date(FechaHoraInicioDescarga)) 
be <- mutate(be, enddate=date(FechaFinFaenahat))
#If FechaFinFaenahat is missing, use maximum FechaFinCala for that trip to get enddate
be <- bind_rows(filter(be, !is.na(enddate)),
                filter(be, is.na(enddate)) %>% dplyr::select(-enddate) %>%
                  left_join(filter(be, is.na(enddate)) %>% 
                              group_by(betripid) %>% 
                              summarise(maxfincala = max(FechaFinCala)) %>% 
                              mutate(enddate = date(maxfincala)) %>% ungroup() %>% dplyr::select(-maxfincala),
                            by='betripid'
                  )
)


#It's great that both datasets have Matricula. I can be sure this is the same vessel
matchvars <- distinct(land, Temporada, Zona, Matricula)

#Given Matricula, season and zone, filter be and land
#Assign landid from land to be
matchTrips <- function(matchrow){
  
  mybe <- filter(be, Temporada==matchvars$Temporada[matchrow] & 
                   Zona==matchvars$Zona[matchrow] & 
                   Matricula==matchvars$Matricula[matchrow]) %>% 
    mutate(landid = as.integer(NA))
  
  if(nrow(mybe)>0){
    
    myland <- filter(land, Temporada==matchvars$Temporada[matchrow] & 
                       Zona==matchvars$Zona[matchrow] & 
                       Matricula==matchvars$Matricula[matchrow])
    
    for(i in 1:nrow(myland)){
      
      landrow <- myland[i,]
      
      #BE that happened before this landing that does not have a landid yet
      berow <- filter(mybe, FechaFinFaenahat <= landrow$FechaHoraInicioDescarga & is.na(landid)) %>%
        #Take be trip that happened closest to landing event
        summarise(FechaFinFaenahat = max(FechaFinFaenahat))
      
      #If end of trip is not more than 24 hours before start of landing, 
      if(as.numeric(difftime(landrow$FechaHoraInicioDescarga,berow$FechaFinFaenahat,units='hours')) < 24){
        
        #give landid to mybe that has this FechaFinFaenahat
        mybe$landid[mybe$FechaFinFaenahat==berow$FechaFinFaenahat] <- landrow$landid
        
      }
    }
    return(mybe)
  }
}

#Apply over rows of matchvars
#Apply over all bins
(myCores <- detectCores())

cl <- makeCluster(18)

clusterExport(cl, "be")
clusterExport(cl, "land")
clusterExport(cl, "matchvars")
clusterExport(cl, "matchTrips")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(lubridate))

matchedbe <- parLapply(cl = cl,
                       1:nrow(matchvars),
                   function(x){
                     
                     matchTrips(x)
                     
                   })

stopCluster(cl)
rm(cl, myCores)

matchedbe <- bind_rows(matchedbe)

#Did not match 6.9% of BE rows to landing events. 
filter(matchedbe, is.na(landid)) %>% nrow() / nrow(be)

#Would be interesting to plot average pj in BE against pj in land
#Want to know 2 degree grid cell by two week group BE are in 
#If multiple groups within trip, use the group for the final set in the trip 
#First create 2 degree grid cell grid over Peru's EEZ
#Load peru EEZ. Downloaded from https://www.marineregions.org/downloads.php on July 2, 2018.
#Version 2 - 2012 (4.96 MB) [Known issues] (created from EEZ version 7)
eez <- st_read("Data/Intersect_IHO_EEZ_v2_2012/eez.shp") %>%
  filter(EEZ == "Peruvian Exclusive Economic Zone (disputed - Peruvian point of view)") %>% 
  st_transform(crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

grid2p <- st_make_grid(
  st_sfc(st_polygon(list(rbind(
    c(st_bbox(eez)["xmin"],st_bbox(eez)["ymax"]), 
    c(st_bbox(eez)["xmax"],st_bbox(eez)["ymax"]), 
    c(st_bbox(eez)["xmax"],st_bbox(eez)["ymin"]), 
    c(st_bbox(eez)["xmin"],st_bbox(eez)["ymin"]),
    c(st_bbox(eez)["xmin"],st_bbox(eez)["ymax"])
  ))),
  crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")),
  cellsize = 2, what = 'polygons') %>%
  st_sf()

#Keep only grid cells that intersect eez
inter <- st_intersects(eez, grid2p) %>% unlist()

grid2p <- grid2p[inter,]

rm(inter)

#Add a cellid column
grid2p <- mutate(grid2p, cellid = 1:nrow(grid2p))

#Save
save(grid2p, file = "Output/Data/grid2p.Rdata")


besf <- st_multipoint(cbind(matchedbe$lon, matchedbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

matchedbe <- st_sf(geometry = besf, matchedbe)

rm(besf)

#Which cell is each centroid inside of
insidecell <- st_intersects(matchedbe, grid2p)

insidecell <- as.numeric(as.character(insidecell))

#Add cellid column onto matched
matchedbe <- mutate(matchedbe, cellid_2p = insidecell)

rm(insidecell, grid2p)

#Create two-week of sample variable. 
twoweek <- as.data.frame(matchedbe) %>% dplyr::select(betripid, enddate)

twoweek <- mutate(twoweek, week = week(enddate), year = year(enddate))

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

#Join onto matchedbe
matchedbe <- left_join(matchedbe, dplyr::select(twoweek, betripid, twowk) %>% distinct(), by = c('betripid'))

#Now grouping/cluster variable
matchedbe <- mutate(matchedbe, twoweek_cellid_2p = paste0(twowk, "_", cellid_2p)) 

#Drop geometry column
matchedbe <- as.data.frame(matchedbe) %>% dplyr::select(-geometry)

#For trips in multiple groups, use group of final cala
#First need to increase set number (Correlativo) for betripids matched to a single landid
updateids <- distinct(matchedbe, betripid, landid) %>% 
  group_by(landid) %>% 
  summarise(count = n()) %>% filter(count > 1 & !is.na(landid)) %>% ungroup()

updateids <- filter(matchedbe, landid %in% updateids$landid) %>% 
  arrange(landid, FechaInicioCala)

for(i in 2:nrow(updateids)){
  if(updateids$landid[i]==updateids$landid[i-1]){
    updateids$Correlativo[i] <- updateids$Correlativo[i-1] + 1 
  }
}

#Bind updateids onto other landids in matchedbe
matchedbe <- filter(matchedbe, landid %not in% updateids$landid) %>% 
  bind_rows(updateids)

#Within landid, make sure Correlativo are sequential 
matchedbe <- arrange(matchedbe, landid, FechaInicioCala)

for(i in 2:nrow(matchedbe)){
  if(!is.na(matchedbe$landid[i]) & !is.na(matchedbe$landid[i-1])){
    if(matchedbe$landid[i]==matchedbe$landid[i-1]){
      matchedbe$Correlativo[i] <- matchedbe$Correlativo[i-1] + 1 
    }
  }
}

#Clean up
rm(twoweek, twowk, counter, i, updateids)

#Now record group of final cala in trip
finalcalgroup <- left_join(dplyr::select(matchedbe, landid, Correlativo, twoweek_cellid_2p),
                           matchedbe %>% group_by(landid) %>% 
                             summarise(maxCor = max(Correlativo)) %>% ungroup(),
                           by = c("landid")
) %>% 
  filter(Correlativo == maxCor & !is.na(landid)) %>% 
  dplyr::select(-maxCor, -Correlativo)

#Now join onto matchedbe (replace twoweek_cellid_2p for observations not missing landid)
matchedbe <- bind_rows(
  filter(matchedbe, is.na(landid)),
  filter(matchedbe, !is.na(landid)) %>% 
    dplyr::select(-twoweek_cellid_2p) %>% 
    left_join(finalcalgroup, by = c("landid"))
)

matchedbe <- rename(matchedbe, bepj = PorcentajeJuveniles, 
                    betons = `ToneladasPesca Declarada`, 
                    betripenddate = enddate, 
                    twoweek_cellid_2p_landidlevel = twoweek_cellid_2p)

#Join on landing information
matchedbe <- left_join(matchedbe, 
               dplyr::select(land, landid, FechaHoraInicioDescarga, FechaHoraTerminoDescarga,
                             TmDescargada, NumeroJuveniles, PorcentajeJuveniles, Moda) %>%
                 rename(landtons = TmDescargada, landNumeroJuveniles = NumeroJuveniles, 
                        landpj = PorcentajeJuveniles, landModa = Moda),
               by='landid')

save(matchedbe, file = "Output/Data/matched_be_landings_belevel.Rdata")

#How do tons differ between two datasets? (within matched observations)
filter(matchedbe, !is.na(landid)) %>% 
  group_by(landid, landtons) %>% 
  summarise(betons = sum(betons)) %>% 
  ungroup() %>% 
  summarise(sum(betons), sum(landtons))
#12348424  / 11372978 #8.6% more tons reported in electronic logbook data
