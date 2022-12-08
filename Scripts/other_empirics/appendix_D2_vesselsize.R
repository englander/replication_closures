rm(list=ls())

library(dplyr); library(ggplot2)
library(sf); library(msm)
library(purrr); library(lubridate)
library(lfe); library(Formula)
library(parallel); library(tidyr); library(cowplot)
library(viridis); library(readr)

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

#Load potential closures created in 4. make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

#Drop outcome variables; that is what I am going to split into fleet categories
rddf <- dplyr::select(rddf, -tons, -numindivids, -numjuv, -nobs, -sdtons, -numadults, -uniquevessels)

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Load vessel information data
owndf <- read_csv("Data/owndf.csv")

#Join length of vessel onto fullbe
fullbe <- left_join(fullbe, 
                    dplyr::select(owndf, Matricula, Temporada, eslora), 
                    by = c("Matricula","Temporada"))

#Create indicator for being above median length
fullbe <- mutate(fullbe, abovemedian = if_else(eslora > median(eslora),1,0))

##Calculate tons and individuals caught in each element of rddf
#Make bedat an sf object
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

besf <- st_sf(geometry = besf, fullbe)

besf <- rename(besf, calatime = FechaInicioCala)

#Given row of rddf, BE data filtered to same time period, and length type
#calculate outcomes
outLength <- function(row, mybesf, lengthtype){
  
  #Filter mybesf to given fleet type
  fleetrows <- filter(mybesf, abovemedian == lengthtype)
  
  if(nrow(fleetrows)){
    
    #BE observations that are spatially inside element of rddf
    inter <- st_intersects(fleetrows, row)
    
    inter <- as.numeric(as.character(inter))  
    
    if(length(which(!is.na(inter))>0)){
      tons <- sum(fleetrows$betons[which(!is.na(inter))],na.rm=T)
      sdtons <- sd(fleetrows$betons[which(!is.na(inter))],na.rm=T)
      numindivids <- sum(fleetrows$numindivids[which(!is.na(inter))],na.rm=T)
      numjuv <- sum(fleetrows$numjuv[which(!is.na(inter))],na.rm=T)
      nobs <- nrow(fleetrows[which(!is.na(inter)),])
      uniquevessels <- length(unique(fleetrows$Matricula[which(!is.na(inter))]))
    } else{
      tons <- 0; numindivids <- 0; numjuv <- 0; nobs <- 0; sdtons = as.numeric(NA); uniquevessels <- 0
    }
  } else{
    tons <- 0; numindivids <- 0; numjuv <- 0; nobs <- 0; sdtons = as.numeric(NA); uniquevessels <- 0
  }
  
  out <- data.frame(abovemedian = lengthtype, tons=tons,
                    #Round number of individuals to integer
                    numindivids=round(numindivids), numjuv=round(numjuv),
                    nobs = nobs, sdtons = sdtons, uniquevessels = uniquevessels)
  
  return(out)
}

outcomesFun <- function(rdrow){
  
  row <- rddf[rdrow,]
  
  #Filter BE observations to within same time period
  mybesf <- filter(besf, row$start<=calatime & calatime<row$end)
  
  if(nrow(mybesf)>0){
    
    #Apply over length types
    out <- map_df(unique(besf$abovemedian), function(x){
      outLength(row, mybesf, x)
    }) %>% 
      #Add bin and rid
      mutate(bin = row$bin, rid = row$rid)
    
  } else{
    
    out <- data.frame(bin=row$bin, rid=row$rid,
                      abovemedian = unique(besf$abovemedian),
                      tons=0, numindivids=0, numjuv=0,
                      nobs = 0, sdtons = as.numeric(NA), uniquevessels = 0)
  }
  
  return(out)
}

#Apply over all bins
(myCores <- detectCores())

cl <- makeCluster(4)

clusterExport(cl, "rddf")
clusterExport(cl, "besf")
clusterExport(cl, "outcomesFun")
clusterExport(cl, "outLength")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lubridate))
clusterEvalQ(cl, library(purrr))


#Apply over rows of rddf
myoutcomes <- parLapply(cl = cl,
                        1:nrow(rddf),
                        function(x){
                          
                          outcomesFun(x)
                          
                        })

stopCluster(cl)
rm(cl, myCores)

myoutcomes <- bind_rows(myoutcomes)

#Join rddf onto myoutcomes
heterodf <- left_join(myoutcomes,
                      as.data.frame(rddf) %>% dplyr::select(-geometry),
                      by = c("bin",'rid'))


heterodf <- arrange(heterodf, tvar, bdist)

#Millions of juveniles
heterodf <- mutate(heterodf, nummjuv = numjuv/10^6) %>%
  #Cluster tons/set and tons/area
  mutate(clusttonsperset = clusttons/clustnobs, clusttonsperarea = clusttons/clustarea_km2) %>%
  #clust millions of juveniles, individuals, and adults
  mutate(clustmjuv = clustnumjuv/10^6, clustmindivids = clustnumindivids/10^6,
         clustmadults = clustnumadults/10^6)

heterodf <- mutate(heterodf, asinhnummjuv = asinh(nummjuv), 
                   asinhtons = asinh(tons))

heterodf$bin <- as.factor(heterodf$bin)
heterodf$bin <- relevel(heterodf$bin, ref="active_in")

heterodf$twoweek_cellid_2p <- as.factor(heterodf$twoweek_cellid_2p)
heterodf$twowk <- as.factor(heterodf$twowk)
heterodf$cellid_2p <- as.factor(heterodf$cellid_2p)

heterodf$abovemedian <- as.factor(heterodf$abovemedian)

#Drop potential closures that have NA for size distribution
heterodf <- filter(heterodf, !is.na(prop12hat))


#Given variable, interact it with bin indicators, giving
interVars <- function(var){
  
  mydf <- heterodf
  names(mydf)[names(mydf)==var] <- "myvar"
  
  #Want to manually interact var with bin indicator so I can look at each bin's coefficient
  #relative to 0 (rather than relative to omitted category)
  bininds <- model.matrix(~bin,data=heterodf) %>% as.data.frame()
  
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

heterodf <- bind_cols(
  heterodf,
  interVars("treatfrac")
)


#Realized juveniles caught
juvcatch <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(
      paste0(grep("_treatfrac",names(heterodf),value=T), ":abovemedian"),
      collapse="+"),
    " + clustnobs + clusttons + clustarea_km2 + 
    kmtocoast + clusttonsperset + clusttonsperarea + ", 
      paste0(grep("prop",names(heterodf),value=T), collapse="+"),
    "| bin:abovemedian + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = heterodf)

jvtab <- summary(juvcatch)[["coefficients"]]

jvtab <- mutate(as.data.frame(jvtab), bin = rownames(jvtab))

jvtab <- jvtab[grep("treatfrac",jvtab$bin),]

jvtab$bin <- gsub("_treatfrac","",jvtab$bin)

jvtab$abovemedian <- 0
jvtab$abovemedian[grep("abovemedian1",jvtab$bin)] <- 1

jvtab$bin <- gsub("abovemedian1","",jvtab$bin)
jvtab$bin <- gsub("abovemedian0","",jvtab$bin)
jvtab$bin <- gsub(":","",jvtab$bin)


#Account for reallocation in tons caught
tonscaught <- felm(
  as.Formula(paste0(
    "asinhtons", "~ ", 
    "treatfrac:abovemedian ",
    " + clustnobs + clusttons + clustarea_km2 + 
    kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(heterodf),value=T),collapse="+"),
    "| bin:abovemedian + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = heterodf)

tonscoef <- summary(tonscaught)[["coefficients"]]
tonscoef <- tonscoef[grep("treatfrac",rownames(tonscoef)),]
tonscoef <- as.data.frame(tonscoef) %>% 
  mutate(abovemedian = gsub("treatfrac:abovemedian","",rownames(tonscoef)))


#Don't reallocate tons from 2017 second season or 2019 second season because they were both 
#shut down well before TAC was reached. 

#Function of fleet type
effectLength <- function(lengthtype){
  
  #Tons caught in state of world I observe in seasons where TAC hit
  tons1 <- sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II" & 
                               fullbe$abovemedian==lengthtype])
  
  #Change in tons because of policy
  ctons <- tons1 - tons1/exp(tonscoef$Estimate[tonscoef$abovemedian==lengthtype]) 
  
  #Average pj outside of treatment window
  avgpjoutside <- filter(fullbe, lead_0==0 & lead_10==0 & lead_20==0 & lead_30==0 & lead_40==0 & lead_50==0 & 
                           active_0==0 & active_10==0 & active_20==0 & active_30==0 & active_40==0 & active_50==0 & 
                           lag1_0==0 & lag1_10==0 & lag1_20==0 & lag1_30==0 & lag1_40==0 & lag1_50==0 & 
                           lag2_0==0 & lag2_10==0 & lag2_20==0 & lag2_30==0 & lag2_40==0 & lag2_50==0 & 
                           lag3_0==0 & lag3_10==0 & lag3_20==0 & lag3_30==0 & lag3_40==0 & lag3_50==0 & 
                           lag4_0==0 & lag4_10==0 & lag4_20==0 & lag4_30==0 & lag4_40==0 & lag4_50==0 & 
                           !is.na(numindivids) & !is.na(bepjhat) & 
                           Temporada!="2017-II" & Temporada!="2019-II" & abovemedian==lengthtype) %>%
    #Weight by tons
    mutate(pjweighted = bepjhat*numindivids) %>%  
    summarise(perjuv = sum(pjweighted)/sum(numindivids)) %>% as.numeric() / 100
  
  
  #Avg weight of individual caught outside of treatment window
  avgweightoutside <- filter(fullbe, lead_0==0 & lead_10==0 & lead_20==0 & lead_30==0 & lead_40==0 & lead_50==0 & 
                               active_0==0 & active_10==0 & active_20==0 & active_30==0 & active_40==0 & active_50==0 & 
                               lag1_0==0 & lag1_10==0 & lag1_20==0 & lag1_30==0 & lag1_40==0 & lag1_50==0 & 
                               lag2_0==0 & lag2_10==0 & lag2_20==0 & lag2_30==0 & lag2_40==0 & lag2_50==0 & 
                               lag3_0==0 & lag3_10==0 & lag3_20==0 & lag3_30==0 & lag3_40==0 & lag3_50==0 & 
                               lag4_0==0 & lag4_10==0 & lag4_20==0 & lag4_30==0 & lag4_40==0 & lag4_50==0 & 
                               !is.na(numindivids) & !is.na(avgweightg) & 
                               Temporada!="2017-II" & Temporada!="2019-II" & abovemedian==lengthtype) %>%
    #Weight by tons
    mutate(weightweighted = avgweightg*numindivids) %>% 
    summarise(avgweightg = sum(weightweighted)/sum(numindivids)) %>% as.numeric()
  
  
  #Decrease in individuals caught outside of treatment window in millions
  #(converting tons to g cancels out conversion to millions)
  chindividsoutside <- -ctons/avgweightoutside
  
  chjuvsoutside <- chindividsoutside*avgpjoutside 
  
  ##Calculate total effect, not accounting for reallocation of tons caught. (in millions)
  toteffect_juv <- filter(heterodf, abovemedian==lengthtype) %>%
    group_by(bin, tvar, bdist) %>% 
    summarise(juv1 = sum(nummjuv)) %>% ungroup()
  
  toteffect_juv <- left_join(toteffect_juv, 
                             filter(jvtab, abovemedian==lengthtype) %>% 
                               dplyr::select(Estimate, bin),
                             by = 'bin') %>% 
    mutate(juv0 = juv1/(exp(Estimate))) %>% 
    mutate(chmjuv = juv1 - juv0) %>%
    arrange(tvar, bdist) %>% ungroup()
  
  #Scale number of juv in fullbe by number of juv in heterodf
  changejuv <- (sum(fullbe$numjuv,na.rm=T) / sum(heterodf$numjuv, na.rm=T)) * 
    sum(toteffect_juv$chmjuv)
  
  #Now can calculate change in juvenile catch due to policy, accounting for reallocation
  chmjuvsstart <- changejuv + chjuvsoutside
  
  #How many juveniles are caught during my sample period in total?
  juv1 <- sum(fullbe$numjuv[fullbe$abovemedian==lengthtype], na.rm=T) / 10^6
  
  #chmjuvsstart = juv(1) - juv(0)
  juv0 <- juv1 - chmjuvsstart
  
  #Change in juvenile catch as percentage
  totper <- chmjuvsstart / juv0
  
  out <- data.frame(abovemedian=lengthtype, chmjuvsstart = chmjuvsstart, totper = totper)
  
  return(out)
  
}


effectLength(1)
# abovemedian chmjuvsstart    totper
# 1           1     43537.04 0.5916097

effectLength(0)
# abovemedian chmjuvsstart    totper
# 1           0     4419.721 0.2303643

#So above median vessels account for 91% of effect
abovemedeffect <- effectLength(1) %>% dplyr::select(chmjuvsstart) %>% as.matrix() %>% as.numeric() / 
  (
    effectLength(1) %>% dplyr::select(chmjuvsstart) %>% as.matrix() %>% as.numeric() + 
      effectLength(0) %>% dplyr::select(chmjuvsstart) %>% as.matrix() %>% as.numeric()
  ); save(abovemedeffect, file = 'Output/TempData/appendix_D2_above_med_length_vessel_frac_effect.Rdata')

#What percent of juveniles are caught by above median vessels? 
above_med_length_vessel_juv_frac <- (sum(heterodf$numjuv[heterodf$abovemedian==1], na.rm=T) / 
    sum(heterodf$numjuv, na.rm=T)); save(above_med_length_vessel_juv_frac, file = 'Output/TempData/appendix_D2_above_med_length_vessel_juv_frac.Rdata')  #0.8340564

#Cannot separate vessel size from fleet size because large vessels are owned by large firms
#e.g. 96% of vessels owned by top 7 firms are above median length
filter(owndf, Temporada=="2019-I") %>% 
  filter(numowned > 11) %>% 
  filter(eslora > 27.7) %>% nrow() / 
  filter(owndf, Temporada=="2019-I") %>% 
  filter(numowned > 11) %>% nrow()