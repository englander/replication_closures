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
rddf <- dplyr::select(rddf, -tons, -numindivids, -numjuv, -nobs, -sdtons, -numadults)

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Load vessel information
owndf <- read_csv("Data/owndf.csv")

#Medium fleets are not top7 firm and have more than one vessel
owndf$fleettype[owndf$Armador %not in% c("TECNOLOGICA DE ALIMENTOS S.A.","PESQUERA DIAMANTE S.A.","CORPORACION PESQUERA INCA S.A.C.",
                                         "PESQUERA EXALMAR S.A.A.", "CFG INVESTMENT S.A.C.", "AUSTRAL GROUP S.A.A", "PESQUERA HAYDUK S.A.") & 
                  owndf$numowned > 1] <- "medium" 

#Join length of vessel onto fullbe
fullbe <- left_join(fullbe, 
                    dplyr::select(owndf, Matricula, Temporada, eslora, fleettype), 
                    by = c("Matricula","Temporada"))

#Create indicator for being above median length
fullbe <- mutate(fullbe, abovemedian = if_else(eslora > median(eslora),1,0))

#What % of medium vessels are above median length?
filter(fullbe, fleettype=="medium" & abovemedian==1) %>% 
  distinct(Matricula) %>% nrow() / 
  filter(fullbe, fleettype=="medium") %>% 
  distinct(Matricula) %>% nrow()

#Filter to medium firms only
fullbe <- filter(fullbe, fleettype=="medium")

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
    } else{
      tons <- 0; numindivids <- 0; numjuv <- 0; nobs <- 0; sdtons = as.numeric(NA)
    }
  } else{
    tons <- 0; numindivids <- 0; numjuv <- 0; nobs <- 0; sdtons = as.numeric(NA)
  }
  
  out <- data.frame(abovemedian = lengthtype, tons=tons,
                    #Round number of individuals to integer
                    numindivids=round(numindivids), numjuv=round(numjuv),
                    nobs = nobs, sdtons = sdtons)
  
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
                      nobs = 0, sdtons = as.numeric(NA))
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

juvcatch$N
jvdf <- summary(juvcatch)$df %>% unique()

jvtab <- summary(juvcatch)[["coefficients"]]

jvtab <- mutate(as.data.frame(jvtab), bin = rownames(jvtab))

jvtab <- jvtab[grep("treatfrac",jvtab$bin),]

jvtab$bin <- gsub("_treatfrac","",jvtab$bin)

jvtab$abovemedian <- 0
jvtab$abovemedian[grep("abovemedian1",jvtab$bin)] <- 1

jvtab$bin <- gsub("abovemedian1","",jvtab$bin)
jvtab$bin <- gsub("abovemedian0","",jvtab$bin)
jvtab$bin <- gsub(":","",jvtab$bin)

#Account for reallocation in tons
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
                               dplyr::select(Estimate, `Cluster s.e.`, bin),
                             by = 'bin') %>% 
    mutate(juv0 = juv1/(exp(Estimate))) %>% 
    mutate(chmjuv = juv1 - juv0) %>%
    arrange(tvar, bdist) %>% ungroup()
  
  #Scale number of juv in fullbe by number of juv in heterodf
  scaleconstant <- (sum(fullbe$numjuv,na.rm=T) / sum(heterodf$numjuv, na.rm=T))
  
  changejuv <- sum(toteffect_juv$chmjuv) * scaleconstant
  
  #Now can calculate change in juvenile catch due to policy, accounting for reallocation
  chmjuvsstart <- changejuv + chjuvsoutside
  
  #How many juveniles are caught during my sample period in total?
  juv1 <- sum(fullbe$numjuv[fullbe$abovemedian==lengthtype], na.rm=T) / 10^6
  
  #chmjuvsstart = juv(1) - juv(0)
  juv0 <- juv1 - chmjuvsstart
  
  #Change in juvenile catch as percentage
  totper <- chmjuvsstart / juv0
  
  out <- list(data.frame(abovemedian=lengthtype, chmjuvsstart = chmjuvsstart, totper = totper),
              list(toteffect_juv, scaleconstant, chjuvsoutside, juv0)
  ) 
  
  return(out)
  
}


effectLength(0)[[1]]
# abovemedian chmjuvsstart    totper
# 1           0     4033.461 0.4454436

effectLength(1)[[1]]
# abovemedian chmjuvsstart    totper
# 1           1     3730.975 0.2448035


##From below, standard error on difference in percent change is 0.05732406
#So the p-value on the difference is 0.0004653779
(1 - pt((0.4454436 - 0.2448035) / 0.05732406, df = jvdf))*2



#Calculate standard error on total percent effect for below median area
secondlist <- effectLength(0)[[2]]

#Get SE for chmjuv for each bin
toteffect_juv <- mutate(secondlist[[1]], chmjuvse = as.numeric(NA))

#Get other values I need
scaleconstant <- secondlist[[2]]; chjuvsoutside <- secondlist[[3]]; juv0 <- secondlist[[4]]

#Calculate delta method standard errors for change in millions of juveniles caught
for(i in 1:nrow(toteffect_juv)){
  
  #Filter toteffect_juv to mybin
  mydf <- toteffect_juv[i,]
  
  myjuv1 <- mydf$juv1
  mycoef <- mydf$Estimate
  myvcov <- mydf$`Cluster s.e.`^2
  
  #Delta se for chmjuvse. Random variable being transformed is myjuvcoef
  chmjuv_delta <- deltamethod( ~ myjuv1 - (myjuv1/(exp(x1))),
                               mycoef,
                               myvcov,
                               ses=T
  )
  
  #Plug this value into toteffect_juv
  toteffect_juv$chmjuvse[i] <- chmjuv_delta
  
}

#Now calculate standard error on total percent change
mycoefs <- toteffect_juv$chmjuv
mybigvcov <- diag(toteffect_juv$chmjuvse^2)

(totperse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                              x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                              x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                              x31 + x32 + x33 + x34 + x35 + x36)*scaleconstant + chjuvsoutside) / 
                          juv0, mycoefs, mybigvcov, ses=T))

#Create a larger coefdf and vcovdf so I can calculate delta se on difference in total percentage change
togcoefs <- mycoefs; togvcov <- mybigvcov 
chjuvsoutside_0 <- chjuvsoutside; juv0_0 <- juv0


#Calculate standard error on total percent effect for above median Length
secondlist <- effectLength(1)[[2]]

#Get SE for chmjuv for each bin
toteffect_juv <- mutate(secondlist[[1]], chmjuvse = as.numeric(NA))

#Get other values I need
scaleconstant <- secondlist[[2]]; chjuvsoutside <- secondlist[[3]]; juv0 <- secondlist[[4]]

#Calculate delta method standard errors for change in millions of juveniles caught
for(i in 1:nrow(toteffect_juv)){
  
  #Filter toteffect_juv to mybin
  mydf <- toteffect_juv[i,]
  
  myjuv1 <- mydf$juv1
  mycoef <- mydf$Estimate
  myvcov <- mydf$`Cluster s.e.`^2
  
  #Delta se for chmjuvse. Random variable being transformed is myjuvcoef
  chmjuv_delta <- deltamethod( ~ myjuv1 - (myjuv1/(exp(x1))),
                               mycoef,
                               myvcov,
                               ses=T
  )
  
  #Plug this value into toteffect_juv
  toteffect_juv$chmjuvse[i] <- chmjuv_delta
  
}

#Now calculate standard error on total percent change
mycoefs <- toteffect_juv$chmjuv
mybigvcov <- diag(toteffect_juv$chmjuvse^2)

(totperse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                              x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                              x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                              x31 + x32 + x33 + x34 + x35 + x36)*scaleconstant + chjuvsoutside) / 
                          juv0, mycoefs, mybigvcov, ses=T))

#Add these coefs and vcov to big so can calculate se on difference in total percent change
togcoefs <- c(togcoefs, mycoefs)
togvcov <- diag(c(diag(togvcov), diag(mybigvcov)))

(difftotperse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                                  x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                                  x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                                  x31 + x32 + x33 + x34 + x35 + x36)*scaleconstant + chjuvsoutside_0) / 
                              juv0_0 - ((x37+x38+x39+x40+x41+x42+x43+x44+x45+x46+
                                           x47+x48+x49+x50+x51+x52+x53+x54+x55+x56+x57+
                                           x58+x59+x60+x61+x62+x63+x64+x65+x66+x67+
                                           x68+x69+x70+x71+x72)*scaleconstant + chjuvsoutside) / 
                              juv0, togcoefs, togvcov, ses=T))
#0.05732406