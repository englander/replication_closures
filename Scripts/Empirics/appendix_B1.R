rm(list=ls())
setwd("C:/Users/englander/Documents/replication_closures")
library(dplyr); library(ggplot2)
library(sf); library(msm)
library(purrr); library(lubridate)
library(lfe); library(Formula)
library(parallel); library(tidyr); library(cowplot)
library(viridis); library(latex2exp)

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

#Drop treatfrac
rddf <- dplyr::select(rddf, -treatfrac)

#Load closed areas created in make_closures_df.R
load("Output/Data/closed.Rdata")

#Only want to use closures declared by PRODUCE from BE data
closed <- filter(closed, days <= 5)

#Filter to closures during season
closed <- filter(closed, !is.na(season) & 
                   (season=="s1_2017" | season=="s2_2017" | 
                      season=="s1_2018" | season=="s2_2018" | season=="s1_2019" | 
                      season=="s2_2019")
)

#Assign the 1 day closure to be a 3 day closure for treatfrac category purposes
closed$days[closed$days==1] <- 3

#Crop closures to North-Central zone
nczone <- c(-85,-15.99999,-70,0)
names(nczone) <- names(st_bbox(closed))
closed <- st_crop(closed, nczone)

rm(nczone)

#Add indicator for closure being above and below median (within-season)
closed <- mutate(closed, area = st_area(geometry)) %>% 
  mutate(area = as.numeric(area))

#Average size of closure is 1,328 km^2
mean(closed$area)/10^6

closed <- left_join(
  closed, 
  as.data.frame(closed) %>% dplyr::select(-geometry) %>% 
  group_by(season) %>% 
  summarise(medarea = median(area)),
  by='season'
)

closed <- mutate(closed, abovemedarea = if_else(area > medarea, 1, 0))

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
  
  out <- do.call("rbind",out)
  
  return(out)
}

#Apply over all buffers
(myCores <- detectCores())

cl <- makeCluster(10)

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

treatType <- function(myrect, actclosed, heterotype){
  
  #Filter closures to given type
  if(heterotype >= 3){
    actclosed <- filter(actclosed, days==heterotype)
  } else{
    #If less than 3, then I want treatfrac for above or below median area 
    actclosed <- filter(actclosed, abovemedarea == heterotype)
  }
  
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
 
  myrect <- as.data.frame(myrect) %>% 
    dplyr::select(rid, bin, treatfrac)
  
  #Rename treatfrac according to heterotype
  names(myrect)[names(myrect)=="treatfrac"] <- paste0("treatfrac",heterotype)
  
  return(myrect) 
}

#Calculate rectangle overlap with actual closures
#Function of start and bin
treatVar <- function(mystart, mybin){
  
  #Rectangles
  myrect <- filter(rddf, start==mystart & bin==mybin) %>% 
    mutate(treatfrac = 0) #Space-time fraction overlapping with actual closed area
  
  #When these rectangles end
  myend <- unique(myrect$end)
  
  #Filter to active closures
  #Only care about closures declared by PRODUCE
  actclosed <- filter(closed, bin==mybin & 
                        start <= myend & 
                        mystart <= end)
  
  #Calculate treatment fraction for each of the treatfrac types
  #(above or below median (0 or 1) and 3-day or 5-day (3 or 5))
  myrect <- lapply(c(0,1,3,5), function(x){
    treatType(myrect, actclosed, x)
  })
  
  #Iteratively join each element
  out <- myrect[[1]]

  for(i in 2:length(myrect)){
    out <- left_join(out, myrect[[i]], by = c("rid","bin"))
  }
  
  return(out)
}

#Apply over all start values for a bin
treatBin <- function(mybin){
  mylist <- lapply(unique(rddf$start[rddf$bin==mybin]), function(x){
    
    treatVar(x, mybin)
  })
  
  out <- do.call("rbind",mylist)
  
  return(out)
}

#Apply over all bins
(myCores <- detectCores())

cl <- makeCluster(10)

clusterExport(cl, "rddf")
clusterExport(cl, "closed")
clusterExport(cl, "treatType")
clusterExport(cl, "treatVar")
clusterExport(cl, "treatBin")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lubridate))

rdf <- parLapply(cl = cl,
                 unique(rddf$bin),
                 function(x){
                   
                   treatBin(x)
                   
                 })

stopCluster(cl)
rm(cl, myCores)

rdf <- do.call("rbind",rdf)

#Join onto rddf
heterodf <- left_join(rddf, rdf, by = c("rid","bin")) %>% 
  as.data.frame() %>% dplyr::select(-geometry)

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
  interVars("treatfrac0"),
  interVars("treatfrac1"),
  interVars("treatfrac3"),
  interVars("treatfrac5")
)

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")


##Hetero by area
{
juvcatch <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(grep("_treatfrac0",names(heterodf),value=T),collapse="+"), " + ",
    paste0(grep("_treatfrac1",names(heterodf),value=T),collapse="+"),
    " + clustnobs + clusttons + clustarea_km2 + 
    kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(heterodf),value=T), collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = heterodf)

  juvcatch$N 
  jvdf <- summary(juvcatch)$df %>% unique()
  
jvtab <- summary(juvcatch)[["coefficients"]]

jvtab <- mutate(as.data.frame(jvtab), bin = rownames(jvtab))

jvtab <- jvtab[grep("treatfrac",jvtab$bin),]

#Area category
jvtab$areatype <- 0
jvtab$areatype[grep("_treatfrac1",jvtab$bin)] <- 1 

jvtab$bin <- gsub("_treatfrac0","",jvtab$bin)
jvtab$bin <- gsub("_treatfrac1","",jvtab$bin)

#Account for reallocation in tons caught
tonscaught <- felm(
  as.Formula(paste0(
    "asinhtons", "~ ", 
    "treatfrac0 + treatfrac1 ",
    " + clustnobs + clusttons + clustarea_km2 + 
    kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(heterodf),value=T),collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = heterodf)

tonscoef <- summary(tonscaught)[["coefficients"]]
tonscoef <- tonscoef[grep("treatfrac",rownames(tonscoef)),]
tonscoef <- as.data.frame(tonscoef) %>% 
  mutate(areatype = gsub("treatfrac","",rownames(tonscoef)))

#Don't reallocate tons from 2017 second season or 2019 second season because they were both 
#shut down well before TAC was reached. 

#Function of area type
effectArea <- function(myarea){
  
  #Tons caught in state of world I observe in seasons where TAC hit
  tons1 <- sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"])
  
  #Change in tons because of policy
  ctons <- tons1 - tons1/exp(tonscoef$Estimate[tonscoef$areatype==myarea]) 
  
  #Average pj outside of treatment window
  avgpjoutside <- filter(fullbe, lead_0==0 & lead_10==0 & lead_20==0 & lead_30==0 & lead_40==0 & lead_50==0 & 
                           active_0==0 & active_10==0 & active_20==0 & active_30==0 & active_40==0 & active_50==0 & 
                           lag1_0==0 & lag1_10==0 & lag1_20==0 & lag1_30==0 & lag1_40==0 & lag1_50==0 & 
                           lag2_0==0 & lag2_10==0 & lag2_20==0 & lag2_30==0 & lag2_40==0 & lag2_50==0 & 
                           lag3_0==0 & lag3_10==0 & lag3_20==0 & lag3_30==0 & lag3_40==0 & lag3_50==0 & 
                           lag4_0==0 & lag4_10==0 & lag4_20==0 & lag4_30==0 & lag4_40==0 & lag4_50==0 & 
                           !is.na(numindivids) & !is.na(bepjhat) & 
                           Temporada!="2017-II" & Temporada!="2019-II") %>%
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
                               Temporada!="2017-II" & Temporada!="2019-II") %>%
    #Weight by tons
    mutate(weightweighted = avgweightg*numindivids) %>% 
    summarise(avgweightg = sum(weightweighted)/sum(numindivids)) %>% as.numeric()
  
  
  #Decrease in individuals caught outside of treatment window in millions
  #(converting tons to g cancels out conversion to millions)
  chindividsoutside <- -ctons/avgweightoutside
  
  chjuvsoutside <- chindividsoutside*avgpjoutside 
  
  ##Calculate total effect, not accounting for reallocation of tons caught. (in millions)
  toteffect_juv <- heterodf %>%
    group_by(bin, tvar, bdist) %>% 
    summarise(juv1 = sum(nummjuv)) %>% ungroup()
  
  toteffect_juv <- left_join(toteffect_juv, 
                             filter(jvtab, areatype==myarea) %>% 
                               dplyr::select(Estimate, bin, `Cluster s.e.`),
                             by = 'bin') %>% 
    mutate(juv0 = juv1/(exp(Estimate))) %>% 
    mutate(chmjuv = juv1 - juv0) %>%
    arrange(tvar, bdist) %>% ungroup()
  
  #Scale number of juv in fullbe by number of juv in heterodf
  scaleconstant <- (sum(fullbe$numjuv,na.rm=T) / sum(heterodf$numjuv, na.rm=T))
  
  changejuv <-  sum(toteffect_juv$chmjuv) * scaleconstant
  
  #Now can calculate change in juvenile catch due to policy, accounting for reallocation
  chmjuvsstart <- changejuv + chjuvsoutside
  
  #How many juveniles are caught during my sample period in total?
  juv1 <- sum(fullbe$numjuv, na.rm=T) / 10^6
  
  #chmjuvsstart = juv(1) - juv(0)
  juv0 <- juv1 - chmjuvsstart
  
  #Change in juvenile catch as percentage
  totper <- chmjuvsstart / juv0
  
  out <- list(data.frame(areatype=myarea, chmjuvsstart = chmjuvsstart, totper = totper),
              list(toteffect_juv, scaleconstant, chjuvsoutside, juv0)
  )  
  return(out)
  
}

effectArea(0)[[1]]
# areatype chmjuvsstart   totper
# 1        0     35373.93 0.335745
#From below SE on totper is 0.06478523

effectArea(1)[[1]]
# areatype chmjuvsstart    totper
# 1        1      51687.9 0.5804659
#From below SE on totper is 0.1053633

##From below, standard error on difference in percent change is 0.1236873
#So the p-value on the difference is  0.04787547
(1 - pt((0.5804659 - 0.335745) / 0.1236873, df = jvdf))*2




#Calculate standard error on total percent effect for below median area
secondlist <- effectArea(0)[[2]]

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
#0.06478523

#Create a larger coefdf and vcovdf so I can calculate delta se on difference in total percentage change
togcoefs <- mycoefs; togvcov <- mybigvcov 
chjuvsoutside_0 <- chjuvsoutside; juv0_0 <- juv0


#Calculate standard error on total percent effect for above median area
secondlist <- effectArea(1)[[2]]

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
#0.1053633

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
#0.1236873

#Clean up
rm(totperse, mycoefs, mybigvcov, myjuv1, juv0, i, chjuvsoutside, 
   chmjuv_delta, mycoef, myvcov, scaleconstant, toteffect_juv, secondlist,
   tonscaught, tonscoef, jvtab, juvcatch, mydf, juv0_0, chjuvsoutside_0,
   togvcov, difftotperse, jvdf, togcoefs)
}

##Hetero by days
{
  
  #Don't include heterogeneity by 4 day closures because not enough of them
  #Also don't include 4 day lag because those estimates are unstable for 5 day closure
  juvcatch <- felm(
    as.Formula(paste0(
      "asinhnummjuv", "~ ", 
      paste0(grep("_treatfrac3",names(heterodf),value=T)[-grep("lag4",grep("_treatfrac3",names(heterodf),value=T))],collapse="+"), " + ",
      paste0(grep("_treatfrac5",names(heterodf),value=T)[-grep("lag4",grep("_treatfrac5",names(heterodf),value=T))],collapse="+"),
      " + clustnobs + clusttons + clustarea_km2 + 
    kmtocoast + clusttonsperset + clusttonsperarea + ", 
      paste0(grep("prop",names(heterodf),value=T), collapse="+"),
      "| bin + twowk:cellid_2p + startdate",
      " | 0 | twoweek_cellid_2p")),
    data = filter(heterodf, tvar!=4))
  
  juvcatch$N 
  jvdf <- summary(juvcatch)$df %>% unique()
  
  jvtab <- summary(juvcatch)[["coefficients"]]
  
  jvtab <- mutate(as.data.frame(jvtab), bin = rownames(jvtab))
  
  jvtab <- jvtab[grep("treatfrac",jvtab$bin),]
  
  #Day category
  jvtab$daytype <- 3
  jvtab$daytype[grep("_treatfrac5",jvtab$bin)] <- 5
  
  
  jvtab$bin <- gsub("_treatfrac3","",jvtab$bin)
  jvtab$bin <- gsub("_treatfrac5","",jvtab$bin)
  
  #Account for reallocation in tons caught
  tonscaught <- felm(
    as.Formula(paste0(
      "asinhtons", "~ ", 
      "treatfrac3 + treatfrac5 ",
      " + clustnobs + clusttons + clustarea_km2 + 
    kmtocoast + clusttonsperset + clusttonsperarea + ", 
      paste0(grep("prop",names(heterodf),value=T),collapse="+"),
      "| bin + twowk:cellid_2p + startdate",
      " | 0 | twoweek_cellid_2p")),
    data =  filter(heterodf, tvar!=4))
  
  tonscoef <- summary(tonscaught)[["coefficients"]]
  tonscoef <- tonscoef[grep("treatfrac",rownames(tonscoef)),]
  tonscoef <- as.data.frame(tonscoef) %>% 
    mutate(daytype = gsub("treatfrac","",rownames(tonscoef)))
  
  #Don't reallocate tons from 2017 second season or 2019 second season because they were both 
  #shut down well before TAC was reached. 
  
  #Function of area type
  effectDays <- function(myday){
    
    #Tons caught in state of world I observe in seasons where TAC hit
    tons1 <- sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"])
    
    #Change in tons because of policy
    ctons <- tons1 - tons1/exp(tonscoef$Estimate[tonscoef$daytype==myday]) 
    
    #Average pj outside of treatment window
    avgpjoutside <- filter(fullbe, lead_0==0 & lead_10==0 & lead_20==0 & lead_30==0 & lead_40==0 & lead_50==0 & 
                             active_0==0 & active_10==0 & active_20==0 & active_30==0 & active_40==0 & active_50==0 & 
                             lag1_0==0 & lag1_10==0 & lag1_20==0 & lag1_30==0 & lag1_40==0 & lag1_50==0 & 
                             lag2_0==0 & lag2_10==0 & lag2_20==0 & lag2_30==0 & lag2_40==0 & lag2_50==0 & 
                             lag3_0==0 & lag3_10==0 & lag3_20==0 & lag3_30==0 & lag3_40==0 & lag3_50==0 & 
                             lag4_0==0 & lag4_10==0 & lag4_20==0 & lag4_30==0 & lag4_40==0 & lag4_50==0 & 
                             !is.na(numindivids) & !is.na(bepjhat) & 
                             Temporada!="2017-II" & Temporada!="2019-II") %>%
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
                                 Temporada!="2017-II" & Temporada!="2019-II") %>%
      #Weight by tons
      mutate(weightweighted = avgweightg*numindivids) %>% 
      summarise(avgweightg = sum(weightweighted)/sum(numindivids)) %>% as.numeric()
    
    
    #Decrease in individuals caught outside of treatment window in millions
    #(converting tons to g cancels out conversion to millions)
    chindividsoutside <- -ctons/avgweightoutside
    
    chjuvsoutside <- chindividsoutside*avgpjoutside 
    
    ##Calculate total effect, not accounting for reallocation of tons caught. (in millions)
    toteffect_juv <- heterodf %>%
      group_by(bin, tvar, bdist) %>% 
      summarise(juv1 = sum(nummjuv)) %>% ungroup() %>%
      filter(tvar != 4)
    
    toteffect_juv <- left_join(toteffect_juv, 
                               filter(jvtab, daytype==myday) %>% 
                                 dplyr::select(Estimate, bin, `Cluster s.e.`),
                               by = 'bin') %>% 
      mutate(juv0 = juv1/(exp(Estimate))) %>% 
      mutate(chmjuv = juv1 - juv0) %>%
      arrange(tvar, bdist) %>% ungroup()
    
    #Scale number of juv in fullbe by number of juv in heterodf
    scaleconstant <- (sum(fullbe$numjuv,na.rm=T) / sum(heterodf$numjuv, na.rm=T))
    
    changejuv <-  sum(toteffect_juv$chmjuv) * scaleconstant
    
    #Now can calculate change in juvenile catch due to policy, accounting for reallocation
    chmjuvsstart <- changejuv + chjuvsoutside
    
    #How many juveniles are caught during my sample period in total?
    juv1 <- sum(fullbe$numjuv, na.rm=T) / 10^6
    
    #chmjuvsstart = juv(1) - juv(0)
    juv0 <- juv1 - chmjuvsstart
    
    #Change in juvenile catch as percentage
    totper <- chmjuvsstart / juv0
    
    out <- list(data.frame(daytype=myday, chmjuvsstart = chmjuvsstart, totper = totper),
                list(toteffect_juv, scaleconstant, chjuvsoutside, juv0)
    )  
    return(out)
    
  }
  
  effectDays(3)[[1]]
  # daytype chmjuvsstart    totper
  # 1       3     45684.71 0.4806452
  #From below SE on totper is 0.05481505
  
  effectDays(5)[[1]]
  # daytype chmjuvsstart    totper
  # 1       5     46657.05 0.4959486
  #From below SE on totper is 0.2244611
  
  ##From below, difference in percent change is 0.2310573
  #So the p-value on the difference is 0.9471936
  (1 - pt((0.4959486 - 0.4806452) / 0.2310573, df = jvdf))*2
  
  

  
  #Calculate standard error on total percent effect for 3 day closures
  secondlist <- effectDays(3)[[2]]
  
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
                                x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + 
                                x30)*scaleconstant + chjuvsoutside) / 
                            juv0, mycoefs, mybigvcov, ses=T))
  #0.05481505
  
  #Create a larger coefdf and vcovdf so I can calculate delta se on difference in total percentage change
  togcoefs <- mycoefs; togvcov <- mybigvcov 
  chjuvsoutside_0 <- chjuvsoutside; juv0_0 <- juv0
  
  
  #Calculate standard error on total percent effect for above median area
  secondlist <- effectDays(5)[[2]]
  
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
                                x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + 
                                x30)*scaleconstant + chjuvsoutside) / 
                            juv0, mycoefs, mybigvcov, ses=T))
  #0.2244611
  
  
  #Add these coefs and vcov to big so can calculate se on difference in total percent change
  togcoefs <- c(togcoefs, mycoefs)
  togvcov <- diag(c(diag(togvcov), diag(mybigvcov)))
  
  (difftotperse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                                    x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                                    x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30)*scaleconstant + chjuvsoutside_0) / 
                                juv0_0 - ((x31 + x32 + x33 + x34 + x35 + x36 + 
                                             x37+x38+x39+x40+x41+x42+x43+x44+x45+x46+
                                             x47+x48+x49+x50+x51+x52+x53+x54+x55+x56+x57+
                                             x58+x59+x60)*scaleconstant + chjuvsoutside) / 
                                juv0, togcoefs, togvcov, ses=T))
  
  #0.2310573
  
  
  #Clean up
  rm(totperse, mycoefs, mybigvcov, myjuv1, juv0, i, chjuvsoutside, 
     chmjuv_delta, mycoef, myvcov, scaleconstant, toteffect_juv, secondlist,
     tonscaught, tonscoef, jvtab, juvcatch, mydf)
}

sessionInfo()
