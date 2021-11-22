#What would happen if instead of closures, regulator highlighted areas with high percentage of adult catch
#1. Subtract change in tons, juveniles and adults caught due to closures policy
#2. Choose 410 potential closures with highest adult percentage (from sets that generate potential closure)
#3. Increase tons caught by 35% within treatment window of these 410 potential closures and calculate
#resulting change in juveniles and adults
#4. Reallocate tons in 4 of 6 seasons where TAC binding and calculate offseting change in juveniles and adults caught
#5. Effect of eliminating closures policy and replacing it with "adult policy" on juvenile catch 
# is 1 + 3 + 4

#Adult policy and no closures policy would decrease juvenile catch by 56% relative to status quo closures policy

rm(list=ls())
setwd("C:/Users/englander/Documents/replication_closures")
library(dplyr); library(ggplot2)
library(sf);library(lubridate); library(glmnet)
library(lfe); library(Formula); library(xtable)
library(purrr); library(cowplot); library(tidyr)
library(latex2exp)

options(scipen=999)
options(lfe.threads=24)

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

#Peru time
Sys.setenv(TZ='America/Lima')

#Load rddf created in 4. make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

rddf <- as.data.frame(rddf) %>% dplyr::select(-geometry)

rddf <- arrange(rddf, tvar, bdist)

#1. From make_figure8.R, I know that 
ctons <- 2722546
chmjuvsstart <- 46829.39
#And percent change in tons caught is 
pctons <- exp(0.29933719)-1

#2. 410 potential closures with highest adult percentage among sets generating potential closures
top410 <- filter(rddf, bin=="active_in" & !is.na(meanpj_weighted)) %>% #note meanpj_weighted was calculated among sets generating potential closures in 4. make_rddf.R
  dplyr::select(rid, meanpj_weighted) %>% 
  arrange(meanpj_weighted)

top410 <- top410[1:410,] %>% dplyr::select(rid) %>% as.matrix() %>% as.integer()
  
#3. Increase tons caught by 35% within treatment window of these potential closures
tons410 <- filter(rddf, rid %in% top410) %>% 
  summarise(tons = sum(tons)) %>% as.matrix() %>% as.numeric()

#Adjust for double-counting
#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

tons410 <- tons410 * (sum(fullbe$betons,na.rm=T) / sum(rddf$tons, na.rm=T))

#Make tons caught 35% larger and calculate difference
chtons3 <- tons410*pctons

#Average weight of individual caught within treatment window of 410 potcl
avgweightg3 <- filter(rddf, rid %in% top410) %>% 
  mutate(weighted = avgweightg * numindivids) %>% 
  summarise(avgweightg = sum(weighted, na.rm=T) / sum(numindivids, na.rm=T)) %>% as.matrix() %>% as.numeric()

#Change in number of individuals caught in millions
#(converting tons to g cancels out conversion to millions)
chindivids3 <- chtons3/avgweightg3

#Average percentage juvenile caught within treatment window of 410 potcl
avgpj3 <- filter(rddf, rid %in% top410) %>% 
  summarise(avgpj3 = sum(numjuv, na.rm=T) / sum(numindivids, na.rm=T)) %>% as.matrix() %>% as.numeric()

chjuv3 <- chindivids3 * avgpj3 #in millions



#4. In 4 of 6 seasons where TAC binding, reallocate chtons3 to outside treatment window of top410
#To do this, need to know avgweight of individual caught outside window and percentage juvenile
chtons4 <- chtons3 * 
  (sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"], na.rm=T) / sum(fullbe$betons, na.rm=T)) #Scale by fraction of tons caught in these four seasons

#tons caught outside treatment window
tonsout <- sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"], na.rm=T) - 
  (
    filter(rddf, rid %in% top410 & season != 's2_2017' & season != 's2_2019') %>%
  summarise(tons = sum(tons)) %>% as.matrix() %>% as.numeric() * 
  #Re-scale to avoid double-counting
    sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"],na.rm=T) / 
    sum(rddf$tons[rddf$season!='s2_2017' & rddf$season!="s2_2019"], na.rm=T)
  )
  
#Individuals caught outside treatment window
nindividsout <- (sum(fullbe$numindivids[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"], na.rm=T) - 
  (
    filter(rddf, rid %in% top410 & season != 's2_2017' & season != 's2_2019') %>%
      summarise(numindivids = sum(numindivids)) %>% as.matrix() %>% as.numeric() * 
      #Re-scale to avoid double-counting
      sum(fullbe$numindivids[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"],na.rm=T) / 
      sum(rddf$numindivids[rddf$season!='s2_2017' & rddf$season!="s2_2019"], na.rm=T)
  ) ) / 10^6 #in millions


#average weight of an individual outside treatment window
avgweightoutside <- tonsout / nindividsout

#So change in individuals is 
chindivids4 <- chtons4 / avgweightoutside

#Similarly, calculate percentage juvenile outside treatment window in 4 of 6 seasons
juvout <- (sum(fullbe$numjuv[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"], na.rm=T) - 
  (
    filter(rddf, rid %in% top410 & season != 's2_2017' & season != 's2_2019') %>%
      summarise(numjuv = sum(numjuv)) %>% as.matrix() %>% as.numeric() * 
      #Re-scale to avoid double-counting
      sum(fullbe$numjuv[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"],na.rm=T) / 
      sum(rddf$numjuv[rddf$season!='s2_2017' & rddf$season!="s2_2019"], na.rm=T)
  )) / 10^6 #in millions


pjout <- juvout / nindividsout
  
#Compensating decrease in juvenile catch
chjuv4 <- chindivids4*pjout


#5. Effect of eliminating closures policy and replacing it with "adult policy" on juvenile catch 
alteffect <- -chmjuvsstart  + chjuv3 - chjuv4

juv1 <- (sum(fullbe$numjuv, na.rm=T) / 10^6)

#counterfactual juvenile with adult policy and no closures policy
juv0 <- juv1 + alteffect #closures policy + (no closures policy + adult policy = alteffect)

#Adult policy and no closures policy would decrease juvenile catch by 56% relative to status quo closures policy
(juv0 - juv1) / juv1


