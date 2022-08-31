#For each potential closure-treatment bin, filter to sets during period. 
#Calculate distance to each of these sets and record five quantiles and mean
#Run main regression, except dependent variable is mean distance to sets

rm(list=ls())

library(dplyr); library(ggplot2); library(lfe)
library(lubridate); library(xtable); library(car)
library(purrr); library(readxl); library(readr)
library(sf); library(rworldmap)
library(collapse); library(furrr); library(Formula)

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

options(scipen=999)

#Peru time
Sys.setenv(TZ='America/Lima')

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

#Created in 3. correct_size_be.R
load("Output/Data/pbe_imp.Rdata")

#Created in make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

#Make sets spatial points
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

besf <- st_sf(geometry = besf, dplyr::select(fullbe, FechaInicioCala) %>% 
                rename(calatime = FechaInicioCala))

#Given potential closure-treatment bin, filter to sets during period. 
#Calculate distance to each of these sets and record five quantiles and mean
distFun <- function(myrowind){
  
  row <- rddf[myrowind,] 
  
  #Filter sets to those during this period
  mysets <- fsubset(besf, calatime >= row$start & calatime <= row$end)
  
  #Distance of each of these sets to potential closure-treatment bin
  mydist <- st_distance(mysets, row) %>% as.matrix() %>% as.numeric()
  
  #Convert to km
  mydist <- mydist / 1000
  
  #Make data frame for column binding
  mydist <- summary(mydist) %>% as.matrix() %>% t() %>% as.data.frame()
  
  names(mydist) <- c("dist_km_min", "dist_km_first_quartile", "dist_km_median", "dist_km_mean", 
                     "dist_km_third_quartile", "dist_km_max")
  
  #Bind to row
  out <- bind_cols(row, mydist) %>% 
    #Also record number of sets
    mutate(dist_km_nobs = nrow(mysets)) %>% 
    #Drop geometry column
    as.data.frame() %>% 
    dplyr::select(-geometry)
  
  
  return(out)
}


plan(multisession, workers = 6)


distlist <- future_map(
  1:nrow(rddf),
  function(x){
    
    try(distFun(x))
    
  })

#Regression data frame
regdf <- bind_rows(distlist)

regdf <- arrange(regdf, tvar, bdist)

#Create control variables
regdf <- regdf %>% #Cluster tons/set and tons/area
  mutate(clusttonsperset = clusttons/clustnobs, clusttonsperarea = clusttons/clustarea_km2) %>%
  #clust millions of juveniles, individuals, and adults
  mutate(clustmjuv = clustnumjuv/10^6, clustmindivids = clustnumindivids/10^6,
         clustmadults = clustnumadults/10^6)

regdf$bin <- as.factor(regdf$bin)
regdf$bin <- relevel(regdf$bin, ref="active_in")

regdf$twoweek_cellid_2p <- as.factor(regdf$twoweek_cellid_2p)
regdf$twowk <- as.factor(regdf$twowk)
regdf$cellid_2p <- as.factor(regdf$cellid_2p)

#Drop potential closures that have NA for size distribution
regdf <- filter(regdf, !is.na(prop12hat))

#Drop 1,134 observations when no sets during period of potential closure-treatment bin (so dist_km_mean is NA)
regdf <- filter(regdf, dist_km_nobs > 0)

#Remaining potential closure-treatment bins have 825 sets on average
mean(regdf$dist_km_nobs)

#Given dependent variable (quantile of distance), re-estimate Eq 1 with this as dependent variable
regFun <- function(mydepvar){
  
  distreg <- felm(
    as.Formula(paste0(
      mydepvar, "~ ", 
      #paste0(grep("_treatfrac",names(regdf),value=T),collapse="+"), #don't care about bin-specific effect; just want single average effect across bins
      "treatfrac ",
      " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
      paste0(grep("prop",names(regdf),value=T),collapse="+"),
      "| bin + twowk:cellid_2p + startdate",
      " | 0 | twoweek_cellid_2p")),
    data =  regdf
  )
  
  return(distreg)
  
}

#Apply over quantiles
reglist <- lapply(c("dist_km_min", "dist_km_first_quartile", "dist_km_median", "dist_km_mean", 
                    "dist_km_third_quartile", "dist_km_max"), function(x){
                      regFun(x)
                    })



#Format coefficient
formCoef <- function(reg, coef, dig){
  
  #Get coefficients from felm object
  mycoefs <- summary(reg)["coefficients"]$coefficients
  
  mycoef <- mycoefs[coef,"Estimate"]
  
  #Round
  roundcoef <- round(mycoef, dig) %>% as.character()
  
  #If rounded to integer, need to add "." to end
  if(length(grep("\\.",roundcoef))==0){
    roundcoef <- paste0(roundcoef, ".")
  }
  
  #Add an extra zero beyond the decimal point if needed to get same length
  #Do coef first
  roundcoef <- sapply(seq_len(length(roundcoef)), function(x){
    if(gsub(".*\\.","",roundcoef[x]) %>% nchar() < dig){
      #Needed length
      zerosneeded <- dig - gsub(".*\\.","",roundcoef[x]) %>% nchar()
      roundcoef[x] <- paste0(roundcoef[x],paste0(rep(0,zerosneeded),collapse=""))
    } else{
      roundcoef[x]
    }
  })
  
  #Add commas if necessary
  roundcoef <- prettyNum(roundcoef, ",")
  
  return(roundcoef)
}

#Given regression object, coefficient of interest, and number of digits to round to, 
#return formatted SE
formSE <- function(reg, coef, dig){
  
  #Get coefficients from felm object
  mycoefs <- summary(reg)["coefficients"]$coefficients
  
  #Get se
  se <- mycoefs[coef,"Cluster s.e."]
  
  #Round
  roundse <- round(se, dig) %>% as.character()
  
  #If rounded to integer, need to add "." to end
  if(length(grep("\\.",roundse))==0){
    roundse <- paste0(roundse, ".")
  }
  
  #Add zeros if necessary
  roundse <- sapply(seq_len(length(roundse)), function(x){
    if(gsub(".*\\.","",roundse[x]) %>% nchar() < dig){
      #Needed length (could need one extra zero or two)
      zerosneeded <- dig - gsub(".*\\.","",roundse[x]) %>% nchar()
      roundse[x] <- paste0(roundse[x],paste0(rep(0,zerosneeded),collapse=""))
    } else{
      roundse[x]
    }
  })
  
  #Add commas if necessary
  roundse <- prettyNum(roundse, ",")
  
  #Add parentheses 
  roundse <- paste0("(",roundse,")")
  
  return(roundse)
}

#Given number and number of digits, round and add zero after decimal if necessary
formNum <- function(num, dig){
  
  #Round
  roundnum <- round(num, dig) %>% as.character()
  
  #If rounded to integer, need to add "." to end
  if(length(grep("\\.",roundnum))==0){
    roundnum <- paste0(roundnum, ".")
  }
  
  #Add an extra zero beyond the decimal point if needed to get same length
  #Do num first
  roundnum <- sapply(seq_len(length(roundnum)), function(x){
    if(gsub(".*\\.","",roundnum[x]) %>% nchar() < dig){
      #Needed length
      zerosneeded <- dig - gsub(".*\\.","",roundnum[x]) %>% nchar()
      roundnum[x] <- paste0(roundnum[x],paste0(rep(0,zerosneeded),collapse=""))
    } else{
      roundnum[x]
    }
  })
  
  #Add commas if necessary
  roundnum <- prettyNum(roundnum, ",")
  
  return(roundnum)
}

table <- matrix(NA, nrow=5, ncol=7)

table[1,] <- c("", "0\\%", "25\\%", "50\\%", "Mean", "75\\%", "100\\%")

table[2,] <- c("","(1)","(2)","(3)","(4)", "(5)", "(6)")

table[3,] <- c("Treatment fraction",
               formCoef(reglist[[1]],"treatfrac",1),
               formCoef(reglist[[2]],"treatfrac",1),
               formCoef(reglist[[3]],"treatfrac",1),
               formCoef(reglist[[4]],"treatfrac",1),
               formCoef(reglist[[5]],"treatfrac",1),
               formCoef(reglist[[6]],"treatfrac",1))


table[4,] <- c("",
               formSE(reglist[[1]],"treatfrac",1),
               formSE(reglist[[2]],"treatfrac",1),
               formSE(reglist[[3]],"treatfrac",1),
               formSE(reglist[[4]],"treatfrac",1),
               formSE(reglist[[5]],"treatfrac",1),
               formSE(reglist[[6]],"treatfrac",1))


table[5,] <- c("Mean dep. var.", 
               sapply(c("dist_km_min", "dist_km_first_quartile", "dist_km_median", "dist_km_mean", 
                        "dist_km_third_quartile", "dist_km_max"), function(x){
                          dplyr::select(regdf, all_of(x)) %>% 
                            as.matrix() %>% as.numeric() %>% mean() %>% 
                            formNum(1)
                        }) %>% unname()
)

myxtable <- xtable(table)

caption(myxtable) <- c("Sets move toward closures")

align(myxtable) <- c("l","l",rep("c",6))

label(myxtable) <- "change_distance_closure"

print(myxtable, floating = TRUE, caption.placement="top",sanitize.text.function = identity,
      include.colnames=FALSE,include.rownames=FALSE,table.placement="tb",
      hline.after=NULL,
      add.to.row=list(
        pos = list(0,0,2,nrow(table) - 1,nrow(table)),
        command = c(
          paste0("\\toprule & \\multicolumn{6}{c}{Dependent variable: Distance quantile (km)}  \\\\ "),
          "\\midrule ",
          "\\midrule ",
          "\\midrule ",
          "\\bottomrule "
        )),
      type = "latex",file="Output/Tables/table3.tex")