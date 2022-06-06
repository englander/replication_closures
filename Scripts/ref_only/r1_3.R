usethis::use_git_config(user.name = "englander", user.email = "genglander@ucsb.edu")
credentials::set_github_pat("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")

#Identify sets within treatment window of a potential closure
#Calculate distance of each set to potential closure
#If set within treatment window of multiple potential closures, record minimum distance
#First minimum distance in time, then minimum distance in km
#If set inside potential closure, distance is negative
#Finally regress distance on treatment fraction (not interacted with bin dummies), 
#six control variables, length distribution, and FE


rm(list=ls())
setwd("C:/Users/englander/Documents/replication_closures")


library(dplyr); library(ggplot2); library(lfe)
library(lubridate); library(xtable); library(car)
library(purrr); library(readxl); library(readr)
library(parallel); library(sf); library(rworldmap)

options(scipen=999)
#Peru time
Sys.setenv(TZ='America/Lima')

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
                      plot.margin = unit(c(0.01,.01,0.01,.01),"in"),
                      plot.tag = element_text(family = "sans", size = 9)
)

#Created in 3. correct_size_be.R
load("Output/Data/pbe_imp.Rdata")

names(fullbe)

#Create indicator for whether set occurred inside treatment bin with statistically significant change in juvenile catch
fullbe <- mutate(fullbe, sigbin = if_else(lead_0==1 |active_10==1 | active_20==1 | active_30==1 | active_40==1 |
                                            active_50==1 | lag1_0==1 | lag1_10==1 | lag1_20==1 |
                                            lag2_0==1 | lag2_10==1 | lag2_20==1,1,0))

fullbe <- rename(fullbe, calatime = FechaInicioCala)

#Day of sample
fullbe <- mutate(fullbe, date = date(calatime) %>% as.factor())

#sf from full be
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

besf <- st_sf(besf)

dist2coast <- st_distance(besf, peru)

kmtocoast <- as.numeric(dist2coast) / 10^3

fullbe <- mutate(fullbe, kmtocoast = kmtocoast)

rm(dist2coast, kmtocoast, besf)

#Make a simple regression table.

fullbe$Matricula <- as.factor(fullbe$Matricula)
fullbe$Temporada <- as.factor(fullbe$Temporada)
fullbe$twoweek_cellid_2p <- as.factor(fullbe$twoweek_cellid_2p)
fullbe$cellid_2p <- as.factor(fullbe$cellid_2p)

#No FE, actual closures
reg1 <- felm(asinhtons ~ sigbin + kmtocoast | 0 | 0 | twoweek_cellid_2p, data = fullbe)

#FE, actual closures
reg2 <- felm(asinhtons ~ sigbin + kmtocoast | Matricula:Temporada + date + 
               Temporada:cellid_2p | 0 | twoweek_cellid_2p, data = fullbe)

#Mean of sigbin
mean(fullbe$sigbin) #0.3912658

##Calculate whether sets occurred inside significant treatment bin of potential closure
#Created in make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

#Only keep relevant columns of potential closures
rddf <- dplyr::select(rddf, rid, start, end, bin, tvar, bdist)

#Re-make besf
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

besf <- st_sf(geometry = besf, fullbe)

#Drop previously defined sigbin, since going to redefine
besf <- dplyr::select(besf, -sigbin)

#Calculate whether each set is inside significant treatment bin of potential closure
inBuf <- function(myrowind){
  
  row <- besf[myrowind,] 
  
  #Filter to potential closure rows and exclude active_in bins
  actpot <- filter(rddf, start<=row$calatime & row$calatime<=end & bin!="active_in")
  
  inter <- st_intersects(row, actpot)
  
  #Filter to potential closure rows that row is spatially inside
  insidepot <- actpot[unlist(inter),] %>%
    as.data.frame() %>% dplyr::select(-geometry) %>%
    distinct(bin)
  
  if(nrow(insidepot)>0){
    
    #Significant treatment bins
    if(
      filter(insidepot, bin=="lead_0" | bin=="active_10" |  bin=="active_20" |  bin=="active_30" |  bin=="active_40" | 
             bin=="active_50" |  bin=="lag1_0" |  bin=="lag1_10" |  bin=="lag1_20" |
             bin=="lag2_0" |  bin=="lag2_10" |  bin=="lag2_20") %>% nrow() > 0){
      sigbin <- 1
    } else{
      sigbin <- 0
    }
    
  } else{
    sigbin <- 0
  }
  
  out <- as.data.frame(row) %>% dplyr::select(-geometry) %>% 
    mutate(sigbin = sigbin)
  
  return(out)
}

(myCores <- detectCores())

cl <- makeCluster(24)

clusterExport(cl, "besf")
clusterExport(cl, "rddf")
clusterExport(cl, "inBuf")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lubridate))


belist <- parLapply(cl = cl,
                    1:nrow(besf),
                    function(x){
                      
                      try(inBuf(x))
                      
                    })

#Regression data frame for regressions 3 and 4
regdf <- bind_rows(belist)

stopCluster(cl)
rm(cl, myCores)

#No FE, potential closures
reg3 <- felm(asinhtons ~ sigbin + kmtocoast | 0 | 0 | twoweek_cellid_2p, data = regdf)

#FE, potential closures
reg4 <- felm(asinhtons ~ sigbin + kmtocoast | Matricula:Temporada + date + 
               Temporada:cellid_2p | 0 | twoweek_cellid_2p, data = regdf)

#Mean of sigbin relative to potential closures
mean(regdf$sigbin) #0.7985614

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

table <- matrix(NA, nrow=7, ncol=5)

table[1,] <- c("","(1)","(2)","(3)","(4)")

table[2,] <- c("$\\mathbb{1}$\\{Near\\}",
               formCoef(reg1,"sigbin",3),
               formCoef(reg2,"sigbin",3),
               formCoef(reg3,"sigbin",3),
               formCoef(reg4,"sigbin",3))

table[3,] <- c("",
               formSE(reg1,"sigbin",3),
               formSE(reg2,"sigbin",3),
               formSE(reg3,"sigbin",3),
               formSE(reg4,"sigbin",3))


table[4,] <- c("Distance to shore (km)",
               formCoef(reg1, "kmtocoast",3),
               formCoef(reg2, "kmtocoast",3),
               formCoef(reg3, "kmtocoast",3),
               formCoef(reg4, "kmtocoast",3))

table[5,] <- c("",
               formSE(reg1, "kmtocoast",3),
               formSE(reg2, "kmtocoast",3),
               formSE(reg3, "kmtocoast",3),
               formSE(reg4, "kmtocoast",3))

table[6,] <- c("Constant",
               formCoef(reg1,"(Intercept)",3),
               "",
               formCoef(reg3,"(Intercept)",3),
               "")

table[7,] <- c("",
               formSE(reg1,"(Intercept)",3),
               "",
               formSE(reg3,"(Intercept)",3),
               "")

myxtable <- xtable(table)

caption(myxtable) <- c("Closures provide valuable information, but the value of this information is competed away")

align(myxtable) <- c("l","l",rep("c",4))

label(myxtable) <- "information_valuable"

print(myxtable, floating = TRUE, caption.placement="top",sanitize.text.function = identity,
      include.colnames=F,include.rownames=F,table.placement="tb",
      hline.after=NULL,
      add.to.row=list(
        pos = list(0,0,1,nrow(table),nrow(table)),
        command = c(
          paste0("\\toprule & \\multicolumn{4}{c}{Dependent variable: asinh(tons)} \\\\ "),
          "\\midrule & \\multicolumn{2}{c}{Actual closures} & \\multicolumn{2}{c}{Potential closures} \\\\",
          "\\midrule ",
          "\\midrule Fixed effects & & X & & X\\\\ ",
          "\\bottomrule \\multicolumn{5}{l}{\\multirow{2}{12cm}{All regressions have 246,914 observations. $\\mathbb{1}$\\{Near\\} is an indicator for whether the set occurred inside a treatment bin in which there is a significant change in juvenile catch because of the temporary spatial closures policy. In Columns 1 and 2, Near is defined relative to actual closures declared by the regulator (mean of this indicator equals .391). In Columns 3 and 4, Near is defined relative to potential closures (mean of this indicator is .799). Electronic logbook data is for all vessels from April 2017 to January 2020. Regressions in Columns 2 and 4 include vessel by season fixed effects, day-of-sample fixed effects, and two-degree grid cell by season fixed effects. Standard errors clustered at level of two-week-of-sample by two-degree grid cell.}} \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ "
        )),
      type = "latex",file="Output/Tables/table1.tex")


sessionInfo()


