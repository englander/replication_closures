rm(list=ls())

library(dplyr); library(ggplot2)
library(sf); library(msm); library(xtable)
library(purrr); library(lubridate)
library(lfe); library(Formula)
library(parallel); library(tidyr); library(cowplot)
library(viridis); library(lwgeom)

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

#Create 48-72 hour leads for potential closures
#This is the period before sets would generate a potential closure
#Wouldn't want to test pre-trends among sets that could generate a potential closure (24-48 hours before),
#because there would be mechanical correlation
predf <- filter(rddf, tvar==0 & bdist==0) %>% 
  mutate(end = start - 48*3600-1) %>% 
  mutate(start = start - 72*3600)

predf$tvar <- -3

predf$bin <- gsub("active","lead3",predf$bin)

#Drop outcome variables and treatment fraction; going to redefine them for 48-72 hours before bins
predf <- dplyr::select(predf, -tons, -numindivids, -numjuv, -nobs, -sdtons, -numadults, -treatfrac)

#Load full BE from PRODUCE where I have imputed size distribution of non-SNP observations
#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

fullbe <- rename(fullbe, calatime = FechaInicioCala)

#Make fullbe an sf object
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

besf <- st_sf(besf, fullbe)

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

closed$days <- as.character(closed$days) %>% as.integer()

closed <- mutate(closed, bin = "lead3_in",bdist=0,tvar=-3)  %>% 
  mutate(end = start - 48*3600-1) %>% 
  mutate(start = start - 72*3600)


#Calculate rectangle overlap with actual closures
#Function of start and bin
treatVar <- function(mystart, mybin){
  
  #Rectangles
  myrect <- filter(predf, start==mystart & bin==mybin) %>% 
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

#Apply over all bins
(myCores <- detectCores())

cl <- makeCluster(6)

clusterExport(cl, "predf")
clusterExport(cl, "closed")
clusterExport(cl, "treatVar")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lubridate))

rdf <- parLapply(cl = cl,
                 unique(predf$start),
                 function(x){
                   
                   treatVar(x,"lead3_in")
                   
                 })

stopCluster(cl)
rm(cl, myCores)

rdf <- do.call("rbind",rdf)

predf <- rdf

outcomesFun <- function(rdrow){
  
  row <- predf[rdrow,]
  
  #Filter BE observations to within same time period
  mybesf <- filter(besf, row$start<=calatime & calatime<row$end)
  
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
  out <- data.frame(bin=row$bin, rid=row$rid,tons=tons,
                    #Round number of individuals to integer
                    numindivids=round(numindivids), numjuv=round(numjuv),
                    nobs = nobs, sdtons = sdtons)
  
  return(out)
}

#Apply over all bins
(myCores <- detectCores())

cl <- makeCluster(6)

clusterExport(cl, "predf")
clusterExport(cl, "besf")
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

myoutcomes <- do.call("rbind",myoutcomes)

#Join onto predf
predf <- left_join(predf, myoutcomes, by = c("bin",'rid'))

predf <- as.data.frame(predf) %>% dplyr::select(-geometry)
rddf <- as.data.frame(rddf) %>% dplyr::select(-geometry)

predf <- mutate(predf, numadults = numindivids - numjuv)

#Row bind with rddf
rddf <- bind_rows(rddf, predf)

#Millions of juveniles
rddf <- mutate(rddf, nummjuv = numjuv/10^6) %>%
  #Cluster tons/set and tons/area
  mutate(clusttonsperset = clusttons/clustnobs, clusttonsperarea = clusttons/clustarea_km2) %>%
  #clust millions of juveniles, individuals, and adults
  mutate(clustmjuv = clustnumjuv/10^6, clustmindivids = clustnumindivids/10^6,
         clustmadults = clustnumadults/10^6)

rddf <- mutate(rddf, asinhnummjuv = asinh(nummjuv))

rddf$bin <- as.factor(rddf$bin)
rddf$bin <- relevel(rddf$bin, ref="active_in")

rddf$twoweek_cellid_2p <- as.factor(rddf$twoweek_cellid_2p)
rddf$cellid_2p <- as.factor(rddf$cellid_2p)

#Drop potential closures that have NA for size distribution
rddf <- filter(rddf, !is.na(prop12hat))

#Need to redefine start date for new bins and twowk
rddf$startdate <- date(rddf$start) %>% as.factor()

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
rddf <- left_join(dplyr::select(rddf, -twowk), dplyr::select(twoweek, rid, bin, twowk), by = c('rid','bin'))

rddf$twowk <- as.factor(rddf$twowk)



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


#Regression table for lead3_in, from no controls until full control specification
reg1 <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(grep("_treatfrac",names(rddf),value=T),collapse="+"),
    "| bin ",
    " | 0 | twoweek_cellid_2p")),
  data = rddf)

#Fixed effects only
reg2 <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(grep("_treatfrac",names(rddf),value=T),collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = rddf)

#Length distribution controls only, no fixed effects
reg3 <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(grep("_treatfrac",names(rddf),value=T),collapse="+"),
    " + ", 
    paste0(grep("prop",names(rddf),value=T),collapse="+"),
    "| bin",
    " | 0 | twoweek_cellid_2p")),
  data = rddf)


#All controls, no fixed effects
reg4 <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(grep("_treatfrac",names(rddf),value=T),collapse="+"),
    " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(rddf),value=T),collapse="+"),
    "| bin",
    " | 0 | twoweek_cellid_2p")),
  data = rddf)

#All controls and fixed effects
reg5 <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(grep("_treatfrac",names(rddf),value=T),collapse="+"),
    " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(rddf),value=T),collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = rddf)


#Make a regression table
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

table <- matrix(NA, nrow=3, ncol=6)

table[1,] <- c("","(1)", "(2)","(3)", "(4)","(5)")

table[2,] <- c("Treatment fraction",formCoef(reg1, "lead3_in_treatfrac", 3),formCoef(reg2, "lead3_in_treatfrac", 3),
               formCoef(reg3, "lead3_in_treatfrac", 3), formCoef(reg4, "lead3_in_treatfrac", 3),
               formCoef(reg5, "lead3_in_treatfrac", 3))

table[3,] <- c("",formSE(reg1, "lead3_in_treatfrac", 3),formSE(reg2, "lead3_in_treatfrac", 3),
               formSE(reg3, "lead3_in_treatfrac", 3), formSE(reg4, "lead3_in_treatfrac", 3),
               formSE(reg5, "lead3_in_treatfrac", 3))

myxtable <- xtable(table)

align(myxtable) <- c("l","l",rep("c",5))

caption(myxtable) <- "Test for difference in pre-period juvenile catch"#: Conditional on controls, no relationship between juvenile catch and closures before closures are announced")

label(myxtable) <- "preperiod_juvenilecatch"

print(myxtable, floating = TRUE, caption.placement="top",sanitize.text.function = identity,
      include.colnames=F,include.rownames=F,table.placement="tb",
      hline.after=NULL,
      add.to.row=list(
        pos = list(0,1,nrow(table), nrow(table),nrow(table),nrow(table)),
        command = c(
          paste0("\\toprule \\multicolumn{6}{c}{Dependent variable: asinh(juvenile catch)} \\\\ \\midrule "),
          " \\midrule ",
          " \\midrule Fixed effects & & X & & & X\\\\ ",
          "Length distribution & & & X & X & X \\\\",
          "Other controls & & & & X & X \\\\",
          "\\bottomrule \\multicolumn{6}{l}{\\multirow{2}{12cm}{Notes: All regressions have 35,113 observations. Dependent variable is the inverse hyperbolic sine of millions of juveniles caught. All regressions estimate treatment effects for all 37 treatment bins, but only the coefficent on treatment fraction for the inside, day-before bin is displayed in this table. Standard errors clustered at level of two-week-of-sample by two-degree grid cell.}} \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ "
        )),
      type = "latex",file="Output/Tables/tableA1.tex")

sessionInfo()
