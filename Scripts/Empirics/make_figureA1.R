rm(list=ls())
setwd("C:/Users/englander/Documents/replication_closures/")
library(dplyr); library(ggplot2)
library(sf); library(msm); library(xtable)
library(purrr); library(lubridate)
library(lfe); library(Formula)
library(parallel); library(tidyr); library(cowplot)
library(viridis); library(lwgeom); library(latex2exp)

options(scipen=999)
options(lfe.threads=24)

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

myThemeStuff <- theme(panel.background = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_line(color = 'black'),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.ticks = element_line(color = "gray5",size=.35),
                      axis.text = element_text(color = "black", size = 8, family="sans"),
                      axis.title = element_text(color = "black", size = 9, family = "sans"),
                      #axis.title.y.right = element_text(angle = 90,hjust=0),
                      axis.title.y = element_text(hjust = .5),
                      legend.key = element_blank(),
                      plot.title = element_text(hjust = 0.5, size = 10), 
                      legend.text=element_text(size=6.5, family = "sans"),
                      legend.title = element_text(size=6.5, family = "sans"),
                      plot.margin = unit(c(0,0,0,0),"in"),
                      plot.tag = element_text(family = "sans", size = 9)
)


#Peru time
Sys.setenv(TZ='America/Lima')

#Load potential closures created in 4. make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

#Create pre-periods
predf <- rbind(
  filter(rddf, tvar==0 & bdist==0) %>% 
    mutate(end = start - 48*3600-1) %>% 
    mutate(start = start - 72*3600, tvar = -3, bin = "lead3_in"),
  filter(rddf, tvar==0 & bdist==0) %>% 
    mutate(end = start - 72*3600-1) %>% 
    mutate(start = start - 96*3600, tvar=-4, bin="lead4_in"),
  filter(rddf, tvar==0 & bdist==0) %>% 
    mutate(end = start - 96*3600-1) %>% 
    mutate(start = start - 120*3600, tvar=-5, bin="lead5_in"),
  filter(rddf, tvar==0 & bdist==0) %>% 
    mutate(end = start - 120*3600-1) %>% 
    mutate(start = start - 144*3600, tvar=-6, bin="lead6_in"),
  filter(rddf, tvar==0 & bdist==0) %>% 
    mutate(end = start - 144*3600-1) %>% 
    mutate(start = start - 168*3600, tvar=-7, bin="lead7_in"),
  filter(rddf, tvar==0 & bdist==0) %>% 
    mutate(end = start - 168*3600-1) %>% 
    mutate(start = start - 192*3600, tvar=-8, bin="lead8_in")
)

#Drop outcome variables and treatment fraction; going to redefine them for pre-bins
predf <- dplyr::select(predf, -tons, -numindivids, -numjuv, -nobs, -sdtons, -numadults, -treatfrac)

#Load full BE from PRODUCE where I have imputed size distribution of non-SNP observations
#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Make fullbe an sf object
besf <- st_multipoint(cbind(fullbe$lon, fullbe$lat))

besf <- st_sfc(besf) %>% st_cast("POINT")

st_crs(besf) <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

besf <- st_sf(besf, fullbe)

besf <- rename(besf, calatime = FechaInicioCala)

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

closed <- rbind(
  mutate(closed, bin = "lead3_in",bdist=0,tvar=-3)  %>% 
    mutate(end = start - 48*3600-1) %>% 
    mutate(start = start - 72*3600),
  mutate(closed, bin = "lead4_in",bdist=0,tvar=-4)  %>% 
    mutate(end = start - 72*3600-1) %>% 
    mutate(start = start - 96*3600),
  mutate(closed, bin = "lead5_in",bdist=0,tvar=-5)  %>% 
    mutate(end = start - 96*3600-1) %>% 
    mutate(start = start - 120*3600),
  mutate(closed, bin = "lead6_in",bdist=0,tvar=-6)  %>% 
    mutate(end = start - 120*3600-1) %>% 
    mutate(start = start - 144*3600),
  mutate(closed, bin = "lead7_in",bdist=0,tvar=-7)  %>% 
    mutate(end = start - 144*3600-1) %>% 
    mutate(start = start - 168*3600),
  mutate(closed, bin = "lead8_in",bdist=0,tvar=-8)  %>% 
    mutate(end = start - 168*3600-1) %>% 
    mutate(start = start - 192*3600))

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

#Apply treatVar over start dates
applyTreat <- function(mybin){
  
  #start dates
  mystarts <- filter(predf, bin==mybin) %>% 
    as.data.frame() %>% dplyr::select(start) %>%
    as.matrix() %>% as.POSIXct() %>% unique()
  
  out <- lapply(mystarts, function(x){
    treatVar(x, mybin)
  })
  
  out <- do.call("rbind", out)
  
  out <- as.data.frame(out) %>%
    dplyr::select(rid, bin, treatfrac)
  
  return(out)
}

#Apply over bins I need treatfrac
(myCores <- detectCores())

cl <- makeCluster(12)

clusterExport(cl, "predf")
clusterExport(cl, "closed")
clusterExport(cl, "treatVar")
clusterExport(cl, "applyTreat")
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(lubridate))

#Apply over rows of predf
treatlist <- parLapply(cl = cl,
                        c("lead3_in","lead4_in","lead5_in","lead6_in","lead7_in","lead8_in"),
                        function(x){
                          
                          applyTreat(x)
                          
                        })

stopCluster(cl)
rm(cl, myCores)

treatdf <- bind_rows(treatlist)

#Join onto predf
predf <- left_join(predf, treatdf, by = c("bin","rid"))

rm(treatlist, treatdf, treatVar, applyTreat)

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

cl <- makeCluster(myCores - 10)

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

myoutcomes <- bind_rows(myoutcomes)

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

#Clean up
rm(i, counter, twoweek, twowk)

rddf$twowk <- as.factor(rddf$twowk)

#Drop potential closures that have NA for size distribution
rddf <- filter(rddf, !is.na(prop12hat))


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



#Unadjusted coefficients
p1 <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(grep("_treatfrac",names(rddf),value=T),collapse="+"),
    "| bin ",
    " | 0 | twoweek_cellid_2p")),
  data = rddf)

#Treatment coefficients with full FE and controls
p2 <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(grep("_treatfrac",names(rddf),value=T),collapse="+"),
    " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(rddf),value=T),collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = rddf)

#Make felm object into a data frame for plotting
plotdf <- function(felmobj){
  
  coefs <- summary(felmobj)[["coefficients"]]
  coefs <- coefs[grep("_in_treatfrac",rownames(coefs)),]
  coefs <- coefs %>% 
    as.data.frame() %>% 
    mutate(bin = rownames(coefs))
  
  #95% confidence intervals
  coefs <- mutate(coefs, ub = Estimate + `Cluster s.e.`*qnorm(.975),
                  lb = Estimate - `Cluster s.e.`*qnorm(.975))
  
  #x-axis values
  coefs$tvar <- "NA"
  
  #For pre-period bins
  coefs$tvar[grep("lead",coefs$bin)] <- paste0("-",
                                               gsub("\\D+","",coefs$bin[grep("lead",coefs$bin)]))
  
  #Want to leave two spaces to represent 9-48 hour before period when sets can generate closure (x=-2,x=-1)
  #Then plot lead9hours at x=0
  coefs$tvar[coefs$bin=="lead9hours_in_treatfrac"] <- "0"
  
  #And active at 2 because closures last three to five days (one space on either side)
  coefs$tvar[coefs$bin=="active_in_treatfrac"] <- "2"
  
  #And lags at their numeric value + 3
  coefs$tvar[grep("lag",coefs$bin)] <- paste0(gsub("\\D+","",coefs$bin[grep("lag",coefs$bin)]))
  
  coefs$tvar <- as.numeric(coefs$tvar)
  
  coefs$tvar[grep("lag",coefs$bin)] <- coefs$tvar[grep("lag",coefs$bin)] + 3
  
  return(coefs)
}

df1 <- plotdf(p1)
df2 <- plotdf(p2)

#Record max and min y-val so axis is consistent between two plots
maxy <- max(
  c(max(df1$ub), max(df2$ub))
)

miny <- min(
  c(min(df1$lb), min(df2$lb))
)

#Make plot given df and title
makeplot <- function(df, tit){

#Text label for period when sets can generate closures
textdf <- data.frame(x=-1.5,y=-.6,lab="Period when sets\ncan generate\nclosure")

#Text label for preperiod
pretext <- data.frame(x=-5.5,y=-.6,lab="Pre-period")

#Text label for after closures end
posttext <- data.frame(x=5,y=-.6,lab="Post-closure")

#Bind together
textdf <- bind_rows(textdf, pretext, posttext)

plot <- ggplot() + 
  geom_point(data=df, aes(x=tvar,y=Estimate)) + 
  geom_errorbar(data=df, aes(x=tvar,ymin=lb,ymax=ub),width=0) + 
  geom_smooth(data=filter(df, tvar <= -3), aes(y=Estimate,x=tvar),
              formula = y~x,method='lm', se=FALSE, col='red') + 
  ggtitle(tit) + 
  geom_text(data=textdf, aes(x=x,y=y,label=lab),size=2.75) + 
  scale_y_continuous("Treatment coefficient and 95% confidence interval",
                     limits = c(miny,maxy),
                     breaks = seq(from=-1,to=3,by=.5)) + 
  scale_x_continuous("",breaks=c(-8,-7,-6,-5,-4,-3,0,2,4,5,6,7),
                     labels = c("-6","-5","-4","-3","-2","-1","Announcement","Closure\nperiod","1","2","3","4")) + 
  myThemeStuff + 
  #geom_vline(aes(xintercept=-2.75),linetype='dashed') + 
  # geom_vline(aes(xintercept=-.25),linetype='dashed') + 
  geom_hline(aes(yintercept=0))


return(plot)
}

#Make two plots
plot1 <- makeplot(df1, "Unadjusted treatment coefficients") + 
  labs(tag = "a") + theme(plot.tag.position = c(.05, 1), 
                          plot.margin = unit(c(.02,0,.1,0),"in"))

plot2 <- makeplot(df2, "Adjusted treatment coefficients") + 
  labs(tag = "b") + theme(plot.tag.position = c(.05, 1), 
                          plot.margin = unit(c(.1,0,.02,0),"in"))

#Stack them
outplot <- plot_grid(plot1, plot2, nrow=2)

ggsave(outplot, file=paste0("Output/Figures/figureA1.png"),
       w=7,h=(7/1.69)*2, units = "in", dpi=1200)

sessionInfo()
