#Estimate effect of policy on (corrected) BE juveniles caught

rm(list=ls())
setwd("C:/Users/englander/Documents/replication_closures")
library(dplyr); library(ggplot2)
library(sf);library(lubridate); library(glmnet)
library(lfe); library(Formula); library(xtable)
library(purrr); library(cowplot); library(tidyr)
library(latex2exp); library(scales)

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

#Load rddf created in make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

rddf <- as.data.frame(rddf) %>% dplyr::select(-geometry)

rddf <- arrange(rddf, tvar, bdist)

#Millions of juveniles
rddf <- mutate(rddf, nummjuv = numjuv/10^6) %>%
  #Tons caught per set
  mutate(tonsperset = tons/nobs) %>%
  #Cluster tons/set and tons/area
  mutate(clusttonsperset = clusttons/clustnobs, clusttonsperarea = clusttons/clustarea_km2) %>%
  #clust millions of juveniles, individuals, and adults
  mutate(clustmjuv = clustnumjuv/10^6, clustmindivids = clustnumindivids/10^6,
         clustmadults = clustnumadults/10^6)

rddf <- mutate(rddf, asinhnummjuv = asinh(nummjuv), 
               asinhtons = asinh(tons))

rddf$bin <- as.factor(rddf$bin)
rddf$bin <- relevel(rddf$bin, ref="active_in")

rddf$twoweek_cellid_2p <- as.factor(rddf$twoweek_cellid_2p)
rddf$twowk <- as.factor(rddf$twowk)
rddf$cellid_2p <- as.factor(rddf$cellid_2p)

#Drop potential closures that have NA for size distribution
rddf <- filter(rddf, !is.na(prop12hat))

#How many clusters are there
unique(rddf$twoweek_cellid_2p) %>% length() #255

#How many observations
nrow(rddf) #34,164 (so dropped 21 potential closures)

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


#Realized juveniles caught
juvcatch <- felm(
  as.Formula(paste0(
  "asinhnummjuv", "~ ", 
  paste0(grep("_treatfrac",names(rddf),value=T),collapse="+"),
  " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
  paste0(grep("prop",names(rddf),value=T),collapse="+"),
  "| bin + twowk:cellid_2p + startdate",
  " | 0 | twoweek_cellid_2p")),
  data =rddf)

jvtab <- summary(juvcatch)[["coefficients"]]

jvtab <- mutate(as.data.frame(jvtab), bin = rownames(jvtab))

jvtab <- jvtab[grep("treatfrac",jvtab$bin),]

jvtab$bin <- gsub("_treatfrac","",jvtab$bin)

#Calculate juv1 in potential closures, 
#then calculate juv0. Then scale up both by ratio of nummjuv in potential closures 
#to total numjuv.
toteffect_juv <- group_by(rddf, bin, tvar, bdist) %>% 
  summarise(juv1 = sum(nummjuv)) %>%
  left_join(jvtab, by = 'bin') %>% 
  mutate(juv0 = juv1/(exp(Estimate))) %>% 
  mutate(chmjuv = juv1 - juv0) %>%
  arrange(tvar, bdist) %>% ungroup()

#Percent effect, not accounting for reallocation
sum(toteffect_juv$chmjuv) / sum(toteffect_juv$juv0) #0.746729

#Scaled total effect see above comment. 
#Created in correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Total effect, not accounting for reallocation
#Rescale (potential closures nummjuv larger because of double-counting.)
(sum(fullbe$numjuv,na.rm=T) / sum(rddf$numjuv, na.rm=T)) * sum(toteffect_juv$chmjuv) #60162.01

toteffect_juv <- dplyr::select(toteffect_juv, -`t value`, -`Pr(>|t|)`)

toteffect_juv <- rename(toteffect_juv, juvcoef = Estimate, juvse = `Cluster s.e.`)

#Calculate delta method standard errors for change in millions of juveniles caught
library(msm)

toteffect_juv <- mutate(toteffect_juv, chmjuvse = as.numeric(NA), 
                        chmjuv_scaled = chmjuv * (sum(fullbe$numjuv,na.rm=T) / sum(rddf$numjuv, na.rm=T)),
                        chmjuvse_scaled = as.numeric(NA), 
                        juv0se = as.numeric(NA))

scaleconstant <- (sum(fullbe$numjuv,na.rm=T) / sum(rddf$numjuv, na.rm=T)) #.394

for(mybin in toteffect_juv$bin){
  
  #Filter toteffect_juv to mybin
  mydf <- filter(toteffect_juv, bin==mybin)
  
  myjuv1 <- mydf$juv1
  myjuvcoef <- mydf$juvcoef
  myjuvvcov <- mydf$juvse^2

  #Delta se for juv0. Random variable being transformed is myjuvcoef
  juv0_delta <- deltamethod(~ (myjuv1 / (exp(x1))), myjuvcoef, myjuvvcov, ses=T)
  
  #Plug this value into toteffect_juv
  toteffect_juv$juv0se[toteffect_juv$bin==mybin] <- juv0_delta
  
  #Delta se for chmjuvse. Random variable being transformed is myjuvcoef
  chmjuv_delta <- deltamethod( ~ myjuv1 - (myjuv1/(exp(x1))),
                            myjuvcoef,
                            myjuvvcov,
                            ses=T
  )
  
  #Plug this value into toteffect_juv
  toteffect_juv$chmjuvse[toteffect_juv$bin==mybin] <- chmjuv_delta
  
  #Delta se for scaled chmjuvse. Random variable being transformed is myjuvcoef
  chmjuv_scaled_delta <- deltamethod( ~ (myjuv1 - (myjuv1/(exp(x1))))*scaleconstant,
                               myjuvcoef,
                               myjuvvcov,
                               ses=T
  )
  
  #Plug this value into toteffect_juv
  toteffect_juv$chmjuvse_scaled[toteffect_juv$bin==mybin] <- chmjuv_scaled_delta
  
}


##Normalize by area
keeprddf <- rddf

#Re-load rddf created in make_rddf.R to get back geometry column
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

#Use same potential closures that I used in estimation
rddf <- filter(rddf, !is.na(prop12hat))

rddf$area_km2 <- st_area(rddf)
rddf$area_km2 <- as.numeric(rddf$area_km2)/ 10^6

#Now can calculate total area for each bdist
annulus_areas <- group_by(as.data.frame(rddf) %>% dplyr::select(-geometry),
                          bdist) %>% 
  summarise(area_km2 = sum(area_km2))

toteffect_juv <- left_join(toteffect_juv, annulus_areas, by = 'bdist')

#Now re-define rddf to be original without geometry
rddf <- keeprddf; rm(keeprddf)

#Normalized scaled chmjuv
toteffect_juv <- mutate(toteffect_juv, chmjuv_scaled_relarea = (chmjuv_scaled/area_km2)*annulus_areas$area_km2[annulus_areas$bdist==0])

##Calculate total effect, not accounting for reallocation of tons caught. (in millions)
changejuv <- (sum(fullbe$numjuv,na.rm=T) / sum(rddf$numjuv, na.rm=T)) * sum(toteffect_juv$chmjuv)

#Decomposing effects:
#Direct effect
toteffect_juv$chmjuv_scaled[toteffect_juv$bin=="active_in"]  #-2176.79

#Inside, post announcement
toteffect_juv$chmjuv_scaled[toteffect_juv$bin=="lead9hours_in"] #1182.788


#Contemporaneous spatial spillovers
toteffect_juv$chmjuv_scaled[toteffect_juv$bin=="active_10" | toteffect_juv$bin=="active_20" | toteffect_juv$bin=="active_30" | 
                              toteffect_juv$bin=="active_40" | toteffect_juv$bin=="active_50"] %>% sum() #44478.46

#Inside, day after
toteffect_juv$chmjuv_scaled[toteffect_juv$bin=="lag1_in"]  #2262.906

#Inside, two days after
toteffect_juv$chmjuv_scaled[toteffect_juv$bin=="lag2_in"]  #2258.567

#Total effect, not accounting for reallocation
sum(toteffect_juv$chmjuv_scaled) #60163.7

#Closures cannot increases tons caught because of TAC, so account for this reallocation
#by estimating how policy affects tons caught on average across treatment bins in treatment window
tonscaught <- felm(
  as.Formula(paste0(
    "asinhtons", "~ ", 
    "treatfrac ",
    " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(rddf),value=T),collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data =rddf)

#Increase tons by 35% (exp( 0.29933719 )-1)
summary(tonscaught)[["coefficients"]]["treatfrac",]

tonscoef <- summary(tonscaught)[["coefficients"]]["treatfrac","Estimate"]

#Don't reallocate tons from 2017 second season or 2019 second season because they were both 
#shut down well before TAC was reached. 

#Tons caught in state of world I observe in seasons where TAC hit
tons1 <- sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"])

#Change in tons because of policy
(ctons <- tons1 - tons1/exp(tonscoef))

#Average pj outside of treatment window
(avgpjoutside <- filter(fullbe, lead_0==0 & lead_10==0 & lead_20==0 & lead_30==0 & lead_40==0 & lead_50==0 & 
                         active_0==0 & active_10==0 & active_20==0 & active_30==0 & active_40==0 & active_50==0 & 
                         lag1_0==0 & lag1_10==0 & lag1_20==0 & lag1_30==0 & lag1_40==0 & lag1_50==0 & 
                         lag2_0==0 & lag2_10==0 & lag2_20==0 & lag2_30==0 & lag2_40==0 & lag2_50==0 & 
                         lag3_0==0 & lag3_10==0 & lag3_20==0 & lag3_30==0 & lag3_40==0 & lag3_50==0 & 
                         lag4_0==0 & lag4_10==0 & lag4_20==0 & lag4_30==0 & lag4_40==0 & lag4_50==0 & 
                         !is.na(numindivids) & !is.na(bepjhat) & Temporada!="2017-II" & Temporada!="2019-II") %>%
  #Weight by number of individuals
  mutate(pjweighted = bepjhat*numindivids) %>%  
  summarise(perjuv = sum(pjweighted)/sum(numindivids)) %>% as.numeric() / 100) #0.09045436


#Avg weight of individual caught outside of treatment window
(avgweightoutside <- filter(fullbe, lead_0==0 & lead_10==0 & lead_20==0 & lead_30==0 & lead_40==0 & lead_50==0 & 
                             active_0==0 & active_10==0 & active_20==0 & active_30==0 & active_40==0 & active_50==0 & 
                             lag1_0==0 & lag1_10==0 & lag1_20==0 & lag1_30==0 & lag1_40==0 & lag1_50==0 & 
                             lag2_0==0 & lag2_10==0 & lag2_20==0 & lag2_30==0 & lag2_40==0 & lag2_50==0 & 
                             lag3_0==0 & lag3_10==0 & lag3_20==0 & lag3_30==0 & lag3_40==0 & lag3_50==0 & 
                             lag4_0==0 & lag4_10==0 & lag4_20==0 & lag4_30==0 & lag4_40==0 & lag4_50==0 & 
                             !is.na(numindivids) & !is.na(avgweightg) & Temporada!="2017-II" & Temporada!="2019-II") %>%
  #Weight by number of individuals
  mutate(weightweighted = avgweightg*numindivids) %>% 
  summarise(avgweightg = sum(weightweighted)/sum(numindivids)) %>% as.numeric())


#Decrease in individuals caught outside of treatment window in millions
#(converting tons to g cancels out conversion to millions)
chindividsoutside <- -ctons/avgweightoutside

chjuvsoutside <- chindividsoutside*avgpjoutside

#Now can calculate change in juvenile catch due to policy, accounting for reallocation
(chmjuvsstart <- changejuv + chjuvsoutside) #46829.39

#How many juveniles are caught during my sample period in total?
#F(1)*pj*individuals/VMS fishing obs
juv1 <- sum(fullbe$numjuv, na.rm=T) / 10^6 

#chmjuvsstart = juv(1) - juv(0)
juv0 <- juv1 - chmjuvsstart 

#Then increase in juvenile catch as a percentage is 
chmjuvsstart / juv0 # 0.4986941

#Calculate standard error on total change in juvenile catch and in total percentage change
mycoefs <- toteffect_juv$chmjuv_scaled
mybigvcov <- diag(toteffect_juv$chmjuvse_scaled^2)

#This includes reallocation, so this is what I want: 
(changebillionsse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                          x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                          x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                          x31 + x32 + x33 + x34 + x35 + x36) + chjuvsoutside) / 
                      1000, mycoefs, mybigvcov, ses=T))
#5.147492

#Now get SE on total percentage change
(totperse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                          x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                          x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                          x31 + x32 + x33 + x34 + x35 + x36) + chjuvsoutside) / 
                      juv0, mycoefs, mybigvcov, ses=T))
#0.0548165

#What is average percentage juvenile among sets within treatment window? (for text of paper)
filter(fullbe, (lead_0==1 | lead_10==1 | lead_20==1 | lead_30==1 | lead_40==1 | lead_50==1 | 
         active_0==1 | active_10==1 | active_20==1 | active_30==1 | active_40==1 | active_50==1 | 
         lag1_0==1 | lag1_10==1 | lag1_20==1 | lag1_30==1 | lag1_40==1 | lag1_50==1 | 
         lag2_0==1 | lag2_10==1 | lag2_20==1 | lag2_30==1 | lag2_40==1 | lag2_50==1 | 
         lag3_0==1 | lag3_10==1 | lag3_20==1 | lag3_30==1 | lag3_40==1 | lag3_50==1 | 
         lag4_0==1 | lag4_10==1 | lag4_20==1 | lag4_30==1 | lag4_40==1 | lag4_50==1) & 
  !is.na(numindivids) & !is.na(bepjhat)) %>%
  #Weight by tons
  mutate(pjweighted = bepjhat*numindivids) %>%  
  summarise(perjuv = sum(pjweighted)/sum(numindivids)) %>% as.numeric() / 100  #0.2537528

#Average percentage juvenile inside closure between announcement and beginning
filter(fullbe, lead_0==1 &
         !is.na(numindivids) & !is.na(bepjhat)) %>%
  #Weight by tons
  mutate(pjweighted = bepjhat*numindivids) %>%  
  summarise(perjuv = sum(pjweighted)/sum(numindivids)) %>% as.numeric() / 100 #0.4677222


##Make main results figures
#Output results in billions rather than millions
toteffect_juv <- mutate(toteffect_juv, chbjuv_scaled = chmjuv_scaled / 10^3,
                        chbjuvse_scaled  = chmjuvse_scaled / 10^3,
                        chbjuv_scaled_relarea = chmjuv_scaled_relarea / 10^3)

#Confidence intervals for chbjuv_scaled
toteffect_juv <- mutate(toteffect_juv, 
                        chbjuv_scaled_ub = chbjuv_scaled + chbjuvse_scaled*qnorm(.975),
                        chbjuv_scaled_lb = chbjuv_scaled - chbjuvse_scaled*qnorm(.975))

finaldf <- toteffect_juv

finaldf$tvar[finaldf$tvar==-1] <- "lead"
finaldf$tvar[finaldf$tvar==0] <- "active"
finaldf$tvar[finaldf$tvar%in%c(1,2,3,4)] <- paste0("lag",finaldf$tvar[finaldf$tvar%in%c(1,2,3,4)])

#Output plot for change in billions of juveniles caught

#Function of one event day and dependent variable
singlePlot <- function(myvar, mytvar, ylab){
  
  #Title
  if(mytvar=="lead"){
    tit <- "After announcement, before closure"
  } else if(mytvar=="active"){
    tit <- "Closure period"
  } else if(mytvar=="lag1"){
    tit <- "1 day after"
  } else{
    tit <- paste0(gsub("lag","",mytvar), " days after")
  }
  
  #Rename desired variable
  usedf <- finaldf
  names(usedf)[names(usedf)==myvar] <- "plotvar"
  
  #Rename lb and ub so can refer to directly
  names(usedf)[names(usedf)==paste0(myvar,"_lb")] <- "lb"
  names(usedf)[names(usedf)==paste0(myvar,"_ub")] <- "ub"
  
  #Want consistent y range across given variable
  myymin <- min(usedf$lb)
  myymax <- max(usedf$ub)
  
  plot <- ggplot(data=filter(usedf, tvar==mytvar),aes(x=bdist)) + 
    geom_hline(aes(yintercept=0)) + 
    geom_point(aes(y=plotvar)) + 
    scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                       labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) +
    scale_y_continuous("",limits = c(myymin,myymax),
                       breaks = breaks_pretty(n=10),
                       labels = label_number(accuracy=1)) + 
    myThemeStuff + 
    ggtitle(tit) + 
    geom_errorbar(data=filter(usedf, tvar==mytvar), aes(x=bdist,ymin=lb,ymax=ub),width=0)
  
  
  
  #If first plot of panel, plot relarea as hollow point but no legend
  if((mytvar!="lead"&mytvar!="lag2") & myvar=="chbjuv_scaled"){
    
    #Rename relarea so can refer to directly
    names(usedf)[names(usedf)==paste0(myvar,"_relarea")] <- "relarea"
    
    plot <- plot + 
      geom_point(data=filter(usedf, tvar==mytvar), aes(x=bdist,y=relarea),shape=2,col='red')
  }
  if((mytvar=="lead"|mytvar=="lag2") & myvar=="chbjuv_scaled"){
    
    #Rename relarea so can refer to directly
    names(usedf)[names(usedf)==paste0(myvar,"_relarea")] <- "relarea"
    
    #Gather so can have legend
    legdf <- gather(
      dplyr::select(usedf, plotvar, relarea, tvar, bdist), `Area-normalized`, pointvar,
      -tvar, -bdist
    )
    
    legdf$`Area-normalized`[legdf$`Area-normalized`=="plotvar"] <- "No"
    legdf$`Area-normalized`[legdf$`Area-normalized`=="relarea"] <- "Yes"
    legdf$`Area-normalized` <- as.factor(legdf$`Area-normalized`)
    
    plot <- ggplot() + 
      geom_errorbar(data=filter(usedf, tvar==mytvar), aes(x=bdist,ymin=lb,ymax=ub),width=0) + 
      geom_point(data=filter(legdf, tvar==mytvar),aes(x=bdist,y=pointvar,
                                                      shape=`Area-normalized`,col=`Area-normalized`)) + 
      scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                         labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) + 
      scale_y_continuous(ylab,limits = c(myymin,myymax),
                         breaks = breaks_pretty(n=10),
                         labels = label_number(accuracy=1)) + 
      scale_shape_manual(values = c(16,2)) + 
      scale_color_manual(values = c("black","red")) + 
      myThemeStuff + 
      geom_hline(aes(yintercept=0)) + 
      ggtitle(tit) + 
      #Put legend in top right corner
      theme(legend.position = c(0.5,.88), 
            legend.margin = margin(0,0,0,0,unit="cm"), 
            legend.key.height = unit(.25, unit = "cm"),
            legend.key.width=unit(0,unit="cm"), 
            legend.key.size=unit(0,unit="cm")
      )
    
  }
 
  return(plot)
}

#Make a 2x3 plot for paper with labeled panels
paperFig <- function(myvar, ylab){
  leadplot <- singlePlot(myvar, "lead", ylab) + 
    labs(tag = "a") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  activeplot <- singlePlot(myvar, "active", ylab)+ 
    labs(tag = "b") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  lag1plot <- singlePlot(myvar, "lag1", ylab)+ 
    labs(tag = "c") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  lag2plot <- singlePlot(myvar, "lag2", ylab)+ 
    labs(tag = "d") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  lag3plot <- singlePlot(myvar, "lag3", ylab)+ 
    labs(tag = "e") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  lag4plot <- singlePlot(myvar, "lag4", ylab)+ 
    labs(tag = "f") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  
  tbt <- plot_grid(leadplot, activeplot, lag1plot, 
                   lag2plot, lag3plot, lag4plot, nrow=2, ncol=3, 
                   rel_widths = c(1.01,1,1))
  
  ggsave(tbt, file=paste0("Output/Figures/fig8.png"),
         w=7,h=(7/1.69)*2, units = "in", dpi=1200)
}


paperFig("chbjuv_scaled", "Change in billions of juveniles caught")

sessionInfo()
