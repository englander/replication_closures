#Make figures 9 and 10

rm(list=ls())

library(dplyr); library(ggplot2); library(lfe)
library(lubridate); library(sf); library(readxl)
library(purrr); library(parallel); 
library(xtable); library(lwgeom); library(tidyr)
library(Formula); library(msm); library(cowplot)

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

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


`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

options(lfe.threads=24)

#Peru time
Sys.setenv(TZ='America/Lima')

#Created in 4. make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

#Created in make_fleetthere_selfthere
load("Output/Data/fleetthere_selfthere.Rdata")

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
  interVars("treatfrac")
)

heterodf$value <- as.factor(heterodf$value)

##Make Figure 10 first

#Fleet there
f3 <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(
      paste0(grep("_treatfrac",names(heterodf),value=T), ":value"),
      collapse="+"),
    " + clustnobs + clusttons + clustarea_km2 + 
    kmtocoast + clusttonsperset + clusttonsperarea + ", 
      paste0(grep("prop",names(heterodf),value=T),collapse="+"),
    "| bin:value + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = filter(heterodf, type=="fleet_already3"))

summary(f3)

f3_df <- f3$df

fleettab <- summary(f3)[["coefficients"]]

fleettab <- mutate(as.data.frame(fleettab), bin = rownames(fleettab))

fleettab <- fleettab[grep("treatfrac",fleettab$bin),]

fleettab$bin <- gsub("_treatfrac","",fleettab$bin)

fleettab$fleetthere <- 0
fleettab$fleetthere[grep("value1",fleettab$bin)] <- 1

fleettab$bin <- gsub("value0","",fleettab$bin)
fleettab$bin <- gsub("value1","",fleettab$bin)
fleettab$bin <- gsub(":","",fleettab$bin)

#Calculate percent effect
fleettab <- mutate(fleettab, perchange = (exp(Estimate)-1)*100, 
                perchangese = as.numeric(NA))

#Calculate delta method standard errors for change in millions of juveniles caught
for(i in 1:nrow(fleettab)){
  
  #Filter toteffect_juv to mybin
  mydf <- fleettab[i,]
  
  mycoef <- mydf$Estimate
  myvcov <- mydf$`Cluster s.e.`^2
  
  #Delta se. Random variable being transformed is mycoef
  mypcse <- deltamethod(~ (exp(x1)-1)*100, mycoef, myvcov, ses=T)
  
  #Plug this value into fleettab
  fleettab$perchangese[i] <- mypcse
  
}

rm(mydf, mycoef, myvcov, i, mypcse)

#Confidence intervals
fleettab <- mutate(fleettab, pc_ub = perchange + perchangese*qnorm(.975),
                pc_lb = perchange - perchangese*qnorm(.975))

#Rename lead bins to match rddf naming convention
fleettab$bin <- gsub("lead1_","lead9hours_",fleettab$bin) 

#Get tvar and bdist from rddf
fleettab <- left_join(fleettab, 
                   distinct(as.data.frame(rddf) %>% dplyr::select(-geometry), bin, tvar, bdist),
                   by = c("bin")
)


#Plot types for same bdist next to each other but slightly offset
fleettab$bdist[fleettab$fleetthere==0] <- fleettab$bdist[fleettab$fleetthere==0] - 1
fleettab$bdist[fleettab$fleetthere==1] <- fleettab$bdist[fleettab$fleetthere==1] + 1

fleettab$fleetthere <- as.factor(fleettab$fleetthere)
fleettab$fleetthere <- relevel(fleettab$fleetthere, ref = "0")


#Make a subfigure (for one of six time periods)
singlePlot <- function(mytvar){
  
  #Title
  if(mytvar==-1){
    tit <- "After announcement, before closure"
  } else if(mytvar==0){
    tit <- "Closure period"
  } else if(mytvar==1){
    tit <- "1 day after"
  } else{
    tit <- paste0(mytvar, " days after")
  }
  
  #Want consistent y range across given variable
  myymin <- min(fleettab$pc_lb)
  myymax <- max(fleettab$pc_ub)
  
  #If not first plot of panel, no legend
  if((mytvar!= -1 & mytvar!= 2)){
    plot <- ggplot(data=filter(fleettab, tvar==mytvar),aes(x=bdist,col=fleetthere,shape=fleetthere)) + 
      geom_hline(aes(yintercept=0)) + 
      geom_point(aes(y=perchange)) + 
      scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                         labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) +
      scale_y_continuous("",limits = c(myymin,myymax),
                         breaks = seq(from = -100, to=800, by = 100),
                         labels = paste0(seq(from = -100, to=800, by = 100), "%")) + 
      myThemeStuff + 
      ggtitle(tit) + 
      scale_color_manual("Firm fished\nthere the\nday before",labels=c("No","Yes"),
                         values = c("dodgerblue4","orange1")) + 
      scale_shape_manual("Firm fished\nthere the\nday before", labels=c("No","Yes"),
                         values = c(17,19)) + 
      geom_errorbar(data=filter(fleettab, tvar==mytvar), aes(ymin=pc_lb,ymax=pc_ub),width=0) + 
      theme(legend.position='none')
  } else{
    
    plot <- ggplot(data=filter(fleettab, tvar==mytvar),aes(x=bdist,col=fleetthere,shape=fleetthere)) + 
      geom_hline(aes(yintercept=0)) + 
      geom_point(aes(y=perchange)) + 
      scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                         labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) +
      scale_y_continuous("Percent change in juveniles caught",limits = c(myymin,myymax),
                         breaks = seq(from = -100, to=800, by = 100),
                         labels = paste0(seq(from = -100, to=800, by = 100), "%")) + 
      myThemeStuff + 
      ggtitle(tit) + 
      scale_color_manual("Firm fished\nthere the\nday before",labels=c("No","Yes"),
                         values = c("dodgerblue4","orange1")) + 
      scale_shape_manual("Firm fished\nthere the\nday before", labels=c("No","Yes"),
                         values = c(17,19)) + 
      geom_errorbar(data=filter(fleettab, tvar==mytvar), aes(ymin=pc_lb,ymax=pc_ub),width=0) + 
      #Legend placement
      theme(legend.position = c(0.8,.85), 
            legend.margin = margin(0,0,0,0,unit="cm"), 
            legend.key.height = unit(.25, unit = "cm"),
            legend.key.width=unit(0,unit="cm"), 
            legend.key.size=unit(0,unit="cm"),
            legend.text=element_text(size=9, family = "sans"),
            legend.title = element_text(size=9, family = "sans"))
  }
  
  
  
  return(plot)
  
}

#Make a 2x3 plot for paper with labeled panels
paperFig <- function(myvar){
  leadplot <- singlePlot(-1) + 
    labs(tag = "a") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  activeplot <- singlePlot(0) + 
    labs(tag = "b") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  lag1plot <- singlePlot(1) + 
    labs(tag = "c") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  lag2plot <- singlePlot(2) + 
    labs(tag = "d") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  lag3plot <- singlePlot(3) + 
    labs(tag = "e") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  lag4plot <- singlePlot(4) + 
    labs(tag = "f") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  
  tbt <- plot_grid(leadplot, activeplot, lag1plot, 
                   lag2plot, lag3plot, lag4plot, nrow=2, ncol=3, 
                   rel_widths = c(1.01,1,1))
  
  ggsave(tbt, file=paste0("Output/Figures/figure10.png"),
         w=7,h=(7/1.69)*2, units = "in", dpi=1200)
}

#Make Figure 10
paperFig("hetero_fleetthere")


##Now calculate total change % for each fleet type


tonscaught <- felm(
  as.Formula(paste0(
    "asinhtons", "~ ", 
    "treatfrac:value",
    " + clustnobs + clusttons + clustarea_km2 + 
    kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(heterodf),value=T),collapse="+"),
    "| bin:value + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = filter(heterodf, type=='fleet_already3'))


tonscoef <- summary(tonscaught)[["coefficients"]]
tonscoef <- tonscoef[grep("treatfrac",rownames(tonscoef)),]
tonscoef <- as.data.frame(tonscoef) %>% 
  mutate(fleetthere = gsub("treatfrac:value","",rownames(tonscoef)))

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Rename lead bins to match rddf naming convention
heterodf$bin <- as.character(heterodf$bin)
heterodf$bin <- gsub("lead1_","lead9hours_",heterodf$bin) 
#Remake as factor
heterodf$bin <- as.factor(heterodf$bin)
heterodf$bin <- relevel(heterodf$bin, ref = "active_in")


#Don't reallocate tons from 2017 second season or 2019 second season because they were both 
#shut down well before TAC was reached. 

#Function of fleet type
effectFleetThere <- function(myvalue){
  
  #Tons caught in state of world I observe in seasons where TAC hit
  tons1 <- filter(heterodf, type=="fleet_already3" & value==myvalue) %>% 
    summarise(sum(tons)) %>% as.matrix() %>% as.numeric()
    
  #Scale it by tons in heterodf compared to total tons in data
  tons1 <- tons1 * 
    (sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"]) / 
       filter(heterodf, type=="fleet_already3") %>% 
       summarise(sum(tons)) %>% as.matrix() %>% as.numeric()
    )
  
  #Change in tons because of policy
  ctons <- tons1 - tons1/exp(tonscoef$Estimate[tonscoef$fleetthere==myvalue]) 
  
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
  toteffect_juv <- filter(heterodf, type=="fleet_already3" & value==myvalue) %>%
    group_by(bin, tvar, bdist) %>% 
    summarise(juv1 = sum(nummjuv)) %>% ungroup()
  
  toteffect_juv <- left_join(toteffect_juv, 
                             filter(fleettab, fleetthere==myvalue) %>% 
                               dplyr::select(Estimate, bin, `Cluster s.e.`),
                             by = 'bin') %>% 
    mutate(juv0 = juv1/(exp(Estimate))) %>% 
    mutate(chmjuv = juv1 - juv0) %>%
    arrange(tvar, bdist) %>% ungroup()
  
  #Scale number of juv in fullbe by number of juv in heterodf
  scaleconstant <- (sum(fullbe$numjuv,na.rm=T) / sum(heterodf$numjuv[heterodf$type=="fleet_already3"], na.rm=T))
  
  changejuv <- sum(toteffect_juv$chmjuv) * scaleconstant
  
  #Now can calculate change in juvenile catch due to policy, accounting for reallocation
  chmjuvsstart <- changejuv + chjuvsoutside
  
  #How many juveniles are caught during my sample period in total?
  #Scale by percentage of juveniles caught in value bin
  juv1 <- (sum(fullbe$numjuv, na.rm=T) / 10^6) * 
    (sum(heterodf$numjuv[heterodf$type=="fleet_already3" & heterodf$value==myvalue], na.rm=T) / 
    sum(heterodf$numjuv[heterodf$type=="fleet_already3"], na.rm=T))

  #chmjuvsstart = juv(1) - juv(0)
  juv0 <- juv1 - chmjuvsstart
  
  #Change in juvenile catch as percentage
  totper <- chmjuvsstart / juv0
  
  out <- list(data.frame(fleetthere=myvalue, chmjuvsstart = chmjuvsstart, totper = totper),
              list(toteffect_juv, scaleconstant, chjuvsoutside, juv0)
  )
  
  
  return(out)
  
}

#Fleet wasn't there
effectFleetThere(0)[[1]]
# fleetthere chmjuvsstart    totper
# 1          0     41336.33 0.7794094
#From below, standard error on totper is  0.05137866

#Fleet was already there
effectFleetThere(1)[[1]]
# fleetthere chmjuvsstart    totper
# 1          1     7439.953 0.1911517
#From below, standard error on totper is 0.03707811


#From below, the standard error in differences is 0.06336049
#So p-value is 
(1 - pt((0.7794094 - 0.1911517) / 0.06336049, df = f3_df))*2

#2/3 of juveniles caught are caught by vessels whose fleet was not already there
(sum(heterodf$numjuv[heterodf$type=="fleet_already3" & heterodf$value==0], na.rm=T) / 
    sum(heterodf$numjuv[heterodf$type=="fleet_already3"], na.rm=T))

#Vessels without fleet there account for % of total effect
effectFleetThere(0)[[1]] %>% dplyr::select(chmjuvsstart) %>% as.matrix() %>% as.numeric() / 
  (
    effectFleetThere(0)[[1]] %>% dplyr::select(chmjuvsstart) %>% as.matrix() %>% as.numeric() + 
      effectFleetThere(1)[[1]] %>% dplyr::select(chmjuvsstart) %>% as.matrix() %>% as.numeric()
  )



#Calculate standard error on total percent effect for fleet there = 0
secondlist <- effectFleetThere(0)[[2]]

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

#So standard error on percent effect is 
(totperse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                              x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                              x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                              x31 + x32 + x33 + x34 + x35 + x36)*scaleconstant + chjuvsoutside) / 
                          juv0, mycoefs, mybigvcov, ses=T)) #0.05137866



#Create a larger coefdf and vcovdf so I can calculate delta se on difference in total percentage change
togcoefs <- mycoefs; togvcov <- mybigvcov
chjuvsoutside_0 <- chjuvsoutside; juv0_0 <- juv0



#Calculate standard error on total percent effect for fleet there = 1
secondlist <- effectFleetThere(1)[[2]]

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

#So standard error on percent effect is 
(totperse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                              x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                              x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                              x31 + x32 + x33 + x34 + x35 + x36)*scaleconstant + chjuvsoutside) / 
                          juv0, mycoefs, mybigvcov, ses=T)) #0.03707811

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
                          juv0, togcoefs, togvcov, ses=T)) #0.06336049


#Clean up
rm(totperse, mycoefs, mybigvcov, myjuv1, juv0, i, chjuvsoutside, 
   chmjuv_delta, mycoef, myvcov, scaleconstant, toteffect_juv, secondlist,
   tonscaught, tonscoef, chjuvsoutside_0, juv0_0)






##Now do effect of self already being there

#Self there
s3 <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(
      paste0(grep("_treatfrac",names(heterodf),value=T), ":value"),
      collapse="+"),
    " + clustnobs + clusttons + clustarea_km2 + 
    kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(heterodf),value=T),collapse="+"),
    "| bin:value + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = filter(heterodf, type=="self_already3"))

summary(s3)

#Save degrees of freedom 
s3_df <- s3$df

selftab <- summary(s3)[["coefficients"]]

selftab <- mutate(as.data.frame(selftab), bin = rownames(selftab))

selftab <- selftab[grep("treatfrac",selftab$bin),]

selftab$bin <- gsub("_treatfrac","",selftab$bin)

selftab$selfthere <- 0
selftab$selfthere[grep("value1",selftab$bin)] <- 1

selftab$bin <- gsub("value0","",selftab$bin)
selftab$bin <- gsub("value1","",selftab$bin)
selftab$bin <- gsub(":","",selftab$bin)

#Calculate percent effect
selftab <- mutate(selftab, perchange = (exp(Estimate)-1)*100, 
                   perchangese = as.numeric(NA))

#Calculate delta method standard errors for change in millions of juveniles caught
for(i in 1:nrow(selftab)){
  
  #Filter toteffect_juv to mybin
  mydf <- selftab[i,]
  
  mycoef <- mydf$Estimate
  myvcov <- mydf$`Cluster s.e.`^2
  
  #Delta se. Random variable being transformed is mycoef
  mypcse <- deltamethod(~ (exp(x1)-1)*100, mycoef, myvcov, ses=T)
  
  #Plug this value into selftab
  selftab$perchangese[i] <- mypcse
  
}

rm(mydf, mycoef, myvcov, i, mypcse)

##Make figure of percent change

#Confidence intervals
selftab <- mutate(selftab, pc_ub = perchange + perchangese*qnorm(.975),
                   pc_lb = perchange - perchangese*qnorm(.975))

#Rename lead bins to match rddf naming convention
selftab$bin <- gsub("lead1_","lead9hours_",selftab$bin) 

#Get tvar and bdist from rddf
selftab <- left_join(selftab, 
                      distinct(as.data.frame(rddf) %>% dplyr::select(-geometry), bin, tvar, bdist),
                      by = c("bin")
)


#Plot types for same bdist next to each other but slightly offset
selftab$bdist[selftab$selfthere==0] <- selftab$bdist[selftab$selfthere==0] - 1
selftab$bdist[selftab$selfthere==1] <- selftab$bdist[selftab$selfthere==1] + 1

selftab$selfthere <- as.factor(selftab$selfthere)
selftab$selfthere <- relevel(selftab$selfthere, ref = "0")




#Also make a single figure for paper
singlePlot <- function(mytvar){
  
  #Title
  if(mytvar==-1){
    tit <- "After announcement, before closure"
  } else if(mytvar==0){
    tit <- "Closure period"
  } else if(mytvar==1){
    tit <- "1 day after"
  } else{
    tit <- paste0(mytvar, " days after")
  }
  
  #Want consistent y range across given variable
  myymin <- min(selftab$pc_lb)
  myymax <- max(selftab$pc_ub)
  
  #If not first plot of panel, no legend
  if((mytvar!= -1 & mytvar!= 2)){
    plot <- ggplot(data=filter(selftab, tvar==mytvar),aes(x=bdist,col=selfthere,shape=selfthere)) + 
      geom_hline(aes(yintercept=0)) + 
      geom_point(aes(y=perchange)) + 
      scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                         labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) +
      scale_y_continuous("",limits = c(myymin,myymax),
                         breaks = seq(from = -100, to=1200, by = 100),
                         labels = paste0(seq(from = -100, to=1200, by = 100), "%")) + 
      myThemeStuff + 
      ggtitle(tit) + 
      scale_color_manual("Vessel fished\nthere the\nday before",labels=c("No","Yes"),
                         values = c("dodgerblue4","orange1")) + 
      scale_shape_manual("Vessel fished\nthere the\nday before", labels=c("No","Yes"),
                         values = c(17,19)) + 
      geom_errorbar(data=filter(selftab, tvar==mytvar), aes(ymin=pc_lb,ymax=pc_ub),width=0) + 
      theme(legend.position='none')
  } else{
    
    plot <- ggplot(data=filter(selftab, tvar==mytvar),aes(x=bdist,col=selfthere,shape=selfthere)) + 
      geom_hline(aes(yintercept=0)) + 
      geom_point(aes(y=perchange)) + 
      scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                         labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) +
      scale_y_continuous("Percent change in juveniles caught",limits = c(myymin,myymax),
                         breaks = seq(from = -100, to=1200, by = 100),
                         labels = paste0(seq(from = -100, to=1200, by = 100), "%")) + 
      myThemeStuff + 
      ggtitle(tit) + 
      scale_color_manual("Vessel fished\nthere the\nday before",labels=c("No","Yes"),
                         values = c("dodgerblue4","orange1")) + 
      scale_shape_manual("Vessel fished\nthere the\nday before", labels=c("No","Yes"),
                         values = c(17,19)) + 
      geom_errorbar(data=filter(selftab, tvar==mytvar), aes(ymin=pc_lb,ymax=pc_ub),width=0) + 
      #Legend placement
      theme(legend.position = c(0.8,.85), 
            legend.margin = margin(0,0,0,0,unit="cm"), 
            legend.key.height = unit(.25, unit = "cm"),
            legend.key.width=unit(0,unit="cm"), 
            legend.key.size=unit(0,unit="cm"), 
            legend.text=element_text(size=9, family = "sans"),
            legend.title = element_text(size=9, family = "sans"))
  }
  
  
  
  return(plot)
}

#Make a 2x3 plot for paper with labeled panels
paperFig <- function(myvar){
  leadplot <- singlePlot(-1) + 
    labs(tag = "a") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  activeplot <- singlePlot(0) + 
    labs(tag = "b") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  lag1plot <- singlePlot(1) + 
    labs(tag = "c") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  lag2plot <- singlePlot(2) + 
    labs(tag = "d") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  lag3plot <- singlePlot(3) + 
    labs(tag = "e") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  lag4plot <- singlePlot(4) + 
    labs(tag = "f") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  
  tbt <- plot_grid(leadplot, activeplot, lag1plot, 
                   lag2plot, lag3plot, lag4plot, nrow=2, ncol=3, 
                   rel_widths = c(1.01,1,1))
  
  ggsave(tbt, file=paste0("Output/Figures/figure9.png"),
         w=7,h=(7/1.69)*2, units = "in", dpi=1200)
}

paperFig("hetero_selfthere")



##Then calculate total change % for each selfthere type

tonscaught <- felm(
  as.Formula(paste0(
    "asinhtons", "~ ", 
    "treatfrac:value",
    " + clustnobs + clusttons + clustarea_km2 + 
    kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(heterodf),value=T),collapse="+"),
    "| bin:value + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = filter(heterodf, type=='self_already3'))


tonscoef <- summary(tonscaught)[["coefficients"]]
tonscoef <- tonscoef[grep("treatfrac",rownames(tonscoef)),]
tonscoef <- as.data.frame(tonscoef) %>% 
  mutate(selfthere = gsub("treatfrac:value","",rownames(tonscoef)))


#Don't reallocate tons from 2017 second season or 2019 second season because they were both 
#shut down well before TAC was reached. 

#Function of selfthere type
effectSelfThere <- function(myvalue){
  
  #Tons caught in state of world I observe in seasons where TAC hit
  tons1 <- filter(heterodf, type=="self_already3" & value==myvalue) %>% 
    summarise(sum(tons)) %>% as.matrix() %>% as.numeric()
  
  #Scale it by tons in heterodf compared to total tons in data
  tons1 <- tons1 * 
    (sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"]) / 
       filter(heterodf, type=="self_already3") %>% 
       summarise(sum(tons)) %>% as.matrix() %>% as.numeric()
    )
  
  #Change in tons because of policy
  ctons <- tons1 - tons1/exp(tonscoef$Estimate[tonscoef$selfthere==myvalue]) 
  
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
  toteffect_juv <- filter(heterodf, type=="self_already3" & value==myvalue) %>%
    group_by(bin, tvar, bdist) %>% 
    summarise(juv1 = sum(nummjuv)) %>% ungroup()
  
  toteffect_juv <- left_join(toteffect_juv, 
                             filter(selftab, selfthere==myvalue) %>% 
                               dplyr::select(Estimate, bin, `Cluster s.e.`),
                             by = 'bin') %>% 
    mutate(juv0 = juv1/(exp(Estimate))) %>% 
    mutate(chmjuv = juv1 - juv0) %>%
    arrange(tvar, bdist) %>% ungroup()
  
  #Scale number of juv in fullbe by number of juv in heterodf
  scaleconstant <- (sum(fullbe$numjuv,na.rm=T) / sum(heterodf$numjuv[heterodf$type=="self_already3"], na.rm=T))
  
  changejuv <- scaleconstant*sum(toteffect_juv$chmjuv)
  
  #Now can calculate change in juvenile catch due to policy, accounting for reallocation
  chmjuvsstart <- changejuv + chjuvsoutside
  
  #How many juveniles are caught during my sample period in total?
  #Scale by percentage of juveniles caught in value bin
  juv1 <- (sum(fullbe$numjuv, na.rm=T) / 10^6) * 
    (sum(heterodf$numjuv[heterodf$type=="self_already3" & heterodf$value==myvalue], na.rm=T) / 
       sum(heterodf$numjuv[heterodf$type=="self_already3"], na.rm=T))
  
  #chmjuvsstart = juv(1) - juv(0)
  juv0 <- juv1 - chmjuvsstart
  
  #Change in juvenile catch as percentage
  totper <- chmjuvsstart / juv0
  
  out <- list(data.frame(fleetthere=myvalue, chmjuvsstart = chmjuvsstart, totper = totper),
              list(toteffect_juv, scaleconstant, chjuvsoutside, juv0)
  )  
  return(out)
  
}

#Vessel wasn't there
effectSelfThere(0)[[1]]
# fleetthere chmjuvsstart    totper
# 1          0     57740.61 0.8749544
#Below, SE on totper is 0.05513569


#Vessel was already there
effectSelfThere(1)[[1]]
# fleetthere chmjuvsstart      totper
# 1          1     116.5841 0.006905187
#Below, SE on totper is 0.03198715


#From below, difference in percent change is 0.06374263
#So the p-value on the difference is 
(1 - pt((0.8749544 - 0.006905187) / 0.06374263, df = s3_df))*2
#0


#% of juveniles caught are caught by vessels who were not already there
(sum(heterodf$numjuv[heterodf$type=="self_already3" & heterodf$value==0], na.rm=T) / 
    sum(heterodf$numjuv[heterodf$type=="self_already3"], na.rm=T))

#Vessels that weren't already there account for % of total effect
effectSelfThere(0)[[1]] %>% dplyr::select(chmjuvsstart) %>% as.matrix() %>% as.numeric() / 
  (
    effectSelfThere(0)[[1]] %>% dplyr::select(chmjuvsstart) %>% as.matrix() %>% as.numeric() + 
      effectSelfThere(1)[[1]] %>% dplyr::select(chmjuvsstart) %>% as.matrix() %>% as.numeric()
  )





#Calculate standard error on total percent effect for self there = 0
secondlist <- effectSelfThere(0)[[2]]

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

(s0se <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                              x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                              x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                              x31 + x32 + x33 + x34 + x35 + x36)*scaleconstant + chjuvsoutside) / 
                          juv0, mycoefs, mybigvcov, ses=T))
# 0.05513569


#Create a larger coefdf and vcovdf so I can calculate delta se on difference in total percentage change
togcoefs <- mycoefs; togvcov <- mybigvcov 
chjuvsoutside_0 <- chjuvsoutside; juv0_0 <- juv0




#Calculate standard error on total percent effect for self there = 1
secondlist <- effectSelfThere(1)[[2]]

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

(s1se <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                          x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                          x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                          x31 + x32 + x33 + x34 + x35 + x36)*scaleconstant + chjuvsoutside) / 
                      juv0, mycoefs, mybigvcov, ses=T))
#0.03198715



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
#0.06374263


sessionInfo()
