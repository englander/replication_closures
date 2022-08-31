rm(list=ls())

library(dplyr); library(ggplot2)
library(sf);library(lubridate); library(glmnet)
library(lfe); library(Formula); library(xtable)
library(purrr); library(cowplot); library(tidyr)
library(latex2exp)

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

#Created in make_actualclosure_regressioncontrol.R
load("Output/TempData/actualclosure_regressioncontrol.Rdata")

regdf <- as.data.frame(regdf) %>% dplyr::select(-geometry)

regdf <- arrange(regdf, tvar, bdist)

#Millions of juveniles
regdf <- mutate(regdf, nummjuv = numjuv/10^6)

regdf <- mutate(regdf, asinhnummjuv = asinh(nummjuv), 
               asinhtons = asinh(tons))

regdf$bin <- as.factor(regdf$bin)
regdf$bin <- relevel(regdf$bin, ref="active_in")

regdf$twoweek_cellid_2p <- as.factor(regdf$twoweek_cellid_2p)
regdf$twowk <- as.factor(regdf$twowk)
regdf$cellid_2p <- as.factor(regdf$cellid_2p)
regdf$startdate <- as.factor(regdf$startdate)

#Drop actual or potential closures that have NA for size distribution
regdf <- filter(regdf, !is.na(prop12hat))

#How many clusters are there
unique(regdf$twoweek_cellid_2p) %>% length()

#Given variable, interact it with bin indicators, giving
interVars <- function(var){
  
  mydf <- regdf
  names(mydf)[names(mydf)==var] <- "myvar"
  
  #Want to manually interact var with bin indicator so I can look at each bin's coefficient
  #relative to 0 (rather than relative to omitted category)
  bininds <- model.matrix(~bin,data=regdf) %>% as.data.frame()
  
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

regdf <- bind_cols(
  regdf,
  interVars("treatfrac")
)

#Control group is potential closure-treatment bins with treatfrac = 0
juvcatch <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(grep("_treatfrac",names(regdf),value=T),collapse="+"),
    " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(regdf),value=T),collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = filter(regdf, closuretype=="actual" | treatfrac==0))



#Number of treatment observations included in regression
filter(regdf, closuretype=="actual") %>% nrow() #14472
filter(regdf, treatfrac==0) %>% nrow() #25127
#Total obs in regression
juvcatch$N #39599

jvtab <- summary(juvcatch)[["coefficients"]]

jvtab <- mutate(as.data.frame(jvtab), bin = rownames(jvtab))

jvtab <- jvtab[grep("treatfrac",jvtab$bin),]

jvtab$bin <- gsub("_treatfrac","",jvtab$bin)

jvtab <- dplyr::select(jvtab, -`t value`, -`Pr(>|t|)`)

jvtab <- rename(jvtab, actualclosure_regressioncontrol = Estimate, se = `Cluster s.e.`)

#Confidence intervals
jvtab <- mutate(jvtab, actualclosure_regressioncontrol_ub = actualclosure_regressioncontrol + se*qnorm(.975), 
                actualclosure_regressioncontrol_lb = actualclosure_regressioncontrol - se*qnorm(.975))

#Separate tvar and bdist variables
jvtab$bdist <- gsub(".*_","",jvtab$bin)
jvtab$bdist <- as.numeric(jvtab$bdist)
jvtab$tvar <- gsub("_.*","",jvtab$bin)
jvtab$bdist[grep("_in",jvtab$bin)] <- 0
jvtab$tvar[jvtab$tvar=="lead9hours"] <- "lead"

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
  usedf <- jvtab
  names(usedf)[names(usedf)==myvar] <- "plotvar"
  
  #Rename lb and ub so can refer to directly
  names(usedf)[names(usedf)==paste0(myvar,"_lb")] <- "lb"
  names(usedf)[names(usedf)==paste0(myvar,"_ub")] <- "ub"
  
  #Want consistent y range across given variable
  myymin <- min(usedf$lb)
  myymax <- max(usedf$ub)
  
  if(mytvar=="lead"|mytvar=="lag2"){
  plot <- ggplot(data=filter(usedf, tvar==mytvar),aes(x=bdist)) + 
    geom_hline(aes(yintercept=0)) + 
    geom_point(aes(y=plotvar)) + 
    scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                       labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) +
    scale_y_continuous(ylab,limits = c(myymin,myymax),
                       breaks = scales::pretty_breaks(n=10)) + 
    myThemeStuff + 
    ggtitle(tit) + 
    geom_errorbar(data=filter(usedf, tvar==mytvar), aes(x=bdist,ymin=lb,ymax=ub),width=0)
  } else{
    plot <- ggplot(data=filter(usedf, tvar==mytvar),aes(x=bdist)) + 
      geom_hline(aes(yintercept=0)) + 
      geom_point(aes(y=plotvar)) + 
      scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                         labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) +
      scale_y_continuous("",limits = c(myymin,myymax),
                         breaks = scales::pretty_breaks(n=10)) + 
      myThemeStuff + 
      ggtitle(tit) + 
      geom_errorbar(data=filter(usedf, tvar==mytvar), aes(x=bdist,ymin=lb,ymax=ub),width=0)
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
  
  ggsave(tbt, file=paste0("Output/Figures/figureA12.pdf"),
         w=7,h=(7/1.69)*2, units = "in", dpi=1200)
}


paperFig("actualclosure_regressioncontrol",
         TeX("$\\beta_{st}$ coefficient and 95% confidence interval (Equation 1)"))


#Calculate level effect for text of paper
jvtab <- rename(jvtab, Estimate = actualclosure_regressioncontrol)

#Calculate juv1 inside closures
#then calculate juv0. Then scale up both by ratio of nummjuv in closures 
#to total numjuv.
toteffect_juv <- group_by(regdf, bin, tvar, bdist) %>% 
  summarise(juv1 = sum(nummjuv)) %>%
  left_join(dplyr::select(jvtab, -tvar, -bdist), by = 'bin') %>% 
  mutate(juv0 = juv1/(exp(Estimate))) %>% 
  mutate(chmjuv = juv1 - juv0) %>%
  arrange(tvar, bdist) %>% ungroup()

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Total effect, not accounting for reallocation
changejuv <- (sum(fullbe$numjuv,na.rm=T) / sum(regdf$numjuv, na.rm=T)) * sum(toteffect_juv$chmjuv)

toteffect_juv <- rename(toteffect_juv, juvcoef = Estimate, juvse = se)

#Calculate delta method standard errors for change in millions of juveniles caught
library(msm)

toteffect_juv <- mutate(toteffect_juv, chmjuvse = as.numeric(NA), 
                        chmjuv_scaled = chmjuv * (sum(fullbe$numjuv,na.rm=T) / sum(regdf$numjuv, na.rm=T)),
                        chmjuvse_scaled = as.numeric(NA), 
                        juv0se = as.numeric(NA))

scaleconstant <- (sum(fullbe$numjuv,na.rm=T) / sum(regdf$numjuv, na.rm=T))

for(mybin in toteffect_juv$bin){
  
  #Filter toteffect_juv to mybin
  mydf <- filter(toteffect_juv, bin==mybin)
  
  #mypj <- mydf$meanpj/100 %>% as.character() %>% as.numeric()
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




#Closures cannot increases tons caught because of TAC, so account for this reallocation
#by estimating how policy affects tons caught
tonscaught <- felm(
  as.Formula(paste0(
    "asinhtons", "~ ", 
    "treatfrac ",
    " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(regdf),value=T),collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data =regdf)

tonscoef <- summary(tonscaught)[["coefficients"]]["treatfrac","Estimate"]

#Don't reallocate tons from 2017 second season or 2019 second season because they were both 
#shut down well before TAC was reached. 

#Tons caught in state of world I observe in seasons where TAC hit
tons1 <- sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"])

#Change in tons because of policy
ctons <- tons1 - tons1/exp(tonscoef) 

#Average pj outside of treatment window
avgpjoutside <- filter(fullbe, lead_0==0 & lead_10==0 & lead_20==0 & lead_30==0 & lead_40==0 & lead_50==0 & 
                         active_0==0 & active_10==0 & active_20==0 & active_30==0 & active_40==0 & active_50==0 & 
                         lag1_0==0 & lag1_10==0 & lag1_20==0 & lag1_30==0 & lag1_40==0 & lag1_50==0 & 
                         lag2_0==0 & lag2_10==0 & lag2_20==0 & lag2_30==0 & lag2_40==0 & lag2_50==0 & 
                         lag3_0==0 & lag3_10==0 & lag3_20==0 & lag3_30==0 & lag3_40==0 & lag3_50==0 & 
                         lag4_0==0 & lag4_10==0 & lag4_20==0 & lag4_30==0 & lag4_40==0 & lag4_50==0 & 
                         !is.na(numindivids) & !is.na(bepjhat) & Temporada!="2017-II" & Temporada!="2019-II") %>%
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
                             !is.na(numindivids) & !is.na(avgweightg) & Temporada!="2017-II" & Temporada!="2019-II") %>%
  #Weight by tons
  mutate(weightweighted = avgweightg*numindivids) %>% 
  summarise(avgweightg = sum(weightweighted)/sum(numindivids)) %>% as.numeric()


#Decrease in individuals caught outside of treatment window in millions
#(converting tons to g cancels out conversion to millions)
chindividsoutside <- -ctons/avgweightoutside

chjuvsoutside <- chindividsoutside*avgpjoutside

#Now can calculate change in juvenile catch due to policy, accounting for reallocation
(chmjuvsstart <- changejuv + chjuvsoutside) #35538.04

#How many juveniles are caught during my sample period in total?
#F(1)*pj*individuals/VMS fishing obs
juv1 <- sum(fullbe$numjuv, na.rm=T) / 10^6 

#chmjuvsstart = juv(1) - juv(0)
juv0 <- juv1 - chmjuvsstart 

#Then increase in juvenile catch as a percentage is 
(chmjuvsstart / juv0) #0.3378291


#Calculate standard error on total change in juvenile catch and in total percentage change
mycoefs <- toteffect_juv$chmjuv_scaled
mybigvcov <- diag(toteffect_juv$chmjuvse_scaled^2)

#This includes reallocation, so this is what I want: 
(changebillionsse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                                      x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                                      x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                                      x31 + x32 + x33 + x34 + x35 + x36) + chjuvsoutside) / 
                                  1000, mycoefs, mybigvcov, ses=T))
#2.881216

#Now get SE on total percentage change
(totperse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                              x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                              x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                              x31 + x32 + x33 + x34 + x35 + x36) + chjuvsoutside) / 
                          juv0, mycoefs, mybigvcov, ses=T))
#0.0273892
