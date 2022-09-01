rm(list=ls())

library(dplyr); library(ggplot2)
library(sf);library(lubridate); library(glmnet)
library(lfe); library(Formula); library(xtable)
library(purrr); library(cowplot); library(tidyr)
library(latex2exp)

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

options(scipen=999)

#Peru time
Sys.setenv(TZ='America/Lima')

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



#Created in make_data_figA13.R
load("Output/Data/data_figA13.Rdata")

#Calculate mean outcome in each bin
meanout <- mutate(synthdf, synthdif = asinhnummjuv_treat - asinhnummjuv_control,
                  asinhtonsdif = asinhtons_treat - asinhtons_control, 
                  nummjuvdif = nummjuv_treat - nummjuv_control) %>%
  group_by(bdist, tvar) %>% 
  summarise(synthdif = mean(synthdif), 
            asinhtonsdif = mean(asinhtonsdif), 
            nummjuvdif = mean(nummjuvdif)) %>% ungroup()

#Function that plots mean outcome for treated and control given tvar
#Function of one event day and dependent variable
singlePlot <- function(myvar, mytvar, ylab){
  
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
  
  #Rename desired variable
  usedf <- meanout
  names(usedf)[names(usedf)==myvar] <- "plotvar"
  
  #Want consistent y range across given variable
  myymin <- min(usedf$plotvar)
  myymax <- max(usedf$plotvar)
  
  if(mytvar==-1|mytvar==2){
    plot <- ggplot(data=filter(usedf, tvar==mytvar),aes(x=bdist)) + 
      geom_hline(aes(yintercept=0)) + 
      geom_point(aes(y=plotvar)) + 
      scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                         labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) +
      scale_y_continuous(ylab,limits = c(myymin,myymax),
                         breaks = scales::pretty_breaks(n=10)) + 
      myThemeStuff + 
      ggtitle(tit)
  } else{
    plot <- ggplot(data=filter(usedf, tvar==mytvar),aes(x=bdist)) + 
      geom_hline(aes(yintercept=0)) + 
      geom_point(aes(y=plotvar)) + 
      scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                         labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) +
      scale_y_continuous("",limits = c(myymin,myymax),
                         breaks = scales::pretty_breaks(n=10)) + 
      myThemeStuff + 
      ggtitle(tit)
  }
  
  return(plot)
}


#Make a 2x3 plot for paper with labeled panels
paperFig <- function(myvar, ylab){
  leadplot <- singlePlot(myvar, -1, ylab) + 
    labs(tag = "a") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  activeplot <- singlePlot(myvar, 0, ylab)+ 
    labs(tag = "b") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  lag1plot <- singlePlot(myvar, 1, ylab)+ 
    labs(tag = "c") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  lag2plot <- singlePlot(myvar, 2, ylab)+ 
    labs(tag = "d") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  lag3plot <- singlePlot(myvar, 3, ylab)+ 
    labs(tag = "e") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  lag4plot <- singlePlot(myvar, 4, ylab)+ 
    labs(tag = "f") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  
  tbt <- plot_grid(leadplot, activeplot, lag1plot, 
                   lag2plot, lag3plot, lag4plot, nrow=2, ncol=3, 
                   rel_widths = c(1.01,1,1))
  
  ggsave(tbt, file=paste0("Output/Figures/figureA13.pdf"),
         w=7,h=(7/1.69)*2, units = "in", dpi=1200)
}

paperFig("synthdif",
         "Treatment juvenile catch minus control juvenile catch")

##Calculate total treatment effect
#Calculate juv1 then calculate juv0. Then scale up both by ratio of nummjuv in data 
#to total numjuv.
toteffect_juv <- group_by(synthdf, tvar, bdist) %>% 
  summarise(juv1_treat = sum(nummjuv_treat), juv1_control=sum(nummjuv_control)) %>%
  mutate(juv1 = juv1_treat + juv1_control) %>% 
  dplyr::select(-juv1_treat, -juv1_control) %>% ungroup() %>% 
  left_join(dplyr::select(meanout, tvar, bdist, synthdif), by = c('tvar','bdist')) %>% 
  mutate(juv0 = juv1/(exp(synthdif))) %>% 
  mutate(chmjuv = juv1 - juv0) %>%
  arrange(tvar, bdist) %>% ungroup()

#Re-scale to avoid double-counting
#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Calculate effect of policy, not accounting for TAC reallocation
#All juveniles caught during data divided by juveniles caught during regression data
changejuv <- sum(toteffect_juv$chmjuv[toteffect_juv$tvar >= -1]) * (sum(fullbe$numjuv,na.rm=T)/10^6) / (toteffect_juv$juv1[toteffect_juv$tvar >= -1] %>% sum())

#Account for reallocation in tons caught
tonscoef <- mean(meanout$asinhtonsdif[meanout$tvar >= -1])

#Tons caught in seasons with binding TAC
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
    summarise(perjuv = sum(pjweighted)/sum(numindivids)) %>% as.numeric() / 100)


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
(chmjuvsstart <- changejuv + chjuvsoutside) #35401.9

#How many juveniles are caught during my sample period in total?
#F(1)*pj*individuals/VMS fishing obs
juv1 <- sum(fullbe$numjuv, na.rm=T) / 10^6 

#chmjuvsstart = juv(1) - juv(0)
juv0 <- juv1 - chmjuvsstart 

#Then increase in juvenile catch as a percentage is 
chmjuvsstart / juv0 #0.3360998