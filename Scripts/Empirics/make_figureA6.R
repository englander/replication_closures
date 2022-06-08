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
                      plot.title = element_text(hjust = 0.5, size = 7), 
                      legend.text=element_text(size=6.5, family = "sans"),
                      legend.title = element_text(size=6.5, family = "sans"),
                      plot.margin = unit(c(0,0.04,0,0),"in"),
                      plot.tag = element_text(family = "sans", size = 9)
)


#Peru time
Sys.setenv(TZ='America/Lima')

#Load rddf created in 4. make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

rddf <- as.data.frame(rddf) %>% dplyr::select(-geometry)

rddf <- arrange(rddf, tvar, bdist)

#Tons caught per set
rddf <- mutate(rddf, tonsperset = tons/nobs) %>%
  #Cluster tons/set and tons/area
  mutate(clusttonsperset = clusttons/clustnobs, clusttonsperarea = clusttons/clustarea_km2) %>%
  #clust millions of juveniles, individuals, and adults
  mutate(clustmjuv = clustnumjuv/10^6, clustmindivids = clustnumindivids/10^6,
         clustmadults = clustnumadults/10^6)

rddf <- mutate(rddf, asinhtonsjuv = asinh(tonsjuv), 
               asinhtons = asinh(tons))

rddf$bin <- as.factor(rddf$bin)
rddf$bin <- relevel(rddf$bin, ref="active_in")

rddf$twoweek_cellid_2p <- as.factor(rddf$twoweek_cellid_2p)
rddf$twowk <- as.factor(rddf$twowk)
rddf$cellid_2p <- as.factor(rddf$cellid_2p)

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

#Realized juvenile catch
juvcatch <- felm(
  as.Formula(paste0(
    "asinhtonsjuv", "~ ", 
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

#Plot treatment coefficients
jvtab <- dplyr::select(jvtab, -`t value`, -`Pr(>|t|)`)

jvtab <- rename(jvtab, juvtons = Estimate, se = `Cluster s.e.`)

#Confidence intervals
jvtab <- mutate(jvtab, juvtons_ub = juvtons + se*qnorm(.975), 
                juvtons_lb = juvtons - se*qnorm(.975))

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
  
  ggsave(tbt, file=paste0("Output/Figures/figureA6.png"),
         w=7,h=(7/1.69)*2, units = "in", dpi=1200)
}

paperFig("juvtons",
         TeX("$\\beta_{st}$ coefficient and 95% confidence interval (Equation 1)"))

##Now calculate level effect for text of paper
jvtab <- rename(jvtab, Estimate = juvtons)

#Calculate juv1 in potential closures, 
#then calculate juv0. Then scale up both by ratio of nummjuv in potential closures 
#to total numjuv.
toteffect_juv <- group_by(rddf, bin, tvar, bdist) %>% 
  summarise(juv1 = sum(tonsjuv)) %>%
  left_join(dplyr::select(jvtab, -tvar, -bdist), by = 'bin') %>% 
  mutate(juv0 = juv1/(exp(Estimate))) %>% 
  mutate(chtonsjuv = juv1 - juv0) %>%
  arrange(tvar, bdist) %>% ungroup()

#Percent effect, not accounting for reallocation
sum(toteffect_juv$chtonsjuv) / sum(toteffect_juv$juv0) 

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Total effect, not accounting for reallocation
(sum(fullbe$tonsjuv,na.rm=T) / sum(rddf$tonsjuv, na.rm=T)) * sum(toteffect_juv$chtonsjuv) 

toteffect_juv <- rename(toteffect_juv, juvcoef = Estimate, juvse = se)

#Calculate delta method standard errors for change in millions of juveniles caught
library(msm)

toteffect_juv <- mutate(toteffect_juv, chtonsjuvse = as.numeric(NA), 
                        chtonsjuv_scaled = chtonsjuv * (sum(fullbe$tonsjuv,na.rm=T) / sum(rddf$tonsjuv, na.rm=T)),
                        chtonsjuvse_scaled = as.numeric(NA), 
                        juv0se = as.numeric(NA))

scaleconstant <- (sum(fullbe$tonsjuv,na.rm=T) / sum(rddf$tonsjuv, na.rm=T))

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
  
  #Delta se for chtonsjuvse. Random variable being transformed is myjuvcoef
  chtonsjuv_delta <- deltamethod( ~ myjuv1 - (myjuv1/(exp(x1))),
                                  myjuvcoef,
                                  myjuvvcov,
                                  ses=T
  )
  
  #Plug this value into toteffect_juv
  toteffect_juv$chtonsjuvse[toteffect_juv$bin==mybin] <- chtonsjuv_delta
  
  #Delta se for scaled chtonsjuvse. Random variable being transformed is myjuvcoef
  chtonsjuv_scaled_delta <- deltamethod( ~ (myjuv1 - (myjuv1/(exp(x1))))*scaleconstant,
                                         myjuvcoef,
                                         myjuvvcov,
                                         ses=T
  )
  
  #Plug this value into toteffect_juv
  toteffect_juv$chtonsjuvse_scaled[toteffect_juv$bin==mybin] <- chtonsjuv_scaled_delta
  
}



##Calculate total effect, not accounting for reallocation of tons caught. (in millions)
changejuv <- (sum(fullbe$tonsjuv,na.rm=T) / sum(rddf$tonsjuv, na.rm=T)) * sum(toteffect_juv$chtonsjuv)


#Closures cannot increases tons caught because of TAC, so account for this reallocation
#by estimating how policy affects tons caught
tonscaught <- felm(
  as.Formula(paste0(
    "asinhtons", "~ ", 
    "treatfrac ",
    " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(rddf),value=T),collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data =rddf)

tonscoef <- summary(tonscaught)[["coefficients"]]["treatfrac","Estimate"]

#Don't reallocate tons from 2017 second season or 2019 second season because they were both 
#shut down well before TAC was reached. 

#Tons caught in state of world I observe in seasons where TAC hit
tons1 <- sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"])

#Change in tons because of policy
ctons <- tons1 - tons1/exp(tonscoef)

#What % of these tons are juvenile tons?
avgjuvtonsoutside <- filter(fullbe, lead_0==0 & lead_10==0 & lead_20==0 & lead_30==0 & lead_40==0 & lead_50==0 & 
                              active_0==0 & active_10==0 & active_20==0 & active_30==0 & active_40==0 & active_50==0 & 
                              lag1_0==0 & lag1_10==0 & lag1_20==0 & lag1_30==0 & lag1_40==0 & lag1_50==0 & 
                              lag2_0==0 & lag2_10==0 & lag2_20==0 & lag2_30==0 & lag2_40==0 & lag2_50==0 & 
                              lag3_0==0 & lag3_10==0 & lag3_20==0 & lag3_30==0 & lag3_40==0 & lag3_50==0 & 
                              lag4_0==0 & lag4_10==0 & lag4_20==0 & lag4_30==0 & lag4_40==0 & lag4_50==0 & 
                              !is.na(numindivids) & !is.na(bepjhat) & Temporada!="2017-II" & Temporada!="2019-II") %>%
  #% of tons that are juveniles
  summarise(juvtons = sum(tonsjuv,na.rm=T), adulttons = sum(tonsadult,na.rm=T)) %>% 
  mutate(perjuvtons = juvtons / (juvtons + adulttons)) %>% 
  dplyr::select(perjuvtons) %>% as.matrix() %>% as.numeric()

#Decrease in tons of juvenile caught outside treatment window because of policy
chtonsjuvoutside <- -ctons*avgjuvtonsoutside

#Now can calculate change in juvenile catch due to policy, accounting for reallocation
(chtonsjuvsstart <- changejuv + chtonsjuvoutside) #331706.1

#How many tons of juveniles are caught during my sample period in total?
juv1 <- sum(fullbe$tonsjuv, na.rm=T) 

#chtonsjuvsstart = juv(1) - juv(0)
juv0 <- juv1 - chtonsjuvsstart 

#Then increase in juvenile catch as a percentage is 
chtonsjuvsstart / juv0 # 0.43847


#Calculate standard error on total change in juvenile catch and in total percentage change
mycoefs <- toteffect_juv$chtonsjuv_scaled
mybigvcov <- diag(toteffect_juv$chtonsjuvse_scaled^2)

#This includes reallocation, so this is what I want
(changetonsjuvse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                                      x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                                      x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                                      x31 + x32 + x33 + x34 + x35 + x36) + chtonsjuvoutside),
                                mycoefs, mybigvcov, ses=T))
# 70144.73

#t-stat
chtonsjuvsstart / changetonsjuvse


#Now get SE on total percentage change
(totperse <- deltamethod(~ ((x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
                               x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + 
                               x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 +x30 + 
                               x31 + x32 + x33 + x34 + x35 + x36) + chtonsjuvoutside) / 
                           juv0, mycoefs, mybigvcov, ses=T))
#0.09272171

sessionInfo()
