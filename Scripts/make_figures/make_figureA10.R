rm(list=ls())

library(dplyr); library(ggplot2); library(msm)
library(sf);library(lubridate); library(glmnet)
library(lfe); library(Formula); library(xtable)
library(purrr); library(cowplot); library(tidyr)
library(readxl); library(latex2exp)

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

#Load rddf created in 4. make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

rddf <- as.data.frame(rddf) %>% dplyr::select(-geometry)

rddf <- arrange(rddf, tvar, bdist)

#Calculate median startdate within season
rddf <- left_join(rddf, 
                  filter(rddf, bin=="active_in") %>% 
  group_by(season) %>% 
  summarise(medstart = median(start)),
  by="season")

#Make sure same rid is in same group 
rddf <- left_join(rddf, 
            filter(rddf, bin=="active_in") %>%
  mutate(secondhalf = if_else(start > medstart, 1, 0)) %>% 
  dplyr::select(secondhalf, rid),
  by='rid'
  )

#Create same indicator for juvenile catch data
#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

fullbe <- left_join(
  fullbe,
  filter(rddf, bin=="active_in") %>% 
    group_by(season) %>% 
    summarise(medstart = median(start)) %>% 
    mutate(Temporada = paste0(substr(season, 4,7),"-",
                              ifelse(substr(season,2,2)=="1","I","II"))),
  by='Temporada')

fullbe <- rename(fullbe, calatime = FechaInicioCala)

fullbe <- mutate(fullbe, secondhalf = if_else(calatime > medstart, 1, 0))

#Also create same indicator for landings data
#Full 2017 to 2019 landings data
land <- read_excel("Data/landings_2017to2019.xlsx")

land <- left_join(
  land,
  filter(rddf, bin=="active_in") %>% 
    group_by(season) %>% 
    summarise(medstart = median(start)) %>% 
    mutate(Temporada = paste0(substr(season, 4,7),"-",
                              ifelse(substr(season,2,2)=="1","I","II"))),
  by='Temporada')

land <- mutate(land, secondhalf = if_else(FechaHoraInicioDescarga > medstart, 1, 0))

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

rddf$secondhalf <- as.factor(rddf$secondhalf)
rddf$secondhalf <- relevel(rddf$secondhalf, ref='0')


juvcatch <- felm(
  as.Formula(paste0(
    "asinhnummjuv", "~ ", 
    paste0(
      paste0(grep("_treatfrac",names(rddf),value=T), ":secondhalf"),
      collapse="+"),
    " + clustnobs + clusttons + clustarea_km2 + 
    kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(rddf),value=T),collapse="+"),
    "| bin:secondhalf + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = rddf)

jvtab <- summary(juvcatch)[["coefficients"]]

jvtab <- mutate(as.data.frame(jvtab), bin = rownames(jvtab))

jvtab <- jvtab[grep("treatfrac",jvtab$bin),]

jvtab$secondhalf <- 0
jvtab$secondhalf[grep("secondhalf1",jvtab$bin)] <- 1

jvtab$bin <- gsub("secondhalf0","",jvtab$bin)
jvtab$bin <- gsub("secondhalf1","",jvtab$bin)
jvtab$bin <- gsub(":","",jvtab$bin)
jvtab$bin <- gsub("_treatfrac","",jvtab$bin)

#Calculate percent effect
jvtab <- mutate(jvtab, perchange = (exp(Estimate)-1)*100, 
                perchangese = as.numeric(NA))

#Calculate delta method standard errors for change in millions of juveniles caught
for(i in 1:nrow(jvtab)){
  
  #Filter toteffect_juv to mybin
  mydf <- jvtab[i,]
  
  mycoef <- mydf$Estimate
  myvcov <- mydf$`Cluster s.e.`^2
  
  #Delta se. Random variable being transformed is mycoef
  mypcse <- deltamethod(~ (exp(x1)-1)*100, mycoef, myvcov, ses=T)
  
  #Plug this value into jvtab
  jvtab$perchangese[i] <- mypcse
  
}

rm(mydf, mycoef, myvcov, i, mypcse)


##Then calculate total change (in levels and %) for each day type 
tonscaught <- felm(
  as.Formula(paste0(
    "asinhtons", "~ ", 
    "treatfrac:secondhalf",
    " + clustnobs + clusttons + clustarea_km2 + 
    kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(rddf),value=T),collapse="+"),
    "| bin:secondhalf + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data = rddf)


tonscoef <- summary(tonscaught)[["coefficients"]]
tonscoef <- tonscoef[grep("treatfrac",rownames(tonscoef)),]
tonscoef <- as.data.frame(tonscoef) %>% 
  mutate(secondhalf = gsub("treatfrac:secondhalf","",rownames(tonscoef)))



#Don't reallocate tons from 2017 second season or 2019 second season because they were both 
#shut down well before TAC was reached. 

#Tons caught in state of world I observe in seasons where TAC hit

#Function of secondhalf
effectHalf <- function(mysecondhalf){
  
  #Tons caught in state of world I observe in seasons where TAC hit
  tons1 <- filter(rddf, secondhalf==mysecondhalf) %>% 
    summarise(sum(tons)) %>% as.matrix() %>% as.numeric()
  
  #Scale it by tons in rddf compared to total tons in data
  tons1 <- tons1 * 
    (sum(fullbe$betons[fullbe$Temporada!="2017-II" & fullbe$Temporada!="2019-II"]) / 
       rddf %>% summarise(sum(tons)) %>% as.matrix() %>% as.numeric()
    )
  
  #Change in tons because of policy
  ctons <- tons1 - tons1/exp(tonscoef$Estimate[tonscoef$secondhalf==mysecondhalf]) 
  
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
  toteffect_juv <- filter(rddf, secondhalf==mysecondhalf) %>%
    group_by(bin, tvar, bdist) %>% 
    summarise(juv1 = sum(nummjuv)) %>% ungroup()
  
  toteffect_juv <- left_join(toteffect_juv, 
                             filter(jvtab, secondhalf==mysecondhalf) %>% 
                               dplyr::select(Estimate, bin, `Cluster s.e.`),
                             by = 'bin') %>% 
    mutate(juv0 = juv1/(exp(Estimate))) %>% 
    mutate(chmjuv = juv1 - juv0) %>%
    arrange(tvar, bdist) %>% ungroup()
  
  #Scale number of juv in fullbe by number of juv in rddf
  scaleconstant <- (sum(fullbe$numjuv,na.rm=T) / sum(rddf$numjuv, na.rm=T))
  
  changejuv <- sum(toteffect_juv$chmjuv) * scaleconstant
  
  #Now can calculate change in juvenile catch due to policy, accounting for reallocation
  chmjuvsstart <- changejuv + chjuvsoutside
  
  #How many juveniles are caught during my sample period in total?
  #Scale by percentage of juveniles caught in secondhalf bin
  juv1 <- (sum(fullbe$numjuv, na.rm=T) / 10^6) * 
    (sum(rddf$numjuv[rddf$secondhalf==mysecondhalf], na.rm=T) / 
       sum(rddf$numjuv, na.rm=T))
  
  #chmjuvsstart = juv(1) - juv(0)
  juv0 <- juv1 - chmjuvsstart
  
  #Change in juvenile catch as percentage
  totper <- chmjuvsstart / juv0
  
  out <- list(data.frame(secondhalf=mysecondhalf, chmjuvsstart = chmjuvsstart, totper = totper),
              list(toteffect_juv, scaleconstant, chjuvsoutside, juv0)
  )
  
  
  return(out)
  
}

effectHalf(0)[[1]]
# secondhalf chmjuvsstart    totper
# 1          0      25744.1 0.4201693

effectHalf(1)[[1]]
# secondhalf chmjuvsstart   totper
# 1          1     15315.86 0.398823


#63% of effect occurs in first half of season
25744.1/ (25744.1 + 15315.86)

#What % of tons caught in first half of season and what % caught in second half?
sum(land$TmDescargada[land$secondhalf==0]) / sum(land$TmDescargada)
#0.5848285

#So not a big difference; 63% of effect occurs in first half of season compared to 58% of tons
#landed are in first half of season

#Get tvar and bdist from rddf
jvtab <- left_join(
  jvtab, 
  as.data.frame(rddf) %>% dplyr::select(tvar, bin, bdist) %>% distinct(),
  "bin"
)

jvtab$tvar[jvtab$tvar==-1] <- "lead"
jvtab$tvar[jvtab$tvar==0] <- "active"
jvtab$tvar[jvtab$tvar%in%c(1,2,3,4)] <- paste0("lag",jvtab$tvar[jvtab$tvar%in%c(1,2,3,4)])

#Plot length types for same bdist next to each other but slightly offset
jvtab$bdist[jvtab$secondhalf==0] <- jvtab$bdist[jvtab$secondhalf==0] - 1
jvtab$bdist[jvtab$secondhalf==1] <- jvtab$bdist[jvtab$secondhalf==1] + 1

jvtab$secondhalf <- as.factor(jvtab$secondhalf)
jvtab$secondhalf <- relevel(jvtab$secondhalf, ref = "0")

#Keep variables of interest and create confidence intervals
jvtab <- dplyr::select(jvtab, -`t value`, -`Pr(>|t|)`) 

names(jvtab)[1:2] <- c("juvcoefsecondhalf", "juvcoefsecondhalfse")

jvtab <- mutate(jvtab, juvcoefsecondhalf_ub = juvcoefsecondhalf + juvcoefsecondhalfse*qnorm(.975), 
                  juvcoefsecondhalf_lb = juvcoefsecondhalf - juvcoefsecondhalfse*qnorm(.975))





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
  
  
  plot <- ggplot(data=filter(usedf, tvar==mytvar),aes(x=bdist,col=secondhalf)) + 
    geom_hline(aes(yintercept=0)) + 
    geom_point(aes(y=plotvar)) + 
    scale_x_continuous("",breaks=seq(from=0,to=50,by=10),
                       labels = c("Inside","10 km","20 km","30 km","40 km","50 km")) +
    scale_y_continuous("",limits = c(myymin,myymax),
                       breaks = scales::pretty_breaks(n=10),
                       labels = scales::comma) + 
    scale_color_manual("Season",labels=c("First half","Second half"),
                       values = c("skyblue","darkorange1")) + 
    myThemeStuff + 
    ggtitle(tit) + 
    geom_errorbar(data=filter(usedf, tvar==mytvar), aes(x=bdist,ymin=lb,ymax=ub),width=0)
  
  #Add legend and y-axis label for left plots
  if(mytvar=="lead"|mytvar=="lag2"){
    plot <- plot + scale_y_continuous(ylab,limits = c(myymin,myymax),
                                      breaks = scales::pretty_breaks(n=10)) + 
      theme(legend.position = c(0.8,.85), 
            legend.margin = margin(0,0,0,0,unit="cm"), 
            legend.key.height = unit(.25, unit = "cm"),
            legend.key.width=unit(0,unit="cm"), 
            legend.key.size=unit(0,unit="cm"))
  }
  if(mytvar!="lead" & mytvar != "lag2"){
    plot <- plot + scale_y_continuous("",limits = c(myymin,myymax),
                                      breaks = scales::pretty_breaks(n=10)) + 
      theme(legend.position='none')
  }
  return(plot)
}



#Also make a 2x3 plot for paper with labeled panels
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
  
  ggsave(tbt, file=paste0("Output/Figures/figureA10.pdf"),
         w=7,h=(7/1.69)*2, units = "in", dpi=1200)
}


paperFig("juvcoefsecondhalf",TeX("$\\beta_{st}$ coefficient and 95% confidence interval (Equation 1)"))