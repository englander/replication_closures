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

##Make Figure A2 first
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
finaldf <- distinct(rddf, bin, tvar, bdist) %>% 
  left_join(jvtab, by = 'bin') %>% 
  arrange(tvar, bdist) %>% ungroup()

finaldf <- dplyr::select(finaldf, -`t value`, -`Pr(>|t|)`)

finaldf <- rename(finaldf, juvcoef = Estimate, juvse = `Cluster s.e.`)

#95% confidence interval
finaldf <- mutate(finaldf, juvcoef_ub = juvcoef + juvse*qnorm(.975),
                  juvcoef_lb = juvcoef - juvse*qnorm(.975))

finaldf$tvar[finaldf$tvar==-1] <- "lead"
finaldf$tvar[finaldf$tvar==0] <- "active"
finaldf$tvar[finaldf$tvar%in%c(1,2,3,4)] <- paste0("lag",finaldf$tvar[finaldf$tvar%in%c(1,2,3,4)])

#Function of yvariable, one event day, dependent variable, figure number, and y-axis ticks 
singlePlot <- function(myvar, mytvar, ylab, fignum, yticks){
  
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
                       labels = label_number(accuracy=yticks)) + 
    myThemeStuff + 
    ggtitle(tit) + 
    geom_errorbar(data=filter(usedf, tvar==mytvar), aes(x=bdist,ymin=lb,ymax=ub),width=0)

  #Add ylabel for leftmost plots
  if(mytvar=="lead"|mytvar=="lag2"){
    plot <- plot + scale_y_continuous(ylab,limits = c(myymin,myymax),
                                      breaks = scales::pretty_breaks(n=10))
  }
  
  return(plot)
}

#Make a 2x3 plot for paper with labeled panels
paperFig <- function(myvar, ylab, fignum, yticks){
  leadplot <- singlePlot(myvar, "lead", ylab, fignum, yticks) + 
    labs(tag = "a") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  activeplot <- singlePlot(myvar, "active", ylab, fignum, yticks)+ 
    labs(tag = "b") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  lag1plot <- singlePlot(myvar, "lag1", ylab, fignum, yticks)+ 
    labs(tag = "c") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in"))
  lag2plot <- singlePlot(myvar, "lag2", ylab, fignum, yticks)+ 
    labs(tag = "d") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  lag3plot <- singlePlot(myvar, "lag3", ylab, fignum, yticks)+ 
    labs(tag = "e") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  lag4plot <- singlePlot(myvar, "lag4", ylab, fignum, yticks)+ 
    labs(tag = "f") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.15,0.04,0,0),"in"))
  
  tbt <- plot_grid(leadplot, activeplot, lag1plot, 
                   lag2plot, lag3plot, lag4plot, nrow=2, ncol=3, 
                   rel_widths = c(1.01,1,1))
  
  ggsave(tbt, file=paste0("Output/Figures/figureA", fignum, ".png"),
         w=7,h=(7/1.69)*2, units = "in", dpi=1200)
}

paperFig("juvcoef", TeX("$\\beta_{st}$ coefficient and 95% confidence interval (Equation 1)"), 
         2, .5)


##Now make Figure A3
rddf <- mutate(rddf, logjuv = if_else(nummjuv > 0, log(nummjuv), as.numeric(NA), as.numeric(NA)))

logreg <- felm(
  as.Formula(paste0(
    "logjuv", "~ ", 
    paste0(grep("_treatfrac",names(rddf),value=T),collapse="+"),
    " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(rddf),value=T),collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data =rddf)

logtab <- summary(logreg)[["coefficients"]]

logtab <- mutate(as.data.frame(logtab), bin = rownames(logtab))

logtab <- logtab[grep("treatfrac",logtab$bin),]

logtab$bin <- gsub("_treatfrac","",logtab$bin)

#Keep variables of interest and create confidence intervals
logtab <- dplyr::select(logtab, -`t value`, -`Pr(>|t|)`) 

names(logtab) <- c("logjuv", "logjuvse", "bin")

logtab <- mutate(logtab, logjuv_ub = logjuv + logjuvse*qnorm(.975), 
                 logjuv_lb = logjuv - logjuvse*qnorm(.975))

#Join onto finaldf
finaldf <- left_join(finaldf, logtab, by = 'bin')

paperFig("logjuv","Log points change in juvenile catch, conditional on positive juvenile catch",
         3, .5)


##Now make Figure A4
rddf <- mutate(rddf, posjuv = if_else(numjuv > 0, 1, 0))

indreg <- felm(
  as.Formula(paste0(
    "posjuv", "~ ", 
    paste0(grep("_treatfrac",names(rddf),value=T),collapse="+"),
    " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(rddf),value=T),collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data =rddf)

indtab <- summary(indreg)[["coefficients"]]

indtab <- mutate(as.data.frame(indtab), bin = rownames(indtab))

indtab <- indtab[grep("treatfrac",indtab$bin),]

indtab$bin <- gsub("_treatfrac","",indtab$bin)

#Keep variables of interest and create confidence intervals
indtab <- dplyr::select(indtab, -`t value`, -`Pr(>|t|)`) 

names(indtab) <- c("indjuv", "indjuvse", "bin")

indtab <- mutate(indtab, indjuv_ub = indjuv + indjuvse*qnorm(.975), 
                 indjuv_lb = indjuv - indjuvse*qnorm(.975))

#Join onto finaldf
finaldf <- left_join(finaldf, indtab, by = 'bin')

paperFig("indjuv","Change in probability juvenile catch exceeds 0", 4, .1)


##Now make Figure A5
rddf <- mutate(rddf, postreat = if_else(treatfrac > 0, 1, 0))

rddf <- bind_cols(
  rddf, 
  interVars("postreat")
)

postreatreg <- felm(
  as.Formula(paste0(
    "posjuv", "~ ", 
    paste0(grep("_postreat",names(rddf),value=T),collapse="+"),
    " + clustnobs + clusttons + clustarea_km2 + kmtocoast + clusttonsperset + clusttonsperarea + ", 
    paste0(grep("prop",names(rddf),value=T),collapse="+"),
    "| bin + twowk:cellid_2p + startdate",
    " | 0 | twoweek_cellid_2p")),
  data =rddf)

postreattab <- summary(postreatreg)[["coefficients"]]

postreattab <- mutate(as.data.frame(postreattab), bin = rownames(postreattab))

postreattab <- postreattab[grep("postreat",postreattab$bin),]

postreattab$bin <- gsub("_postreat","",postreattab$bin)

#Keep variables of interest and create confidence intervals
postreattab <- dplyr::select(postreattab, -`t value`, -`Pr(>|t|)`) 

names(postreattab) <- c("postreatjuv", "postreatjuvse", "bin")

postreattab <- mutate(postreattab, postreatjuv_ub = postreatjuv + postreatjuvse*qnorm(.975), 
                      postreatjuv_lb = postreatjuv - postreatjuvse*qnorm(.975))

#Join onto finaldf
finaldf <- left_join(finaldf, postreattab, by = 'bin')

paperFig("postreatjuv",TeX("$\\beta_{st}$ coefficient and 95% confidence interval (Equation 1)"),
         5, .05)

sessionInfo()
