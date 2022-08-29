#Figure C1. A two by two plot of a) uncorrected length distribution,
#b) corrected length distribution, 
#c) distribution of reported and d) corrected percentage juvenile values
options(scipen=999)

rm(list=ls())

library(sf); library(dplyr); library(ggplot2); library(lubridate)
library(purrr); library(readxl); library(tidyr); library(cowplot)
library(xtable); library(scales)

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

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
                      plot.margin = unit(c(0.1,0,0.02,0),"in"),
                      plot.tag = element_text(family = "sans", size = 9),
                      plot.tag.position = c(.05, 1)
)

Sys.setenv(TZ='America/Lima')

#Plot size distribution of individual anchoveta caught during my study period
#Load BE data. Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

fullbe_corrected <- fullbe

#Total number of individuals caught
numcaught <- sum(fullbe_corrected$numindivids, na.rm=T)

lengthbins <- grep("hat",names(fullbe_corrected),value=T)
lengthbins <- grep("prop",lengthbins,value=T)

#Select these columns as well as tons so can weight distribution by tons
lengthdf <- dplyr::select(fullbe_corrected, lengthbins, numindivids)

#Drop NA values
lengthdf <- filter(lengthdf, !is.na(prop12hat) & !is.na(numindivids))

lengthdf <- gather(lengthdf, lengthbin, prop, -numindivids)

lengthdf$lengthbin <- gsub("prop","",lengthdf$lengthbin)
lengthdf$lengthbin <- gsub("hat","",lengthdf$lengthbin)
lengthdf$lengthbin <- gsub("p",".",lengthdf$lengthbin)
lengthdf$lengthbin <- as.numeric(lengthdf$lengthbin)

lengthdf <- mutate(lengthdf, numlength = numindivids*prop) %>% 
  group_by(lengthbin) %>% 
  summarise(prop = sum(numlength)/numcaught) %>% ungroup()

#Median caught fish is between 13 and 13.5 cm
sum(lengthdf$prop[lengthdf$lengthbin>=13])
sum(lengthdf$prop[lengthdf$lengthbin<=13])

#Median caught juvenile is between 11 and 11.5 cm
sum(lengthdf$prop[lengthdf$lengthbin>=11 & lengthdf$lengthbin<12]) / sum(lengthdf$prop[lengthdf$lengthbin<12])
sum(lengthdf$prop[lengthdf$lengthbin<=11 & lengthdf$lengthbin<12]) / sum(lengthdf$prop[lengthdf$lengthbin<12])

#Median caught adult is between 13.5 and 14 cm cm
sum(lengthdf$prop[lengthdf$lengthbin>=13.5 & lengthdf$lengthbin>=12]) / sum(lengthdf$prop[lengthdf$lengthbin>=12])
sum(lengthdf$prop[lengthdf$lengthbin<=13.5 & lengthdf$lengthbin>=12]) / sum(lengthdf$prop[lengthdf$lengthbin>=12])

#What % are juveniles?
filter(lengthdf, lengthbin < 12) %>% summarise(sum(prop)) #.183



#12 cm bin means individuals between 12 and 12.5 cm, 
#And geom_col takes x coordinate as right interval
#so shift bins so they reflect correct interval in plot
lengthdf <- mutate(lengthdf, lengthbin = lengthbin + .5)

lengthdf_corrected <- lengthdf


##Now created the same df, but for uncorrected length distribution
#Created in 2. impute_size_be.R
load("Output/Data/pbe_imp_uncorrected.Rdata")

#Total number of individuals caught
numcaught <- sum(fullbe$numindivids, na.rm=T)

lengthbins <- grep("hat",names(fullbe),value=T)
lengthbins <- grep("prop",lengthbins,value=T)

#Select these columns as well as tons so can weight distribution by tons
lengthdf <- dplyr::select(fullbe, lengthbins, numindivids)

#Drop NA values
lengthdf <- filter(lengthdf, !is.na(prop12hat) & !is.na(numindivids))

lengthdf <- gather(lengthdf, lengthbin, prop, -numindivids)

lengthdf$lengthbin <- gsub("prop","",lengthdf$lengthbin)
lengthdf$lengthbin <- gsub("hat","",lengthdf$lengthbin)
lengthdf$lengthbin <- gsub("p",".",lengthdf$lengthbin)
lengthdf$lengthbin <- as.numeric(lengthdf$lengthbin)

lengthdf <- mutate(lengthdf, numlength = numindivids*prop) %>% 
  group_by(lengthbin) %>% 
  summarise(prop = sum(numlength)/numcaught) %>% ungroup()

#Make the two plots have same y-axis scale
maxpropval <- max(
  max(lengthdf$prop, lengthdf_corrected$prop)
)


#b
(lengthdist_corrected <- ggplot() + 
  geom_col(data=lengthdf_corrected, aes(x=lengthbin,y=prop)) + 
  scale_x_continuous("Half-cm length interval (cm)",breaks=c(seq(from=3.25,to=17.25,by=2),18.75),
                     labels = c(seq(from=3,to=17,by=2),18.5)) + 
  scale_y_continuous("Mean proportion of individuals in length interval",
                     breaks=seq(from=0,to=.12,by=.02), 
                     limits = c(0, maxpropval)) + 
  myThemeStuff + 
  geom_vline(aes(xintercept=12.25),col='red') + 
  theme(axis.text = element_text(color = "black", size = 8, family="sans"),
        axis.title = element_text(color = "black", size = 8, family = "sans"),
        plot.margin = unit(c(.05, .01, .1, .1), 'in'),
        plot.tag.position = c(.08, .985)) + 
  ggtitle("Corrected length distribution") + 
  labs(tag = "b")
)


##Now make a similar plot but with uncorrected length distribution

#Label juveniles as less than 12 cm
labtext <- data.frame(x=8.5,y=.1,text="Juveniles are < 12 cm")

#Median caught fish is between 13 and 13.5 cm
sum(lengthdf$prop[lengthdf$lengthbin>=13])
sum(lengthdf$prop[lengthdf$lengthbin<=13])

#What % are juveniles?
filter(lengthdf, lengthbin < 12) %>% summarise(sum(prop)) #.103

#12 cm bin means individuals between 12 and 12.5 cm, 
#And geom_col takes x coordinate as right interval
#so shift bins so they reflect correct interval in plot
lengthdf <- mutate(lengthdf, lengthbin = lengthbin + .5)

(uncorrected_lengthdist <- ggplot() + 
  geom_col(data=lengthdf, aes(x=lengthbin,y=prop)) + 
  scale_x_continuous("Half-cm length interval (cm)",breaks=c(seq(from=3.25,to=17.25,by=2),18.75),
                     labels = c(seq(from=3,to=17,by=2),18.5)) + 
    scale_y_continuous("Mean proportion of individuals in length interval",
                       breaks=seq(from=0,to=.12,by=.02), 
                       limits = c(0, maxpropval)) +
    myThemeStuff + 
  geom_vline(aes(xintercept=12.25),col='red') + 
  geom_text(data=labtext,aes(x=x,y=y,label=text),col='red', size = 3) + 
  theme(axis.text = element_text(color = "black", size = 8, family="sans"),
        axis.title = element_text(color = "black", size = 8, family = "sans"),
        plot.margin = unit(c(.05, .1, .1, .01), 'in'),
        plot.tag.position = c(.08, .985)) + 
  ggtitle("Uncorrected length distribution") + 
  labs(tag = "a")
)



##Finally plot distribution of reported (c) and corrected percentage juvenile values (d)
#Load BE data. Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

(reported_pj <- ggplot(data = fullbe, aes(x = bepj)) + 
  geom_histogram(binwidth = 1, closed = "left") + 
    myThemeStuff +
    scale_x_continuous("",breaks = seq(from = 0, to = 100, by = 20),
                       labels=paste0(seq(from=0,to=100,by=20),"%"),
                       limits = c(-1,100)) +
    scale_y_continuous("Count", breaks = breaks_extended(n = 5), 
                       labels = label_comma(), 
                       limits = c(0, 150000)) + 
    theme(axis.text = element_text(color = "black", size = 8, family="sans"),
          axis.title = element_text(color = "black", size = 8, family = "sans"),
          plot.margin = unit(c(.1, .1, .05, .01), 'in'),
          plot.tag.position = c(.08, .98)) + 
    ggtitle("Percentage juvenile reported to regulator") + 
    labs(tag = "c")
  
)

(corrected_pj <- ggplot(data = fullbe, aes(x = bepjhat)) + 
    geom_histogram(binwidth = 1, closed = "left") + 
    myThemeStuff +
    scale_x_continuous("",breaks = seq(from = 0, to = 100, by = 20),
                       labels=paste0(seq(from=0,to=100,by=20),"%"),
                       limits = c(-1,100)) +
    scale_y_continuous("Count", breaks = breaks_extended(n = 5), 
                       labels = label_comma(),
                       limits = c(0, 150000)) + 
    theme(axis.text = element_text(color = "black", size = 8, family="sans"),
          axis.title = element_text(color = "black", size = 8, family = "sans"),
          plot.margin = unit(c(.1, .01, .05, .1), 'in'),
          plot.tag.position = c(.08, .98)) + 
    ggtitle("Corrected percentage juvenile") + 
    labs(tag = "d")
  
)



tbt <- plot_grid(uncorrected_lengthdist, lengthdist_corrected, 
                 reported_pj, corrected_pj, 
                 nrow = 2, ncol=2, 
                 rel_heights = c(1, 1))

ggsave(tbt, file = "Output/Figures/figureC1.png",
       width=7,height=(7/1.69) * 4/3,units='in',dpi=1200)
