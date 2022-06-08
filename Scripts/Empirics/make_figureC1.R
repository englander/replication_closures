options(scipen=999)

rm(list=ls())

library(sf); library(dplyr); library(ggplot2); library(lubridate)
library(purrr); library(readxl); library(tidyr); library(cowplot)
library(xtable)

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


#Label juveniles as less than 12 cm
labtext <- data.frame(x=9,y=.1,text="Juveniles are < 12 cm")

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

lengthdist <- ggplot() + 
  geom_col(data=lengthdf, aes(x=lengthbin,y=prop)) + 
  scale_x_continuous("Half-cm length interval (cm)",breaks=c(seq(from=3.25,to=17.25,by=2),18.75),
                     labels = c(seq(from=3,to=17,by=2),18.5)) + 
  scale_y_continuous("Average proportion of individuals in length interval",breaks=seq(from=0,to=.12,by=.02)) + 
  myThemeStuff + 
  geom_vline(aes(xintercept=12.25),col='red') + 
  geom_text(data=labtext,aes(x=x,y=y,label=text),col='red') + 
  theme(axis.text = element_text(color = "black", size = 8, family="sans"),
        axis.title = element_text(color = "black", size = 9, family = "sans"))

ggsave("Output/Figures/figureC1.png",
       lengthdist, width=7,height=7/1.69,units='in',dpi=1200)

sessionInfo()
