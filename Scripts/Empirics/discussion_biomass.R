#Digitize number of individuals in each length bin at start of first season 2017
#This is population distribution (in the water), not distribution of catch

rm(list=ls())
setwd("C:/Users/gabee/Documents/replication_closures")

library(dplyr); library(sf); library(ggplot2); library(lubridate)
library(readr); library(purrr); library(readxl); library(imager)

#Peru time
Sys.setenv(TZ='America/Lima')

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
                      plot.margin = unit(c(0,.1,0,.01),"in"),
                      plot.tag = element_text(family = "sans", size = 9)
)


#length-weight equation from IMARPE (2019)
lengthweight <- function(length){
  #length in cm; weight in g
  weight <- .0036*length^3.238
  return(weight)
}

#Load Figure 6 from IMARPE pdf as an image
lfd <- load.image("Data/EvaluacionHidroacusticaRecPelagicosCrucero170304_march2017.jpg")

lfd <- imrotate(lfd, 0.4)
#Rotate image slightly

plot(lfd)
abline(v = 550)

glfd <- grayscale(lfd)

glfd <- as.data.frame(glfd)

#density is 0 beyond 17.5 cm so can crop out lengths > 17.5 (x < 1900)
glfd <- filter(glfd, y < 1540 & y > 1000 & x > 550 & x < 1900)

#Reverse y-axis ordering so plot is right-side-up
glfd <- mutate(glfd, y=-y)

#I need the x-values corresponding to length bin tick marks
ggplot() + geom_raster(data=glfd,aes(x=x,y=y,fill=value)) + scale_fill_gradient(low='black',high='white') +
  myThemeStuff + 
  geom_hline(aes(yintercept=-1535), col = 'red') +
  geom_vline(aes(xintercept=645), col = 'red')

ggplot() + geom_raster(data=glfd,aes(x=x,y=y,fill=value)) + scale_fill_gradient(low='black',high='white') +
  myThemeStuff + 
  geom_vline(aes(xintercept=975), col = 'red')

filter(glfd, y == -1535 & x > 961)

#Only going to integrate from 6 cm to 17 cm since these are min and max lengths 
#975 is 6 cm. 1845 is 17 cm
majticks <- data.frame(x = c(seq(from = 975, to = 1858, by = (1858 - 975) / (17 - 6) ))) %>% 
  mutate(y = -1535, majtick=1, x = round(x))

#Check that I got the major ticks right.
ggplot() + geom_raster(data=glfd,aes(x=x,y=y,fill=value)) + scale_fill_gradient(low='black',high='white') + 
  geom_point(data=majticks, aes(x=x,y=y),col='red')

#Give majticks its corresponding length value
majticks <- dplyr::select(majticks, x) %>% mutate(length = 6:17)

#Also calculate minor ticks (6.5, 7.5 etc.)
minticks <- majticks %>% mutate(length = length + .5) %>%
  mutate(leadx = lead(x, 1)) %>% 
  mutate(halfdist = (x+leadx)/2) %>% 
  filter(length!=17.5) %>% 
  mutate(halfdist = round(halfdist))

#Plot minor ticks 
ggplot() + geom_raster(data=glfd,aes(x=x,y=y,fill=value)) + scale_fill_gradient(low='black',high='white') + 
  geom_point(data=mutate(minticks, y=-1535), aes(x=halfdist,y=y),col='red')

#Row bind minticks with majticks
majticks <- bind_rows(majticks, 
                      minticks %>% dplyr::select(halfdist, length) %>% rename(x = halfdist))

majticks <- arrange(majticks, x)

#highest and lowest y values 
bottomy <- -1531 #y value of 0 % of individuals
topy <- -1016 #y value corresponding to 14% of individuals

ggplot() + geom_raster(data=glfd,aes(x=x,y=y,fill=value)) + scale_fill_gradient(low='black',high='white') + 
  geom_hline(aes(yintercept=-1016), col = 'red')

#Function to find distribution line, applying over length values
#Also watch out for horizontal border line at top of figure
findprop <- function(lengths){
  
  #Apply over lengths
  out <- map_df(lengths, function(x){
    
    #if length >= 7 or <= 16.5, don't search over the x-axis for darkest point
    #Proportion is above the x-axis, but the blackest point might be the x-axis itself 
    #(x-axis might be darker than the point on the density curve)
    if(x >= 7 & x <= 16.5){
      
      #if length >= 15, don't look higher than half-way up y-axis because true density is never higher than this
      if(x >= 15){
        blackestval <- min(glfd$value[glfd$y>= (bottomy + 5) & glfd$y<=ceiling(topy - (topy-bottomy)/2) & 
                                      glfd$x==majticks$x[majticks$length==x]])
      }else{
        #Highest proportion is about .13, so never search above this value (topy - 10 instead of topy)
        blackestval <- min(glfd$value[glfd$y>= (bottomy + 5) & glfd$y<= (topy - 10) & 
                                        glfd$x==majticks$x[majticks$length==x]])
      }
    }else{
      
      #if length < 7 or > 16.5, can search over x-axis, but don't look higher than halfway up y-axis
      blackestval <- min(glfd$value[glfd$y>= bottomy & glfd$y<=ceiling(topy - (topy-bottomy)/2) & 
                                      glfd$x==majticks$x[majticks$length==x]])
      
    }
    
    #Find y-value this pixel value represents
    blackesty <- glfd$y[glfd$y>=bottomy & glfd$y<=topy & glfd$x==majticks$x[majticks$length==x] & 
                          glfd$value==blackestval]
    
    #Vertical line at 12 cm seems to be darker than density line, so code this one manually
    if(x == 12){
      blackesty <- -1430
    }
    
    #Calculate proportion blackesty represents relative to height of plot
    #Take mean in case there are multiple y-values with equally black pixel values
    prop <- mean(((blackesty-bottomy) / (topy-bottomy)))*.14 #.14 is maximum value
    
    out <- data.frame(length = x, prop=prop)
  })
  
  return(out)
}

#Calculate proportion of individuals at each half cm length interval
props <- findprop(majticks$length)

#proportions sum to almost exactly one
sum(props$prop) #1.000117

#Downweight each bin so that they sum to exactly one
props <- mutate(props, prop = prop / sum(props$prop))

#Check how my percentage juvenile compares to 72% displayed on IMARPE plot
sum(props$prop[props$length < 12]) #0.7241098
#Very accurate


