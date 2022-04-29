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


#Load Figure 6 from IMARPE pdf as an image
#Full name of pdf is Informe "Evaluacion Hidroacustica de Recursos Pelagicos" Crucero 1703-04
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

#Drop length bins that represent 0 percent of individuals
props <- filter(props, prop > 0)

#proportions sum to almost exactly one
sum(props$prop) #1.000117

#Downweight each bin so that they sum to exactly one
props <- mutate(props, prop = prop / sum(props$prop))

#Check how my percentage juvenile compares to 72% displayed on IMARPE plot
sum(props$prop[props$length < 12]) #0.7241098
#Very accurate


##Now calculate status quo harvest (harvest with closures policy)
#length-weight equation from IMARPE (2019)
lengthweight <- function(length){
  #length in cm; weight in g
  weight <- .0036*length^3.238
  return(weight)
}


#length-age from Salvatteci and Mendo (2005)
#Length is in cm and age is in years. 
lengthage <- function(length){
  
  #Using average values from Tabla 2
  age <- -.0193 - (1/.96)*log(1 - length/19.35)
  
  return(age)
}
 
#Survival equation from Salvatteci and Mendo (2005)
survival <- function(ntm1, timestep){
  
  #Number of survivors in next period
  nt <- ntm1*exp(-.8*timestep)
  
  return(nt)
}

#Length given age from Salvatteci and Mendo (2005)
agelength <- function(age){
  
  #Using average values from Tabla 2
  length <- 19.35*(1 - exp(-.96*(age + .0193)))
  
  return(length)
}


#Status quo start from this population distribution 
#Harvest what we observe in this data and simulate until equilibrium
props <- mutate(props, age_years = lengthage(length))

#Weight of fish in each length bin
props <- mutate(props, weight = lengthweight(length))

#Normalized biomass in each length interval
props <- mutate(props, biomass = weight*prop)

#Average weight of a fish in the population
sum(props$biomass)
#multiplying by number of individuals in population would give total biomass (non-normalized) in grams

#Total biomass is 7.78 million tons (Informe "Evaluacion Hidroacustica de Recursos Pelagicos" Crucero 1703-04)

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Catch about 4 million tons per year
sum(fullbe$betons) / (3*10^6)

#Normalized harvest per day is 
harvest_per_day <- ((sum(fullbe$betons) / (3*10^6) / 365) / 7.78) * sum(props$biomass)

#Juvenile harvest in weight
harvest_juv_per_day <- (sum(fullbe$tonsjuv, na.rm=T) / (sum(fullbe$tonsjuv, na.rm=T) + 
                                                     sum(fullbe$tonsadult, na.rm=T))) * 
  harvest_per_day

#Recruitment -- number of 7 cm individuals next period
#Perea et al. (2011) finds anchoveta 12 - 14 cm contribute 40% of eggs while 
#anchoveta >= 14 cm contribute 60% of eggs. 
#So recruitment = rc*(.4*N_12-14 + .6*N_>=14)
#Solve for recruitment constant (rc) that will allow me to calculate counterfactual recruitment
#(using different values of N_12-14 and N_>=14 but same alpha)
recruit_constant <- props$prop[props$length == 7] / (
  .4*sum(props$prop[props$length >= 12 & props$length < 14]) + 
    .6*sum(props$prop[props$length >= 14])
)

#Want recruitment to acrue each day
#Takes about 15 days for a 7 cm anchoveta to reach 7.5 cm
(props$age_years[props$length == 7.5] - props$age_years[props$length == 7]) * 365

#To get per day recruitment, divide recruit constant by this length of time
recruit_constant <- recruit_constant / ((props$age_years[props$length == 7.5] - 
                                          props$age_years[props$length == 7]) * 365)


#Set "new" values as starting values. will update these in simulation
sqdf <- mutate(props,newlength = length, newprop = prop, newage_years = age_years, 
               newweight = weight, newbiomass = biomass)

#Starting value for biomass difference is just biomass values themselves
biomassdif <- sqdf$biomass

#Recruitment accrues each day, but a new class of 7 cm anchoveta only enter the population 
#(the data frame) when the initial class of 7 cm anchoveta have reached 7.5 cm
prop7 <- 0

while(max(abs(biomassdif)) > 1e-5){
  
  #Starting levels of biomass
  startbiomass <- sqdf$newbiomass
  
  #juvenile and adult biomass
  juvbiomass <- sum(sqdf$newbiomass[sqdf$length < 12])
  adultbiomass <- sum(sqdf$newbiomass[sqdf$length >= 12])
  
  #Increase age by one day
  sqdf <- mutate(sqdf, newage_years = newage_years + 1/365) %>% 
    #one day of growth
    mutate(newlength = agelength(newage_years)) %>% 
    mutate(newweight = lengthweight(newlength)) %>%
    #One day of death
    mutate(newprop = survival(newprop, 1/365)) 
  
  
  #Recruitment accrues each day
  prop7 <- prop7 + recruit_constant*(
    .4*sum(sqdf$newprop[sqdf$newlength >= 12 & sqdf$newlength < 14]) + 
      .6*sum(sqdf$newprop[sqdf$newlength >= 14])
  )
  
  sqdf <- sqdf %>% 
    # #One day of harvest. New normalized biomass of fish in each length bin
    # mutate(newbiomass = if_else(newlength < 12, 
    #                             newweight*newprop - 
    #                               harvest_juv_per_day*(newbiomass / juvbiomass),
    #                             newweight*newprop - 
    #                               (harvest_per_day - harvest_juv_per_day)*(newbiomass / adultbiomass)
    # )) %>% 
    mutate(newbiomass = newweight*newprop) %>% 
    #New proportion of individuals
    mutate(newprop = newbiomass / newweight)
  
  #which rows have negative newprop? that means there are no more individuals of this length bin
  whichneg <- which(sqdf$newprop < 0)
  
  if(length(whichneg) > 0){
    
    #Drop these elements from sqdf 
    sqdf <- sqdf[-whichneg,]
    
    #Also drop them from startbiomass so when re-calculating biomassdif the two vectors have the same number of elements
    startbiomass <- startbiomass[-whichneg]
    
  }
  
  #which rows have age > 3? 3 is max age of anchoveta so drop any rows older than 3
  which3 <- which(sqdf$newage_years > 3)
  
  if(length(which3) > 0){
    
    #Drop these elements from sqdf 
    sqdf <- sqdf[-which3,]
    
    #Also drop them from startbiomass so when re-calculating biomassdif the two vectors have the same number of elements
    startbiomass <- startbiomass[-which3]
    
  }
  
  #Re-calculate biomassdif
  biomassdif <- startbiomass - sqdf$newbiomass
  
  #A new class of 7 cm anchoveta only enter the population 
  #(the data frame) when the initial class of 7 cm anchoveta have reached 7.5 cm
  if(min(sqdf$newlength) >= 7.5){
    
  sqdf <- bind_rows(
      data.frame(length = 7, prop = 0.01331884, age_years = 0.4484462, weight = 1.96214, 
                 biomass = 0.02613342, 
                 newlength = 7, newprop = prop7, newage_years = 0.4484462, newweight = 1.96214, 
                 newbiomass = 0.02613342),
      sqdf
    )
  
  #reset prop7
  prop7 <- 0
    
  }
  
  #Maximum length of anchoveta is 20.59 cm: Tabla 4 IMARPE (2019)
  sqdf <- filter(sqdf, newlength <= 20.59)
  
}


#I COULD SIMULATE TO EQUILIBRIUM WITHOUT HARVEST, THEN ADD HARVEST W.R.T
#EQUILIBRIUM BIOMASS RATHER THAN STARTING BIOMASS. PROBABLY I'M HARVESTING TOO MUCH.













#One day of growth
props <- mutate(props, newage_years = age_years + 1/365) %>% 
  mutate(newlength = agelength(newage_years)) 

#One day of death
props <- mutate(props, newprop = survival(prop, 1/365))

#New weight of fish in each length bin (after one day of growth and death)
props <- mutate(props, newweight = lengthweight(newlength))

#New normalized biomass of fish in each length bin
props <- mutate(props, newbiomass = newweight*newprop)

#Average weight of a fish in the population has increased after one day of growth and death
sum(props$newbiomass)

#Harvest fish in each length bin s.t. biomass = newbiomass
#biomass in each length interval constant because quotas set in tons, not number of individuals
#anchoveta >= 15.5 cm die faster than the grow (in weight terms), so set harvest = 0 for these three length bins
props <- mutate(props, harvest_weight = if_else(newbiomass - biomass >= 0, newbiomass - biomass, 0))

#Status quo normalized (numindivids = 1) equilibrium harvest
sq_harvest_weight <- sum(props$harvest_weight)

##Multiplying sq_harvest by number of individuals in population gives one day of harvest in grams
#if biomass is 9 trillion grams (9 million tons), number of individuals in population is 
sq_numindivids <- 9*10^12 / sum(props$biomass)

#harvest 32,000 tons per day
sq_harvest_weight * sq_numindivids / 10^6 #convert tons to g

#since there are two fishing seasons of about three months each year
#status quo harvest is 5.76 million tons; pretty close to what we observe in data
sq_harvest_weight * sq_numindivids / 10^6 * 180

#How many individuals are caught in status quo? (in proportion of individuals)
props <- mutate(props, harvest_individs = harvest_weight / newweight)

#Recruitment -- number of 7 cm individuals next period -- is same as number of 7 cm 
#individuals this period.
#Will need to allow differential recruitment in counterfactual though. 
#Perea et al. (2011) finds anchoveta 12 - 14 cm contribute 40% of eggs while 
#anchoveta >= 14 cm contribute 60% of eggs. 
#So recruitment = rc*(.4*N_12-14 + .6*N_>=14)
#Solve for recruitment constant (rc) that will allow me to calculate counterfactual recruitment
#(using different values of N_12-14 and N_>=14 but same alpha)
recruit_constant <- props$prop[props$length == 7] / (
  .4*sum(props$prop[props$length >= 12 & props$length < 14]) + 
    .6*sum(props$prop[props$length >= 14])
)













#step forward one day, harvest same weight as sq_harvest but allocated differently 
#across length bins (1/1.5 fewer juveniles), add recruitment once initially 7 cm bin reaches 7.5 cm bin.
#calculate difference between period t and t+1 harvest in each length bin
#and stop simulation once harvest in each length bin converges. 

#Closures increase juvenile catch by 50%, so calculate number of juveniles caught 
#in counterfactual where catch 50% fewer
props <- mutate(props, harvest_individs_counter = if_else(length < 12, harvest_individs/1.5, as.numeric(NA)))

#Calculate the weight of juveniles caught in counterfactual
props <- mutate(props, harvest_weight_counter = harvest_individs_counter * newweight)

#Total juvenile harvest in normalized g
juv_harvest_weight_counter <- sum(props$harvest_weight_counter, na.rm = T)

#Catch same weight in counterfactual as in status quo, so need to catch more adults in counterfactual
adult_harvest_weight_counter <- sq_harvest_weight - juv_harvest_weight_counter


#Upweight counterfactual harvest of adults so that total counterfactual harvest equals total sq harvest
props <- mutate(props, harvest_weight_counter = if_else(length >= 12, 
                                                        harvest_weight * 
                                                          (adult_harvest_weight_counter / sum(props$harvest_weight[props$length >= 12])),
                                                        harvest_weight_counter))

#Confirm counterfactual and status quo total harvest are the same
sum(props$harvest_weight) == sum(props$harvest_weight_counter)

#Now fill in missing harvest_individs_counter values
props <- mutate(props, harvest_individs_counter = if_else(length >= 12, 
                                                          harvest_weight_counter / newweight,
                                                          harvest_individs_counter))

#What % more adult individuals are caught in counterfactual?
#137% more in this simulation (very different from reduced form result, but this simulation is just illustrative)
(sum(props$harvest_individs_counter[props$length >= 12]) - 
    sum(props$harvest_individs[props$length >= 12])) / 
  sum(props$harvest_individs[props$length >= 12])


counterdf <- dplyr::select(props, length, prop, age_years, weight, biomass, harvest_individs_counter) %>% 
  #And set "new" values as starting values. will update these in simulation
  mutate(newlength = length, newprop = prop, newage_years = age_years, newweight = weight, newbiomass = biomass)

#Starting value for biomass difference is just biomass values themselves
biomassdif <- counterdf$biomass

#Currently harvest is too high
#What if just find equilibrium without harvest to get simulation working?

#THIS IS WHERE I AM
#Got convergence, but newbiomass is much lower than initial biomass, which seems weird

while(max(abs(biomassdif)) > 1e-5){
  
  #Starting levels of biomass
  startbiomass <- counterdf$newbiomass
  
  #Increase age by one day
  counterdf <- mutate(counterdf, newage_years = newage_years + 1/365) %>% 
    #one day of growth
    mutate(newlength = agelength(newage_years)) %>% 
    mutate(newweight = lengthweight(newlength)) %>%
    #One day of death
    mutate(newprop = survival(newprop, 1/365)) %>% 
    #One day of harvest
    #mutate(newprop = newprop - harvest_individs_counter) %>% 
    #New normalized biomass of fish in each length bin
    mutate(newbiomass = newweight*newprop)
  
  #which rows have negative newprop? that means there are no more individuals of this length bin
  whichneg <- which(counterdf$newprop < 0)
  
  if(length(whichneg) > 0){
    
    #Drop these elements from counterdf 
    counterdf <- counterdf[-whichneg,]
    
    #Also drop them from startbiomass so when re-calculating biomassdif the two vectors have the same number of elements
    startbiomass <- startbiomass[-whichneg]
    
  }

  #which rows have age > 3? 3 is max age of anchoveta so drop any rows older than 3
  which3 <- which(counterdf$newage_years > 3)
  
  if(length(which3) > 0){
    
    #Drop these elements from counterdf 
    counterdf <- counterdf[-which3,]
    
    #Also drop them from startbiomass so when re-calculating biomassdif the two vectors have the same number of elements
    startbiomass <- startbiomass[-which3]
    
  }
  
  #Re-calculate biomassdif
  biomassdif <- startbiomass - counterdf$newbiomass
  
  if(min(counterdf$newlength) >= 7.5){
    
    #add recruitment
    prop7 <- recruit_constant*(
      .4*sum(counterdf$newprop[counterdf$newlength >= 12 & counterdf$newlength < 14]) + 
        .6*sum(counterdf$newprop[counterdf$newlength >= 14])
    )
    
    counterdf <- bind_rows(
      data.frame(length = 7, prop = 0.01331884, age_years = 0.4484462, weight = 1.96214, 
                 biomass = 0.02613342, harvest_individs_counter = 0.0001125001, 
                 newlength = 7, newprop = prop7, newage_years = 0.4484462, newweight = 1.96214, 
                 newbiomass = 0.02613342),
      counterdf
    )
    
  }

  #Maximum length of anchoveta is 20.59 cm: Tabla 4 IMARPE (2019)
  counterdf <- filter(counterdf, newlength <= 20.59)
  
}
