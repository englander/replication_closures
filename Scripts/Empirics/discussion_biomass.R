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

#Distribute proportion uniformly over .01 cm interval within each .5 cm interval
props <- map_df(1:nrow(props), function(x){
  
  myrow <- props[x,]
  
  out <- data.frame(length = seq(from = myrow$length, to = (myrow$length + 0.49), by = 0.01))
  
  out <- mutate(out, prop = myrow$prop / nrow(out))
  
  out
  
})

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

#Fraction of harvest by weight that is juvenile
harvest_juv_frac <- (sum(fullbe$tonsjuv, na.rm=T) / (sum(fullbe$tonsjuv, na.rm=T) + 
                                                       sum(fullbe$tonsadult, na.rm=T)))

#Given harvest per day, simulate growth, natural mortality, recruitment (reproduction), 
#and harvest until biomass converges. Return equilibrium biomass if returnbiomass == T, 
#otherwise return equilibrium population (sqdf)
sq_sim <- function(myhpd, returnbiomass){
  
  #Set "new" values as starting values. will update these in simulation
  sqdf <- mutate(props, newlength = length, newprop = prop, newage_years = age_years, 
                 newweight = weight, newbiomass = biomass)
  
  #Equilibrium is when total biomass doesn't change
  biomassdif <- sqdf$biomass %>% sum()
  
  while(max(abs(biomassdif)) > 1e-5){
    
    #Starting levels of biomass
    startbiomass <- sum(sqdf$newbiomass)
    
    #Increase age by one day
    sqdf <- mutate(sqdf, newage_years = newage_years + 1/365) %>% 
      #one day of growth
      mutate(newlength = agelength(newage_years)) %>% 
      mutate(newweight = lengthweight(newlength)) %>%
      #One day of death
      mutate(newprop = survival(newprop, 1/365))
    
    #Recruitment accrues each day
    prop7 <- recruit_constant*(
      .4*sum(sqdf$newprop[sqdf$newlength >= 12 & sqdf$newlength < 14]) + 
        .6*sum(sqdf$newprop[sqdf$newlength >= 14])
    )
    
    #Add this new recruited 7 cm length class to sqdf 
    sqdf <- bind_rows(
      data.frame(length = 7, prop = as.numeric(NA), age_years = as.numeric(NA), weight = as.numeric(NA),
                 biomass = as.numeric(NA),
                 newlength = 7, newprop = prop7, newage_years = 0.4484462, newweight = 1.96214,
                 newbiomass = 0.02613342),
      sqdf
    )
    
    #juvenile and adult biomass for allocating harvest proportionally
    juvbiomass <- sum(sqdf$newbiomass[sqdf$length < 12])
    adultbiomass <- sum(sqdf$newbiomass[sqdf$length >= 12])
    
    #Harvest
    sqdf <- sqdf %>% 
      #One day of harvest. New normalized biomass of fish in each length bin
      mutate(newbiomass = if_else(newlength < 12,
                                  newweight*newprop -
                                    myhpd*harvest_juv_frac*(newbiomass / juvbiomass),
                                  newweight*newprop -
                                    myhpd*(1 - harvest_juv_frac)*(newbiomass / adultbiomass)
      )) %>%
      #mutate(newbiomass = newweight*newprop) %>% 
      #New proportion of individuals
      mutate(newprop = newbiomass / newweight)
    
    #which rows have negative newprop? that means there are no more individuals of this length bin
    whichneg <- which(sqdf$newprop < 0)
    
    if(length(whichneg) > 0){
      
      #Drop these elements from sqdf 
      sqdf <- sqdf[-whichneg,]
      
    }
    
    #which rows have age > 3? 3 is max age of anchoveta so drop any rows older than 3
    which3 <- which(sqdf$newage_years > 3)
    
    if(length(which3) > 0){
      
      #Drop these elements from sqdf 
      sqdf <- sqdf[-which3,]
      
    }
    
    #Maximum length of anchoveta is 20.59 cm: Tabla 4 IMARPE (2019)
    sqdf <- filter(sqdf, newlength <= 20.59)
    
    #Re-calculate biomassdif
    biomassdif <- startbiomass - sum(sqdf$newbiomass)
    
  }
  
  if(returnbiomass == T){
    #Output equilibrium biomass
    out <- data.frame(harvest_per_day = myhpd, eqbiomass = sum(sqdf$newbiomass))
  } else{
    #Otherwise return equilibrium population
    out <- sqdf
  }
  
  return(out)
  
}

sq_biomass <- map_df(c(seq(from = 0, to = harvest_per_day, by = .001), harvest_per_day), function(x){
  sq_sim(x, T)
})

#The simulation is highly unstable: 
# sq_biomass
# harvest_per_day    eqbiomass
# 1       0.00000000 1.411371e+01
# 2       0.00100000 0.000000e+00
# 3       0.00200000 0.000000e+00
# 4       0.00300000 0.000000e+00
# 5       0.00400000 1.199643e+01
# 6       0.00500000 0.000000e+00
# 7       0.00600000 6.031774e-06
# 8       0.00700000 1.650649e-05
# 9       0.00800000 2.999305e-05
# 10      0.00900000 3.917417e-05
# 11      0.01000000 1.046065e+01
# 12      0.01100000 5.379711e-05
# 13      0.01200000 6.423925e-05
# 14      0.01300000 5.850758e-05
# 15      0.01305276 6.309077e-05

#Use the highest harvest per day that has a biomass size close to starting biomass: 
# 11      0.01000000 1.046065e+01

sqdf <- sq_sim(.01, F)

#Plot biomass by length
ggplot(data = sqdf, aes(x = length, y = newbiomass)) + 
  geom_line() + 
  geom_vline(aes(xintercept = 12))

#Plot number of individuals by length
ggplot(data = sqdf, aes(x = length, y = newprop)) + 
  geom_line() + 
  geom_vline(aes(xintercept = 12))


##Now do counterfactual. Same total harvest per day, but 44% fewer juveniles by weight 
harvest_juv_frac_counter <- harvest_juv_frac / 1.44

counter_sim <- function(myhpd, returnbiomass){
  
  #Set "new" values as starting values. will update these in simulation
  counterdf <- mutate(props, newlength = length, newprop = prop, newage_years = age_years, 
                      newweight = weight, newbiomass = biomass)
  
  #Equilibrium is when total biomass doesn't change
  biomassdif <- counterdf$biomass %>% sum()
  
  while(max(abs(biomassdif)) > 1e-5){
    
    #Starting levels of biomass
    startbiomass <- sum(counterdf$newbiomass)
    
    #Increase age by one day
    counterdf <- mutate(counterdf, newage_years = newage_years + 1/365) %>% 
      #one day of growth
      mutate(newlength = agelength(newage_years)) %>% 
      mutate(newweight = lengthweight(newlength)) %>%
      #One day of death
      mutate(newprop = survival(newprop, 1/365))
    
    #Recruitment accrues each day
    prop7 <- recruit_constant*(
      .4*sum(counterdf$newprop[counterdf$newlength >= 12 & counterdf$newlength < 14]) + 
        .6*sum(counterdf$newprop[counterdf$newlength >= 14])
    )
    
    #Add this new recruited 7 cm length class to counterdf 
    counterdf <- bind_rows(
      data.frame(length = 7, prop = as.numeric(NA), age_years = as.numeric(NA), weight = as.numeric(NA),
                 biomass = as.numeric(NA),
                 newlength = 7, newprop = prop7, newage_years = 0.4484462, newweight = 1.96214,
                 newbiomass = 0.02613342),
      counterdf
    )
    
    #juvenile and adult biomass for allocating harvest proportionally
    juvbiomass <- sum(counterdf$newbiomass[counterdf$length < 12])
    adultbiomass <- sum(counterdf$newbiomass[counterdf$length >= 12])
    
    #Harvest
    counterdf <- counterdf %>% 
      #One day of harvest. New normalized biomass of fish in each length bin
      mutate(newbiomass = if_else(newlength < 12,
                                  newweight*newprop -
                                    myhpd*harvest_juv_frac_counter*(newbiomass / juvbiomass),
                                  newweight*newprop -
                                    myhpd*(1 - harvest_juv_frac_counter)*(newbiomass / adultbiomass)
      )) %>%
      #mutate(newbiomass = newweight*newprop) %>% 
      #New proportion of individuals
      mutate(newprop = newbiomass / newweight)
    
    #which rows have negative newprop? that means there are no more individuals of this length bin
    whichneg <- which(counterdf$newprop < 0)
    
    if(length(whichneg) > 0){
      
      #Drop these elements from counterdf 
      counterdf <- counterdf[-whichneg,]
      
    }
    
    #which rows have age > 3? 3 is max age of anchoveta so drop any rows older than 3
    which3 <- which(counterdf$newage_years > 3)
    
    if(length(which3) > 0){
      
      #Drop these elements from counterdf 
      counterdf <- counterdf[-which3,]
      
    }
    
    #Maximum length of anchoveta is 20.59 cm: Tabla 4 IMARPE (2019)
    counterdf <- filter(counterdf, newlength <= 20.59)
    
    #Re-calculate biomassdif
    biomassdif <- startbiomass - sum(counterdf$newbiomass)
    
  }
  
  if(returnbiomass == T){
    #Output equilibrium biomass
    out <- data.frame(harvest_per_day = myhpd, eqbiomass = sum(counterdf$newbiomass))
  } else{
    #Otherwise return equilibrium population
    out <- counterdf
  }
  
  return(out)
  
}

#Population converges to 0. This is because age-structured models are very unstable
counter_sim(.01, T)

#Try simulating over a range of values
counter_biomass <- map_df(c(seq(from = 0, to = harvest_per_day, by = .001), harvest_per_day), function(x){
  counter_sim(x, T)
})
# harvest_per_day    eqbiomass
# 1       0.00000000 1.411371e+01
# 2       0.00100000 0.000000e+00
# 3       0.00200000 0.000000e+00
# 4       0.00300000 1.239315e+01
# 5       0.00400000 0.000000e+00
# 6       0.00500000 0.000000e+00
# 7       0.00600000 9.680707e-07
# 8       0.00700000 4.711383e-06
# 9       0.00800000 8.514486e-06
# 10      0.00900000 1.247117e-05
# 11      0.01000000 1.770528e-05
# 12      0.01100000 1.922137e-05
# 13      0.01200000 2.268590e-05
# 14      0.01300000 2.062169e-05
# 15      0.01305276 2.579606e-05

#The closest harvest per day to status quo harvest per day is 
# 4       0.00300000 1.239315e+01

counterdf <- counter_sim(0.003, F)

#Plot biomass by length
ggplot(data = counterdf, aes(x = length, y = newbiomass)) + 
  geom_line() + 
  geom_vline(aes(xintercept = 12))

#Plot number of individuals by length
ggplot(data = counterdf, aes(x = length, y = newprop)) + 
  geom_line() + 
  geom_vline(aes(xintercept = 12))


#Plot status quo and counterfactual biomass by length
ggplot() + 
  geom_line(data = sqdf, aes(x = length, y = newbiomass)) + 
  geom_line(data = counterdf, aes(x = length, y = newbiomass), col = 'red')











