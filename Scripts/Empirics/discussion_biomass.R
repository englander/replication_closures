#Digitize number of individuals in each length bin at start of first season 2017
#This is population distribution (in the water), not distribution of catch

usethis::use_git_config(user.name = "englander", user.email = "genglander@ucsb.edu")
credentials::set_github_pat("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")

rm(list=ls())
setwd("C:/Users/gabee/Documents/replication_closures")

library(dplyr); library(sf); library(ggplot2); library(lubridate)
library(readr); library(purrr); library(readxl); library(imager)
library(furrr); library(tidyr); library(cowplot)

#Peru time
Sys.setenv(TZ='America/Lima')

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
survival <- function(ntm1, timestep, decay){
  
  #Number of survivors in next period
  nt <- ntm1*exp(-decay*timestep)
  
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

#Biomass at start of season is 8 or 9 million tons
#Since two fishing seasons, catch about 25% of biomass each season
#Define this "adaptive management" harvest within while loop

#Fraction of harvest by weight that is juvenile
harvest_juv_frac <- (sum(fullbe$tonsjuv, na.rm=T) / (sum(fullbe$tonsjuv, na.rm=T) + 
                                                       sum(fullbe$tonsadult, na.rm=T)))

#Clean up
rm(topy, lfd, bottomy, minticks, majticks, glfd, fullbe)

#Given whether want to return summed biomass or full population size structure data frame,
#whether convergence condition is based on biomass (default) or want to run simulation for certain number of years,
#decay rate for survival function, recruitment constant, and juvenile fraction (in biomass terms),
#simulate recruitment (reproduction), growth, natural mortality, 
#and harvest until biomass converges. Return equilibrium biomass if returnbiomass == T, 
#otherwise return equilibrium population (sqdf)
sq_sim <- function(returnbiomass, convergecondition, decay, myrc, myjuvfrac){
  
  #Set "new" values as starting values. will update these in simulation
  sqdf <- mutate(props, newlength = length, newprop = prop, newage_years = age_years, 
                 newweight = weight, newbiomass = biomass)
  
  harvestdf <- data.frame(year = as.numeric(NULL), 
                          season = as.numeric(NULL),
                          adultbiomass_start = as.numeric(NULL), #adult biomass at start of season
                          juvbiomass_start = as.numeric(NULL),
                          adultbiomass_end = as.numeric(NULL), #adult biomass at end of season
                          juvbiomass_end = as.numeric(NULL),
                          juvharvest = as.numeric(NULL), #total juvenile harvest in biomass terms over course of this season)
                          adultharvest = as.numeric(NULL) #total adult harvest in biomass terms over course of this season)
  )
                        
    
  mytime <- 0
  juvharvestcounter <- 0 #count juv harvest each day, then record total harvest over season after season ends
  adultharvestcounter <- 0 #count juv harvest each day, then record total harvest over season after season ends
  
  #Simulate population for convergecondition years and as long as population is not extinct
  while(mytime < convergecondition*365 & sum(sqdf$newbiomass) > 0){
    
    #Increase age by one day 
    sqdf <- mutate(sqdf, newage_years = newage_years + 1/365) %>% 
      #one day of growth
      mutate(newlength = agelength(newage_years)) %>% 
      mutate(newweight = lengthweight(newlength)) %>% 
      #One day of death
      mutate(newprop = survival(newprop, 1/365, decay))
    
    #Recruitment accrues each day
    prop7 <- myrc*(
      .4*sum(sqdf$newprop[sqdf$newlength >= 12 & sqdf$newlength < 14]) + 
        .6*sum(sqdf$newprop[sqdf$newlength >= 14])
    )
    
    #Add new recruited 7 cm length class to sqdf (new recruited class does not grow or die this period)
    sqdf <- bind_rows(
      data.frame(length = 7, prop = as.numeric(NA), age_years = as.numeric(NA), weight = as.numeric(NA),
                 biomass = as.numeric(NA),
                 newlength = 7, newprop = prop7, newage_years = 0.4484462, newweight = 1.96214) %>% 
        mutate(newbiomass = newweight*newprop),
      sqdf
    )
    
    #Harvest if within season and adult biomass > 4. 
    #Two fishing seasons of 91 days each.
    if(mytime > 31 & #first fishing season starts about one month after biomass measurement
       (
         (mod(mytime, 365) > 31 & mod(mytime, 365) < 123) | 
         (mod(mytime, 365) > 214 & mod(mytime, 365) < 306)
       )){
      
      #Total biomass at start of season
      if(round(mod(mytime, 365)) %in% c(32, 215)){
        
        #adult and juv biomass at start of season
        myadultbiomass_start <- sum(sqdf$newbiomass[sqdf$newlength >= 12])
        myjuvbiomass_start <- sum(sqdf$newbiomass[sqdf$newlength < 12])
        
        #year and season of this season
        myseason <- ifelse(round(mod(mytime, 365)) == 32, 1, 2)
        myyear <- ceiling(mytime / 365)
        
        #If adult biomass < 4, don't allow fishing this season
        if(sum(sqdf$newbiomass[sqdf$newlength >= 12]) < 4){
          
          myhpd <- 0
          
        }else{
          
          #Otherwise harvest 1/4 of biomass
          #total biomass
          totbiomass <- sum(sqdf$newbiomass)
          
          #harvest 25% of total biomass
          totharvest <- totbiomass / 4
          
          #Harvest per day (91 days in season)
          myhpd <- totharvest / 91
        }
        
      }
      
      #Juvenile and adult total biomass this period 
      juvbiomass <- sum(sqdf$newbiomass[sqdf$newlength < 12])
      adultbiomass <- sum(sqdf$newbiomass[sqdf$newlength >= 12])  
      
      #Calculate harvest this period in each length interval
      harvest_todaydf <- mutate(sqdf, 
                                harvest_today = if_else(newlength < 12,
                                                        #myjuvfrac % of total harvest (myhpd) comes from juveniles. 
                                                        #Then allocate that harvest across juveniles in proportion to relative biomass of each length interval ((newbiomass / juvbiomass))
                                                        myhpd*myjuvfrac*(newbiomass / juvbiomass),
                                                        myhpd*(1 - myjuvfrac)*(newbiomass / adultbiomass)
                                ))
      
      #Add today's juvenile harvest to juvenile harvest_counter
      juvharvestcounter <- juvharvestcounter + 
        sum(harvest_todaydf$harvest_today[harvest_todaydf$newlength < 12])
      
      #Add today's juvenile harvest to juvenile harvest_counter
      adultharvestcounter <- adultharvestcounter + 
        sum(harvest_todaydf$harvest_today[harvest_todaydf$newlength >= 12])
      
      #Harvest
      sqdf <- sqdf %>% 
        #One day of harvest. New normalized biomass of fish in each length bin
        mutate(newbiomass = if_else(newlength < 12,
                                    newweight*newprop -
                                      #myjuvfrac % of total harvest (myhpd) comes from juveniles. 
                                      #Then allocate that harvest across juveniles in proportion to relative biomass of each length interval ((newbiomass / juvbiomass))
                                      myhpd*myjuvfrac*(newbiomass / juvbiomass),
                                    newweight*newprop -
                                      myhpd*(1 - myjuvfrac)*(newbiomass / adultbiomass)
        )) %>%
        #If newbiomass has become NaN (this occurs when juvbiomass = 0 or adultbiomass = 0), 
        #reset newbiomass to 0
        mutate(newbiomass = if_else(is.nan(newbiomass), 0, newbiomass)) %>% 
        #New proportion of individuals
        mutate(newprop = newbiomass / newweight)
      
    }else{
      #Otherwise we are not in fishing season and just update biomass to reflect newweight and newprop
      
      #If end of season, record adult and juv biomass
      #and add a row to harvestdf to record values of interest
      if(round(mod(mytime, 365)) %in% c(123, 306)){
        
        #adult and juv biomass at end of season
        myadultbiomass_end <- sum(sqdf$newbiomass[sqdf$newlength >= 12])
        myjuvbiomass_end <- sum(sqdf$newbiomass[sqdf$newlength < 12])
        
        harvestdf <- bind_rows(harvestdf, 
                               data.frame(
                                 year = myyear, season = myseason, 
                                 adultbiomass_start = myadultbiomass_start,
                                 juvbiomass_start = myjuvbiomass_start,
                                 adultbiomass_end = myadultbiomass_end,
                                 juvbiomass_end = myjuvbiomass_end,
                                 juvharvest = juvharvestcounter, 
                                 adultharvest = adultharvestcounter
                               ))
        
        #Reset harvest counters for next season
        juvharvestcounter <- 0
        adultharvestcounter <- 0
        
      }
      
      sqdf <- mutate(sqdf, newbiomass = newweight*newprop)
    }
    
    
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
    
    mytime <- mytime + 1
    
  }
  

  if(returnbiomass == T){
    #Output equilibrium biomass
    out <- data.frame(eqbiomass = sum(sqdf$newbiomass), 
                      sim_days = mytime)
  } else{
    #Otherwise return equilibrium population and harvestdf
    out <- list(
      sqdf %>% 
      mutate(sim_days = mytime),
      harvestdf)
  }
  
  return(out)
  
}


#Stock goes extinct with current recruitment constant, fishing mortality, and natural mortality
#Since recruitment constant is made up, while fishing mortality comes from data and natural mortality
#comes from Salvatteci and Mendo (2005), choose the recruitment constant such that 
#biomass after 100 years is close to initial biomass. 
# plan(multisession, workers = 10)
# 
# rclist <- future_map(seq(from = 6.75, to = 7, by = 0.001), function(x){
#   
#   eqbiomass <- try(sq_sim(returnbiomass = T, convergecondition = 100, 
#                           decay = 0.8, myrc = recruit_constant * x, 
#                           myjuvfrac = harvest_juv_frac)) %>% 
#     mutate(rcfactor = x)
#   
# })
# 
# rcdf <- bind_rows(rclist)

#recruit_constant * 6.853 gives biomass after 20 years very close to initial biomass:

#Status quo simulation over 100 years
g1 <- sq_sim(F, convergecondition = 100, decay = 0.8, myrc = recruit_constant * 6.9, 
            myjuvfrac = harvest_juv_frac)

sum(g1[[1]]$newbiomass) #8.96961
#Compare to status quo initial biomass
sum(props$biomass) #8.967002

#Plot number of individuals by length
ggplot() + 
  geom_line(data = props, aes(x = length, y = prop)) + 
  geom_line(data = g1[[1]], aes(x = newlength, y = newprop), col = 'red') 

#Plot biomass by length
ggplot() + 
  geom_line(data = props, aes(x = length, y = biomass)) + 
  geom_line(data = g1[[1]], aes(x = newlength, y = newbiomass), col = 'red')

#Juvenile and adult harvest over time
ggplot(data = g1[[2]] %>% mutate(yeartime = if_else(season == 2, year + .5, year)) %>% 
         pivot_longer(cols = c(juvharvest, adultharvest), names_to = 'harvesttype', values_to = 'harvest') %>% 
         mutate(harvesttype = as.factor(harvesttype)), 
       aes(x = yeartime, y = harvest, col = harvesttype)
         ) + 
  geom_line()

#How does changing harvest_juv_frac affect final biomass and harvest?
#It seems like it reduces it. biomass and harvest would be higher if caught fewer juveniles.
g0 <- sq_sim(F, convergecondition = 20, decay = 0.8, myrc = recruit_constant * 6.853, 
             myjuvfrac = harvest_juv_frac / 1.44)

sum(g0[[1]]$newbiomass) #13.81285
sum(g0[[2]]$adultharvest) + sum(g0[[2]]$juvharvest) #143.3653

g1 <- sq_sim(F, convergecondition = 20, decay = 0.8, myrc = recruit_constant * 6.853, 
             myjuvfrac = harvest_juv_frac)

sum(g1[[1]]$newbiomass) #8.96961
sum(g1[[2]]$adultharvest) + sum(g1[[2]]$juvharvest) #117.6033

#compare distributions
ggplot() + 
  geom_line(data = g0[[1]], aes(x = newlength, y = newprop), col = 'red') + 
  geom_line(data = g1[[1]], aes(x = newlength, y = newprop)) 

#Plot biomass by length
ggplot() + 
  geom_line(data = g0[[1]], aes(x = newlength, y = newbiomass), col = 'red') + 
  geom_line(data = g1[[1]], aes(x = newlength, y = newbiomass))

#Harvest start out same, and then overharvesting of juveniles causes divergence
ggplot() + 
  geom_line(data = g0[[2]] %>% 
              mutate(yeartime = if_else(season == 2, year + .5, year)),
            aes(x = yeartime, y = harvest), col = 'red') + 
  geom_line(data = g1[[2]] %>% 
              mutate(yeartime = if_else(season == 2, year + .5, year)),
            aes(x = yeartime, y = harvest))



##Make a 2 x 2 plot. 

#First plot is total harvest, adult harvest, and juvenile harvest each season over time, 
#in 
harvestdf <- bind_rows(
  g1[[2]] %>% mutate(yeartime = if_else(season == 2, year + .5, year), 
                     scenario = 'status quo') %>% 
    dplyr::select(yeartime, scenario, adultharvest, juvharvest) %>% 
    rename(adult = adultharvest, juvenile = juvharvest) %>%
    mutate(total = adult + juvenile) %>%
    pivot_longer(cols = c(total, adult, juvenile), names_to = 'harvesttype', values_to = 'harvestquantity'),
  g0[[2]] %>% mutate(yeartime = if_else(season == 2, year + .5, year), 
                     scenario = 'counterfactual') %>% 
    dplyr::select(yeartime, scenario, adultharvest, juvharvest) %>% 
    rename(adult = adultharvest, juvenile = juvharvest) %>%
    mutate(total = adult + juvenile) %>%
    pivot_longer(cols = c(total, adult, juvenile), names_to = 'harvesttype', values_to = 'harvestquantity'),
)

#Scale harvest into millions of tons. initial biomass is 7.78 million tons
harvestdf <- mutate(harvestdf, harvestquantity = 
                      harvestquantity * (7.78 / sum(props$biomass)))

harvestdf$harvesttype <- as.factor(harvestdf$harvesttype)
harvestdf$harvesttype <- relevel(harvestdf$harvesttype, ref = "total")

harvestdf$scenario <- as.factor(harvestdf$scenario)
harvestdf$scenario <- relevel(harvestdf$scenario, ref = 'counterfactual')


(harvest_time <- ggplot(data = harvestdf, 
                        aes(x = yeartime, y = harvestquantity, col = harvesttype, linetype = scenario)) + 
    geom_line() +
    myThemeStuff + 
    scale_x_continuous("Year", breaks = c(1, 6, 11, 16, 20),
                       labels = c(2017, 2022, 2027, 2032, 2036)) + 
    scale_y_continuous("Harvest (millions of tons)") + 
    scale_color_manual("Harvest type", values = c("black", "dodgerblue4","orange1")) + 
    scale_linetype_manual(values = c("dotted", "solid")) + 
    ggtitle("Harvest each season") + 
    labs(tag = "a") + theme(plot.tag.position = c(.05, 1), 
                            plot.margin = unit(c(.05,0.04,.1,0),"in")))



#Second plot is total biomass, adult biomass, and juv biomass in status quo and counterfactual 
#over time. Use biomass recorded at start of each season.
biomassdf <- bind_rows(
  g1[[2]] %>% mutate(yeartime = if_else(season == 2, year + .5, year), 
                     scenario = 'status quo') %>% 
    dplyr::select(yeartime, scenario, adultbiomass_start, juvbiomass_start) %>% #
    rename(adult = adultbiomass_start, juvenile = juvbiomass_start) %>%
    mutate(total = adult + juvenile) %>%
    pivot_longer(cols = c(total, adult, juvenile), names_to = 'biomasstype', values_to = 'biomassquantity'),
  g0[[2]] %>% mutate(yeartime = if_else(season == 2, year + .5, year), 
                     scenario = 'counterfactual') %>% 
    dplyr::select(yeartime, scenario, adultbiomass_start, juvbiomass_start) %>% #
    rename(adult = adultbiomass_start, juvenile = juvbiomass_start) %>%
    mutate(total = adult + juvenile) %>%
    pivot_longer(cols = c(total, adult, juvenile), names_to = 'biomasstype', values_to = 'biomassquantity') 
)

#Scale biomass into millions of tons. initial biomass is 7.78 million tons
biomassdf <- mutate(biomassdf, biomassquantity = 
                      biomassquantity * (7.78 / sum(props$biomass)))

biomassdf$biomasstype <- as.factor(biomassdf$biomasstype)
biomassdf$biomasstype <- relevel(biomassdf$biomasstype, ref = "total")

biomassdf$scenario <- as.factor(biomassdf$scenario)
biomassdf$scenario <- relevel(biomassdf$scenario, ref = 'counterfactual')


(biomass_time <- ggplot(data = biomassdf, 
       aes(x = yeartime, y = biomassquantity, col = biomasstype, linetype = scenario)) + 
  geom_line() +
  myThemeStuff + 
  scale_x_continuous("Year", breaks = c(1, 6, 11, 16, 20),
                     labels = c(2017, 2022, 2027, 2032, 2036)) + 
  scale_y_continuous("Biomass (millions of tons)") + 
  scale_color_manual("biomass value", values = c("black", "dodgerblue4","orange1")) + 
  scale_linetype_manual(values = c("dotted", "solid")) + 
  ggtitle("Biomass at start of season") + 
  labs(tag = "b") + theme(plot.tag.position = c(.05, 1), 
                          plot.margin = unit(c(.05,0,.1,0.04),"in")))


#Third plot proportion of individuals in each length interval
#as measured in March 2017, in status quo at end of simulation,
#and in counterfactual at end of simulation
propdf <- bind_rows(
  dplyr::select(props, length, prop) %>% mutate(scenario = 'measured, 2017'), 
  g1[[1]] %>% dplyr::select(newprop, newlength) %>% mutate(scenario = 'status quo') %>%
    #Normalize proportion to 1 (currently units are normalized number of individuals)
    mutate(newprop = newprop / sum(g1[[1]]$newprop)) %>% 
    rename(prop = newprop, length = newlength),
  g0[[1]] %>% dplyr::select(newprop, newlength) %>% mutate(scenario = 'counterfactual') %>%
    #Normalize proportion to 1 (currently units are normalized number of individuals)
    mutate(newprop = newprop / sum(g0[[1]]$newprop)) %>% 
    rename(prop = newprop, length = newlength)
)

propdf$scenario <- as.factor(propdf$scenario)
propdf$scenario <- relevel(propdf$scenario, ref = 'status quo')
propdf$scenario <- relevel(propdf$scenario, ref = 'counterfactual')

#Group individuals into common length intervals
#status quo and counterfactual are already in common length intervals
unique(propdf$length[propdf$scenario == 'counterfactual'])[unique(propdf$length[propdf$scenario == 'counterfactual']) %in% 
    unique(propdf$length[propdf$scenario == 'status quo'])] %>% length() == 
  unique(propdf$length[propdf$scenario == 'counterfactual']) %>% length()

#Save these intervals
commonintervals <- unique(propdf$length[propdf$scenario == 'counterfactual'])

#But need to redistribute measured 2017 into these same length intervals
#First aggregate measured 2017 values to original interval
m17 <- filter(propdf, scenario == 'measured, 2017')


rm(x, lb, ub, out, nextinter, belowfrac, firstrow, myrow, previnter, secondrow, abovefrac)
#Distribute proportion uniformly over interval defined by status quo and counterfactual intervals
m17 <- map_df(1:100, function(x){
  
  myrow <- m17[x,]
  
  #which intervals is this length below
  lb <- commonintervals[commonintervals < myrow$length] %>% max()
  ub <- commonintervals[commonintervals > myrow$length] %>% min()
  
  if(myrow$length == 7){
    lb <- 7
  }
  
  #is length between a common interval?
  if(round(ub - myrow$length, digits = 4) < 0.01 & 
     min(myrow$length - lb, ub - myrow$length) != 0){ #excluding length = 7
    
    #If myrow$length below ub, 
    #need to output two rows: density between mylength and ub, 
    #and density between ub and next common interval
    if(which(round(c(myrow$length - lb, ub - myrow$length), digits = 4) < 0.01) == 2){
      
      firstrow <- data.frame(length = seq(from = myrow$length, to = ub, by = 0.0001))
      
      #Fraction of original length interval that is below ub
      belowfrac <- (ub - myrow$length) / .01
      
      firstrow <- mutate(firstrow, prop = (myrow$prop / nrow(firstrow)) * belowfrac)
      
      #Now sum to lb
      firstrow <- data.frame(length = lb, prop = sum(firstrow$prop))
      
      #Then remainder goes to next interval
      nextinter <- commonintervals[commonintervals > ub] %>% min()
      
      out <- bind_rows(firstrow, 
                       data.frame(length = nextinter, prop = myrow$prop - firstrow$prop))
    }else{
      
      #Otherwise myrow$length > ub and need to allocate some of mass to below interval
      #First calculate mass above lb
      secondrow <- data.frame(length = seq(from = lb, to = myrow$length, by = 0.0001))
      
      #Fraction of original length interval that is above lb
      abovefrac <- (myrow$length - lb) / .01
      
      secondrow <- mutate(secondrow, prop = (myrow$prop / nrow(secondrow)) * abovefrac)
      
      #Now sum to lb
      secondrow <- data.frame(length = lb, prop = sum(secondrow$prop))
      
      #Then remainder goes to previous
      previnter <- commonintervals[commonintervals < lb] %>% max()
      
      out <- bind_rows(secondrow, 
                       data.frame(length = previnter, prop = myrow$prop - secondrow$prop))
      
      
    }
      
  }else{
    
    #If not allocate entire mass to this common length interval
    out <- data.frame(length = lb, prop = myrow$prop)
    
  }
  
  out
  
}) #%>% 
  #group_by(length) %>% 
  #summarise(prop = sum(prop))

arrange(m17, length) %>% head()

group_by(m17, length) %>% summarise(prop = sum(prop))

#Given row of propdf, save the 

propdf$commonlength <- cut(propdf$length, , include.lowest = T)

#Count number of 

#Calculate number of individuals rather than prop, as measured in 2017
n_2017 <- 7.78 * 10^6 * 10^6 /#biomass in g
  mutate(props, avgweight = prop*weight) %>%   #Average weight of individual in grams, measured in 2017.
  summarise(sum(avgweight)) %>%  #prop already sums to 1 so don't need to divide by sum(prop)
  as.matrix() %>% as.numeric()

propdf <- mutate(propdf, n = prop * n_2017)

group_by(propdf, scenario) %>% 
  summarise(sum(n))

ggplot(data = propdf, aes(x = length, y = n, linetype = scenario)) + 
  geom_line() + 
  myThemeStuff + 
  ggtitle("Length distribution of population") + 
  labs(tag = "c") + theme(plot.tag.position = c(.05, 1), 
                          plot.margin = unit(c(.05,0,.1,0.04),"in")) + 
  scale_linetype_manual(values = c("dotted", "solid", 'dashed')) + 
  ylab("Number of individuals") + 
  scale_x_continuous("Length (cm)", breaks = seq(from = 7, to = 18, by = 1))



#Fourth plot is biomass
beqdf <- bind_rows(
  dplyr::select(props, length, biomass) %>% mutate(scenario = 'measured, 2017'), 
  g1[[1]] %>% dplyr::select(newbiomass, newlength) %>% mutate(scenario = 'status quo') %>%
    rename(prop = newprop, length = newlength),
  g0[[1]] %>% dplyr::select(newprop, newlength) %>% mutate(scenario = 'counterfactual') %>%
    #Normalize proportion to 1 (currently units are normalized number of individuals)
    mutate(newprop = newprop / sum(g0[[1]]$newprop)) %>% 
    rename(prop = newprop, length = newlength)
)

propdf$scenario <- as.factor(propdf$scenario)
propdf$scenario <- relevel(propdf$scenario, ref = 'counterfactual')

ggplot(data = propdf, aes(x = length, y = prop, linetype = scenario)) + 
  geom_line() + 
  myThemeStuff + 
  ggtitle("Length distribution of population") + 
  labs(tag = "c") + theme(plot.tag.position = c(.05, 1), 
                          plot.margin = unit(c(.05,0,.1,0.04),"in")) + 
  scale_linetype_manual(values = c("dotted", "solid")) + 
  ylab("Proportion of individuals") + 
  scale_x_continuous("Length (cm)", breaks = seq(from = 7, to = 18, by = 1))
  



tbt <- plot_grid(leadplot, activeplot, lag1plot, 
                 lag2plot, lag3plot, lag4plot, nrow=2, ncol=3, 
                 rel_widths = c(1.01,1,1))

ggsave(tbt, file=paste0("Output/Figures/figure8.png"),
       w=7,h=(7/1.69)*2, units = "in", dpi=1200)