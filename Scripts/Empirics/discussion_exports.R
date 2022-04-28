#Salvatteci and Mendo (2005) take size distribution of landings in a month,
#let that "population" (the fish landed that month) grow until percent juvenile 
#equals 10%. Difference in tons between this "population" and actual landings
#is the difference in landings. 
#I reproduce this method, butinstead of projecting until juvenile catch equals 10%, 
#project until juvenile catch is lower by closures treatment effect.

rm(list=ls())

setwd("C:/Users/gabee/Documents/replication_closures")

library(dplyr); library(tidyr); library(purrr)

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Want number of individuals caught in each bin
#so multiply numindivids by propbin
propbins <- grep("hat",names(fullbe),value=T)
propbins <- grep("prop",propbins,value=T)

ni <- dplyr::select(fullbe, numindivids, all_of(propbins))

ni <- filter(ni, !is.na(numindivids) & !is.na(prop10hat))

names(ni) <- gsub("hat","",names(ni))
names(ni) <- gsub("prop","",names(ni))
names(ni) <- gsub("p","\\.",names(ni))

ni <- gather(ni, length, proportion, -numindivids)

ni$length <- ni$length %>% as.numeric()

#Shift length so that it corresponds to middle of length interval
ni$length <- ni$length + .25

ni <- mutate(ni, numlength = numindivids*proportion)

#Now sum to "population level"
ni <- group_by(ni, length) %>% 
  summarise(numindivids = sum(numlength)) %>% ungroup()

#Distribute individuals uniformly over interval
distribFun <- function(rowind){
  
  row <- ni[rowind,]
  
  out <- data.frame(length = seq(from=row$length-.25,to=row$length+.24,by=.01),
                    numindivids = row$numindivids/50)
  
  return(out)
}

#Apply over rows of ni
ni <- map_df(1:nrow(ni), function(x){
  distribFun(x)
})


#Projection works by moving forward one month
#individuals grow and they die. 
#So I need to calculate age given length
#length-age from S & M
#Length is in cm and age is in years. 
lengthage <- function(length){
  
  #Using average values from Tabla 2
  age <- -.0193 - (1/.96)*log(1 - length/19.35)
  
  return(age)
}

#Survival equation
survival <- function(ntm1, timestep){
  
  #Number of survivors in next period
  nt <- ntm1*exp(-.8*timestep)
  
  return(nt)
}

#Length given age
agelength <- function(age){
  
  #Using average values from Tabla 2
  length <- 19.35*(1 - exp(-.96*(age + .0193)))
  
  return(length)
}


#Calculate age given length
ni <- mutate(ni, age = lengthage(length))

#How many juveniles are currently caught?
juv1 <- filter(ni, length < 12) %>% summarise(sum(numindivids)) %>% 
  as.matrix() %>% as.numeric() 

#Closures increase juvenile catch by 50%, so project until number of juveniles 
#caught is 1/1.50 of juv1
juv0 <- juv1

projdf <- ni; numdays <- 0

#So move forward one day at a time until number of juveniles caught is <= 1/1.5 of juv1
while (juv0/juv1 >= 1/1.5) {

  numdays <- numdays + 1
  
  #One day of growth
  projdf <- mutate(projdf, newage = age + 1/365) %>% 
    mutate(newlength = agelength(newage)) %>% 
    #One day of death
    mutate(newnumindivids = survival(numindivids, 1/365))
  
  #Update projdf by renaming columns
  projdf <- dplyr::select(projdf, -length, -numindivids, -age)
  
  names(projdf) <- gsub("new","",names(projdf))
  
  #Calculate new juv0
  juv0 <- filter(projdf, length < 12) %>% summarise(sum(numindivids)) %>% 
    as.matrix() %>% as.numeric()
  
}

#So it took 22 days for juv0 to be 1/1.5 lower than juv1
#Now calculate how much this population weighs
#I will use length-weight equation from IMARPE (2019) to be consistent with rest of paper,
#rather than (similar) length-weight equation in S & M
lengthweight <- function(length){
  #length in cm; weight in g
  weight <- .0036*length^3.238
  return(weight)
}

#There are fewer individuals than before because of natural mortality
sum(projdf$numindivids) / sum(ni$numindivids) #0.9529249

projdf <- mutate(projdf, weight = lengthweight(length))

#Convert grams to tons
(projtons <- mutate(projdf, weight = numindivids*weight / 10^6) %>% 
  summarise(sum(weight)) %>% as.matrix() %>% as.numeric()) #13423376 tons

#Tons caught in data
ni <- mutate(ni, weight = lengthweight(length))

(actualtons <- mutate(ni, weight = numindivids*weight / 10^6) %>% 
  summarise(sum(weight)) %>% as.matrix() %>% as.numeric()) #13135589

#So could have caught 287786.6 more tons
actualtons - projtons

#Or 2.1% more tons
(actualtons - projtons) / projtons



#Or 38 million lower export revenues per year (1788500000 is 2017 USD export revenues)
1788500000*((actualtons - projtons) / projtons)
#This is going to be an underestimate because does not account for greater reproduction of stock



sessionInfo()
