rm(list=ls())

#Peru time
Sys.setenv(TZ='America/Lima')


library(ggplot2); library(viridis); library(geosphere); library(purrr)
library(sf); library(dplyr); library(cowplot); library(parallel)

#Turn off spherical geometry since I wrote these scripts before sf v1
sf::sf_use_s2(FALSE) 

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
                      plot.margin = unit(c(0,0,0,0),"in"),
                      plot.tag = element_text(family = "sans", size = 9)
)


#Median-sized closure and 10 km-wide buffers for 2017-2019 in North-Central zone
{
  #Load closed areas. Created in make_closures_df.R
  load("Output/Data/closed.Rdata")
  
  #Filter to closures during study period
  closed <- filter(closed, !is.na(season) & 
                     (season=="s1_2017" | season=="s2_2017" | 
                        season=="s1_2018" | season=="s2_2018" | season=="s1_2019" | 
                        season=="s2_2019")
  )
  
  bbox <- c(-84, -15.9999, -74, -3.5) 
  names(bbox) <- names(st_bbox(closed))
  
  closed <- st_crop(closed, bbox)
  
  #Only plot closures declared by PRODUCE
  closed <- filter(closed, days <= 5)
  
  #To calculate areas, project closures
  medarea <- st_transform(closed, st_crs("+proj=laea +lon_0=-76.5"))
  
  #Area
  medarea <- mutate(medarea, area_km2 = st_area(medarea)/10^6) %>% 
    mutate(area_km2 = as.numeric(area_km2)) %>%
    mutate(areadif = abs(area_km2 - median(area_km2)))
  
  medclosed <- as.data.frame(medarea) %>% dplyr::select(areadif) %>% 
    summarise(min(areadif)) %>% as.matrix() %>% as.numeric()
  
  medclosed <- filter(medarea, areadif==medclosed)[1,]
  
  
  ##Create annuli
  #Function of buffer distance
  buF <- function(myb){
    #dist given in m,  
    buffered <- st_buffer(medclosed, dist = myb*1000) %>%
      #but I want bdist variable to be in km
      mutate(bdist = myb)
    
    return(buffered)
  }
  
  #Apply over buffer distances
  buffers <- lapply(seq(from=10,to=50,by=10), function(x){
    buF(x)
  })
  
  buffers <- do.call("rbind", buffers)
  
  #Bind with medclosed
  medclosed <- rbind(mutate(medclosed, bdist = 0), buffers)
  
  #Iteratively difference to create annuli
  annuli <- dplyr::select(medclosed, bdist)
  
  annuli <- st_transform(annuli, st_crs("+proj=longlat +datum=WGS84 +no_defs"))
  
  for(i in 5:1){
    annuli[i+1,] <- st_difference(annuli[i+1,],annuli[i,]) %>% dplyr::select(bdist)
  }
  
  closure_annuli <- annuli
  
  rm(buffers, closed, medarea, medclosed, i, annuli, bbox, buF)
}

#Created in 3. correct_be.R
load("Output/Data/pbe_imp.Rdata")

#Calculate average percentage juvenile and total juvenile catch in each treatmetn bin
#Function to use in apply that multiplies perjuv by given bin column
#Returns NA if bin column = 0, since that means set did not occur in that bin
#so don't want to use that perjuv value in calculating average perjuv catch in that bin
#Then calculate mean percentage juvenile for that bin
meanFun <- function(mybin){
  
  usepj <- dplyr::select(fullbe, bepjhat, numindivids, mybin, betons)
  
  usepj[,5] <- usepj[,1]*usepj[,2]*usepj[,3]
  
  #Weighted pj is NA if not in treatzone
  usepj[,5][usepj[,3]==0] <- NA 
  #Same for numindivids
  usepj[,2][usepj[,3]==0] <- NA
  #Same for tons
  usepj[,4][usepj[,3]==0] <- NA
  
  
  #Output the weighted mean of this new column
  outmean <- sum(usepj[,5],na.rm=T) / sum(usepj[,2],na.rm=T)
  
  #Output df of this mean, sum of juvenile catch, and name of bin
  out <- data.frame(bin = mybin, meanpj = outmean, juvcatch = sum(usepj[,5],na.rm=T),
                    tons = sum(usepj[,4],na.rm=T))
  
  return(out)
}

#Apply over bins
avgpj <- map_df(grep("lead|active|lag",names(fullbe),value=T),function(x){
  meanFun(x)
})


#Add bdist column
avgpj$bdist <- gsub(".*_","",avgpj$bin)
avgpj$bdist <- as.numeric(avgpj$bdist)

#And tvar column
avgpj$tvar <- gsub("_.*","",avgpj$bin)

#Join with annuli
plotdf <- left_join(avgpj, closure_annuli, by = 'bdist')

plotdf <- st_sf(plotdf)

#Plot meanpj given tvar
plotPJ <- function(mytvar){
  
  usedat <- filter(plotdf, tvar == mytvar)
  
  #Record maximum value of meanpj to set fill limit in plot
  maxval <- max(plotdf$meanpj)
  minval <- min(plotdf$meanpj)
  
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
  
  plot <- ggplot(data=usedat) + geom_sf(aes(fill=meanpj)) + 
    scale_fill_viridis("Percentage\njuvenile", limits = c(15,50),  
                       breaks = c(15, 20, 25, 30, 35, 40, 45, 50)) + 
    myThemeStuff + 
    ggtitle(tit) + 
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),legend.position = 'none',
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  
  return(plot)
}



##Also make a plot for paper that is 3x2
#Determine correct proportions for plot
myheight <- distGeo(c(-79,-10.787231),c(-79,-9.546358))
mywidth <- distGeo(c(-79.121430,-10.16679),c(-77.911960,-10.16679))

leadplot <- plotPJ("lead") + 
  theme(plot.margin = unit(c(.01,0,.3,0),"in"), 
        plot.title = element_text(hjust = 0.5, size = 11, family = 'sans'))
activeplot <- plotPJ("active") + 
  theme(plot.margin = unit(c(.01,0,.3,0),"in"), 
        plot.title = element_text(hjust = 0.5, size = 11, family = 'sans'))
lag1plot <- plotPJ("lag1")+ 
  theme(plot.margin = unit(c(.01,0,.3,0),"in"), 
        plot.title = element_text(hjust = 0.5, size = 11, family = 'sans'))
lag2plot <- plotPJ("lag2")+ 
  theme(plot.margin = unit(c(.01,0,.3,0),"in"), 
        plot.title = element_text(hjust = 0.5, size = 11, family = 'sans'))
lag3plot <- plotPJ("lag3")+ 
  theme(plot.margin = unit(c(.01,0,.3,0),"in"), 
        plot.title = element_text(hjust = 0.5, size = 11, family = 'sans'))
lag4plot <- plotPJ("lag4")+ 
  theme(plot.margin = unit(c(.01,0,.3,0),"in"), 
        plot.title = element_text(hjust = 0.5, size = 11, family = 'sans'))

tbt <- plot_grid(leadplot, activeplot, lag1plot, lag2plot, lag3plot, lag4plot, nrow=3, ncol=2)

legend_tbt <- get_legend(
  leadplot + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", 
          legend.key.width=unit(2,"cm"),
          legend.title = element_text(vjust = 1, size = 10),
          legend.text = element_text(size=10,family='sans'))
)

# add the legend underneath the row we made earlier. 
tbt <- plot_grid(tbt, legend_tbt, ncol = 1, rel_heights = c(1, .05))

ggsave(tbt, file="Output/Figures/figureC2.pdf",
       w=7,h=(7/2)*(myheight/mywidth)*3, units = "in", dpi=1200)