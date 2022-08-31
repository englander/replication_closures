#Make a figure representing Nash equilibria

rm(list=ls())

library(ggplot2); library(dplyr)
library(latex2exp); library(tidyr)

#Set parameters: 
mu_g <- 1; mu_k <- 1; mu_h <- 1.5
alpha <- 1; alpha_h <- 1.75
I <- 1

#Values for the three lines
e_g <- function(I_g){
  
  mu_g - alpha*(I_g)
  
}

e_h <- function(I_g){
  
  mu_h - alpha_h*I_g
  
}


e_k <- function(I_g){
  
  mu_k - alpha * I + alpha * I_g
  
}

df <- seq(from=0,to=1,by=.01) %>% as.data.frame() %>% as_tibble()
names(df) <- "I_g"

df <- mutate(df, k = e_k(I_g), g = e_g(I_g), h = e_h(I_g))

#Don't plot h all the way to x-axis
df$h[df$h < .1] <- as.numeric(NA)

#Gather
df <- gather(df, key = 'location', value = 'e_pi', -I_g)

#Drop NA
df <- filter(df, !is.na(e_pi))

df$location <- as.factor(df$location)


#Plot
(plot <- ggplot() + 
    geom_line(data=df, aes(x=I_g, y=e_pi, col=location)) + 
    theme(panel.background = element_rect(fill=NA),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(color = "black", size = 12, family="sans"),
          axis.title = element_text(color = "black", size = 12, family = "sans"),
          plot.margin = unit(c(0.05,0.15,.05,.1),"in"),
          legend.position='none') + 
    scale_x_continuous(TeX("Number of other vessels choosing $g$ or $h$ $(\\Iota_{-i, g}$ or $\\Iota_{-i, h})$"),
                       breaks = c(0, 1), 
                       labels = c(TeX("$\\Iota_k = \\Iota$"), TeX("$\\Iota_g$ or $\\Iota_h = \\Iota$"))) + 
    scale_y_continuous(TeX("Vessel i's expected profit from fishing in given location"),
                       breaks = c(mu_g, mu_h),
                       labels = c(TeX("$\\tilde{\\mu_g}$"), TeX("$\\tilde{\\mu_{h}}$"))) + 
    annotate(geom='text', x=I-.05, y=mu_k + .08, label=TeX("$k$", output='character'), parse=TRUE) +
    annotate(geom='text', x=0.025, y=mu_g + .08, label=TeX("$g$", output='character'), parse=TRUE) +
    annotate(geom='text', x=0.025, y=mu_h + .1, label=TeX("$h$", output='character'), parse=TRUE, col='red') +
    annotate(geom='point', x=(mu_g - mu_k)/2*alpha + .5*I,
             y=e_g((mu_g - mu_k)/2*alpha + .5*I), color='black') + 
    annotate(geom='text', x=(mu_g - mu_k)/2*alpha + .5*I,
             y=e_g((mu_g - mu_k)/2*alpha + .5*I) - .1, color='black',
             label = TeX("($\\Iota_g^*$, $\\Iota_k^*$)")) + 
    annotate(geom='point', x=(mu_h - mu_k)/(alpha_h + alpha) + (alpha/(alpha_h + alpha))*I,
             y=e_h((mu_h - mu_k)/(alpha_h + alpha) + (alpha/(alpha_h + alpha))*I), color='red') + 
    annotate(geom='text', x=(mu_h - mu_k)/(alpha_h + alpha) + (alpha/(alpha_h + alpha))*I,
             y=e_h((mu_h - mu_k)/(alpha_h + alpha) + (alpha/(alpha_h + alpha))*I) + .13, color='red',
             label = TeX("($\\Iota_h^*$, $\\Iota_k^*$)")) + 
    scale_color_manual("",values=c("black","red","black"))
)




ggsave(plot=plot, filename = "Output/Figures/figureB2.pdf",
       units = "in", width = 7, h=(7/1.5), dpi = 900)


#Also plot vessels choices when C=0 and when C=1
rm(list=ls())

library(ggplot2); library(rworldmap); library(cowplot)
library(sf); library(dplyr); library(geosphere)

#Median-sized potential closure
#Load potential closures created in 4. make_rddf*.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

rddf <- filter(rddf, bin=="active_in")

#To calculate areas, project
medarea <- st_transform(rddf, st_crs("+proj=laea +lon_0=-76.5"))

#Area
medarea <- mutate(medarea, area_km2 = st_area(medarea)/10^6) %>% 
  mutate(area_km2 = as.numeric(area_km2)) %>%
  mutate(areadif = abs(area_km2 - median(area_km2)))

medclosed <- as.data.frame(medarea) %>% dplyr::select(areadif) %>% 
  summarise(min(areadif)) %>% as.matrix() %>% as.numeric()

set.seed(20200422)

medclosed <- filter(medarea, areadif==medclosed) %>% 
  sample_n(1)

#Create 20 km buffers
buffer <- st_buffer(medclosed, dist = 20000)

#Subtract medclosed from buffer
h <- st_difference(buffer, medclosed)

#Transform to lonlat 
h <- st_transform(h, st_crs(rddf))
g <- st_transform(buffer, st_crs(rddf))

#Shift g to right for plotting
k <- st_geometry(g)
k <- k + c(.75, 0)
k <- st_sf(k, crs = st_crs(g))

#Labels
g <- mutate(g, lab = "g")
h <- mutate(h, lab = 'h')
k <- mutate(k, lab = 'k')



rm(medarea, rddf, medclosed, buffer)

(c0 <- ggplot() + 
    geom_sf(data=g,fill='dodgerblue2',col='dodgerblue2',alpha=.25) + 
    geom_sf(data=k,fill='purple',col='purple',alpha=.25) + 
    geom_sf_text(data=g, aes(label=lab),size=10,col='dodgerblue2') + 
    geom_sf_text(data=k, aes(label=lab),size=10,col='purple') + 
    scale_x_continuous("",limits = c(st_bbox(g)["xmin"], st_bbox(k)["xmax"])) +
    ggtitle("When C=0 (no closure), \nvessels fish in g or k") + 
  theme(panel.background = element_blank(),
    axis.title = element_blank(),
    #text = element_text(color = "black", size = 20, family="sans"),
    panel.border = element_rect(color='black',fill=NA),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.margin = unit(c(0,.25,0,0),"in"))
)

(c1 <- ggplot() + 
    geom_sf(data=h,fill='dodgerblue2',col='dodgerblue2',alpha=.25) + 
    geom_sf(data=k,fill='purple',col='purple',alpha=.25) + 
    geom_sf_text(data=h, aes(label=lab),size=10,col='dodgerblue2') + 
    geom_sf_text(data=k, aes(label=lab),size=10,col='purple') + 
    scale_x_continuous("",limits = c(st_bbox(g)["xmin"], st_bbox(k)["xmax"])) + 
    ggtitle("When C=1 (closure), \nvessels fish in h or k") + 
    theme(panel.background = element_blank(),
          axis.title = element_blank(),
          #text = element_text(color = "black", size = 20, family="sans"),
          panel.border = element_rect(color='black',fill=NA),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14),
          plot.margin = unit(c(0,0,0,.25),"in"))
)

#Plot them side by side
tb1 <- plot_grid(c0, c1, ncol=2, nrow=1)

ggsave(tb1, file="Output/Figures/figureB1.pdf",
       w=7,h=2.5, units = "in", dpi=1200)
