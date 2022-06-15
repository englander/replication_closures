rm(list=ls())

library(dplyr); library(xtable); library(readr)
library(readxl)

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

options(scipen=999)

#Load ownership information
load("Data/owndf.Rdata")

#Large fleets are the top7 firms
owndf <- mutate(owndf, fleettype = as.character(NA))

owndf$fleettype[owndf$Armador %in% c("TECNOLOGICA DE ALIMENTOS S.A.","PESQUERA DIAMANTE S.A.","CORPORACION PESQUERA INCA S.A.C.",
                                     "PESQUERA EXALMAR S.A.A.", "CFG INVESTMENT S.A.C.", "AUSTRAL GROUP S.A.A", "PESQUERA HAYDUK S.A.")] <- "large"

#Medium fleets are not top7 firm and have more than one vessel
owndf$fleettype[owndf$Armador %not in% c("TECNOLOGICA DE ALIMENTOS S.A.","PESQUERA DIAMANTE S.A.","CORPORACION PESQUERA INCA S.A.C.",
                                         "PESQUERA EXALMAR S.A.A.", "CFG INVESTMENT S.A.C.", "AUSTRAL GROUP S.A.A", "PESQUERA HAYDUK S.A.") & 
                  owndf$numowned > 1] <- "medium" 

owndf$fleettype[owndf$numowned==1] <- "singleton"

#Get landings for each vessel each season
land <- read_xlsx("Data/landings_2017to2019.xlsx")

#What percent of landings are in South? 0.05197118
sum(land$TmDescargada[land$Zona=="Sur"]) / sum(land$TmDescargada)

#Drop Southern landings
land <- filter(land, Zona!="Sur")

#Calculate tons landed for each vessel each season
land <- group_by(land, Matricula, Temporada) %>% 
  summarise(tons = sum(TmDescargada)) %>% ungroup()

#Join onto owndf
owndf <- left_join(owndf, land, by = c("Matricula","Temporada"))

#Drop observations if no landings observations
owndf <- filter(owndf, !is.na(tons))

#What % of landings come from SUPNEP vessels? (Union contract mentioned in Section 2.1)
#TASA, AUSTRAL GROUP SAA, CANTABRIA SA, DIAMANTE SA, LOS  HALCONES y CARAL SA
#http://www.perupesquero.org/web/pescadores-industriales-recibiran-22-4-del-8-por-participacion-de-pesca/
supnep <- c("TECNOLOGICA DE ALIMENTOS S.A.","AUSTRAL GROUP S.A.A","PESQUERA CANTABRIA S.A.",
           "PESQUERA DIAMANTE S.A.","LOS HALCONES S.A.","PESQUERA CARAL S.A.")

#33%
filter(owndf, Armador %in% supnep) %>% 
  summarise(sum(tons)) %>% as.matrix() %>% as.numeric() / sum(land$tons) #0.3304142

#Vessel statistics independent of season
vesstats <- distinct(owndf, Matricula, casco, tipo_presevacion, eslora, potencia_motor, capbod_m3)

##Summary statistics for Section 2.3
#Percent by material: 40% steel and 60% wood
group_by(vesstats, casco) %>% 
  summarise(count = n()) %>% mutate(prop = count / nrow(vesstats))

vesstats$potencia_motor <- as.numeric(vesstats$potencia_motor)
vesstats$capbod_m3 <- as.numeric(vesstats$capbod_m3)

#Average length, motor and capacity by casco
group_by(vesstats, casco) %>% 
  summarise(eslora = mean(eslora,na.rm=T), 
            potencia_motor = mean(potencia_motor, na.rm=T), 
            capbod_m3 = mean(capbod_m3,na.rm=T))

#What % of steel vessels belong to non-singletons
filter(owndf, fleettype!='singleton' & casco=="ACERO NAVAL") %>% nrow() / 
  filter(owndf, casco=="ACERO NAVAL") %>% nrow()

#What % of wood vessels belong to non-singletons
filter(owndf, fleettype!='singleton' & casco=="MADERA") %>% nrow() / 
  filter(owndf, casco=="MADERA") %>% nrow()

#Make Table D1
tab <- matrix(NA, ncol=5,nrow=13)

tab[1,] <- c("","(1)","(2)","(3)","(4)")

tab[2,] <- c("Minimum", min(owndf$tons), min(owndf$tons[owndf$fleettype=="large"]), 
             min(owndf$tons[owndf$fleettype=="medium"]), 
             min(owndf$tons[owndf$fleettype=="singleton"]))

tab[3,] <- c("Mean", round(mean(owndf$tons),2), round(mean(owndf$tons[owndf$fleettype=="large"]),2), 
             round(mean(owndf$tons[owndf$fleettype=="medium"]), 2),
             round(mean(owndf$tons[owndf$fleettype=="singleton"]),2))

tab[4,] <- c("Median", round(median(owndf$tons),2), round(median(owndf$tons[owndf$fleettype=="large"]),2), 
             round(median(owndf$tons[owndf$fleettype=="medium"]), 2),
             round(median(owndf$tons[owndf$fleettype=="singleton"]),2))

tab[5,] <- c("Max", round(max(owndf$tons),2), round(max(owndf$tons[owndf$fleettype=="large"]),2), 
             round(max(owndf$tons[owndf$fleettype=="medium"]), 2),
             round(max(owndf$tons[owndf$fleettype=="singleton"]),2))

#Number of vessels in each fleettype each season
nves <- group_by(owndf, fleettype, Temporada) %>% 
  summarise(count = n()) %>% ungroup()

#Also add an "all" column: summed over fleettypes within season
nves <- bind_rows(nves, 
                  group_by(nves, Temporada) %>% 
                    summarise(count = sum(count)) %>% ungroup() %>% mutate(fleettype="all"))

tab[6,] <- c("Minimum", min(nves$count[nves$fleettype=="all"]), min(nves$count[nves$fleettype=="large"]),
              min(nves$count[nves$fleettype=="medium"]),
              min(nves$count[nves$fleettype=="singleton"]))

tab[7,] <- c("Mean", round(mean(nves$count[nves$fleettype=="all"]),2), round(mean(nves$count[nves$fleettype=="large"]),2),
              round(mean(nves$count[nves$fleettype=="medium"]), 2),
              round(mean(nves$count[nves$fleettype=="singleton"]),2))

tab[8,] <- c("Median", median(nves$count[nves$fleettype=="all"]), median(nves$count[nves$fleettype=="large"]),
              median(nves$count[nves$fleettype=="medium"]),
              median(nves$count[nves$fleettype=="singleton"]))

tab[9,] <- c("Maximum", max(nves$count[nves$fleettype=="all"]), max(nves$count[nves$fleettype=="large"]),
              max(nves$count[nves$fleettype=="medium"]),
              max(nves$count[nves$fleettype=="singleton"]))

tab[10,] <- c("Minimum", round(min(owndf$eslora),2), round(min(owndf$eslora[owndf$fleettype=="large"]),2), 
             round(min(owndf$eslora[owndf$fleettype=="medium"]), 2),
             round(min(owndf$eslora[owndf$fleettype=="singleton"]),2))

tab[11,] <- c("Mean", round(mean(owndf$eslora),2), round(mean(owndf$eslora[owndf$fleettype=="large"]),2), 
             round(mean(owndf$eslora[owndf$fleettype=="medium"]), 2),
             round(mean(owndf$eslora[owndf$fleettype=="singleton"]),2))

tab[12,] <- c("Median", round(median(owndf$eslora),2), round(median(owndf$eslora[owndf$fleettype=="large"]),2), 
              round(median(owndf$eslora[owndf$fleettype=="medium"]), 2),
              round(median(owndf$eslora[owndf$fleettype=="singleton"]),2))

tab[13,] <- c("Max", round(max(owndf$eslora),2), round(max(owndf$eslora[owndf$fleettype=="large"]),2), 
              round(max(owndf$eslora[owndf$fleettype=="medium"]), 2),
              round(max(owndf$eslora[owndf$fleettype=="singleton"]),2))


myxtable <- xtable(tab)

caption(myxtable) <- "Vessel characteristics in the six fishing seasons of 2017, 2018, and 2019"
  
align(myxtable) <- c("l","l",rep("c",4))

label(myxtable) <- "vesstats"

print(myxtable, floating = TRUE, caption.placement="top",sanitize.text.function = identity,
      include.colnames=F,include.rownames=F,table.placement="tb",
      hline.after=NULL,
      add.to.row=list(
        pos = list(0,1,5,9,nrow(tab)),
        command = c(
          #paste0(" & \\multirow{2}{*}{All vessels} & \\multirow{2}{*}{Large-firm vessels} & \\multirow{2}{*}{Medium-firm vessels} & \\multirow{2}{*}{Singleton vessels} \\\\ & & & & \\\\"),
          "\\toprule \  & All vessels & Large-firm vessels & Medium-firm vessels & Singleton vessels \\\\",
          "\\midrule \\multicolumn{5}{l}{A. Average tons landed per season} \\\\ ",
          "\\midrule \\multicolumn{5}{l}{B. Average number of active vessels per season} \\\\ ",
          "\\midrule \\multicolumn{5}{l}{C. Vessel length (m)} \\\\ ",
          "\\bottomrule  "
          )),
      type = "latex",file="Output/Tables/tableD1.tex")

#What percent of landings do each fleettype account for?
sum(owndf$tons[owndf$fleettype=="large"]) / sum(owndf$tons)

sum(owndf$tons[owndf$fleettype=="medium"]) / sum(owndf$tons)

sum(owndf$tons[owndf$fleettype=="singleton"]) / sum(owndf$tons)

sessionInfo()