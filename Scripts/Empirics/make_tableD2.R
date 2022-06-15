#Calculate bepj- bepjhat for the three types of vessels: 
#singleton, medium firm, and large firm

rm(list=ls())

library(dplyr)
library(lubridate); library(xtable)

options(scipen=999)

#Peru time
Sys.setenv(TZ='America/Lima')

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

#Given number and number of digits, round and add zero after decimal if necessary
formNum <- function(num, dig){
  
  #Round
  roundnum <- round(num, dig) %>% as.character()
  
  #If rounded to integer, need to add "." to end
  if(length(grep("\\.",roundnum))==0){
    roundnum <- paste0(roundnum, ".")
  }
  
  #Add an extra zero beyond the decimal point if needed to get same length
  #Do num first
  roundnum <- sapply(seq_len(length(roundnum)), function(x){
    if(gsub(".*\\.","",roundnum[x]) %>% nchar() < dig){
      #Needed length
      zerosneeded <- dig - gsub(".*\\.","",roundnum[x]) %>% nchar()
      roundnum[x] <- paste0(roundnum[x],paste0(rep(0,zerosneeded),collapse=""))
    } else{
      roundnum[x]
    }
  })
  
  #Add commas if necessary
  roundnum <- prettyNum(roundnum, ",")
  
  return(roundnum)
}

#Created in correct_be*.R
load("Output/Data/pbe_imp.Rdata")

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

#Calculate individuals-weighted reported pj and corrected pj
pjdif <- left_join(fullbe, dplyr::select(owndf, Matricula, Temporada, fleettype)) %>%  
  filter(!is.na(numindivids)) %>%
  mutate(bepj_weighted = bepj*numindivids, 
         bepjhat_weighted = bepjhat * numindivids) %>% 
  group_by(fleettype) %>% 
  summarise(bepj = sum(bepj_weighted) / sum(numindivids), 
            bepjhat = sum(bepjhat_weighted) / sum(numindivids)) %>% 
  ungroup()

tabdf <- mutate(pjdif, pjdif = bepj - bepjhat)

#Make table 
tab <- matrix(NA, ncol = 4, nrow = 4)

tab[1, ] <- c("Vessel type", "Reported \\% juvenile", "Corrected \\% juvenile", 
              "$\\Delta$ \\% juvenile")

tab[2:4, 1] <- c("Large-firm vessels", "Medium-firm vessels", "Singleton vessels") 

tab[2:4, 2] <- sapply(tabdf$bepj, function(x){
  formNum(x, 1)
})

tab[2:4, 3] <- sapply(tabdf$bepjhat, function(x){
  formNum(x, 1)
})

tab[2:4, 4] <- sapply(tabdf$pjdif, function(x){
  formNum(x, 1)
})

myxtable <- xtable(tab)

caption(myxtable) <- c("Large-firm vessels underreport percentage juvenile more")

align(myxtable) <- c("l","l",rep("c",3))

label(myxtable) <- "large_underreport"

print(myxtable, floating = TRUE, caption.placement="top",sanitize.text.function = identity,
      include.colnames=FALSE,include.rownames=FALSE,table.placement="tb",
      hline.after=NULL,
      add.to.row=list(
        pos = list(0,1,nrow(tab)),
        command = c(
          paste0("\\toprule "),
          "\\midrule ",
          "\\bottomrule "
        )),
      type = "latex", file = "Output/Tables/tableD2.tex")


