rm(list=ls())
setwd("C:/Users/englander/Documents/replication_closures")

library(dplyr); library(lubridate); library(sf)
library(lfe); library(xtable); library(Formula)

options(scipen=999)
options(lfe.threads=24)

#Load rddf created in 4. make_rddf.R
load("Output/Data/rddf_10km_lead1tolag4_3dayrect.Rdata")

rddf <- as.data.frame(rddf) %>% dplyr::select(-geometry)

rddf <- arrange(rddf, tvar, bdist)

#Millions of juveniles
rddf <- mutate(rddf, nummjuv = numjuv/10^6) %>%
  #Tons caught per set
  mutate(tonsperset = tons/nobs) %>%
  #Cluster tons/set and tons/area
  mutate(clusttonsperset = clusttons/clustnobs, clusttonsperarea = clusttons/clustarea_km2) %>%
  #clust millions of juveniles, individuals, and adults
  mutate(clustmjuv = clustnumjuv/10^6, clustmindivids = clustnumindivids/10^6,
         clustmadults = clustnumadults/10^6)

rddf <- mutate(rddf, asinhnummjuv = asinh(nummjuv), 
               asinhtons = asinh(tons))

rddf$bin <- as.factor(rddf$bin)
rddf$bin <- relevel(rddf$bin, ref="active_in")

rddf$twoweek_cellid_2p <- as.factor(rddf$twoweek_cellid_2p)
rddf$twowk <- as.factor(rddf$twowk)
rddf$cellid_2p <- as.factor(rddf$cellid_2p)

#Drop potential closures that have NA for size distribution
rddf <- filter(rddf, !is.na(prop12hat))

balancevars <- c("kmtocoast","clusttonsperset","clusttonsperarea")

#Try predicting outcome from balancevars and save fitted values
fitted <- lm(as.Formula(
  paste0("asinhnummjuv ~ ",paste0(balancevars,collapse='+'))),
  data=rddf)

fitted <- predict(fitted, rddf)

#Mutate onto rddf
rddf <- mutate(rddf, fitted = fitted)

#Add fitted onto balancevars
balancevars <- c(balancevars, "fitted")

#Correlation between outcome and covariate
outcome_cov <- function(cov){
  
  #Regress outcome on variable
  reg1 <- as.Formula(paste0("asinhnummjuv ~ ",cov,"| 0 | 0| twoweek_cellid_2p")) %>%
    felm(data=rddf)
  
  return(reg1)
}

#Correlation between covariate and treatment fraction
cov_treatfrac <- function(cov){
  
  #Regress variable on treatment fraction, fixed effects, and length distribution
  reg2 <- as.Formula(paste0(cov,"~ treatfrac + ",
                            paste0(grep("prop",names(rddf),value=T),collapse="+"),
                            "| bin + twowk:cellid_2p + startdate",
                            " | 0 | twoweek_cellid_2p")) %>% 
    felm(data=rddf)
  
  return(reg2)
}

#Format coefficient
formCoef <- function(reg, coef, dig){
  
  #Get coefficients from felm object
  mycoefs <- summary(reg)["coefficients"]$coefficients
  
  mycoef <- mycoefs[coef,"Estimate"]
  
  #Round
  roundcoef <- round(mycoef, dig) %>% as.character()
  
  #If rounded to integer, need to add "." to end
  if(length(grep("\\.",roundcoef))==0){
    roundcoef <- paste0(roundcoef, ".")
  }
  
  #Add an extra zero beyond the decimal point if needed to get same length
  #Do coef first
  roundcoef <- sapply(seq_len(length(roundcoef)), function(x){
    if(gsub(".*\\.","",roundcoef[x]) %>% nchar() < dig){
      #Needed length
      zerosneeded <- dig - gsub(".*\\.","",roundcoef[x]) %>% nchar()
      roundcoef[x] <- paste0(roundcoef[x],paste0(rep(0,zerosneeded),collapse=""))
    } else{
      roundcoef[x]
    }
  })
  
  #Add commas if necessary
  roundcoef <- prettyNum(roundcoef, ",")
  
  return(roundcoef)
}

#Given regression object, coefficient of interest, and number of digits to round to, 
#return formatted SE
formSE <- function(reg, coef, dig){
  
  #Get coefficients from felm object
  mycoefs <- summary(reg)["coefficients"]$coefficients
  
  #Get se
  se <- mycoefs[coef,"Cluster s.e."]
  
  #Round
  roundse <- round(se, dig) %>% as.character()
  
  #If rounded to integer, need to add "." to end
  if(length(grep("\\.",roundse))==0){
    roundse <- paste0(roundse, ".")
  }
  
  #Add zeros if necessary
  roundse <- sapply(seq_len(length(roundse)), function(x){
    if(gsub(".*\\.","",roundse[x]) %>% nchar() < dig){
      #Needed length (could need one extra zero or two)
      zerosneeded <- dig - gsub(".*\\.","",roundse[x]) %>% nchar()
      roundse[x] <- paste0(roundse[x],paste0(rep(0,zerosneeded),collapse=""))
    } else{
      roundse[x]
    }
  })
  
  #Add commas if necessary
  roundse <- prettyNum(roundse, ",")
  
  #Add parentheses 
  roundse <- paste0("(",roundse,")")
  
  # #Add stars if applicable
  # tstat <- mycoefs[coef,"Estimate"] / se
  # pval <- 2*pt(abs(tstat), df = summary(reg)["df"]$df[1], lower.tail = FALSE)
  # 
  # if(pval<.1 & pval>=.05){
  #   roundse <- paste0(roundse,"*")
  # } else if(pval<.05 & pval>=.01){
  #   roundse <- paste0(roundse,"**")
  # } else if(pval < .01){
  #   roundse <- paste0(roundse,"***")
  # }
  
  return(roundse)
}

#Table of correlation between outcome and covariate
covtab <- matrix(NA, nrow=9, ncol=5)

#Regress outcome on all three variables
fitreg <- felm(asinhnummjuv ~ kmtocoast + clusttonsperset + clusttonsperarea | 0 | 0 | twoweek_cellid_2p,
               data=rddf)

covtab[1,] <- c("","(1)","(2)","(3)","(4)")

covtab[2,] <- c("DistToCoast",
                formCoef(outcome_cov("kmtocoast"),"kmtocoast",4),
                "","",
                formCoef(fitreg,"kmtocoast",4)
)

covtab[3,] <- c("",
                formSE(outcome_cov("kmtocoast"),"kmtocoast",4),
                "","",
                formSE(fitreg,"kmtocoast",4)
)


covtab[4,] <- c("TonsPerSet","",
                formCoef(outcome_cov("clusttonsperset"),"clusttonsperset",4),
                "",
                formCoef(fitreg,"clusttonsperset",4)
)

covtab[5,] <- c("","",
                formSE(outcome_cov("clusttonsperset"),"clusttonsperset",4),
                "",
                formSE(fitreg,"clusttonsperset",4)
)

covtab[6,] <- c("TonsPerArea","","",
                formCoef(outcome_cov("clusttonsperarea"),"clusttonsperarea",4),
                formCoef(fitreg,"clusttonsperarea",4)
)

covtab[7,] <- c("","","",
                formSE(outcome_cov("clusttonsperarea"),"clusttonsperarea",4),
                formSE(fitreg,"clusttonsperarea",4)
)

covtab[8,] <- c("Intercept",
                formCoef(outcome_cov("kmtocoast"),"(Intercept)",4),
                formCoef(outcome_cov("clusttonsperset"),"(Intercept)",4),
                formCoef(outcome_cov("clusttonsperarea"),"(Intercept)",4),
                formCoef(fitreg,"(Intercept)",4)
)

covtab[9,] <- c("",
                formSE(outcome_cov("kmtocoast"),"(Intercept)",4),
                formSE(outcome_cov("clusttonsperset"),"(Intercept)",4),
                formSE(outcome_cov("clusttonsperarea"),"(Intercept)",4),
                formSE(fitreg,"(Intercept)",4)
)

myxtable <- xtable(covtab)

caption(myxtable) <- c("Correlation between juvenile catch and measures of fishing productivity")

align(myxtable) <- c("l","l",rep("c",4))

label(myxtable) <- "balance_relevantvariables"

print(myxtable, floating = TRUE, caption.placement="top",sanitize.text.function = identity,
      include.colnames=F,include.rownames=F,table.placement="tb",
      hline.after=NULL,
      add.to.row=list(
        pos = list(0,1,nrow(covtab)),
        command = c(
          paste0("\\toprule & \\multicolumn{4}{c}{Dependent variable: asinh(juvenile catch)} \\\\ "),
          "\\midrule ",
          "\\bottomrule \\multicolumn{5}{l}{\\multirow{2}{11cm}{All regressions have 34,164 observations. Dependent variable is inverse hyperbolic sine of millions of juveniles caught in a potential closure-treatment bin. Standard errors clustered at level of two-week-of-sample by two-degree grid cell.}} \\\\\\\\\\\\\\\\\\\\\\\\ "
        )),
      type = "latex",file="Output/Tables/tableA2.tex")




#Now make balance on observables table
baltab <- matrix(NA,ncol=5,nrow=4)

baltab[1,] <- c("","DistToCoast","TonsPerSet","TonsPerArea","FittedVals")

baltab[2,] <- c("","(1)","(2)","(3)","(4)")

baltab[3,] <- c("Treatment",
                formCoef(cov_treatfrac("kmtocoast"),"treatfrac",3),
                formCoef(cov_treatfrac("clusttonsperset"),"treatfrac",3),
                formCoef(cov_treatfrac("clusttonsperarea"),"treatfrac",3),
                formCoef(cov_treatfrac("fitted"),"treatfrac",3)
)

baltab[4,] <- c("fraction",
                formSE(cov_treatfrac("kmtocoast"),"treatfrac",3),
                formSE(cov_treatfrac("clusttonsperset"),"treatfrac",3),
                formSE(cov_treatfrac("clusttonsperarea"),"treatfrac",3),
                formSE(cov_treatfrac("fitted"),"treatfrac",3)
)

myxtable <- xtable(baltab)

caption(myxtable) <- c("Test for balance on measures of fishing productivity")

align(myxtable) <- c("l","l",rep("c",4))

label(myxtable) <- "balance"

print(myxtable, floating = TRUE, caption.placement="top",sanitize.text.function = identity,
      include.colnames=F,include.rownames=F,table.placement="tb",
      hline.after=NULL,
      add.to.row=list(
        pos = list(0,2,nrow(baltab)),
        command = c(
          paste0("\\toprule  "),
          "\\midrule ",
          "\\bottomrule \\multicolumn{5}{l}{\\multirow{2}{13cm}{All regressions have 34,164 observations and control for two-week-of-sample by two-degree-grid-cell fixed effects, day-of-sample fixed effects, and potential closure-level length distribution. Standard errors clustered at level of two-week-of-sample by two-degree grid cell.}} \\\\\\\\\\\\\\\\\\\\\\\\ "
        )),
      type = "latex",file="Output/Tables/tableA3.tex")

sessionInfo()
