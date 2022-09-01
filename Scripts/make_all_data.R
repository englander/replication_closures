#Remember to run RUN THIS FIRST.R first
#Once you've done that, you can start to run this script

#Create a spatial object with boundaries of actual closures declared by regulator 
#also includes their start and end times
source("Scripts/make_data/0. make_closures_df.R")

#Match PRODUCE Bitacora Electronica to PRODUCE landings data
source("Scripts/make_data/1. match_be_landings.R")

#Impute length distribution and calculate number of juveniles caught by each set
source("Scripts/make_data/2. impute_size_be.R")

#Correct percentage juvenile and length distribution
source("Scripts/make_data/3. correct_be.R")

#Create potential closures. 
#After running this script, can make all main text tables, 
#and Figures 1-7
source("Scripts/make_data/4. make_rddf.R")

#Create data required for creating Figures 8 and 9
source("Scripts/make_data/5. make_fleetthere_selfthere.R")

#Create data required for creating Figures 12 and 13
source("Scripts/make_data/6. make_actualclosure_regressioncontrol.R")

#Create data required for creating Figure 13
#Warning: This script requires approximately 900 CPU hours to run
#Instead of running this script, you may skip to make_figureA13.R
#and use the output data provided by this script
source("Scripts/make_data/7. make_data_figA13.R")