# Replication files for "Information and Spillovers from Targeting Policy in Peru's Anchoveta Fishery"

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7041707.svg)](https://doi.org/10.5281/zenodo.7041707)

## Data and Code Availability Statement

The paper uses public, non-confidential data from Peru's Ministry of Production (PRODUCE). The archive contains the data in the folder "Data/". The files from PRODUCE are BE_2017to2019_allvessels.xlsx, landings_2017to2019.xlsx, closures_2014to2019.xlsx, and owndf.csv. BE_2017to2019_allvessels.xlsx are the Bitácora Electrónica data for all vessels, landings_2017to2019.xlsx are the landings data for all vessels, closures_2014to2019.xslx are the temporary spatial closures that PRODUCE declared between 2014 and 2019, and owndf.csv are ownership information and vessel characteristics for all vessels.

The paper uses public, non-confidential data from Sociedad Nacional de Pesquería (SNP). The archive contains the data in the folder "Data/". The file from SNP is bedat_snp.csv. This file contains Bitácora Electrónica data for SNP vessels.

The paper uses public, non-confidential boundary data of Peru's Exclusive Economic Zone from the Flanders Marine Institute (Flanders Marine Institute, 2012). The archive contains the data in the folder "Data/Intersect_IHO_EEZ_v2_2012".

The paper uses public, non-confidential population-level length distribution data from Peru's scientific agency (IMARPE, 2017). The archive contains the data in the folder "Data/". The file from IMARPE is EvaluacionHidroacusticaRecPelagicosCrucero170304_march2017.jpg. This file is page 23 from their public report (IMARPE, 2017). 


## Computational Requirements

-Software: R. I used Version 4.1.2, but other versions should work too, especially those >= 4.1.0.

-Packages: There are many of them. They are all recorded in renv.lock file. When you run Scripts/RUN THIS FIRST.R, the renv package will automatically install all of them.

-OS: I used Windows 10. Other versions of Windows, as well as Mac and Linux, should work too. 

-CPU: I have Intel(R) Xeon(R) Gold 6132 CPUE @ 2.60GHz 2.60 GHz (2 processors). This is the equivalent of 16 cores. If you have fewer than 16 cores, you will still be able to run the scripts, but your runtime will be longer.

-Memory: 128 GB. You can run most scripts with as little as 4 GB of memory though. A few scripts will require 128 GB of memory though, such as Scripts/make_figures/make_figureA11.R. 

-Necessary disk space: 2 GB

-Wall clock-time: 920 hours. 900 hours for Scripts/make_data/7. make_data_figA13.R, and 20 hours for all other scripts. If you don't want to run Scripts/make_data/7. make_data_figA13.R, you may skip to make_figureA13.R since I provide the data necessary for creating Figure A13 in Output/Data/data_figA13.Rdata. 

## Instructions for Data Preparation and Analysis

### Data preparation

First, run Scripts/RUN THIS FIRST.R. That script will install all R packages you need. It installs the same package versions I used to facilitate reproducibility. 

Then, open Scripts/make_all_data.R. If you skip the last line, source("Scripts/make_data/7. make_data_figA13.R"), it will run much faster. You can see the sub-scripts that make_all_data.R runs in Scripts/make_data folder. 

### Analysis

Scripts/make_all_figures.R creates all figures in the paper. Scripts/make_tables.R creates all tables in the paper. Files in Scripts/other_empirics folder contain calculations that are described in the paper but which do not produce a table or figure.

I provide output from Scripts/make_all_data.R in Output/Data folder. So if you want to skip the data preparation stage, and go right to reproducing tables and figures, you may do so.
