#This script ensures you are using the same package versions that I did

install.packages('devtools')

devtools::install_version("renv", version = "0.15.5", repos = "http://cran.us.r-project.org")

#If using Windows OS, set download method
if(Sys.info()["sysname"] == "Windows"){
  Sys.setenv(RENV_DOWNLOAD_METHOD = "wininet")
}

renv::restore()
#If R asks whether you want to "activate", enter Y into the console