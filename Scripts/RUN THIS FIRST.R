#This script ensures you are using the same package versions that I did

install.packages('devtools')

devtools::install_version("renv", version = "0.15.5", repos = "http://cran.us.r-project.org")

#If using Windows OS, set download method
if(Sys.info()["sysname"] == "Windows"){
  Sys.setenv(RENV_DOWNLOAD_METHOD = "wininet")
}

renv::restore()
#If R asks whether you want to "activate", enter Y into the console

#If renv::restore() fails, check whether Rtools is installed on your computer
#(it's usually installed on the C drive, e.g. C:/rtools...)
#If Rtools is not installed, install it by following the instructions here: 
#https://cran.r-project.org/bin/windows/Rtools/rtools40.html
#Then re-run renv::restore() in the console