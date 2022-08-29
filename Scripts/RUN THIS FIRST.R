#This script ensures you are using the same package versions that I did

#If using Windows OS, set download method
if(Sys.info()["sysname"] == "Windows"){
  Sys.setenv(RENV_DOWNLOAD_METHOD = "wininet")
}
renv::restore()