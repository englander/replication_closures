#This script ensures you are using the same package version that I did

#If using Windows OS, set download method
if(Sys.info()["sysname"] == "Windows"){
  Sys.setenv(RENV_DOWNLOAD_METHOD = "wininet")
}
renv::restore()