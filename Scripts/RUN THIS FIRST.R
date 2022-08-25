#This script ensures you are using the same package version that I did
install.packages('devtools')

devtools::install_packages("renv", version = "0.15.2", 
                           repos = "http://cran.us.r-project.org")

renv::activate()

renv::restore()
