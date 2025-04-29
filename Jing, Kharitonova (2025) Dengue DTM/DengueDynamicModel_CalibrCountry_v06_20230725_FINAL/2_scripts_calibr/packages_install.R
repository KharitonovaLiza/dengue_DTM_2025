##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA
#  DESCRIPTION: This script installs required packages (if they are not installed yet)
###########################################################################################################################################################################################################

.libPaths(.libPaths())                                    #Destination folder for installed packages (R installation)
options(repos="http://cran.us.r-project.org")             #Source repository for package downloads

if (!("devtools" %in% rownames(installed.packages()))){install.packages("devtools")}  
if (!("shiny" %in% rownames(installed.packages()))){install.packages("shiny")}
if (!("shinythemes" %in% rownames(installed.packages()))){install.packages("shinythemes")}
devtools::install_github("dreamRs/shinyWidgets")                                          
if (!("DT" %in% rownames(installed.packages()))){install.packages("DT")}                  
if (!("openxlsx" %in% rownames(installed.packages()))){install.packages("openxlsx")}    
if (!("xlsx" %in% rownames(installed.packages()))){install.packages("xlsx")}              
if (!("deSolve" %in% rownames(installed.packages()))){install.packages("deSolve")}        
if (!("ggplot2" %in% rownames(installed.packages()))){install.packages("ggplot2")}            
if (!("reshape2" %in% rownames(installed.packages()))){install.packages("reshape2")}          
if (!("rhandsontable" %in% rownames(installed.packages()))){install.packages("rhandsontable")}          
if (!("plotly" %in% rownames(installed.packages()))){install.packages("plotly")}                   
if (!("png" %in% rownames(installed.packages()))){install.packages("png")}          
if (!("shinytoastr" %in% rownames(installed.packages()))){install.packages("shinytoastr")}          
if (!("shinyjs" %in% rownames(installed.packages()))){install.packages("shinyjs")}          
if (!("htmltools" %in% rownames(installed.packages()))){install.packages("htmltools")}          
if (!("shinyBS" %in% rownames(installed.packages()))){install.packages("shinyBS")}          
if (!("bsplus" %in% rownames(installed.packages()))){install.packages("bsplus")}          
if (!("extrafont" %in% rownames(installed.packages()))){install.packages("extrafont")}          
if (!("scales" %in% rownames(installed.packages()))){install.packages("scales")}       
if (!("Rmisc" %in% rownames(installed.packages()))){install.packages("Rmisc")}          
if (!("stringr" %in% rownames(installed.packages()))){install.packages("stringr")}          
if (!("dplyr" %in% rownames(installed.packages()))){install.packages("dplyr")}          
if (!("EnvStats" %in% rownames(installed.packages()))){install.packages("EnvStats")}
if (!("shinycssloaders" %in% rownames(installed.packages()))){install.packages("shinycssloaders")}
if (!("knitr" %in% rownames(installed.packages()))){install.packages("knitr")}
if (!("kableExtra" %in% rownames(installed.packages()))){install.packages("kableExtra")}
if (!("whoami" %in% rownames(installed.packages()))){install.packages("whoami")}
if (!("reshape" %in% rownames(installed.packages()))){install.packages("reshape")}
if (!("testthat" %in% rownames(installed.packages()))){install.packages("testthat")}