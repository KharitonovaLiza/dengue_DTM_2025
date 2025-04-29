##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  SCRIPT NAME: compile_C_code.R
#  AUTHORS: Olivier CRISTEAU, Anna Tytula, Aurelien JAMOTTE
#  LIST OF FUNCTIONS:
#  - CompileModelFromCcode
###########################################################################################################################################################################################################

#Notes:
#if you have: In system(cmd) : 'make' not found - before running this function run code:
#	install.packages("devtools")
#	library(devtools)
# Rtools should be installed, check https://cran.r-project.org/bin/windows/Rtools/ and information on your version, for example for Rtools 4.0 it is needed to put Rtools on the PATH 
#	Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) 
#	Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")
# write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)

#This function allows compiling the C script into a shared library (dll) file. This function is only run when the C code is modified
#It does not need to be run at each model run

CompileModelFromCcode<-function(path_src, 
                                script_name){

  #Get the name of the C file
  name_c_file<-paste(script_name, ".c", sep="")
  
  #Create a share library file (.dll on Windows)
  system(paste("R  CMD SHLIB ", file.path(path_src,name_c_file), sep="")) 
  
  #Get the name of the dll file
  name_dll_file<-paste(script_name, .Platform$dynlib.ext, sep="")
  
  #Remove the .o file that was created so that it does not interfer with future compilations
  file.remove(file.path(path_src,paste(script_name,".o", sep=""))) 
  return(0)
}
