##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA
#  DESCRIPTION: This script allows saving the list of run scenarios into a csv file "list_simulations" located in the folder "22_results_aggr"
###########################################################################################################################################################################################################

#Note: Before running the script it's recommended to create a copy of the file "TEMPLATE_list_simulations.csv" and rename it to "list_simulations.csv"
#      This will guarantee that if the script is run multiple times, there will be no left-over data from the previous runs in the csv file

rm(list=ls())

#--- Paths to folders --------------
path_wd<-getwd()
path_data<-file.path(path_wd, "3_data")                             #Folder "3_data"
path_lists<-file.path(path_data, "00_lists")                        #Folder "3_data\00_lists"
path_res_aggr<-file.path(path_data, "22_results_aggr")              #Folder "3_data\22_results_aggr"
path_scripts_app<<-file.path(path_wd, "2_scripts_app")              #Script for Shiny app (R)

#--- Load the list of simulations --
source(file.path(path_scripts_app,"1_common.R"))                   #Source the script with the function Load_ListR
list_sim<-Load_ListR(path_lists, file_name="list_simulations")    #R object containing the list of all run simulations (may contain multiple rows per scenario if multiple simulations were run) 

#--- Save as csv -------------------
write.csv(list_sim, file=file.path(path_res_aggr, "list_simulations.csv"),row.names=FALSE) 


