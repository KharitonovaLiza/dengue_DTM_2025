##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA
#  DESCRIPTION: This script allows exporting into a csv file the results for one selected scenario, including
#               - Nb of infections over time, by selected variable
#               - Population size over time, for selected vac status(es)
###########################################################################################################################################################################################################

rm(list=ls())

#=========================================================================
#--- Parameters to provide -----------------------------------------------
#=========================================================================

#Name of the scenario to export (as Scen_XXXXX)
scen_to_export<-"Scen_00005"       

#---------
#Name for the results file(s)
file_name<-"R11_CU12_31_PopSize_Vac"

#---------
#Export population size by:
#1-overall; 2-by age; 3-by vac status
group_by<-2

#Include population in:
incl_age<-c(1:101)           #Years of age to incl (min 0, max 100); +1 to correctly indicate the position in the vector
incl_vac_status<-c(2,3)    #Vac statuses to include (1-unvaccinated, 2-vaccinated as seroneg, 3-vaccinated as seropos)


#=========================================================================
#--- Code that exports the results ---------------------------------------
#=========================================================================

#-------------------------------------------------------
#--- Paths to folders ----------------------------------
#-------------------------------------------------------
path_wd<-getwd()
path_data<-file.path(path_wd, "3_data")                              #Folder "3_data"
path_lists<-file.path(path_data, "00_lists")                         #Folder "3_data\00_lists"
path_results<-file.path(path_data, "21_results")                     #Folder "3_data\21_results"
path_results_aggr<-file.path(path_data, "22_results_aggr")           #Folder "3_data\22_results_aggr"
path_scripts_model<<-file.path(path_wd, "1_scripts_model")           #Scripts for transmission model (R)
path_scripts_app<-file.path(path_wd, "2_scripts_app")                #Script for Shiny app (R)
path_inputs_epi<<-file.path(path_data, "11_inputs_epi")              #Population & transmission inputs
path_inputs_cost_vac<<-file.path(path_data, "16_inputs_cost_vac")    #Vaccination cost inputs
path_inputs_QoL<<-file.path(path_data, "18_inputs_QoL")              #Quality-of-life inputs

#-------------------------------------------------------
#--- Source the required functions ---------------------
#-------------------------------------------------------
debugSource(file.path(path_scripts_model, "11_parms_common.R"))  
debugSource(file.path(path_scripts_app,"1_common.R"))               
debugSource(file.path(path_scripts_app,"6_results_functions.R"))

#-------------------------------------------------------
#--- Load list of simulation ---------------------------
#-------------------------------------------------------
list_sim<-Load_ListR(path_lists, file_name="list_simulations")     #R object containing the list of all run simulations (may contain multiple rows per scenario if multiple simulations were run) 

#-------------------------------------------------------
#--- Full list of scenarios ----------------------------
#-------------------------------------------------------
#This object is required for the proper functioning of the function ResultsPage_DataEpi and ResultsPage_DataEco below
nb_scen<<-length(unique(list_sim[,1]))                                        #Number of unique scenarios
if(nb_scen==1){
  list_scen<-matrix(NA, nrow=nb_scen, ncol=16)                                #List of all scenarios (one line per scenario); removing the columns for start year and timeframe (simulation-specific)
  list_scen[1,]<-list_sim[match(unique(list_sim[,1]),list_sim[,1]),-c(7,8)] 
  colnames(list_scen)<-colnames(list_sim)[-c(7,8)]
  list_scen<<-list_scen
}else{
  list_scen<<-list_sim[match(unique(list_sim[,1]),list_sim[,1]),-c(7,8)]        
}

#-------------------------------------------------------
#--- Get the list of start years & timeframe -----------
#-------------------------------------------------------
start_years<-as.integer(list_sim[which(list_sim[,"scen_ref"]==scen_to_export), "year_start"]) 
timeframe  <-as.integer(list_sim[which(list_sim[,"scen_ref"]==scen_to_export), "timeframe"])[1] #Only taking the first element (assuming that all simulations for this scenario were run with the same timeframe)
vac_incl   <-as.logical(list_sim[which(list_sim[,"scen_ref"]==scen_to_export), "vac_incl"])[1]

message("Nb start years: ", length(start_years))
message("Timeframe: ", timeframe, " years")

#-------------------------------------------------------
#--- Initialize the objects for the results ------------
#-------------------------------------------------------
#Define the nb of columns in the final object
nb_col<-NA
if(group_by==1){nb_col<-1}                                 #Overall pop counts --> 1 column 
if(group_by==2){nb_col<-101}                               #Pop counts by age --> 101 columns
if(group_by==3){if(vac_incl){nb_col<-3}else{nb_col<-1}}    #Pop counts by vac status --> 3 columns  

#Initialize the final object
pop_host_2D<-matrix(NA, nrow=timeframe, ncol=nb_col)

#Labels column (when needed)
if(group_by==3){
  if(vac_incl){
    colnames(pop_host_2D)<-c("Unvac.", "Vac. as seroneg.", "Vac. as seropos.")  
  }else{
    colnames(pop_host_2D)<-c("Unvac.")
  }
}


#-------------------------------------------------------
#--- Prepare the result data for the scenario ----------
#-------------------------------------------------------

#Epi data contains:
#- A 7D object with the nb of infection
#- A 4D object with the nb of newly vaccinated/screened individuals

epi_data<-ResultsPage_DataEpi(scen_ref=scen_to_export,   #Scenario reference (e.g. "Scen_00001")
                              start_years=start_years,   #Vector with start years (years of vaccine introduction of corresponding periods in the scenarios w/o vaccination)  
                              timeframe_max=timeframe)   #Timeframe

#--- Population size ---------
#Population size is a 4D array with the population counts by year, age, vac status and start year
#This 4D object is generated by aggregating the 3D result objects for each individual simulation (i.e. start year)

#Extract the 4D object
pop_host_4D<-epi_data$pop_host_4D$pop_host

#Initialize a new 3D object which will contain average results across all simulations
pop_host_3D<-array(NA, dim=dim(pop_host_4D)[1:3])

#Calculate the average results across all simulations (same logic as for the nb of infections above)
pop_host_3D<-apply(pop_host_4D, c(1:3), sum)/(dim(pop_host_4D)[4])

#-------------------------------------------------------
#--- Prepare the final object for the pop counts -------
#-------------------------------------------------------
if(group_by==1){pop_host_2D[]<-apply(pop_host_3D[,incl_age,incl_vac_status],c(1)  ,sum)}    #Overall pop counts --> 1 column 
if(group_by==2){pop_host_2D[]<-apply(pop_host_3D[,        ,incl_vac_status],c(1,2),sum)}    #Pop counts by age --> 101 columns
if(group_by==3){pop_host_2D[]<-apply(pop_host_3D[,incl_age,               ],c(1,3),sum)}    #Pop counts by vac status --> 3 columns  


#-------------------------------------------------------
#--- Save results as csv -------------------------------
#-------------------------------------------------------
write.csv(pop_host_2D, file=file.path(path_results_aggr, paste(file_name, ".csv", sep="")),row.names=FALSE) 