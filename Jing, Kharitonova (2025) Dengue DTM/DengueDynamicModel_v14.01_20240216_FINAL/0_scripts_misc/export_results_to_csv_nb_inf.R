##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA
#  DESCRIPTION: This script allows exporting into a csv file the results for one selected scenario, including
#               - Nb of infections over time, by selected variable
###########################################################################################################################################################################################################

rm(list=ls())

#=========================================================================
#--- Parameters to provide -----------------------------------------------
#=========================================================================

#Name of the scenario to export (as Scen_XXXXX)
scen_to_export<-"Scen_00002"       

#---------
#Name for the results file(s)
file_name<-"R11_CU12_31_ZeroEff_NbInfSympt_Unvac6"

#---------
#Export nb of infections by:
#1-overall; 2-by age; 3-by serotype; 4-by inf type (primary, secondary, post-secondary); 5-by vac status; 6-by severity
group_by<-1

#Include infections in:
incl_age<-c(1:101)              #Years of age to incl (min 0, max 100); +1 to correctly indicate the position in the vector
incl_serotype<-c(1,2,3,4)       #Serotypes to include
incl_type<-c(1,2,3)             #Types to include (1-primary, 2-secondary, 3-post-secondary)
incl_vac_status<-c(1,2,3)       #Vac statuses to include (1-unvaccinated, 2-vaccinated as seroneg, 3-vaccinated as seropos)
incl_severity<-c(1,2,3,4,5,6)   #Severity to include:
#                               # 1-Asymptomatic                                  
#                               # 2-Symptomatic, mild, non-hospitalised
#                               # 3-Symptomatic, mild, hospitalised  
#                               # 4-Symptomatic, severe, non-hospitalised
#                               # 5-Symptomatic, severe, hospitalised, non-fatal
#                               # 6-Symptomatic, severe, hospitalised, fatal

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
#Labels for objects columns
labels_serotype<-c("DENV-1", "DENV-2", "DENV-3", "DENV-4")
labels_type<-c("Primary", "Secondary", "Post-secondary")

if(vac_incl){labels_vac_status<-c("Unvac.", "Vac. as seroneg.", "Vac. as seropos.")}else{labels_vac_status<-c("Unvac.")}

labels_severity<-c("Asympt",
                   "Sympt mild non-hosp",
                   "Sympt mild hosp",
                   "Sympt severe non-hosp",
                   "Sympt severe hosp non-fatal",
                   "Sympt severe hosp fatal")

#---------
#Object for the nb of infections
#Define the nb of columns in the final object
nb_col<-NA
if(group_by==1){nb_col<-1}                                  #Overall nb infections --> 1 column
if(group_by==2){nb_col<-101}                                #Nb infections by age --> 101 columns
if(group_by==3){nb_col<-4}                                  #Nb infections by serotype --> 4 columns
if(group_by==4){nb_col<-3}                                  #Nb infections by type --> 3 columns
if(group_by==5){if(vac_incl){nb_col<-3}else{nb_col<-1}}     #Nb infections by vac status --> 3 columns
if(group_by==6){nb_col<-6}                                  #Nb infections by severity --> 6 columns

#Initialize the final object
nb_inf_2D<-matrix(NA, nrow=timeframe, ncol=nb_col)

#Label columns (when needed)
if(group_by==3){colnames(nb_inf_2D)<-labels_serotype}    #If by serotype
if(group_by==4){colnames(nb_inf_2D)<-labels_type}        #If by type
if(group_by==5){colnames(nb_inf_2D)<-labels_vac_status}  #If by vac status
if(group_by==6){colnames(nb_inf_2D)<-labels_severity}    #If by severity

#-------------------------------------------------------
#--- Prepare the result data for the scenario ----------
#-------------------------------------------------------

#Epi data contains:
#- A 7D object with the nb of infection
#- A 4D object with the nb of newly vaccinated/screened individuals

epi_data<-ResultsPage_DataEpi(scen_ref=scen_to_export,   #Scenario reference (e.g. "Scen_00001")
                              start_years=start_years,   #Vector with start years (years of vaccine introduction of corresponding periods in the scenarios w/o vaccination)  
                              timeframe_max=timeframe)   #Timeframe


#--- Nb of infections --------
#Nb of infections is a 7D array with the nb of inf by year, age, serotype, inf type, vac status, severity (see below) and start year (i.e. year of vac introduction)
#This 7D object is generated by aggregating the 6D result objects for each individual simulation (i.e. each start year)

#Extract the 7D object with nb inf
nb_inf_7D<-epi_data$nb_inf_7D

#Initialize a new 6D object which will contain average results across all simulations
nb_inf_6D<-array(NA, dim=dim(nb_inf_7D)[1:6])

#Calculate the average results across all simulations (start years)
nb_inf_6D<-apply(nb_inf_7D,c(1:6),sum)/(dim(nb_inf_7D)[7])   #Summing for each level of each dimension and dividing by the length of the last dimension to obtain average across the simulations, i.e. start years

#-------------------------------------------------------
#--- Prepare the final object for the nb of infections -
#-------------------------------------------------------
if(group_by==1){nb_inf_2D[]<-apply(nb_inf_6D[,incl_age,incl_serotype,incl_type,incl_vac_status,incl_severity]           ,c(1)  ,sum)}      #Overall nb infections
if(group_by==2){nb_inf_2D[]<-apply(nb_inf_6D[,        ,incl_serotype,incl_type,incl_vac_status,incl_severity]           ,c(1,2),sum)}      #Nb infections by age --> 101 columns
if(group_by==3){nb_inf_2D[]<-apply(nb_inf_6D[,incl_age,             ,incl_type,incl_vac_status,incl_severity]           ,c(1,3),sum)}      #Nb infections by serotype --> 4 columns  
if(group_by==4){nb_inf_2D[]<-apply(nb_inf_6D[,incl_age,incl_serotype,         ,incl_vac_status,incl_severity]           ,c(1,4),sum)}      #Nb infections by type --> 3 columns
if(group_by==5){nb_inf_2D[]<-apply(nb_inf_6D[,incl_age,incl_serotype,incl_type,               ,incl_severity]           ,c(1,5),sum)}      #Nb infections by vac status --> 3 columns
if(group_by==6){nb_inf_2D[]<-apply(nb_inf_6D[,incl_age,incl_serotype,incl_type,incl_vac_status,             ,drop=FALSE],c(1,6),sum)}      #Nb infections by severity --> 6 columns

#-------------------------------------------------------
#--- Save results as csv -------------------------------
#-------------------------------------------------------
write.csv(nb_inf_2D, file=file.path(path_results_aggr, paste(file_name, ".csv", sep="")),row.names=FALSE) 



