##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA
#  DESCRIPTION: This script allows exporting into a csv file the results for selected scenarios, including
#               - Nf of infections by type, over time
#               - Nb of newly vaccinated/screened subjects, over time
#               - Costs by type, over time
#               - DALYs, over time
###########################################################################################################################################################################################################

#Note: This code is based on the following assumptions:
#      - All the scenarios that are being exported were run over the timeframe that is equal to or longer than the variable timeframe below
#      - All the scenarios were run with the same number of simulations (i.e. years of vaccine introduction)
#      - Scen_00002 is a scenario with vaccination
#      - Vector control is not included (if it is, it's cost will be ignored)

rm(list=ls())

#=========================================================================
#--- Parameters to provide -----------------------------------------------
#=========================================================================

#-- General parameters ---------------------------------
timeframe<-30            #Simulation timeframe, years (all simulations should be run over this or longer timeframe)
disc_rate_outcomes<-3    #Discount rate for outcomes, as percentage points (e.g. for 3%, the variable should be equal to 3)
disc_rate_cost<-3        #Similarly for the costs


#-- Simulations/scenarios to export --------------------
export_all_scen<-TRUE     #Flag whether all the existing simulations should be exported
export_these_scen<-c()    #If the flag above is FALSE, provide vector with the number of simulations that should be exported (NOTE: there may be multiple simulations per scenario; the values in the vector refer to the SIMULATIONS and not SCENAIROS)

#-- Flag to output confidence intervals ----------------
output_CIs<-FALSE

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
path_results_export<-file.path(path_results_aggr, "export")          #Folder "3_data\22_results_aggr\export"
path_scripts_model<<-file.path(path_wd, "1_scripts_model")           #Scripts for transmission model (R)
path_scripts_app<-file.path(path_wd, "2_scripts_app")                #Script for Shiny app (R)
path_inputs_epi<<-file.path(path_data, "11_inputs_epi")              #Population & transmission inputs
path_inputs_cost_ttt<<-file.path(path_data, "15_inputs_cost_ttt")    #Cost of treatment inputs
path_inputs_cost_vac<<-file.path(path_data, "16_inputs_cost_vac")    #Vaccination cost inputs
path_inputs_QoL<<-file.path(path_data, "18_inputs_QoL")              #Quality-of-life inputs

#-------------------------------------------------------
#--- Source the required functions ---------------------
#-------------------------------------------------------
debugSource(file.path(path_scripts_model, "11_parms_common.R"))
debugSource(file.path(path_scripts_model, "25_parms_cost.R"))  
debugSource(file.path(path_scripts_app,"1_common.R"))               
debugSource(file.path(path_scripts_app,"6_results_functions.R"))
debugSource(file.path(path_scripts_app,"7_results_QoL_functions.R"))
debugSource(file.path(path_scripts_app,"7_results_costs_functions.R"))

#-------------------------------------------------------
#--- Load list of simulation ---------------------------
#-------------------------------------------------------
list_sim<-Load_ListR(path_lists, file_name="list_simulations")     #R object containing the list of all run simulations (may contain multiple rows per scenario if multiple simulations were run) 

#-------------------------------------------------------
#--- List of scenarios to export -----------------------
#-------------------------------------------------------
if(export_all_scen){                                                         #If all scenarios are exported
  scen_refs<-unique(list_sim[,1])                                            #...List of scenario references (unique because there may be multiple simulations per scenario)  
  scen_refs_name<-unique(list_sim[,c(1,5)])                                  #...Matrix with scenario reference/name
}else{                                                                       #If some scenarios/simulations are exported, similar logic, but only taking some of the rows from list_sim
  scen_refs<-unique(list_sim[export_these_scen,1])                           
  scen_refs_name<-unique(list_sim[export_these_scen,c(1,5)])                 
}

nb_scen_to_export<-length(scen_refs)                                         #Total nb of unique scenarios to export

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
#--- Initialize the objects for the results ------------
#-------------------------------------------------------
#Nb infections, by category, over time
nb_inf_2D_asympt                     <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)       #Nb of inf of each severity (see description below), by year (row) and scenario (column)
nb_inf_2D_sympt_mild_non_hosp        <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)       
nb_inf_2D_sympt_mild_hosp            <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)       
nb_inf_2D_sympt_severe_non_hosp      <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)
nb_inf_2D_sympt_severe_hosp_non_fatal<-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)
nb_inf_2D_sympt_severe_hosp_fatal    <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)

#Nb newly vaccinated/screened, over time
nb_new_vac_2D   <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #Nb newly vaccinated, by year (row) and scenario (column)
nb_new_screen_2D<-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #Same for nb newly screened

#Costs, by category, over time
cost_2D_direct_med       <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #Direct medical cost by year (row) and scenario (column)
cost_2D_direct_non_med   <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #Direct non-medical cost
cost_2D_income_lost      <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #Income lost
cost_2D_school_lost      <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #Cost of school absenteeism
cost_2D_screening        <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #Cost of serological screening
cost_2D_vaccine_and_admin<-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #Cost of vaccine acquisition & administration
cost_2D_vaccination_misc <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #Misc cost of vaccination campaign

#DALYs, over time
DALY_2D_mild_non_hosp        <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #DALYs for mild non-hospitalized dengue, by year (row) and scenario (column)
DALY_2D_mild_hosp      		   <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #DALYs for mild hospitalized dengue
DALY_2D_severe_non_hosp   	 <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #DALYs for severe non-hospitalized dengue
DALY_2D_severe_hosp_non_fatal<-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #DALYs for severe hospitalized non-fatal dengue
DALY_2D_severe_hosp_fatal    <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #DALYs for fatal dengue 
DALY_2D_long_term       	   <-matrix(NA,nrow=timeframe,ncol=nb_scen_to_export)  #DALYs for persistent (long-term) dengue

#-------------------------------------------------------
#--- Get the start years used for each scenario --------
#-------------------------------------------------------
#Note: Here it is assumed that all scenarios that are being exported were run with with these start years
#      The list is taken from Scen_00002, assuming that this is a scenario with vaccination
start_years<-as.integer(list_sim[which(list_sim[,"scen_ref"]=="Scen_00002"), "year_start"])    

#-------------------------------------------------------
#--- Get results for each scenario ---------------------
#-------------------------------------------------------

#Looping through all scenarios
for(i in 1:nb_scen_to_export){

  message("Exporting results for scenarios ", i, "/", nb_scen_to_export)
  
  #===================================
  #--- Epi data ----------------------
  #===================================
  
  #--- Load epi data -----------                           
  #Epi data contains:
  #- A 7D object with the nb of infection
  #- A 4D object with the nb of newly vaccinated/screened individuals
  
  epi_data<-ResultsPage_DataEpi(scen_ref=scen_refs[i],     #Scenario reference (e.g. "Scen_00001")
                                start_years=start_years,   #Vector with start years (years of vaccine introduction of corresponding periods in the scenarios w/o vaccination)  
                                timeframe_max=timeframe)   #Timeframe
  
  
  #--- Get nb of infections ----
  #Nb of infections is a 7D array with the nb of inf by year, age, serotype, inf type, vac status, severity (see below) and start year (i.e. year of vac introduction)
  #This 7D object is generated by aggregating the 6D result objects for each individual simulation (i.e. each start year)
  #The severity dimension has 6 levels:
  #1/ Asymptomatic                                  
  #2/ Symptomatic, mild, non-hospitalised
  #3/ Symptomatic, mild, hospitalised  
  #4/ Symptomatic, severe, non-hospitalised --> assumed to be zero (all severe are hospitalized)           
  #5/ Symptomatic, severe, hospitalised, non-fatal
  #6/ Symptomatic, severe, hospitalised, fatal
  
  #Take the object with nb inf
  nb_inf_7D<-epi_data$nb_inf_7D
  
  #Get the nb inf by category
  nb_inf_2D_asympt[,i]                     <-apply(nb_inf_7D[(1:timeframe),,,,,1,],1,sum)/(dim(nb_inf_7D)[7])   #Nb of asymptomatic infections (dividing by the length of the last dimension to obtain average across the simulations, i.e. start years)
  nb_inf_2D_sympt_mild_non_hosp[,i]        <-apply(nb_inf_7D[(1:timeframe),,,,,2,],1,sum)/(dim(nb_inf_7D)[7])
  nb_inf_2D_sympt_mild_hosp[,i]            <-apply(nb_inf_7D[(1:timeframe),,,,,3,],1,sum)/(dim(nb_inf_7D)[7])
  nb_inf_2D_sympt_severe_non_hosp[,i]      <-apply(nb_inf_7D[(1:timeframe),,,,,4,],1,sum)/(dim(nb_inf_7D)[7])
  nb_inf_2D_sympt_severe_hosp_non_fatal[,i]<-apply(nb_inf_7D[(1:timeframe),,,,,5,],1,sum)/(dim(nb_inf_7D)[7])
  nb_inf_2D_sympt_severe_hosp_fatal[,i]    <-apply(nb_inf_7D[(1:timeframe),,,,,6,],1,sum)/(dim(nb_inf_7D)[7])
  
  #Apply discounting
  disc_vector<-ResultsPage_DiscountVector(disc_rate_outcomes, timeframe)      #Create a vector with multiplier for each year in the timeframe
  
  nb_inf_2D_asympt[,i]                     <-nb_inf_2D_asympt[,i]*disc_vector
  nb_inf_2D_sympt_mild_non_hosp[,i]        <-nb_inf_2D_sympt_mild_non_hosp[,i]*disc_vector
  nb_inf_2D_sympt_mild_hosp[,i]            <-nb_inf_2D_sympt_mild_hosp[,i]*disc_vector
  nb_inf_2D_sympt_severe_non_hosp[,i]      <-nb_inf_2D_sympt_severe_non_hosp[,i]*disc_vector
  nb_inf_2D_sympt_severe_hosp_non_fatal[,i]<-nb_inf_2D_sympt_severe_hosp_non_fatal[,i]*disc_vector
  nb_inf_2D_sympt_severe_hosp_fatal[,i]    <-nb_inf_2D_sympt_severe_hosp_fatal[,i]*disc_vector
  
  #--- Get nb vaccinated -------
  #New vac is a 4D array with the nb of newly vaccinated screened by year, age, type (vaccinated or screened) and start year (i.e. year of vac introduction)
  new_vac_4D<-epi_data$pop_host_4D$new_vac
  
  nb_new_vac_2D[,i]   <-apply(new_vac_4D[(1:timeframe),,1,],1,sum)/(dim(new_vac_4D)[4])     #Newly vaccinated, by year (average across simulations)
  nb_new_screen_2D[,i]<-apply(new_vac_4D[(1:timeframe),,2,],1,sum)/(dim(new_vac_4D)[4])     #Newly screened, by year (average across simulations)
  
  #===================================
  #--- Cost data ---------------------
  #===================================
  #--- Load eco data -----------
  #Epi data contain:
  #- A 4D object with costs
  #- A 4D object with DALYs
  
  eco_data<-ResultsPage_DataEco(scen_ref=scen_refs[i],    #Scenario reference (e.g. "Scen_00001")
                                start_years=start_years,  #Vector with start years (years of vaccine introduction of corresponding periods in the scenarios w/o vaccination)  
                                timeframe_max=timeframe)  #Timeframe  

  #--- Get costs ---------------
  #Costs is a 4D array with the costs by year, age, category (see below) and start year (i.e. year of vac introduction)
  #This 4D object is generated by aggregating the 3D result objects for each individual simulation (i.e. each start year)
  #The category dimension has 6 levels:
  #1/ Direct medical costs 
  #2/ Direct non-medical costs (societal perspective only)
  #3/ Income lost (societal perspective only)
  #4/ Cost of serological screening 
  #5/ Cost of vaccine acquistion & administration
  #6/ Misc. vaccination campaign cost
  
  #Take the object with costs
  costs_4D<-eco_data$costs_4D$cost_ttt_vac   
  
  #--------
  #Get the objects by cost category
  cost_2D_direct_med       [,i]<-apply(costs_4D[(1:timeframe),,1,],1,sum)/(dim(costs_4D)[4])  #Direct medical costs (dividing by the length of the last dimension, as we need average across simulations)
  cost_2D_direct_non_med   [,i]<-apply(costs_4D[(1:timeframe),,2,],1,sum)/(dim(costs_4D)[4])  #Similarly for other categories
  cost_2D_income_lost      [,i]<-apply(costs_4D[(1:timeframe),,3,],1,sum)/(dim(costs_4D)[4])
  cost_2D_school_lost      [,i]<-apply(costs_4D[(1:timeframe),,4,],1,sum)/(dim(costs_4D)[4])
  cost_2D_screening        [,i]<-apply(costs_4D[(1:timeframe),,5,],1,sum)/(dim(costs_4D)[4])
  cost_2D_vaccine_and_admin[,i]<-apply(costs_4D[(1:timeframe),,6,],1,sum)/(dim(costs_4D)[4])
  cost_2D_vaccination_misc [,i]<-apply(costs_4D[(1:timeframe),,7,],1,sum)/(dim(costs_4D)[4])
  
  #--------
  #Estimate mortality-related productivity (if it is included)
  if(eco_data$costs_4D$parms_cost_ttt$mortality_related_productivity){
    
    #Extract the array with fatal cases for the timeframe of interest
    fatal_cases_3D<-eco_data$costs_4D$fatal_cases_3D[(1:timeframe),,,drop=F]
    
    #Create a vector with the value of income lost, by age of death
    income_lost_to_death_by_age<-IncomeLostToDeath(life_expectancy=eco_data$costs_4D$life_expectancy,
                                                   parms_cost_ttt=eco_data$costs_4D$parms_cost_ttt,
                                                   disc_rate_costs=disc_rate_cost)  
    
    #Convert in an array with the same dimensions as fatal_cases_3D (timeframe, 101, nb start years)
    income_lost_to_death_by_age<-array(rep(income_lost_to_death_by_age, each=timeframe),dim=dim(fatal_cases_3D))
    
    #Estimate the cost of fatal cases by multiplying the array with nb of cases and income lost
    fatal_cases_cost_3D<-fatal_cases_3D*income_lost_to_death_by_age
    
    #Calculate average over start years
    fatal_cases_cost<-apply(fatal_cases_cost_3D[(1:timeframe),,],1,sum)/(dim(costs_4D)[4])
    
    #Add the final estimations to the previous estimate of income lost
    cost_2D_income_lost[,i]<-cost_2D_income_lost[,i]+fatal_cases_cost
  }
  
  #--------
  #Apply discounting
  disc_vector<-ResultsPage_DiscountVector(disc_rate_cost, timeframe)              #Re-calculate the discounting vector with the discount rate for costs
  
  cost_2D_direct_med       [,i]<-cost_2D_direct_med       [,i]*disc_vector        #Apply this vector to each vector with costs, for a given scenario
  cost_2D_direct_non_med   [,i]<-cost_2D_direct_non_med   [,i]*disc_vector
  cost_2D_income_lost      [,i]<-cost_2D_income_lost      [,i]*disc_vector
  cost_2D_school_lost      [,i]<-cost_2D_school_lost      [,i]*disc_vector
  cost_2D_screening        [,i]<-cost_2D_screening        [,i]*disc_vector
  cost_2D_vaccine_and_admin[,i]<-cost_2D_vaccine_and_admin[,i]*disc_vector
  cost_2D_vaccination_misc [,i]<-cost_2D_vaccination_misc [,i]*disc_vector
  
  #--- Get DALYs ---------------Q
  #List QoL_4D contains the following objects:
  #- 4D object QoL_cases_4D with the nb of cases by year, age, severity (mild, severe or fatal) and start year
  #- parms_QoL - a list of all parameters related to QoL (disability weight, disability duration etc)
  #- life_expectancy - vector with life expectancy by age
  
  QoL_4D<-eco_data$QoL_4D
  
  #Transform 4D data into 3D data (using a function previously written for the interface)
  QoL_3D<-ResultsPage_QoL_Data3D(scen_ref=scen_refs[i],            #Scenario reference
                                 QoL_4D_101y=QoL_4D,               #QoL data (list of several objects, see above)
                                 incl_crit_age=c(1:101),           #Ages to be included (here, all)
                                 timeframe=c(1,timeframe),         #First and last year of the timeframe of interest
                                 disc_rate=disc_rate_outcomes)	   #Discount rate for outcomes
       
  #Transform 3D data into 2D data (also using a previously written function; discounting already included in the functions)
  
  QoL_2D_mild_non_hosp        <-ResultsPage_QoL_Data2D(QoL_3D=QoL_3D$QoL_3D_by_series_mild_non_hosp,    
                                                       timeframe=c(1,timeframe),  
                                                       disc_rate=disc_rate_outcomes)
  
  QoL_2D_mild_hosp            <-ResultsPage_QoL_Data2D(QoL_3D=QoL_3D$QoL_3D_by_series_mild_hosp,    
                                                       timeframe=c(1,timeframe),  
                                                       disc_rate=disc_rate_outcomes)
  
  QoL_2D_severe_non_hosp      <-ResultsPage_QoL_Data2D(QoL_3D=QoL_3D$QoL_3D_by_series_severe_non_hosp,    
                                                       timeframe=c(1,timeframe),  
                                                       disc_rate=disc_rate_outcomes)
  
  QoL_2D_severe_hosp_non_fatal<-ResultsPage_QoL_Data2D(QoL_3D=QoL_3D$QoL_3D_by_series_severe_hosp_non_fatal,    
                                                       timeframe=c(1,timeframe),  
                                                       disc_rate=disc_rate_outcomes) 
  
  QoL_2D_severe_hosp_fatal    <-ResultsPage_QoL_Data2D(QoL_3D=QoL_3D$QoL_3D_by_series_severe_hosp_fatal,    
                                                       timeframe=c(1,timeframe),  
                                                       disc_rate=disc_rate_outcomes) 									 
  
  QoL_2D_long_term			      <-ResultsPage_QoL_Data2D(QoL_3D=QoL_3D$QoL_3D_by_series_long_term,    
                                                       timeframe=c(1,timeframe),  
                                                       disc_rate=disc_rate_outcomes) 
 								 
	  
	DALY_2D_mild_non_hosp[,i]     	 <-QoL_2D_mild_non_hosp$over_time
	DALY_2D_mild_hosp[,i]     		   <-QoL_2D_mild_hosp$over_time
	DALY_2D_severe_non_hosp[,i]    	 <-QoL_2D_severe_non_hosp$over_time
	DALY_2D_severe_hosp_non_fatal[,i]<-QoL_2D_severe_hosp_non_fatal$over_time
	DALY_2D_severe_hosp_fatal[,i]    <-QoL_2D_severe_hosp_fatal$over_time
	DALY_2D_long_term[,i]			       <-QoL_2D_long_term$over_time

}

#-------------------------------------------------------
#--- Prepare results for saving ------------------------
#-------------------------------------------------------
#Add column headers with scenario references
colnames(nb_inf_2D_asympt)                     <-scen_refs
colnames(nb_inf_2D_sympt_mild_non_hosp)        <-scen_refs
colnames(nb_inf_2D_sympt_mild_hosp)            <-scen_refs
colnames(nb_inf_2D_sympt_severe_non_hosp)      <-scen_refs
colnames(nb_inf_2D_sympt_severe_hosp_non_fatal)<-scen_refs
colnames(nb_inf_2D_sympt_severe_hosp_fatal)    <-scen_refs
colnames(nb_new_vac_2D)                        <-scen_refs
colnames(nb_new_screen_2D)                     <-scen_refs
colnames(cost_2D_direct_med)                   <-scen_refs
colnames(cost_2D_direct_non_med)               <-scen_refs
colnames(cost_2D_income_lost)                  <-scen_refs
colnames(cost_2D_school_lost)                  <-scen_refs
colnames(cost_2D_screening)                    <-scen_refs
colnames(cost_2D_vaccine_and_admin)            <-scen_refs
colnames(cost_2D_vaccination_misc)             <-scen_refs
colnames(DALY_2D_mild_non_hosp)                <-scen_refs
colnames(DALY_2D_mild_hosp)   	               <-scen_refs
colnames(DALY_2D_severe_non_hosp)              <-scen_refs
colnames(DALY_2D_severe_hosp_non_fatal)        <-scen_refs
colnames(DALY_2D_severe_hosp_fatal)            <-scen_refs
colnames(DALY_2D_long_term)                    <-scen_refs



#-------------------------------------------------------
#--- Save results as csv -------------------------------
#-------------------------------------------------------
write.csv(nb_inf_2D_asympt,                      file=file.path(path_results_export, "res_nb_inf_asympt.csv"                     ),row.names=TRUE) 
write.csv(nb_inf_2D_sympt_mild_non_hosp,         file=file.path(path_results_export, "res_nb_inf_sympt_mild_non_hosp.csv"        ),row.names=TRUE)
write.csv(nb_inf_2D_sympt_mild_hosp,             file=file.path(path_results_export, "res_nb_inf_sympt_mild_hosp.csv"            ),row.names=TRUE)
write.csv(nb_inf_2D_sympt_severe_non_hosp,       file=file.path(path_results_export, "res_nb_inf_sympt_severe_non_hosp.csv"      ),row.names=TRUE)
write.csv(nb_inf_2D_sympt_severe_hosp_non_fatal, file=file.path(path_results_export, "res_nb_inf_sympt_severe_hosp_non_fatal.csv"),row.names=TRUE)
write.csv(nb_inf_2D_sympt_severe_hosp_fatal,     file=file.path(path_results_export, "res_nb_inf_sympt_severe_hosp_fatal.csv"    ),row.names=TRUE)
write.csv(nb_new_vac_2D,                         file=file.path(path_results_export, "res_nb_new_vac.csv"                        ),row.names=TRUE)
write.csv(nb_new_screen_2D,                      file=file.path(path_results_export, "res_nb_new_screen.csv"                     ),row.names=TRUE)
write.csv(cost_2D_direct_med,                    file=file.path(path_results_export, "res_cost_direct_med.csv"                   ),row.names=TRUE)
write.csv(cost_2D_direct_non_med,                file=file.path(path_results_export, "res_cost_direct_non_med.csv"               ),row.names=TRUE) 
write.csv(cost_2D_income_lost,                   file=file.path(path_results_export, "res_cost_income_lost.csv"                  ),row.names=TRUE) 
write.csv(cost_2D_school_lost,                   file=file.path(path_results_export, "res_cost_school_lost.csv"                  ),row.names=TRUE)
write.csv(cost_2D_screening,                     file=file.path(path_results_export, "res_cost_screening.csv"                    ),row.names=TRUE) 
write.csv(cost_2D_vaccine_and_admin,             file=file.path(path_results_export, "res_cost_vaccine_and_admin.csv"            ),row.names=TRUE) 
write.csv(cost_2D_vaccination_misc,              file=file.path(path_results_export, "res_cost_vaccination_misc.csv"             ),row.names=TRUE) 
write.csv(DALY_2D_mild_non_hosp,                 file=file.path(path_results_export, "res_DALYs_mild_non_hosp.csv"               ),row.names=TRUE)
write.csv(DALY_2D_mild_hosp,                     file=file.path(path_results_export, "res_DALYs_mild_hosp.csv"                   ),row.names=TRUE)
write.csv(DALY_2D_severe_non_hosp,               file=file.path(path_results_export, "res_DALYs_severe_non_hosp.csv"             ),row.names=TRUE)
write.csv(DALY_2D_severe_hosp_non_fatal,         file=file.path(path_results_export, "res_DALYs_severe_hosp_non_fatal.csv"       ),row.names=TRUE)
write.csv(DALY_2D_severe_hosp_fatal,             file=file.path(path_results_export, "res_DALYs_severe_hosp_fatal.csv"           ),row.names=TRUE)
write.csv(DALY_2D_long_term,                     file=file.path(path_results_export, "res_DALYs_long_term.csv"                   ),row.names=TRUE)

#Also saving a csv with scenario names
write.csv(scen_refs_name, file=file.path(path_results_export, "res_scen_refs_name.csv"),row.names=TRUE) 

#-------------------------------------------------------
#--- Lower and upper CI  -------------------------------
#-------------------------------------------------------
if(output_CIs){
  
  if (!("EnvStats" %in% rownames(installed.packages()))){install.packages("EnvStats")}
  require('EnvStats')
  
  ########################################################################################
  #This function calculates 95% confidence interval for a log-normal distribution
  CI_lognormal<-function(x){
    if(length(x)==length(x[which(x==0)])){
      mean=0
      lower=0
      upper=0
    }else{
      res<-elnorm(x, method = "mvue", ci = TRUE, ci.type = "two-sided", ci.method = "exact", conf.level = 0.95)
      quantile_lognorm<-qlnorm(p=c(0.025, 0.975), meanlog = res$parameters[1], sdlog = res$parameters[2], lower.tail = TRUE, log.p = FALSE)
      mean=exp(res$parameters[1])
      lower=quantile_lognorm[1]
      upper=quantile_lognorm[2]
    }
    
    return(list(mean=unname(mean),lower=unname(lower),upper=unname(upper)))
  }
  
  
  #Nb infections, by category, over time
  nb_inf_3D_asympt                     <-array(NA,dim=c(timeframe,dim(nb_inf_7D)[7],nb_scen_to_export))      
  nb_inf_3D_sympt_mild_non_hosp        <-array(NA,dim=c(timeframe,dim(nb_inf_7D)[7],nb_scen_to_export))        
  nb_inf_3D_sympt_mild_hosp            <-array(NA,dim=c(timeframe,dim(nb_inf_7D)[7],nb_scen_to_export))        
  nb_inf_3D_sympt_severe_non_hosp      <-array(NA,dim=c(timeframe,dim(nb_inf_7D)[7],nb_scen_to_export))  
  nb_inf_3D_sympt_severe_hosp_non_fatal<-array(NA,dim=c(timeframe,dim(nb_inf_7D)[7],nb_scen_to_export))  
  nb_inf_3D_sympt_severe_hosp_fatal    <-array(NA,dim=c(timeframe,dim(nb_inf_7D)[7],nb_scen_to_export))  
  nb_inf_3D_total                      <-array(NA,dim=c(timeframe,dim(nb_inf_7D)[7],nb_scen_to_export)) 
  nb_inf_3D_sympt                      <-array(NA,dim=c(timeframe,dim(nb_inf_7D)[7],nb_scen_to_export)) 
  nb_inf_3D_sympt_hosp                 <-array(NA,dim=c(timeframe,dim(nb_inf_7D)[7],nb_scen_to_export)) 
  nb_inf_3D_sympt_non_hosp             <-array(NA,dim=c(timeframe,dim(nb_inf_7D)[7],nb_scen_to_export)) 
  nb_inf_3D_mild                       <-array(NA,dim=c(timeframe,dim(nb_inf_7D)[7],nb_scen_to_export)) 
  nb_inf_3D_severe                     <-array(NA,dim=c(timeframe,dim(nb_inf_7D)[7],nb_scen_to_export)) 
  #Nb infections avoided, by category, and start year and scenario
  nb_inf_avoided_total           <-array(0,dim=c(dim(nb_inf_7D)[7],nb_scen_to_export-1))    
  nb_inf_avoided_asympt          <-array(0,dim=c(dim(nb_inf_7D)[7],nb_scen_to_export-1))      
  nb_inf_avoided_sympt           <-array(0,dim=c(dim(nb_inf_7D)[7],nb_scen_to_export-1))       
  nb_inf_avoided_sympt_hospit    <-array(0,dim=c(dim(nb_inf_7D)[7],nb_scen_to_export-1)) 
  nb_inf_avoided_sympt_non_hospit<-array(0,dim=c(dim(nb_inf_7D)[7],nb_scen_to_export-1)) 
  nb_inf_avoided_mild            <-array(0,dim=c(dim(nb_inf_7D)[7],nb_scen_to_export-1)) 
  nb_inf_avoided_severe          <-array(0,dim=c(dim(nb_inf_7D)[7],nb_scen_to_export-1)) 
  nb_inf_avoided_fatal           <-array(0,dim=c(dim(nb_inf_7D)[7],nb_scen_to_export-1)) 
  
  #Looping through all scenarios
  for(i in 1:nb_scen_to_export){
    
    message("Exporting results for scenarios ", i, "/", nb_scen_to_export)
    
    #--- Load epi data -----------                           
    epi_data<-ResultsPage_DataEpi(scen_ref=scen_refs[i],     #Scenario reference (e.g. "Scen_00001")
                                  start_years=start_years,   #Vector with start years (years of vaccine introduction of corresponding periods in the scenarios w/o vaccination)  
                                  timeframe_max=timeframe)   #Timeframe
    
    
    #Take the object with nb inf
    nb_inf_7D<-epi_data$nb_inf_7D
    #Get the nb inf by category
    nb_inf_3D_asympt[,,i]                     <-apply(nb_inf_7D[(1:timeframe),,,,,1,,drop=FALSE],c(1,7),sum)   #Nb of asymptomatic infections (dividing by the length of the last dimension to obtain average across the simulations, i.e. start years)
    nb_inf_3D_sympt_mild_non_hosp[,,i]        <-apply(nb_inf_7D[(1:timeframe),,,,,2,,drop=FALSE],c(1,7),sum)
    nb_inf_3D_sympt_mild_hosp[,,i]            <-apply(nb_inf_7D[(1:timeframe),,,,,3,,drop=FALSE],c(1,7),sum)
    nb_inf_3D_sympt_severe_non_hosp[,,i]      <-apply(nb_inf_7D[(1:timeframe),,,,,4,,drop=FALSE],c(1,7),sum)
    nb_inf_3D_sympt_severe_hosp_non_fatal[,,i]<-apply(nb_inf_7D[(1:timeframe),,,,,5,,drop=FALSE],c(1,7),sum)
    nb_inf_3D_sympt_severe_hosp_fatal[,,i]    <-apply(nb_inf_7D[(1:timeframe),,,,,6,,drop=FALSE],c(1,7),sum)
    
    #Apply discounting
    disc_vector<-ResultsPage_DiscountVector(disc_rate_outcomes, timeframe)      #Create a vector with multiplier for each year in the timeframe
    disc_vector<-array(rep(disc_vector,each=dim(nb_inf_7D)[7]),dim(nb_inf_3D_asympt[,,i]))
    
    nb_inf_3D_asympt[,,i]                     <-nb_inf_3D_asympt[,,i]*disc_vector
    nb_inf_3D_sympt_mild_non_hosp[,,i]        <-nb_inf_3D_sympt_mild_non_hosp[,,i]*disc_vector
    nb_inf_3D_sympt_mild_hosp[,,i]            <-nb_inf_3D_sympt_mild_hosp[,,i]*disc_vector
    nb_inf_3D_sympt_severe_non_hosp[,,i]      <-nb_inf_3D_sympt_severe_non_hosp[,,i]*disc_vector
    nb_inf_3D_sympt_severe_hosp_non_fatal[,,i]<-nb_inf_3D_sympt_severe_hosp_non_fatal[,,i]*disc_vector
    nb_inf_3D_sympt_severe_hosp_fatal[,,i]    <-nb_inf_3D_sympt_severe_hosp_fatal[,,i]*disc_vector
    
    nb_inf_3D_sympt_non_hosp[,,i]<-nb_inf_3D_sympt_mild_non_hosp[,,i]+nb_inf_3D_sympt_severe_non_hosp[,,i]
    nb_inf_3D_sympt_hosp[,,i]<-nb_inf_3D_sympt_mild_hosp[,,i]+nb_inf_3D_sympt_severe_hosp_non_fatal[,,i]+nb_inf_3D_sympt_severe_hosp_fatal[,,i] 
    nb_inf_3D_sympt[,,i]     <-nb_inf_3D_sympt_hosp[,,i]+nb_inf_3D_sympt_non_hosp[,,i]    
    nb_inf_3D_total[,,i]     <-nb_inf_3D_asympt[,,i]+nb_inf_3D_sympt[,,i]
    nb_inf_3D_mild[,,i]      <-nb_inf_3D_sympt_mild_non_hosp[,,i] +nb_inf_3D_sympt_mild_hosp[,,i] 
    nb_inf_3D_severe[,,i]    <-nb_inf_3D_sympt_severe_non_hosp[,,i]+nb_inf_3D_sympt_severe_hosp_non_fatal[,,i]+nb_inf_3D_sympt_severe_hosp_fatal[,,i] 
    
    if(i>1){
      
      nb_inf_avoided_total[,i-1]           <-(apply(nb_inf_3D_total[,,1],2,sum)-apply(nb_inf_3D_total[,,i], 2,sum))/apply(nb_inf_3D_total[,,1],2,sum)  
      nb_inf_avoided_asympt[,i-1]          <-(apply(nb_inf_3D_asympt[,,1],2,sum)-apply(nb_inf_3D_asympt[,,i],2,sum))/apply(nb_inf_3D_asympt[,,1],2,sum) 
      nb_inf_avoided_sympt[,i-1]           <-(apply(nb_inf_3D_sympt[,,1],2,sum)-apply(nb_inf_3D_sympt[,,i],2,sum))/apply(nb_inf_3D_sympt[,,1],2,sum)      
      nb_inf_avoided_sympt_hospit[,i-1]    <-(apply(nb_inf_3D_sympt_hosp[,,1],2,sum)-apply(nb_inf_3D_sympt_hosp[,,i],2,sum))/apply(nb_inf_3D_sympt_hosp[,,1],2,sum)     
      nb_inf_avoided_sympt_non_hospit[,i-1]<-(apply(nb_inf_3D_sympt_non_hosp[,,1],2,sum)-apply(nb_inf_3D_sympt_non_hosp[,,i],2,sum))/apply(nb_inf_3D_sympt_non_hosp[,,1],2,sum)     
      nb_inf_avoided_mild[,i-1]            <-(apply(nb_inf_3D_mild[,,1],2,sum)-apply(nb_inf_3D_mild[,,i],2,sum))/apply(nb_inf_3D_mild[,,1],2,sum)     
      nb_inf_avoided_severe[,i-1]          <-(apply(nb_inf_3D_severe[,,1],2,sum)-apply(nb_inf_3D_severe[,,i],2,sum))/apply(nb_inf_3D_severe[,,1],2,sum)     
      nb_inf_avoided_fatal[,i-1]           <-(apply(nb_inf_3D_sympt_severe_hosp_fatal[,,1],2,sum)-apply(nb_inf_3D_sympt_severe_hosp_fatal[,,i],2,sum))/apply(nb_inf_3D_sympt_severe_hosp_fatal[,,1],2,sum)     
    }
  }
  
  nb_inf_avoided_lower_CI<-matrix(NA,nrow=nb_scen_to_export-1,ncol=8)
  nb_inf_avoided_upper_CI<-matrix(NA,nrow=nb_scen_to_export-1,ncol=8)
  
  for(i in 1:(nb_scen_to_export-1)){
    CI<-CI_lognormal(nb_inf_avoided_total[,i])
    nb_inf_avoided_lower_CI[i,1]<-CI$lower
    nb_inf_avoided_upper_CI[i,1]<-CI$upper
    CI<-CI_lognormal(nb_inf_avoided_asympt[,i])
    nb_inf_avoided_lower_CI[i,2]<-CI$lower
    nb_inf_avoided_upper_CI[i,2]<-CI$upper
    CI<-CI_lognormal(nb_inf_avoided_sympt[,i])
    nb_inf_avoided_lower_CI[i,3]<-CI$lower
    nb_inf_avoided_upper_CI[i,3]<-CI$upper
    CI<-CI_lognormal(nb_inf_avoided_sympt_hospit[,i])
    nb_inf_avoided_lower_CI[i,4]<-CI$lower
    nb_inf_avoided_upper_CI[i,4]<-CI$upper
    CI<-CI_lognormal(nb_inf_avoided_sympt_non_hospit[,i])
    nb_inf_avoided_lower_CI[i,5]<-CI$lower
    nb_inf_avoided_upper_CI[i,5]<-CI$upper
    CI<-CI_lognormal(nb_inf_avoided_mild[,i])
    nb_inf_avoided_lower_CI[i,6]<-CI$lower
    nb_inf_avoided_upper_CI[i,6]<-CI$upper
    CI<-CI_lognormal(nb_inf_avoided_severe[,i])
    nb_inf_avoided_lower_CI[i,7]<-CI$lower
    nb_inf_avoided_upper_CI[i,7]<-CI$upper
    CI<-CI_lognormal(nb_inf_avoided_fatal[,i])
    nb_inf_avoided_lower_CI[i,8]<-CI$lower
    nb_inf_avoided_upper_CI[i,8]<-CI$upper
  }
  
  write.csv(nb_inf_avoided_lower_CI, file=file.path(path_results_export, "nb_inf_avoided_lower_CI.csv"),row.names=TRUE) 
  write.csv(nb_inf_avoided_upper_CI, file=file.path(path_results_export, "nb_inf_avoided_upper_CI.csv"),row.names=TRUE) 
}





