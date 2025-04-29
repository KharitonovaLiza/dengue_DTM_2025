
##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA
#  DESCRIPTION: This script allows running a batch of scenarios without the interface, but by listing the scenarios to run in a CSV file
###########################################################################################################################################################################################################

#Note: Prior to running this script a list of scenarios to run has to be provided in the csv file named "run_scenarios_list.csv" located in the same folder as this script

rm(list=ls())

via_interface<<-FALSE

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------- Declare paths to model sub-folders ------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Working directory
path_wd<<-getwd()

#Model codes
path_scripts_misc<<-file.path(path_wd, "0_scripts_misc")                  #Miscellaneous helper scripts
path_scripts_model<<-file.path(path_wd, "1_scripts_model")                #Scripts for transmission model (R)
path_scripts_app<<-file.path(path_wd, "2_scripts_app")                    #Script for Shiny app (R)
path_scripts_src<<-file.path(path_scripts_model, "src")                   #Script for infection process (C)


#Model inputs
path_data<<-file.path(path_wd, "3_data")
path_lists<<-file.path(path_data, "00_lists")
path_inputs_epi<<-file.path(path_data, "11_inputs_epi")                   #Population & transmission inputs
path_inputs_vac_eff<<-file.path(path_data, "12_inputs_vac_eff")           #Vaccine effects inputs
path_inputs_vac_strat<<-file.path(path_data, "13_inputs_vac_strat")       #Vaccination strategy inputs
path_inputs_vec_eff<<-file.path(path_data, "14_inputs_vec_eff")           #Vector control effects inputs
path_inputs_cost_ttt<<-file.path(path_data, "15_inputs_cost_ttt")         #Cost of treatment inputs
path_inputs_cost_vac<<-file.path(path_data, "16_inputs_cost_vac")         #Vaccination cost inputs
path_inputs_cost_vec<<-file.path(path_data, "17_inputs_cost_vec")         #Vector control cost inputs
path_inputs_QoL<<-file.path(path_data, "18_inputs_QoL")                   #Quality-of-life inputs
path_results<<-file.path(path_data, "21_results")  
path_user_data<<-file.path(path_data, "31_user_data")
path_log<<-file.path(path_data,"41_model_logs")

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------- Source required functions -------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Loading packages
debugSource(file.path(path_scripts_app, "0_packages_load.R"))

#Interface codes
debugSource(file.path(path_scripts_app,"1_common.R"))
debugSource(file.path(path_scripts_app,"2_cover_page.R"))
debugSource(file.path(path_scripts_app,"3_inputs_page.R"))
debugSource(file.path(path_scripts_app,"3_inputs_page_cost_ttt.R"))
debugSource(file.path(path_scripts_app,"3_inputs_page_cost_vac.R"))
debugSource(file.path(path_scripts_app,"3_inputs_page_epi_pop.R"))
debugSource(file.path(path_scripts_app,"3_inputs_page_epi_trans.R"))
debugSource(file.path(path_scripts_app,"3_inputs_page_gen.R"))
debugSource(file.path(path_scripts_app,"3_inputs_page_QoL.R"))
debugSource(file.path(path_scripts_app,"3_inputs_page_vac_eff.R"))
debugSource(file.path(path_scripts_app,"3_inputs_page_vac_strat.R"))
debugSource(file.path(path_scripts_app,"3_inputs_page_vec.R"))
debugSource(file.path(path_scripts_app,"4_inputs_page_run_queue.R"))
debugSource(file.path(path_scripts_app,"5_scen_selection_page.R"))
debugSource(file.path(path_scripts_app,"6_results_functions.R"))
debugSource(file.path(path_scripts_app,"6_results_page.R"))
debugSource(file.path(path_scripts_app,"7_results_costs_functions.R"))
debugSource(file.path(path_scripts_app,"7_results_costs_page.R"))
debugSource(file.path(path_scripts_app,"7_results_nb_inf_functions.R"))
debugSource(file.path(path_scripts_app,"7_results_nb_inf_page.R"))
debugSource(file.path(path_scripts_app,"7_results_pop_host_functions.R"))
debugSource(file.path(path_scripts_app,"7_results_pop_host_page.R"))
debugSource(file.path(path_scripts_app,"7_results_QoL_functions.R"))
debugSource(file.path(path_scripts_app,"7_results_QoL_page.R"))
debugSource(file.path(path_scripts_app,"7_results_CE_functions.R"))
debugSource(file.path(path_scripts_app,"7_results_CE_page.R"))
debugSource(file.path(path_scripts_app,"7_results_threshold_functions.R"))
debugSource(file.path(path_scripts_app,"7_results_threshold_page.R"))
debugSource(file.path(path_scripts_app,"7_results_BI_page.R"))

#Model codes
debugSource(file.path(path_scripts_model, "01_model_run.R"))
debugSource(file.path(path_scripts_model, "11_parms_common.R"))
debugSource(file.path(path_scripts_model, "12_parms_theta.R"))
debugSource(file.path(path_scripts_model, "13_parms_psy.R"))
debugSource(file.path(path_scripts_model, "21_parms_epi.R"))
debugSource(file.path(path_scripts_model, "22_parms_vac_eff.R"))
debugSource(file.path(path_scripts_model, "23_parms_vac_strat.R")) 
debugSource(file.path(path_scripts_model, "24_parms_vec_eff.R"))
debugSource(file.path(path_scripts_model, "25_parms_cost.R"))
debugSource(file.path(path_scripts_model, "26_parms_QoL.R"))
debugSource(file.path(path_scripts_model, "31_age_vac_wrapper.R"))
debugSource(file.path(path_scripts_model, "32_ageing.R"))
debugSource(file.path(path_scripts_model, "33_vaccination.R"))
debugSource(file.path(path_scripts_model, "41_results_epi.R"))
debugSource(file.path(path_scripts_model, "51_module_eco.R"))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------- Add shared objects --------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#List with the sequence of start years for automatically generated simulations
list_start_year_int<<-read.csv(file=file.path(path_lists, "list_start_year_int.csv"), stringsAsFactors=FALSE)
list_start_year_int<<-list_start_year_int$list_start_year_int

#List with timeframes to introduce interventions & of maximum timeframe after intervention introduction
list_model_timeframes<<-read.csv(file=file.path(path_lists, "list_model_timeframes.csv"), stringsAsFactors=FALSE)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------- Run all scenarios in the list ---------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#--- Read the list of scenarios to run --------------------------
scen_to_run<<-as.matrix(read.csv(file=file.path(path_scripts_misc,  "run_scenarios_list.csv")))

#--- Loop through each scenario in the list ---------------------
for(scen_nb in 1:dim(scen_to_run)[1]){ 
  
  #--- Get scenario settings --------------------------------------
  scen_ref<<-scen_to_run[scen_nb,'scen_ref']                 #Scenario reference in the format Scen_XXXXX; unique identified for each combination of input sets; results are saved under this name
  set_epi<<-scen_to_run[scen_nb,'set_epi']                   #Reference of the input set with population & transmission inputs
  set_vac_eff<<-scen_to_run[scen_nb,'set_vac_eff']           #Reference of the input set with vaccine effects inputs (always provided, even if a scenario w/o vaccination)
  set_vac_strat<<-scen_to_run[scen_nb,'set_vac_strat']       #Reference of the input set with vaccination strategy inputs (always provided, even if a scenario w/o vaccination)
  set_vec_eff<<-scen_to_run[scen_nb,'set_vec_eff']           #Reference of the input set with vector control inputs (always provided, even if a scenario w/o vector control)
  set_cost_ttt<<-scen_to_run[scen_nb,'set_cost_ttt']         #Reference of the input set with cost of treatment inputs
  set_cost_vac<<-scen_to_run[scen_nb,'set_cost_vac']         #Reference of the input set with vaccination cost inputs (always provided, even if a scenario w/o vaccination)
  set_cost_vec<<-scen_to_run[scen_nb,'set_cost_vec']         #Reference of the input set with vector control cost (always provided, even if a scenario w/o vector control)
  set_QoL<<-scen_to_run[scen_nb,'set_QoL']                   #Reference of the input set with quality-of-life inputs
  scen_ref_res_epi<<-scen_to_run[scen_nb,'scen_ref_res_epi'] #Reference of a scenario with the same combination of population & transmission, vaccination and vector control inputs (if this combination was already run, the results are copied and only the eco module is run; may remain empty)
  vac_switch<<-scen_to_run[scen_nb,'vac_incl']               #Flag for vaccination inclusion (TRUE or FALSE)
  vec_switch<<-scen_to_run[scen_nb,'vec_incl']               #Flag for vector control inclusion (TRUE or FALSE)
  run_epi<<-scen_to_run[scen_nb,'run_epi']                   #Flag to run the transmission module (TRUE or FALSE)
  save_init_state<<-scen_to_run[scen_nb,'save_init_state']   #Flag to save initial states (if a model w/o vaccination with an epi set that hasn't been run before)  
  timeframe<<-as.integer(scen_to_run[scen_nb,'timeframe'])   #Simulation time frame
  
  #--- Loop through each simulation (start year) ------------------  #For each scenario we may want to run multiple simulations, introducing vaccine at different years
  nb_start_years<-scen_to_run[scen_nb,'nb_start_years']              #Get the nb of start years to run (variable "nb_start_years" in the csv file)
  
  for(year_start_id in 1:nb_start_years){                            
    
    print(paste("Running ", scen_ref, ", start year ", year_start_id))
      
    year_start<<-list_start_year_int[year_start_id]                   #Draw the start year from the pre-defined list "list_start_year_int"
    Model_Run(input=NA, output=NA, session=NA, queue=NA)             #Run the model
    
    #--- Update the list of simulations -----------------------------
    list_sim<-Load_ListR(path_lists, file_name="list_simulations")   #List of all simulations (may include multiple lines for the same scenario, with different years of vaccination introduction)
    list<-list_sim[,, drop=F]                                        #Renaming to "list" 

    sims_to_add<-matrix(scen_to_run[scen_nb,1:18],nrow=1,ncol=18)                     #Preparing a new row for the list of simulations (same information as in the csv file)
    sims_to_add[1,7]<-year_start                                                      #Replacing the variable "nb_start_years" by the specific start year used in this simulation
    sims_to_add[1,2]<-as.Date(Sys.Date(), format="%Y-%m-%d")                          #Update the date (so that it's the date the simulation was run and not when it was created)
    list<-rbind(list, sims_to_add[,, drop=F])                                         #Add the new simulation to the list
    list<-list[order(list[,"scen_ref"], as.numeric(list[,"year_start"])),, drop=F]    #Sort by scenario reference and start year
    list_sim<-save(list, file=file.path(path_lists, "list_simulations.Rdata"))        #Save the updated list of simulations before running the next simulation
  }
}

