##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA
#  DESCRIPTION: This script allows re-running the eco part of the model without the interface (to obtain results for an existing scenario but with a new set of eco inputs)
###########################################################################################################################################################################################################

#Note: The previously obtained results will be replaced
#      Before running this script the csv files with eco inputs (costs, QoL) should be modified as needed
#      If the objective is to run new scenarios (rather than replacing old ones), the script "run_scenarios.R" is more appropriate

rm(list=ls())

via_interface<-FALSE

#--- Paths to folders --------------
path_wd<-getwd()
path_data<-file.path(path_wd, "3_data")                             #Folder "3_data"
path_lists<-file.path(path_data, "00_lists")                        #Folder "3_data\00_lists"
path_res_aggr<-file.path(path_data, "22_results_aggr")              #Folder "3_data\22_results_aggr"
path_scripts_app<<-file.path(path_wd, "2_scripts_app")              #Script for Shiny app (R)

#--- Load the list of simulations --
source(file.path(path_scripts_app,"1_common.R"))                   #Source the script with the function Load_ListR
list_sim<-Load_ListR(path_lists, file_name="list_simulations")     #R object containing the list of all run simulations (may contain multiple rows per scenario if multiple simulations were run) 

#Loop through all the simulations
for(scen_id in 1:length(list_sim[,1])){
  
  #--- Get the simulation parameters ---------------------------------
  scen_ref<-list_sim[scen_id,1]                              #Scenario reference in the format Scen_XXXXX; unique identified for each combination of input sets; results are saved under this name
  set_epi<-list_sim[scen_id,'set_epi']                       #Reference of the input set with population & transmission inputs
  set_vac_eff<-list_sim[scen_id,'set_vac_eff']               #Reference of the input set with vaccine effects inputs (always provided, even if a scenario w/o vaccination)
  set_vac_strat<-list_sim[scen_id,'set_vac_strat']           #Reference of the input set with vaccination strategy inputs (always provided, even if a scenario w/o vaccination)
  set_vec_eff<-list_sim[scen_id,'set_vec_eff']               #Reference of the input set with vector control inputs (always provided, even if a scenario w/o vector control)
  set_cost_ttt<-list_sim[scen_id,'set_cost_ttt']             #Reference of the input set with cost of treatment inputs
  set_cost_vac<-list_sim[scen_id,'set_cost_vac']             #Reference of the input set with vaccination cost inputs (always provided, even if a scenario w/o vaccination)
  set_cost_vec<-list_sim[scen_id,'set_cost_vec']             #Reference of the input set with vector control cost (always provided, even if a scenario w/o vector control)
  set_QoL<-list_sim[scen_id,'set_QoL']                       #Reference of the input set with quality-of-life inputs
  scen_ref_res_epi<-list_sim[scen_id,'scen_ref']             #Reference of a scenario with the same combination of population & transmission, vaccination and vector control inputs (if this combination was already run, the results are copied and only the eco module is run; may remain empty)
  vac_switch<-as.logical(list_sim[scen_id,'vac_incl'])       #Flag for vaccination inclusion (TRUE or FALSE)
  vec_switch<-as.logical(list_sim[scen_id,'vec_incl'])       #Flag for vector control inclusion (TRUE or FALSE)
  run_epi<-FALSE                                             #Flag to run the transmission module (TRUE or FALSE)
  save_init_state<-FALSE                                     #Flag to save initial states (if a model w/o vaccination with an epi set that hasn't been run before)  
  year_start<-as.numeric(list_sim[scen_id,'year_start'])     #Year of simulation start (year 1 for scenarios without interventions; year of intervention start otherwise)
  timeframe<-as.numeric(list_sim[scen_id,'timeframe'])       #Simulation time frame
  
  
  #--- Declare paths to model folders --------------------------------
  #Model codes
  path_scripts_model<<-file.path(path_wd, "1_scripts_model")                #Scripts for transmission model (R)
  path_scripts_app<<-file.path(path_wd, "2_scripts_app")                    #Script for Shiny app (R)
  path_scripts_src<<-file.path(path_scripts_model, "src")                   #Script for infection process (C)
  
  #Model inputs
  path_data<<-file.path(path_wd, "3_data")
  path_lists<<-file.path(path_data, "00_lists")
  path_inputs_epi<<-file.path(path_data, "11_inputs_epi")                   #Population & transmission inputs
  path_inputs_vac_eff<<-file.path(path_data, "12_inputs_vac_eff")           #Vacine effects inputs
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
  #Model version & last update date
  model_parms<<-read.csv(file=file.path(path_wd, "model_version.csv"), stringsAsFactors=FALSE)
  
  #List of users
  user_list<<-c(read.csv(file=file.path(path_lists, "list_users.csv"), stringsAsFactors=FALSE)$name)
  
  #List of simulations & scenarios
  list_sim<<-Load_ListR(path_lists, file_name="list_simulations")               #List of all simulations (may include multiple lines for the same scenario)
  nb_sim<<-length(list_sim[,1])                                                 #Total nb of simulations
  nb_scen<<-length(unique(list_sim[,1]))                                        #Number of unique scenarios
  if(nb_scen==1){
    list_scen<-matrix(NA, nrow=nb_scen, ncol=16)                                #List of all scenarios (one line per scenario); removing the columns for start year and timeframe (simulation-specific)
    list_scen[1,]<-list_sim[match(unique(list_sim[,1]),list_sim[,1]),-c(7,8)] 
    colnames(list_scen)<-colnames(list_sim)[-c(7,8)]
    list_scen<<-list_scen
  }else{
    list_scen<<-list_sim[match(unique(list_sim[,1]),list_sim[,1]),-c(7,8)]        
  }
  
  #Lists of inputs sets
  list_sets_epi<<-Load_ListR(path_lists, file_name="list_input_sets_epi")                  #Population & transmission inputs
  list_sets_vac_eff<<-Load_ListR(path_lists, file_name="list_input_sets_vac_eff")          #Vaccine effects inputs
  list_sets_vac_strat<<-Load_ListR(path_lists, file_name="list_input_sets_vac_strat")      #Vaccination strategy inputs
  list_sets_vec_eff<<-Load_ListR(path_lists, file_name="list_input_sets_vec_eff")          #Vector control effects inputs
  list_sets_cost_ttt<<-Load_ListR(path_lists, file_name="list_input_sets_cost_ttt")        #Treatment cost inputs
  list_sets_cost_vac<<-Load_ListR(path_lists, file_name="list_input_sets_cost_vac")        #Vaccination cost inputs
  list_sets_cost_vec<<-Load_ListR(path_lists, file_name="list_input_sets_cost_vec")        #Vector control cost inputs
  list_sets_QoL<<-Load_ListR(path_lists, file_name="list_input_sets_QoL")                  #Quality of life inputs
  
  #Counter of input sets and scenarios
  list_scen_sets_counter<<-Load_ListR(path_lists, file_name="list_scen_sets_counter")
  
  #List with the sequence of start years for automatically generated simulations
  list_start_year_int<<-read.csv(file=file.path(path_lists, "list_start_year_int.csv"), stringsAsFactors=FALSE)
  list_start_year_int<<-list_start_year_int$list_start_year_int
  
  #List with timeframes to introduce interventions & of maximum timeframe after intervention introduction
  list_model_timeframes<<-read.csv(file=file.path(path_lists, "list_model_timeframes.csv"), stringsAsFactors=FALSE)
  
  #List of labels for results display controls 
  list_labels<<-Parms_Import(path=path_lists, filename="list_labels")
  
  #List of default age groups
  age_gr_default<<-read.csv(file=file.path(path_lists, "list_age_groups_default.csv"), stringsAsFactors=FALSE)
  
  
  Model_Run(input=NA,
            output=NA,
            session=NA,
            queue=NA)
  
}

