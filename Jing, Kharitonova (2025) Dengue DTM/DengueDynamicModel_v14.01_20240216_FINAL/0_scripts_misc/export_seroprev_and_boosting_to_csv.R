##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA
#  DESCRIPTION: This script allow analyzing the occurence of boosting
###########################################################################################################################################################################################################

#This script allow analyzing the indirect effects and the impact of boosting by creating the several objects:
#- Seroprevalence without vaccination (using object seroprev)
#- Seroprevalence with vaccination (using object seroprev)
#- Seroprevalence without vaccination, only in cohorts who would be vaccinated (using object seroprev)
#- Seroprevalence with vaccination, only in cohorts who are vaccinated (using object seroprev)
#- Seroprevalence with vaccination in cohorts who are vaccinated (using object trace)
#- Distribution of the vaccinated subject by episode post-vaccination, by year since last boosting (using object trace)

#Results are exported to csv files located in the folder "3_data\22_results_aggr\export_impact_of_vac"
#Note: Seroprevalence with vaccination is affected by two factors - vaccine efficacy against infection itself and indirect effects
#      Since the object seroprev does not differentiate between vaccinated and unvaccinated individuals, we can only see the overall reduction in seroprevalence
#      If the vaccine does not protect against infection, this is appropriate (both vaccinated and unvaccinated will experience the same indirect effects)

#Note: When comparing the seroprevalence estimated from the objects "seroprev" and "trace", one should keep in mind that:
#      - In the object "seroprev" the values for a specific age are at the END of the year; in the object "trace" they are at the BEGINNING of the year; so there may be discrepancies
#      - The object "seroprev" is estimated directly in the transmission module, while the object "trace" is approximated by counting the nb of infections
#      - The results obtained from the object "seroprev" for vaccinated cohorts only do not take coverage into account (i.e. these are not the vaccinated individuals, but all individuals of eligible age)

rm(list=ls())

#=========================================================================
#--- Scenarios to export -------------------------------------------------
#=========================================================================

ref_scen_wo_vac<-"Scen_00009"         #Reference of a scenario without vaccination
ref_scen_with_vac<-"Scen_00012"       #Reference of a scenario with vaccination 

#=========================================================================
#--- Code that exports the results ---------------------------------------
#=========================================================================

#----------------------------------------------------------------
#---- Prepare general parameter ---------------------------------
#----------------------------------------------------------------
#-- Paths to folders -------------
path_wd<-getwd()                                                            #Working directory

path_data<-file.path(path_wd, "3_data")                                     #Folder "3_data"
path_inputs_epi<<-file.path(path_data, "11_inputs_epi")                     #Population & transmission inputs
path_inputs_vac_eff<<-file.path(path_data, "12_inputs_vac_eff")             #Vaccine effects inputs
path_inputs_vac_strat<<-file.path(path_data, "13_inputs_vac_strat")         #Vaccination strategy inputs

path_results<-file.path(path_data, "21_results")                            #Folder "3_data\21_results" with results of the simulations
path_results_aggr<-file.path(path_data, "22_results_aggr")                  #Folder "3_data\22_results_aggr"
path_results_export<-file.path(path_results_aggr, "export_impact_of_vac")   #Folder "3_data\22_results_aggr\export_impact_of_vac"
path_lists<-file.path(path_data, "00_lists")                                #Folder "3_data\00_lists"
path_scripts_model<<-file.path(path_wd, "1_scripts_model")                  #Scripts for transmission model (R)
path_scripts_app<<-file.path(path_wd, "2_scripts_app")                      #Script for Shiny app (R)

#-- Source functions -------------
debugSource(file.path(path_scripts_app,"1_common.R")) 
debugSource(file.path(path_scripts_model, "11_parms_common.R"))
debugSource(file.path(path_scripts_model, "12_parms_theta.R"))
debugSource(file.path(path_scripts_model, "13_parms_psy.R"))
debugSource(file.path(path_scripts_model, "21_parms_epi.R"))
debugSource(file.path(path_scripts_model, "22_parms_vac_eff.R"))
debugSource(file.path(path_scripts_model, "23_parms_vac_strat.R")) 

#-- Used input sets --------------
#Load the list of simulations
list_sim<-Load_ListR(path_lists, file_name="list_simulations")                                    

#Get settings that were used to run the scenario with vaccination
timeframe<-min(as.integer(list_sim[which(list_sim[,"scen_ref"]==ref_scen_with_vac), "timeframe"]))    #Timeframe (taking the min in case different simulations were run with different timeframe)

start_years<-as.integer(list_sim[which(list_sim[,"scen_ref"]==ref_scen_with_vac), "year_start"])      #List of start years (when vaccination was introduced)
nb_start_years<-length(start_years)

set_epi<-list_sim[which(list_sim[,"scen_ref"]==ref_scen_with_vac), "set_epi"][1]                      #Set of epi (population & transmission parameters); only taking the first element, as it's necessarily the same for all simulations
set_vac_eff<-list_sim[which(list_sim[,"scen_ref"]==ref_scen_with_vac), "set_vac_eff"][1]              #Set of vaccine efficacy parameters
set_vac_strat<-list_sim[which(list_sim[,"scen_ref"]==ref_scen_with_vac), "set_vac_strat"][1]          #Set of vaccination strategy parameters

#-- Load epi and vac parms -------
parms_epi<-ParmsEpi_Load(path_inputs_epi, set_epi, calibr=FALSE) 
parms_vac_strat<-ParmsVacStrat_Load(path_inputs_vac_strat, set_vac_strat, parms_epi)
parms_vac_eff<-ParmsVacEff_Load(path_inputs_vac_eff, set_vac_eff, parms_vac_strat, timeframe)

#-- Get age of vac ---------------
routine_age<-parms_vac_strat$parms_vac_strat$age_routine                  #Age of routine vaccination, years
catch_up_incl<-parms_vac_strat$parms_vac_strat$catchup_switch             #Inclusion of catch-up
if(catch_up_incl){
  catch_up_age_low<-parms_vac_strat$parms_vac_strat$age_catchup_low         #Youngest age in catch-up campaign, years
  catch_up_age_high<-parms_vac_strat$parms_vac_strat$age_catchup_high       #Oldest age in catch-up campaign, years
}else{
  catch_up_age_low<-NA
  catch_up_age_high<-NA
}

#---------------------------------------------------------------- #Seroprevalence will be different in each simulation, as each simulation starts with a different year 
#---- Get average seroprevalence across simulations ------------- #(i.e. there may be a peak or a drop in incidence etc)
#---------------------------------------------------------------- #Aggregating average seroprevalence across simulations

#-- Create new objects ------------------
#At each simulation an object "seroprev" is saved with seroprevalence by simulation year, age and serostatus (5 levels - history of 0 infections, 1 infection, ..., 4 infection)
#Creating new objects with an additional dimension for each start year
seroprev_wo_vac_new  <-array(NA, dim=c(timeframe, 101, 5, nb_start_years))   #Dimensions are year, age, serostatus, year of vac introduction
seroprev_with_vac_new<-array(NA, dim=c(timeframe, 101, 5, nb_start_years))

#-- Get seroprev wo vaccination ---------
#Load predicted seroprev
#For a scenario without vaccination there's just one object to load (seroprevalence from calendar year 1 to year 101)
load(file.path(path_results, paste(ref_scen_wo_vac,'_y1_Epi.Rdata', sep="")))  #Loading epi results for a scenario w/o vac
seroprev_wo_vac_old<-results_epi$seroprev                                      #Taking the object for seroprevalence 
rm(results_epi); gc()                                                          #Removing the larger object results_epi

#Fill the object for each start year
for(y in 1:nb_start_years){                     #Loop through each start year (simulation)
  year_first<-start_years[y]                    #First year to be taken (for a scenario w/o vaccination it's the calendar year of when the vaccine would be introduced)
  year_last <-year_first+timeframe-1            #Last year to be taken
  seroprev_wo_vac_new[,,,y]<-seroprev_wo_vac_old[(year_first:year_last),,]   #Fill the new object by taking the corresponding years
}
rm(seroprev_wo_vac_old); gc()

#-- Get seroprev with vaccination -------
#Load predicted serorevalence.
#For scenarios with vaccination there are multiple result objects (one for each start year), so loading them in a loop
for(y in 1:nb_start_years){
  load(file.path(path_results, paste(ref_scen_with_vac, "_y", start_years[y], "_Epi.Rdata", sep="")))   #Loading epi results
  seroprev_with_vac_old<-results_epi$seroprev                                                           #Taking the object for seroprevalence
  rm(results_epi); gc()                                                                                 #Removing the larger object results_epi
  seroprev_with_vac_new[,,,y]<-seroprev_with_vac_old[(1:timeframe),,]                                   #Saving results into a bigger object (for a scenario with vaccination taking results from year 1, i.e. year when vaccine was introduced)
  rm(seroprev_with_vac_old); gc()
}

#-- Get average across simulations ------
seroprev_wo_vac_avg  <-apply(seroprev_wo_vac_new  ,c(1,2,3), mean)       #Average nb of individuals at each year, each age, with each serostatus, across all simulations
seroprev_with_vac_avg<-apply(seroprev_with_vac_new,c(1,2,3), mean) 
rm(seroprev_wo_vac_new, seroprev_with_vac_new); gc()

#----------------------------------------------------------------
#---- Seroprevalence wo and with vaccination (all cohorts) ------
#----------------------------------------------------------------

#-- Initialise final objects ------------
seroprev_in_all_cohorts_wo_vac_final  <-matrix(NA, nrow=timeframe, ncol=5)    #One row per year, one column per serostatus
seroprev_in_all_cohorts_with_vac_final<-matrix(NA, nrow=timeframe, ncol=5)

#-- Fill the final objects --------------
seroprev_in_all_cohorts_wo_vac_final  <-apply(seroprev_wo_vac_avg  [,,,drop=F],c(1,3),sum)   #Sum across all ages, for each simulation year and each serostatus
seroprev_in_all_cohorts_with_vac_final<-apply(seroprev_with_vac_avg[,,,drop=F],c(1,3),sum)

#-- Save the final objects --------------
colnames(seroprev_in_all_cohorts_wo_vac_final)  <-c("naive", "history_of_1inf", "history_of_2inf", "history_of_3inf", "history_of_4inf")
colnames(seroprev_in_all_cohorts_with_vac_final)<-c("naive", "history_of_1inf", "history_of_2inf", "history_of_3inf", "history_of_4inf")

write.csv(seroprev_in_all_cohorts_wo_vac_final,   file=file.path(path_results_export, paste(ref_scen_wo_vac,   "_seroprev_in_all_cohorts.csv", sep="")), row.names=TRUE)
write.csv(seroprev_in_all_cohorts_with_vac_final, file=file.path(path_results_export, paste(ref_scen_with_vac, "_seroprev_in_all_cohorts.csv", sep="")), row.names=TRUE) 

#----------------------------------------------------------------
#---- Seroprevalence wo and with vaccination (vac cohorts) ------
#----------------------------------------------------------------
#Only keeping the parts of the object that correspond to the vaccinated age cohorts, by year and serostatus
#Note: here we are ignoring the coverage (i.e. we look at the cohorts that are vaccinated, but not specifically at vaccinated individuals)
#      If we assume that the vaccine does not protect against infection, the seroprevalence is only affected by the indirect effects, which are the same in vaccinated and unvaccinated

#-- Initialise final objects ------------
seroprev_in_vac_cohorts_wo_vac_final  <-matrix(NA, nrow=timeframe, ncol=5)    #One row per year, one column per serostatus
seroprev_in_vac_cohorts_with_vac_final<-matrix(NA, nrow=timeframe, ncol=5)

#-- Fill the final objects --------------
for(year in 1:timeframe){                      #Loop through each year of the simulation
  
  #Determine youngest and oldest vaccinated cohorts
  age_low_years<-routine_age                   #At any year the youngest vaccinated cohort are those vaccinated as part of the routine
  if(!catch_up_incl){                          #If catch-up is not included...
    age_high_years<-routine_age+year-1         #... the oldest cohort are those who were vaccinated routinely at year 1
  }else{                                       #If catch-up is included...
    age_high_years<-catch_up_age_high+year-1   #... the oldest cohort are those who were the oldest at the moment of catch-up vaccination
  }
  
  age_low_ind<-age_low_years+1                 #Converting age in years into an index (position in the vector) as age starts with 0 and indexing starts with 1
  age_high_ind<-age_high_years+1               
  
  #Fill the final object
  seroprev_in_vac_cohorts_wo_vac_final  [year,]<-apply(seroprev_wo_vac_avg  [year,(age_low_ind:age_high_ind),,drop=F],3,sum)   #Sum for a given year, given age groups, for each serostatus
  seroprev_in_vac_cohorts_with_vac_final[year,]<-apply(seroprev_with_vac_avg[year,(age_low_ind:age_high_ind),,drop=F],3,sum)
  
}

rm(seroprev_wo_vac_avg, seroprev_with_vac_avg); gc()

#-- Save the final objects --------------
colnames(seroprev_in_vac_cohorts_wo_vac_final)  <-c("naive", "history_of_1inf", "history_of_2inf", "history_of_3inf", "history_of_4inf")
colnames(seroprev_in_vac_cohorts_with_vac_final)<-c("naive", "history_of_1inf", "history_of_2inf", "history_of_3inf", "history_of_4inf")

write.csv(seroprev_in_vac_cohorts_wo_vac_final,   file=file.path(path_results_export, paste(ref_scen_wo_vac,   "_seroprev_in_vac_cohorts.csv", sep="")), row.names=TRUE)
write.csv(seroprev_in_vac_cohorts_with_vac_final, file=file.path(path_results_export, paste(ref_scen_with_vac, "_seroprev_in_vac_cohorts.csv", sep="")), row.names=TRUE) 

#----------------------------------------------------------------
#---- Get an average trace for all simulations ------------------
#----------------------------------------------------------------

#-- Create a new object -----------------
#A 5D object "trace" is saved for each simulation with vaccination (if boosting is included)
#Creating a 6D object with an additional dimension for each simulation (year of vaccine introduction)

#Estimating nb of vaccinated cohorts
vac_cohorts_nb<-NA                                                #Total nb of cohorts that will be vaccinated during the simulation
vac_cohorts_nb<-timeframe                                         #At a minimum, the number of vaccinated cohorts will be equal to the nb of simulated years (as one cohort will be routinely vaccinated every year)
if(catch_up_incl){                                                #If catch-up vaccination in the 1st year is included, estimating the additional nb of vaccinated cohorts
  nb_catchup_cohorts<-catch_up_age_high-catch_up_age_low+1
  vac_cohorts_nb<-vac_cohorts_nb+nb_catchup_cohorts               #... adding them to the previoulsy estimated nb
}

#Creating an object
trace_new<-array(NA, dim=c(vac_cohorts_nb,              #Nb vaccinated cohorts
                           5,                           #Nb infections before vaccination (0 to 4)
                           timeframe,                   #Simulation year
                           5,                           #Current nb of infections (0 to 4)
                           timeframe,                   #Nb of years since the last infection ("tunnel" state)
                           nb_start_years))             #Nb simulation years (years of vaccine introduction)
                           

#-- Fill the new object -----------------
for(y in 1:nb_start_years){
  load(file.path(path_results, paste(ref_scen_with_vac, "_y", start_years[y], "_Trace.Rdata", sep="")))   #Loading the object "trace" for the current start year
  trace_new[,,,,,y]<-trace[,,(1:timeframe),,(1:timeframe)]                                                #Saving the object trace to the last dimension of the object trace_new (only taking the "slice" that corresponds to the timeframe of interest)
  rm(trace); gc()
}

#-- Getting average values --------------
trace_avg<-array(NA, dim=c(vac_cohorts_nb, 5, timeframe, 5, timeframe))   
trace_avg<-apply(trace_new, c(1:5), mean)   #Average value for each dimension in the original object trace, for each start year  
rm(trace_new); gc()

#----------------------------------------------------------------
#---- Get seroprevalence in vac cohorts from the trace ----------
#----------------------------------------------------------------
#Estimating the seroprevalence from the trace based on the distribution of subjects by their current number of infections (4th dimension of the object trace_new)

#Initialize and fill a new object
trace_seroprev<-matrix(NA, nrow=timeframe, ncol=5)            #One row per year, one column per serostatus
trace_seroprev<-apply(trace_avg, c(3,4), sum)

#Save
colnames(trace_seroprev)<-c("naive", "history_of_1inf", "history_of_2inf", "history_of_3inf", "history_of_4inf")
write.csv(trace_seroprev,         file=file.path(path_results_export, paste(ref_scen_with_vac, "_seroprev_vac_cohorts_from_trace.csv", sep="")), row.names=TRUE)


#----------------------------------------------------------------
#---- Get distribution by episode from the trace ----------------
#----------------------------------------------------------------

#Initialize a new object
trace_distr_by_episode<-matrix(0, nrow=timeframe, ncol=4*timeframe+1)     #One row per year, one column per episode the individuals are at risk for (1 to 4) and per tunnel year (i.e. how long they've been at risk for this episode); 5th column for those no longer at risk

#Fill the new object
for(nb_inf_before_vac in 0:3){                                #Loop through the "strata" (i.e. categories of vaccinees with different nb of infections before vaccination; ignoring those who had 4 infections before vaccination as they are no longer at risk of a new infection)
  for(nb_inf_current in nb_inf_before_vac:3){                 #Loop through the number of current infections starting with the nb of infections before vac and until a maximum of 3 (those with 4 are no longer at risk)
    nb_episode<-nb_inf_current-nb_inf_before_vac+1            #Calculate the episode for which the individual is at risk
    
    nb_inf_before_vac_ind<-nb_inf_before_vac+1                #Index in the array (as the nb of inf starts with 0 and indexing starts with 1)
    nb_inf_current_ind<-nb_inf_current+1
    
    col_first<-(nb_episode-1)*timeframe+1
    col_last<-col_first+timeframe-1
    
    trace_distr_by_episode[,(col_first:col_last)]<-trace_distr_by_episode[,(col_first:col_last)]+apply(trace_avg[,nb_inf_before_vac_ind,,nb_inf_current_ind, ,drop=FALSE],c(3,5),sum)   #Estimating the nb of individuals with the current combinations of nb of inf before vaccination and now, for each simulation year; adding this to the previously estimating values. Here, making a sum across "tunnel" state to facilitate the interpretation of the results
  }
}

trace_distr_by_episode[,(4*timeframe+1)]<-apply(trace_avg[,,,5, ,drop=FALSE],3,sum)    #Adding into the last column those who currently have 4 infections (no matter how many they had before being vaccinated)

#Save
write.csv(trace_distr_by_episode, file=file.path(path_results_export, paste(ref_scen_with_vac, "_distr_by_episode_from_trace.csv", sep="")), row.names=TRUE)

#----------------------------------------------------------------
#---- Get prob sympt from the trace -----------------------------
#----------------------------------------------------------------
#Distribution by episode only allows estimating the effectiveness in all vaccinated cohorts over time
#It does not take into account that the underlying probability will vary as well (e.g. as a vaccinated cohort becomes older, their chance of symptomatic/hospitalized disease will go down)
#Estimating the probability of symptomatic (hospitalized and non-hospitalized) in vaccinated cohorts over time, taking into account the underlying probability and the efficacy

#Initialize new objects
prob_sympt_vac_cohorts<-matrix(0, nrow=timeframe, ncol=3)                  #Probability of symptomatic infection, one row per simulation year, one column per category (non-hospitalized, hospitalized, total)
colnames(prob_sympt_vac_cohorts)<-c("non_hosp", "hosp", "total")

#Estimate the probability of sympt non-hospit and sympt hosp, by infection nb
prob_sympt_non_hosp_by_inf_nb<-rep(NA, 4)
prob_sympt_hosp_by_inf_nb    <-rep(NA, 4)

prob_sympt_non_hosp_by_inf_nb[1]<-parms_epi$parms_epi_other$prop_sympt_DENV1_unvac[1]*(1-parms_epi$parms_epi_other$prop_hospit_sympt_DENV1_unvac[1])
prob_sympt_non_hosp_by_inf_nb[2]<-parms_epi$parms_epi_other$prop_sympt_DENV1_unvac[2]*(1-parms_epi$parms_epi_other$prop_hospit_sympt_DENV1_unvac[2])
prob_sympt_non_hosp_by_inf_nb[3]<-parms_epi$parms_epi_other$prop_sympt_DENV1_unvac[3]*(1-parms_epi$parms_epi_other$prop_hospit_sympt_DENV1_unvac[3])
prob_sympt_non_hosp_by_inf_nb[4]<-parms_epi$parms_epi_other$prop_sympt_DENV1_unvac[3]*(1-parms_epi$parms_epi_other$prop_hospit_sympt_DENV1_unvac[3])

prob_sympt_hosp_by_inf_nb[1]<-parms_epi$parms_epi_other$prop_sympt_DENV1_unvac[1]*parms_epi$parms_epi_other$prop_hospit_sympt_DENV1_unvac[1]
prob_sympt_hosp_by_inf_nb[2]<-parms_epi$parms_epi_other$prop_sympt_DENV1_unvac[2]*parms_epi$parms_epi_other$prop_hospit_sympt_DENV1_unvac[2]
prob_sympt_hosp_by_inf_nb[3]<-parms_epi$parms_epi_other$prop_sympt_DENV1_unvac[3]*parms_epi$parms_epi_other$prop_hospit_sympt_DENV1_unvac[3]
prob_sympt_hosp_by_inf_nb[4]<-parms_epi$parms_epi_other$prop_sympt_DENV1_unvac[3]*parms_epi$parms_epi_other$prop_hospit_sympt_DENV1_unvac[3]


for(year in 1:timeframe){                              #Loop through every year in the timeframe

  #Calculate the population at risk
  pop_at_risk<-sum(trace_avg[,,year,(1:4),])           #At risk = currently have 0 to 3 infections --> levels 1 to 4 in the 4th dimensions); ignoring the nb of infections before vac
  
  #Loop through every combination of inf before vac and now
  for(nb_inf_before_vac in 0:3){                                #Loop through the "strata" (i.e. categories of vaccinees with different nb of infections before vaccination; ignoring those who had 4 infections before vaccination as they are no longer at risk of a new infection)
    for(nb_inf_current in nb_inf_before_vac:3){                 #Loop through the number of current infections starting with the nb of infections before vac and until a maximum of 3 (those with 4 are no longer at risk)
      
      #Determine positionin the array
      nb_inf_before_vac_ind<-nb_inf_before_vac+1                #Index in the array (as the nb of inf starts with 0 and indexing starts with 1)
      nb_inf_current_ind<-nb_inf_current+1
      
      #Determine episode nb
      nb_episode<-nb_inf_current-nb_inf_before_vac+1
      
      #Determine current eff (depends on the serostatus, but ignoring the serotypes)
      if(nb_inf_before_vac==0){        #If seronegative at vaccination
        eff_non_hosp<-parms_vac_eff$parms_vac_eff$eff_sympt_non_hosp[1,1,nb_episode,]   #Taking 1st level of 1st dimension (DENV-1), 1st level of 2nd dimension (seronegative), current episode nb --> vector with as many elements as years in the timeframe (tunnel states)
        eff_hosp    <-parms_vac_eff$parms_vac_eff$eff_sympt_hosp    [1,1,nb_episode,]
      }else{                           #If seropositive at vaccination
        eff_non_hosp<-parms_vac_eff$parms_vac_eff$eff_sympt_non_hosp[1,2,nb_episode,]   #Same as above, but 2nd level of 2nd dimension (seropositive)
        eff_hosp    <-parms_vac_eff$parms_vac_eff$eff_sympt_hosp    [1,2,nb_episode,]
      }
      
      #Determine underlying prob of sympt (depends on the current nb of infections)
      prob_sympt_non_hosp<-prob_sympt_non_hosp_by_inf_nb[nb_inf_current+1]              #+1 because nb of infections starts with zero
      prob_sympt_hosp    <-prob_sympt_hosp_by_inf_nb[nb_inf_current+1]                  #+1 because nb of infections starts with zero
      
      #Determine prob of sympt, after eff is applied
      prob_sympt_non_hosp_final<-rep(NA, timeframe)
      prob_sympt_hosp_final    <-rep(NA, timeframe)
      
      prob_sympt_non_hosp_final<-prob_sympt_non_hosp*(1-eff_non_hosp)
      prob_sympt_hosp_final    <-prob_sympt_hosp*(1-eff_hosp)
      
      #Get the proportion of individuals with the current combination of infections before and now (as % of the total vaccinated population at risk)
      pop_prop<-apply(trace_avg[,nb_inf_before_vac_ind,year,nb_inf_current_ind, ,drop=FALSE], 5, sum)/pop_at_risk                                   #Sum across all vaccinated cohorts, for a given combination of infections before and now, for each tunnel year --> vector with as many elements as years in the timeframe (year since last boosting)
      
      #Calculate final weighted effectiveness and add it to the previously calculated values
      prob_sympt_vac_cohorts[year,1]<-prob_sympt_vac_cohorts[year,1]+sum(pop_prop*prob_sympt_non_hosp_final)
      prob_sympt_vac_cohorts[year,2]<-prob_sympt_vac_cohorts[year,2]+sum(pop_prop*prob_sympt_hosp_final)
      
      #Final column is the sum of non-hosp and hosp
      prob_sympt_vac_cohorts[year,3]<-prob_sympt_vac_cohorts[year,1]+prob_sympt_vac_cohorts[year,2]
    }
  }
}

#Save
write.csv(prob_sympt_vac_cohorts, file=file.path(path_results_export, paste(ref_scen_with_vac, "_prob_sympt_vac_cohorts.csv", sep="")), row.names=TRUE)
