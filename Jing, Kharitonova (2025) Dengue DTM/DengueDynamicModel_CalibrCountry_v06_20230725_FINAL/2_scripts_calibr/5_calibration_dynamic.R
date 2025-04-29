##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  SCRIPT NAME: 5_calibration_dynamic.R
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA,   Olivier CRISTEAU , Aurélien JAMOTTE
###########################################################################################################################################################################################################

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Function that runs dynamic calibration and saves its results

CalibrDyn_Run<-function(){
  
  #--------- Prepare the population & transmission inputs ------------
  parms_epi<-ParmsEpi_Initialise(parms_epi_selected=parms_epi_selected, calibration=TRUE)         #Calculate and format the final population & transmission inputs 
     
  #-------------------------------------------------------------------
  #--------- Run dynamic calibration ---------------------------------
  #-------------------------------------------------------------------
  #Set the iteration counter
  optim_dyn_iteration<<-1
  
  #Initialize an empty vector for the parameters being calibrated
  parms_to_calibrate<-c()
  
  #Set the flags for whether betaVH and seasonality are being calibrated (required to position the corresponding elements in the vector parms_to_calibrate)
  calibr_betaVH_incl<<-0
  calibr_season_incl<<-0
  
  #-- IF MULTIPLE PARMS --------
  #If seasonality parameters are being calibrated (i.e. multiple parameters being calibrated at the same time - betaVH AND seasonality parms)
  if(calibr_seasonality){
    
    #Initialize a vector with the parameters to be calibrated
    if(calibr_betaVH){
      calibr_betaVH_incl<<-1  
      parms_to_calibrate<-c(parms_to_calibrate, 
                            parms_epi$parms_epi_dyn_model$betaVH[1])  #Add initial value of betaVH to the vector with parameters to calibrate
    }
    if(calibr_seasonality){
      calibr_season_incl<<-1
      parms_to_calibrate<-c(parms_to_calibrate, 
                            parms_epi$parms_epi_dyn_model$season_p1,  #Same logic as for betaVH above 
                            parms_epi$parms_epi_dyn_model$season_p2)  
    }

    #Run optimisation (using function "optim" for multi-dimensional optimization)
    runtime_calibr_dyn<-system.time(calibr_dyn_out<-optim(parms_to_calibrate,                 #Vector containing the values of the parameters to be calibrated
                                                          CalibrDyn_Target_Multiple,          #Target function to be minimized during the optimization 
                                                          parms_epi=parms_epi,                #Population & transmission inputs
                                                          control=list(trace=0,               #Optimization parameters
                                                                       reltol=optim_reltol)))
  #-- IF ONLY BETAVH -----------
  #If only betaVH is being calibrated
  }else{
    if(save_SP_data_only){
      calibr_dyn_out_xxx<-CalibrDyn_Target_Single(betaVH=parms_epi$parms_epi_dyn_model$betaVH[1],                    
                                                  parms_epi=parms_epi)
    }else{
      runtime_calibr_dyn<-system.time(calibr_dyn_out<-optimize(CalibrDyn_Target_Single,         #Target function to be minimized during the optimization
                                                               interval=c(0,1),                 #Range of possible values of betaVH  
                                                               parms_epi=parms_epi)) 
    }
  }
  
  #-------------------------------------------------------------------
  #--- Once optimization is over, save the final results -------------
  #-------------------------------------------------------------------
  #Nb parameters that were being calibrated
  nb_calibr_parms<-0
  if(calibr_betaVH){nb_calibr_parms<-nb_calibr_parms+1}
  if(calibr_seasonality){nb_calibr_parms<-nb_calibr_parms+2}
   
  #------------------------------------------------------------------- Adding the best-fitting values to the input file is useful for the next steps of the calibration (e.g. when betaVH is calibrated and we want to calibrate seasonality given best-fitting beta)
  #--------- Add best-fitting values to the input file --------------- It's also required to do the final run to output the initial state and later to run the model
  #-------------------------------------------------------------------
  #Import inputs from the file "inputs_calibr"
  parms<-read.csv(file=file.path(path_data_calibr, paste(epi_file_name, ".csv", sep="")), stringsAsFactors=FALSE)
  parms_list<-as.list(parms)
  
  #Replace the values for betaVH and seasonality if they were being calibrated
  if(nb_calibr_parms==1){
    parms_list$betaVH<-rep(calibr_dyn_out$minimum,4)                    #Repeating four times (same betaVH for each serotype)
  }else{
    if(calibr_betaVH){
      parms_list$betaVH<-rep(calibr_dyn_out$par[1],4)                   #Repeating four times (same betaVH for each serotype)
    }
    if(calibr_seasonality){
      parms_list$season_p1<-calibr_dyn_out$par[calibr_betaVH_incl+1]
      parms_list$season_p2<-calibr_dyn_out$par[calibr_betaVH_incl+2]
    }
  }
  
  #Save the modified file
  write.csv(data.frame(makePaddedDataFrame(parms_list)), file=file.path(path_data_calibr, paste(epi_file_name, ".csv", sep="")),row.names=FALSE) 

}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Target function that is minimized during optimisation (if multiple parameters are being calibrated)

CalibrDyn_Target_Multiple<-function(parms_to_calibrate,        #Vector with the parameters to be calibrated (four values of betaVH, two values of seasonal parameters)
                                    parms_epi){                #Population & transmission inputs

  #----- Prepare the parameters for the current run ------------------         #Plug the current values of the calibrated parameters into the list with population & transmission inputs

  #-- betaVH ---------
  if(calibr_betaVH){
    parms_epi$parms_epi_dyn_model$betaVH[]<-rep(parms_to_calibrate[1],4)         #Repeating the value of betaVH four times (for each serotype)
  }
  #-- Seasonality ----
  if(calibr_seasonality){
    parms_epi$parms_epi_dyn_model$season_p1<-parms_to_calibrate[1+calibr_betaVH_incl]
    parms_epi$parms_epi_dyn_model$season_p2<-parms_to_calibrate[2+calibr_betaVH_incl]
  }
  
  #--------- Create an initial state for vector population -----------       
  Vini_asymmetric<-c(rep(0,9))                                                                        #Vector with nine zeroes (for nine compartments used for mosquito population)
  Vini_asymmetric[1]<-(parms_epi$parms_epi_other$NH*parms_epi$parms_epi_dyn_model$ratio_VH-1)         #All (but one) mosquitoes are assumed to be susceptible
  Vini_asymmetric[6:9]<-prop_serotype                                                                 #Introducing a total of one infectious mosquito (for all four serotypes) and distributing it across serotypes based on the proportions specified in the script "0_calibration_launch.R"
  
  coeff_Jan1=(1+parms_epi$parms_epi_dyn_model$season_p1*sin(2*3.14159265358979*(1/365)+parms_epi$parms_epi_dyn_model$season_p2))    #Applying seasonality. Proportion of vectors present on January 1st with the given values of seasonal parameters. Note: In the base case the calibration for archetypes does not include seasonality, but this option is kept in case it is needed
  Vini_Jan1=Vini_asymmetric*coeff_Jan1 
  
  #--------- Create an initial state for host population -------------
  Hini_fully_naive<-array(rep(0,N_age_groups*5*5*5*5), dim=c(N_age_groups,5,5,5,5))       #Create an array where each element is the nb of host in a specific age group with a specific status for each serotype; populate with zeros
  Hini_fully_naive[,1,1,1,1]<-parms_epi$parms_epi_other$NHi                               #For each age group, put all the hosts of this age into the "cell" with status (1,1,1,1), which means fully susceptible to all serotypes. This indicates that the entire human population is dengue-naive
  
  #----- Prepare the initial state for the current run ---------------
  yini<-c(as.vector(Hini_fully_naive),                                   #Final vector with the initial state includes all the compartments for hosts (transformed from array form to vector form),
          Vini_Jan1,                                                     #Vector with 9 elements representing the number of mosquitoes on January 1
          rep(0, N_age_groups*4*3))                                      #Additional vector used to record incidence (one element per age group, serotype and infection type - primary, secondary, post-secondary)
  
  #----- Run the simulation -------------------------------------------------------------------- #The model is being run first for the burn-in period, then for the analysis period
  out_2nd_period<-Simulation_Run_TwoPeriods(yini=yini,
                                            parms=parms_epi) 
  
  #----- Split the simulation results ---------------------------------------------------------- #The ODE solver output is split into categories (incidence, seroprevalence, etc.)
  out_2nd_period_split<-Simulation_SplitResults(out=out_2nd_period$out_wo_hosts,
                                                parms=parms_epi,
                                                hosts_included=FALSE,
                                                keep_first=FALSE)         
  
  #----- Aggregate the results required for the calculation of likelihood ----------------------
  data_for_likelihood<-CalibrDyn_FormatResults(out_split=out_2nd_period_split,
                                               parms_epi=parms_epi)
  
  
  #----- Calculate log-likelihood --------------------------------------------------------------
  loglike<-CalibrDyn_CalcLogLike(incidence_by_age_predicted=data_for_likelihood$incidence_rate_by_age,
                                 incidence_by_age_observed=parms_epi$parms_epi_other$incidence_by_age_observed,
                                 incidence_by_month_predicted=data_for_likelihood$incidence_rate_by_month,
                                 incidence_by_month_observed=parms_epi$parms_epi_other$incidence_by_month_observed)
  if(is.na(loglike)){
    loglike=100000          #If could not calculate log-likelihood, setting it to a high positive value
  }
  
  #----- Save the results (if the best run so far) ---------------------------------------------
  if((optim_dyn_iteration==1)||                                      #If the first run (calibration iteration) or...
     (optim_dyn_iteration!=1 && loglike<best_loglike)){              #... not the first run but better fit than before
    
    best_loglike<<-loglike
    write.csv(data_for_likelihood$seroprevalence_2states,    file=file.path(path_data_calibr,paste("Best_fit_seroprevalence.csv",           sep="")), row.names=FALSE)
    write.csv(data_for_likelihood$incidence_rate_by_age,     file=file.path(path_data_calibr,paste("Best_fit_incidence_rate_by_age.csv",    sep="")), row.names=FALSE)
    write.csv(data_for_likelihood$incidence_rate_by_month,   file=file.path(path_data_calibr,paste("Best_fit_incidence_rate_by_month.csv",  sep="")), row.names=FALSE)
    write.csv(data_for_likelihood$nb_inf_yearly_by_serotype, file=file.path(path_data_calibr,paste("Best_fit_nb_inf_by_serotype.csv",       sep="")), row.names=FALSE)
    
    init_state<-list(Hini=out_2nd_period$out_hosts_array, Vini=out_2nd_period$out_vectors)       #Saving the initial state file
    save(init_state, file=file.path(path_data_calibr, paste(epi_file_name, "_init_state_y1.Rdata", sep="")))
    
  } 
  
  #----- Print the information on the current run ----------------------------------------------
  print("")
  print(paste("ITERATION ", optim_dyn_iteration, sep=""))
  if(calibr_betaVH){
    print(paste("betaVH: ", parms_to_calibrate [1]))
  }else{
    print(paste("betaVH: ", parms_epi$parms_epi_dyn_model$betaVH[1]))
  }
  if(calibr_seasonality){
    print(paste("season_p1: ", parms_to_calibrate [1+calibr_betaVH_incl]))
    print(paste("season_p2: ", parms_to_calibrate [2+calibr_betaVH_incl]))   
  }
  
  print(paste("Likelihood: ", loglike, "(Best so far: ", best_loglike, ")", sep="" ))
  
  optim_dyn_iteration<<-optim_dyn_iteration+1
  
  #------ Return positive log-likelihood (value to be minimized) ------------------------------- 
  return(loglike)
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Target function that is minimised during optimisation (if only betaVH is being calibrated)

CalibrDyn_Target_Single<-function(betaVH,                    #Vector with the parameters to be calibrated (four values of betaVH, two values of seasonal parameters)
                                  parms_epi){                #Population & transmission inputs
  
   
  #----- Prepare the parameters for the current run ------------------         #Plug the current values of the calibrated parameters into the list with population & transmission inputs

  #BetaVH
  parms_epi$parms_epi_dyn_model$betaVH[]<-rep(betaVH,4)                        #Repeating the value of betaVH four times (for each serotype)
  # message("betaVH DENV1: ", parms_epi$parms_epi_dyn_model$betaVH[1])
  # message("betaVH DENV2: ", parms_epi$parms_epi_dyn_model$betaVH[2])
  # message("betaVH DENV3: ", parms_epi$parms_epi_dyn_model$betaVH[3])
  # message("betaVH DENV4: ", parms_epi$parms_epi_dyn_model$betaVH[4])
  
  #--------- Create an initial state for vector population -----------       
  Vini_asymmetric<-c(rep(0,9))                                                                        #Vector with nine zeroes (for nine compartments used for mosquito population)
  Vini_asymmetric[1]<-(parms_epi$parms_epi_other$NH*parms_epi$parms_epi_dyn_model$ratio_VH-1)         #All (but one) mosquitoes are assumed to be susceptible
  Vini_asymmetric[6:9]<-prop_serotype                                                                 #Introducing a total of one infectious mosquito (for all four serotypes) and distributing it across serotypes based on the proportions specified in the script "0_calibration_launch.R"
  
  
  coeff_Jan1=(1+parms_epi$parms_epi_dyn_model$season_p1*sin(2*3.14159265358979*(1/365)+parms_epi$parms_epi_dyn_model$season_p2))    #Applying seasonality. Proportion of vectors present on January 1st with the given values of seasonal parameters. Note: In the base case the calibraiton for archetypes does not include seasonality, but this option is kept in case it is needed
  Vini_Jan1=Vini_asymmetric*coeff_Jan1 
  
  #--------- Create an initial state for host population -------------
  Hini_fully_naive<-array(rep(0,N_age_groups*5*5*5*5), dim=c(N_age_groups,5,5,5,5))       #Create an array where each element is the nb of host in a specific age group with a specific status for each serotype; populate with zeros
  Hini_fully_naive[,1,1,1,1]<-parms_epi$parms_epi_other$NHi                               #For each age group, put all the hosts of this age into the "cell" with status (1,1,1,1), which means fully susceptible to all serotypes. This indicates that the entire human population is dengue-naive
  
  
  #----- Prepare the initial state for the current run ---------------
  yini<-c(as.vector(Hini_fully_naive),                                   #Final vector with the initial state includes all the compartments for hosts (transformed from array form to vector form),
          Vini_Jan1,                                                     #Vector with 9 elements representing the number of mosquitoes on January 1
          rep(0, N_age_groups*4*3))                                      #Additional vector used to record incidence (one element per age group, serotype and infection type - primary, secondary, post-secondary)
  
  #----- Run the simulation -------------------------------------------------------------------- #The model is being run first for the burn-in period, then for the analysis period
  out_2nd_period<-Simulation_Run_TwoPeriods(yini=yini,
                                            parms=parms_epi) 
  
  #----- Split the simulation results ---------------------------------------------------------- #The ODE solver output is split into categories (incidence, seroprevalence, etc.)
  out_2nd_period_split<-Simulation_SplitResults(out=out_2nd_period$out_wo_hosts,
                                                parms=parms_epi,
                                                hosts_included=FALSE,
                                                keep_first=FALSE)         
  
  #----- Aggregate the results required for the calculation of likelihood ----------------------
  data_for_likelihood<-CalibrDyn_FormatResults(out_split=out_2nd_period_split,
                                               parms_epi=parms_epi)
  
  
  #----- Calculate log-likelihood --------------------------------------------------------------
  loglike=NA
  if(!save_SP_data_only){
  loglike<-CalibrDyn_CalcLogLike(incidence_by_age_predicted=data_for_likelihood$incidence_rate_by_age,
                                 incidence_by_age_observed=parms_epi$parms_epi_other$incidence_by_age_observed,
                                 incidence_by_month_predicted=data_for_likelihood$incidence_rate_by_month,
                                 incidence_by_month_observed=parms_epi$parms_epi_other$incidence_by_month_observed)
  }
  
  if(is.na(loglike)){
    loglike=100000          #If could not calculate log-likelihood, setting it to a high value
  }
  
  #----- Save the results (if the best run so far) ---------------------------------------------
  if((optim_dyn_iteration==1)||                                      #If the first run (calibration iteration) or...
     (optim_dyn_iteration!=1 && loglike<best_loglike)||
      save_SP_data_only){              #... not the first run but better fit than before
    
    best_loglike<<-loglike
    write.csv(data_for_likelihood$seroprevalence_2states,   file=file.path(path_data_calibr,paste("Best_fit_seroprevalence.csv",         sep="")), row.names=FALSE)
    write.csv(data_for_likelihood$incidence_rate_by_age,    file=file.path(path_data_calibr,paste("Best_fit_incidence_rate_by_age.csv",  sep="")), row.names=FALSE)
    write.csv(data_for_likelihood$incidence_rate_by_month,  file=file.path(path_data_calibr,paste("Best_fit_incidence_rate_by_month.csv",sep="")), row.names=FALSE)
    write.csv(data_for_likelihood$nb_inf_yearly_by_serotype,file=file.path(path_data_calibr,paste("Best_fit_nb_inf_by_serotype.csv",     sep="")), row.names=FALSE)

    init_state<-list(Hini=out_2nd_period$out_hosts_array, Vini=out_2nd_period$out_vectors)       #Saving the initial state  
    save(init_state, file=file.path(path_data_calibr, paste(epi_file_name, "_init_state_y1.Rdata", sep="")))
    
  } 
  
  #----- Print the information on the current run ----------------------------------------------
  print("")
  print(paste("ITERATION ", optim_dyn_iteration, sep=""))
  print(paste("betaVH: ", betaVH))
  print(paste("season_p1: ", parms_epi$parms_epi_dyn_model$season_p1))
  print(paste("season_p2: ", parms_epi$parms_epi_dyn_model$season_p2))   
  print(paste("Likelihood: ", loglike, "(Best so far: ", best_loglike, ")", sep="" ))
  
  optim_dyn_iteration<<-optim_dyn_iteration+1
  
  #------ Return positive log-likelihood (value to be minimised) ------------------------------- 
  return(loglike)
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Function that formats the results obtained at a given calibration iteration. These results are then used for the calculation of log-likelihood

CalibrDyn_FormatResults<-function(out_split,        #Results of the iteration split into categories
                                  parms_epi){     
  
  with(out_split,{
    
    #------------- Prepare time variables ----------------------
    times_year<-floor(((times-1)/365+1))             #The year of the simulation for each day
    times_month<-floor((times-1)%%365/(365/12))+1    #The month of the simulation for each day (from 1 to 12)
    
    times_year_month<-(unique(cbind(times_year, times_month)))[,1]    #The year number of repeated for each month included in this year
    times_months_year<-(unique(cbind(times_year, times_month)))[,2]   #The month number is repeated for each year
    
    #------------- Calculate seroprevalence --------------------
    #Calculate average seroprevalence (as proportion of the age group)
    
	  seroprevalence_mean=apply(seroprevalence, c(2,3), mean)/apply(NHit, 2, mean)  
    
    #Aggregate for two serostatuses
    seroprevalence_2states=cbind(seroprevalence_mean[,1],                           #Fully susceptible hosts
                                 apply(seroprevalence_mean[,2:8],1,sum))            #Hosts with the history of at least one infection
  
    #------------- Calculate incidence, by serotype ------------
    nb_inf_daily_by_serotype<-apply(incidence, c(1,3),sum)                                       #Nb infections, for each day & serotype, sum across all ages & types
    nb_inf_yearly_by_serotype<-aggregate(nb_inf_daily_by_serotype, list(year=times_year),sum)    #Same, but yearly
    colnames(nb_inf_yearly_by_serotype)<-c("Year", "Serotype 1", "Serotype 2", "Serotype 3", "Serotype 4")
    
    #------------- Calculate incidence, by age ------------
    p_sympt_hospit<-parms_epi$parms_epi_other$prop_sympt_unvac*parms_epi$parms_epi_other$prop_hospit_sympt_unvac          #Proportion symptomatic hospitalised (of all infections)
    p_sympt_non_hospit<-parms_epi$parms_epi_other$prop_sympt_unvac*(1-parms_epi$parms_epi_other$prop_hospit_sympt_unvac)  #Proportion symptomatic non hospitalised (of all infections)
    incidence_sympt_hosp   <-array(NA,dim=dim(incidence))     #Array with symptomatic hospitalized incidence by time, age, serotype and type
    incidence_sympt_non_hosp<-array(NA,dim=dim(incidence))    #Array with symptomatic non-hospitalized incidence by time, age, serotype and type
    for(year in 1:(dim(incidence)[1])){                       #Note: loop variable is named "year", but the object incidence actually contains the nb of infectios by day
      incidence_sympt_hosp[year,,,]<-incidence[year,,,]*p_sympt_hospit
      incidence_sympt_non_hosp[year,,,]<-incidence[year,,,]*p_sympt_non_hospit
    }

    #------------- Calculate sympt  reported incidence, by age --
    nb_age_groups<-length(parms_epi$parms_epi_other$incidence_by_age_observed[,1])             #Nb of age groups used to provide observed incidence rate
    incidence_rate_by_age<-rep(NA, nb_age_groups)                                              #Vector that will be used to record the predicted incidence rate of symptomatic dengue
    
    for(age_group in 1:nb_age_groups){
      
      #Determine the start and end age of the age group
      age_low<-parms_epi$parms_epi_other$incidence_by_age_observed[age_group,1]+1              #Lower bound of the current age group (+1 as age starts with zero and indexing starts with 1); first column of the object "observed incidence" and row corresponding to the current age group
      age_high<-parms_epi$parms_epi_other$incidence_by_age_observed[age_group,2]+1             #Upper bound of the current age group
      
      #Nb of symptomatic infections in the current age groups, by day
      nb_inf_daily_sympt_hosp     <-apply(incidence_sympt_hosp[,(age_low:age_high),,], 1, sum)       #Nb of hospitalized dengue infections in the current age group (sum for each day across the ages of interest, all serotypes and types) --> vector with as many elements as days in the 2nd simulation period
      nb_inf_daily_sympt_non_hosp  <-apply(incidence_sympt_non_hosp[,(age_low:age_high),,], 1, sum)  #Same for symptomatic non-hospitalized dengue
      
      #Nb of symptomatic infections, yearly
      nb_inf_yearly_sympt_hosp    <-aggregate(nb_inf_daily_sympt_hosp,    list(year=times_year), sum)[,-1]   #Vector with as many elements as years in the 2nd simulation period
      nb_inf_yearly_sympt_non_hosp <-aggregate(nb_inf_daily_sympt_non_hosp, list(year=times_year), sum)[,-1]

      #Nb of hospitalized cases that would be reported, by type (annual average across the 2nd simulation period)
      nb_sympt_mean_sympt_hosp     <-mean(nb_inf_yearly_sympt_hosp)/parms_epi$parms_epi_other$expansion_factor_hospit   #To obtain the number of cases that would be reported, dividing the predicted nb of cases by the expansion factors (real cases = reported cases * expansion factor)
      
      #Nb of non-hospitalized cases that would be reported (same as above)
      nb_sympt_mean_sympt_non_hosp <-mean(nb_inf_yearly_sympt_non_hosp)/parms_epi$parms_epi_other$expansion_factor_non_hospit         
    
      #Total average annual nb of symptomatic reported cases
      nb_sympt_mean<-sum(nb_sympt_mean_sympt_hosp,
                         nb_sympt_mean_sympt_non_hosp)
      
      #Estimating the population size in the age groups
      age_group_pop<-sum(parms_epi$parms_epi_other$NHi[age_low:age_high])
      
      #Estimating the incidence rate per 100k (or whatever denominator was used for the observed incidence)
      incidence_rate<-nb_sympt_mean/age_group_pop*parms_epi$parms_epi_other$incidence_observed_per
      
      #Recording the value in the final vector
      incidence_rate_by_age[age_group]<-incidence_rate
    }
    
    #------------- Calculate incidence, by month ---------------
    #Nb of infections, by day
    nb_inf_daily_sympt_hosp    <-0   #Reset the variables
    nb_inf_daily_sympt_non_hosp<-0

    nb_inf_daily_sympt_hosp    <-apply(incidence_sympt_hosp,    1, sum)   #Nb of symptomatic hospitalized dengue infections --> vector with as many elements as days in the 2nd simulation period
    nb_inf_daily_sympt_non_hosp<-apply(incidence_sympt_non_hosp, 1, sum)  #Same for symptomatic non hospitalized

    #Nb symptomatic infections, by month
    nb_inf_monthly_sympt_hosp    <-aggregate(nb_inf_daily_sympt_hosp,    list(month=times_month, year=times_year), sum)[,-c(1,2)]
    nb_inf_monthly_sympt_non_hosp<-aggregate(nb_inf_daily_sympt_non_hosp, list(month=times_month, year=times_year), sum)[,-c(1,2)]

    #Average nb of symptomatic infections, by month
    nb_inf_average_monthly_sympt_hosp    <-aggregate(nb_inf_monthly_sympt_hosp,    list(month=times_months_year), mean)[,-1]
    nb_inf_average_monthly_sympt_non_hosp<-aggregate(nb_inf_monthly_sympt_non_hosp, list(month=times_months_year), mean)[,-1]
    
    #Average nb of REPORTED hospitalized symptomatic infections, by month
    nb_inf_monthly_average_hosp_reported    <-nb_inf_average_monthly_sympt_hosp/parms_epi$parms_epi_other$expansion_factor_hospit
    
    #Average nb of REPORTED non-hospitalized symptomatic infections, by month
    nb_inf_monthly_average_non_hosp_reported<-nb_inf_average_monthly_sympt_non_hosp/parms_epi$parms_epi_other$expansion_factor_non_hospit
    
    #Total average nb of REPORTED symptomatic infections (hosp & non-hosp combined)
    nb_inf_monthly_average_total_reported   <-nb_inf_monthly_average_hosp_reported+nb_inf_monthly_average_non_hosp_reported
    
    #Average incidence rate per month
    total_pop<-parms_epi$parms_epi_other$NH
    incidence_rate_by_month<-nb_inf_monthly_average_total_reported/total_pop*parms_epi$parms_epi_other$incidence_observed_per
    
    
    #------------- Return data for the calculation of log-likelihood -
    return(list(times_months_year=times_months_year,
                times_year_month=times_year_month, 
                seroprevalence_2states=seroprevalence_2states, 
                incidence_rate_by_age=incidence_rate_by_age,
                incidence_rate_by_month=incidence_rate_by_month,
                nb_inf_yearly_by_serotype=nb_inf_yearly_by_serotype))
  })
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Function that calculates positive log-likelihood based on the formatted results of the calibration iteration

CalibrDyn_CalcLogLike<-function(incidence_by_age_predicted,
                                incidence_by_age_observed,
                                incidence_by_month_predicted,
                                incidence_by_month_observed){
  
  #---- Calculate log-likelihood for incidence ---------------       
  #The log-likelihood is calculated by comparing observed rate of symptomatic dengue to the average rate predicted in the second simulation period
  
  #Initialize the values
  log_like_by_age_pos<-0
  log_like_by_age_neg<-0
  
  #If betaVH was calibrated, estimate log-likelihood
  if(calibr_betaVH){
    for(age_group in 1:length(incidence_by_age_observed[,1])){
      loglike_age_group<-incidence_by_age_observed[age_group,3]*log(incidence_by_age_predicted[age_group])-incidence_by_age_predicted[age_group]
      log_like_by_age_pos<-log_like_by_age_pos+loglike_age_group
    }
    log_like_by_age_neg<-(-log_like_by_age_pos) #With this method log-likelihood is positive --> converting to negative to that it can be minimized  
  }
  
  #---- Calculate log-likelihood for seasonality -------------       
  log_like_by_month_pos<-0
  log_like_by_month_neg<-0
  
  if(calibr_seasonality){
    log_like_by_month_pos<-sum(incidence_by_month_observed*log(incidence_by_month_predicted)-incidence_by_month_predicted)
    log_like_by_month_neg<-(-log_like_by_month_pos)
  }
  
  #---- Aggregate final log-likelihood -----------------------       
  loglike<-log_like_by_age_neg+log_like_by_month_neg
  
  return(loglike)
}


#####################################################################################
#These functions is used to export inputs. As the columns in the input file have different length, the NAs are added to make them the same length
na.pad<-function(x,len){x[1:len]}

makePaddedDataFrame<-function(l,...){
  l<-l%>% replace(.=="NULL", NA)
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}
