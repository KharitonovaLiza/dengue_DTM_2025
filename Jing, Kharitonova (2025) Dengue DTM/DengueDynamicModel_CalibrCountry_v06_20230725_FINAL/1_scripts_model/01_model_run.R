##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contains functions that run the model and display progress window
###########################################################################################################################################################################################################

############################################################
#This function runs the model
Model_Run<-function(input, output, session,
                    queue){

  #------------------------------------------------------------------
  #--- Set  integration parameters ----------------------------------
  #------------------------------------------------------------------
  ode_method<-"ode23"            #"ode23" - an integration method for ODE systems using 2nd and 3rd order Runge-Kutta-Fehlberg formulas with automatic step size
                                 #"euler" - an integration method with fixed step size
  ode_hini<-0                    #Initial step size to be attempted. If 0, the initial step size is determined automatically by solvers with flexible time step
  ode_atol<-1e-06                #Absolute error tolerance. Value of 1e-06 means an absolute error of 0.000001 is tolerated for each compartment during the integration.
                                 #As the population consist of thousands of people, this level of accuracy was deemed sufficient
  #------------------------------------------------------------------
  #--- Run the model for each simulation in the run list ------------
  #------------------------------------------------------------------
  if(via_interface){                        #Number of simulations to be run (1 if testing the model outside the interface or as many as specified in the interface)
    nb_sim<-length(queue[,"scen_ref"])      
  }else{
    nb_sim<-1 
  }
  
  run_success<-rep(0, nb_sim)               #Vector with as many elements as simulations; if the value is 0 the model run failed; if the value is 1, the run was successfull 

  if(via_interface){progress_matrix<-ModelRun_UpdateProgress(queue=queue, progress_matrix=NULL, current_sim=0, current_status=NULL)}
  
  for(sim_nb in 1:nb_sim){

    if(via_interface){progress_matrix<-ModelRun_UpdateProgress(queue=queue, progress_matrix=progress_matrix, current_sim=sim_nb, current_status="Loading parameters")}
    if(via_interface){GetStackInfo()}
    
    #------------------------------------------------------------------
    #--- Start error handling -----------------------------------------
    #------------------------------------------------------------------
    #It is set to 0 at the simulation start.The value is then changed to 1 if the simulation is successfully completed.
    #If for any reason the simulation could not be completed, the value will remain 0 ( results will not be saved in the model interface).
    run_success[sim_nb]<-tryCatch({  
      
      #------------------------------------------------------------------
      #--- Get scenario parameters --------------------------------------    #All parameters describing a scenario are stored as global variables
      #------------------------------------------------------------------
      if(via_interface){
        parms_gen<-list(                                                   #If running the model via interface, then taking the values from the object "queue"...
          scen_ref=queue[sim_nb,"scen_ref"],                               #...Scenario reference in the format Scen_XXXXX; unique identified for each combination of input sets; results are saved under this name
          set_epi=queue[sim_nb,"set_epi"],                                 #...Reference of the input set with population & transmission inputs
          set_vac_eff=queue[sim_nb,"set_vac_eff"],                         #...Reference of the input set with vaccine effects inputs (always provided, even if a scenario w/o vaccination)
          set_vac_strat=queue[sim_nb,"set_vac_strat"],                     #...Reference of the input set with vaccination strategy inputs (always provided, even if a scenario w/o vaccination)
          set_vec_eff=queue[sim_nb,"set_vec_eff"],                         #...Reference of the input set with vector control inputs (always provided, even if a scenario w/o vector control)
          set_cost_ttt=queue[sim_nb,"set_cost_ttt"],                       #...Reference of the input set with cost of treatment inputs
          set_cost_vac=queue[sim_nb,"set_cost_vac"],                       #...Reference of the input set with vaccination cost inputs (always provided, even if a scenario w/o vaccination)
          set_cost_vec=queue[sim_nb,"set_cost_vec"],                       #...Reference of the input set with vector control cost (always provided, evene if a scenario w/o vector control)
          set_QoL=queue[sim_nb,"set_QoL"],                                 #...Reference of the input set with quality-of-life inputs
          scen_ref_res_epi=queue[sim_nb,"scen_ref_res_epi"],               #...Reference of a scenario with the same combination of population & transmission, vaccination and vector control inputs (if this combination was already run, the results are copied and only the eco module is run; may remain empty)
          vac_switch=as.logical(queue[sim_nb,"vac_switch"]),               #...Flag for vaccination inclusion (TRUE or FALSE)
          vec_switch=as.logical(queue[sim_nb,"vec_switch"]),               #...Flag for vector control inclusion (TRUE or FALSE)
          run_epi=as.logical(queue[sim_nb,"run_epi"]),                     #...Flag to run the transmission module (TRUE or FALSE)
          save_init_state=as.logical(queue[sim_nb,"save_init_state"]),     #...Flag to save initial states (if a model w/o vaccination with an epi set that hasn't been run before)  
          year_start=as.integer(queue[sim_nb,"year_start"]),               #...Year of simulation start (year 1 for scenarios without interventions; year of intervention start otherwise)
          timeframe=as.integer(queue[sim_nb,"timeframe"]))                 #...Simulation time frame
      }else{
        parms_gen<-list(                                                   #Otherwise taking the values specified manuallyin the script "test_wo_interface.R" (folder "0_scripts_misc")
          scen_ref=scen_ref,          
          set_epi=set_epi,           
          set_vac_eff=set_vac_eff,       
          set_vac_strat=set_vac_strat,     
          set_vec_eff=set_vec_eff,       
          set_cost_ttt=set_cost_ttt,
          set_cost_vac=set_cost_vac,
          set_cost_vec=set_cost_vec,
          set_QoL=set_QoL,
          scen_ref_res_epi=scen_ref_res_epi,
          vac_switch=as.logical(vac_switch),
          vec_switch=as.logical(vec_switch),
          run_epi=as.logical(run_epi),
          save_init_state=as.logical(save_init_state),
          year_start=as.integer(year_start),
          timeframe=as.integer(timeframe)
        )
      }
      
      if(via_interface){GetStackInfo()}
      
      #------------------------------------------------------------------
      #--- Load & create model parameters -------------------------------
      #------------------------------------------------------------------

      #Population & transmission parameters
      parms_epi<-ParmsEpi_Load(path_inputs_epi, parms_gen$set_epi, calibr=FALSE)   
      parms_mod<-parms_epi 
      
      if(via_interface){GetStackInfo()}
      
      #Vaccine effects & vaccination strategy parameters (if vaccination is included & transmission module needs to be run)
      if(parms_gen$vac_switch && parms_gen$run_epi){                                                                            
        parms_vac_strat<-ParmsVacStrat_Load(path_inputs_vac_strat, parms_gen$set_vac_strat, parms_epi)
        parms_vac_eff<-ParmsVacEff_Load(path_inputs_vac_eff, parms_gen$set_vac_eff, parms_vac_strat, parms_gen$timeframe)
        parms_mod<-c(parms_mod, parms_vac_eff, parms_vac_strat)
        rm(parms_vac_eff, parms_vac_strat); gc()
      }
      
      if(via_interface){GetStackInfo()}
      
      #Vector control effects (if vector control is included & transmission module needs to be run)
      if(parms_gen$vec_switch && parms_gen$run_epi){    
        parms_vec_eff<-ParmsVecEff_Load(path_inputs_vec_eff, parms_gen$set_vec_eff)
        parms_mod<-c(parms_mod, parms_vec_eff)
        rm(parms_vec_eff); gc()
      }
      rm(parms_epi); gc()
       
      if(via_interface){GetStackInfo()}
      
      #Additional parameters
      N_age_groups<-parms_mod$parms_epi_other$N_age_groups          #Number of age groups used to run the model
      NstatesH<-N_age_groups*5*5*5*5                                #Size of the system for hosts (nb. of compartments per one vaccination status)
      if(!parms_gen$vac_switch){                                    #Nb. of levels for vaccination status
        vac_levels<-1                                               #No vaccination --> 1 status (unvac)
      }else{
        vac_levels<-3                                               #Vaccination --> 3 statuses (unvac, vac seronegative, vac seropositive)
      }      
      parms_gen<-c(parms_gen, N_age_groups=N_age_groups, NstatesH=NstatesH, vac_levels=vac_levels)
      
      #------------------------------------------------------------------
      #--- Determine start & end year -----------------------------------
      #------------------------------------------------------------------
      year_start<-parms_gen$year_start
      year_end<-year_start+parms_gen$timeframe-1
      
      if(via_interface){GetStackInfo()}
      
      #------------------------------------------------------------------
      #--- Run the transmission module (if needed) ----------------------
      #------------------------------------------------------------------    
      if(parms_gen$run_epi){
        
        #------------------------------------------  #The extra compartments are used to record:
        #--- Nb. of extra compartments ------------  #-Nb. of hosts, by vaccination status & by age --> 1 col. per age per included vaccination status
        #------------------------------------------  #-Vector pop. size --> 1 col.
        #                                            #-Nb. of hosts, by serostatus & by age --> 8 col. per age group (susceptible; history of one inf, by serotype; history of 2 inf; history of 3 inf; history of 4 inf.)   
        nout<-(8+parms_gen$vac_levels)*N_age_groups+1  
        
        #------------------------------------------
        #--- Load the C model ---------------------
        #------------------------------------------
        if(!parms_gen$vac_switch){script_name<-"M1_no_vac"}else{script_name<-"M2_vac"}        #Name of the C script (w/o extension)
        dll_name<-paste(script_name,.Platform$dynlib.ext, sep="")                             #Name of the dll file
        dyn.load(file.path(path_scripts_src, dll_name), local=TRUE, now=FALSE)                #Load the model
       	if(via_interface){GetStackInfo()}
       	
       	#------------------------------------------  #Vector times contains the simulation days, starting from the last day of the year immediately before the first simulation year
        #--- List of time points ------------------  #(to account for the initial state) and until the last day of the last simulation year
        #------------------------------------------
        times_sim<-((year_start-1)*365):(year_end*365)         
        
        #------------------------------------------  #Initial state is obtained from the calibration process. It reflects the number of hosts and vectors in each compartment under the given epidemiological settings (i.e. input set epi)
        #--- Prepare initial state ----------------  #Initial states are stored in the folder with population & transmission inputs under the name "Epi_XXXXX_init_state_yYY.Rdata" where Y is the year
        #------------------------------------------  #Each such object contains two sub-objects:
        #                                            # - A 5D array with the size of each compartment for the hosts on the day before the simulation starts
        #                                            # - A vector with the size of each compartment for the vectors on the day before the simulation starts
        initpath<-file.path(path_inputs_epi, paste(parms_gen$set_epi, "_init_state_y", year_start, ".Rdata", sep=""))
        if(via_interface){GetStackInfo()}
        
        #Load the initial state (system size for vectors & for unvaccinated hosts)
        load(initpath)
        if(via_interface){GetStackInfo()}
        
        #Add an extra dimension for vaccination status (if a model with vaccination) & transform the array into a vector
        #Since there are no vaccinated hosts at the simulation start, all new elements in the array are equal to zero
        if(!parms_gen$vac_switch){                                #Model without vaccination --> no need for an extra dimension (no compartments for vaccinated individuals in the dynamic model)
          Hini_new<-array(0,dim=c(N_age_groups,5,5,5,5))
          Hini_new<-init_state$Hini
          Hini_vect<-as.vector(Hini_new)
        } else{                                                   #Model with vaccination --> an extra dimension with three levels (one per vaccination status)
          Hini_new<-array(0,dim=c(N_age_groups,5,5,5,5,vac_levels))
          Hini_new[,,,,,1]<-init_state$Hini
          Hini_vect<-c(as.vector(Hini_new[,,,,,1]), as.vector(Hini_new[,,,,,2:vac_levels]))
        }    
        if(via_interface){GetStackInfo()}
        
        #Create a matrix to keep track of newly vaccinated and screened individuals (by age)
        #In a model w/o vaccination this matrix will remain empty, but it is required by the design of the function AgeingVaccination
        new_vac<-matrix(0, parms_gen$timeframe, (N_age_groups*2))        	        
        colnames(new_vac)<-c(apply(cbind(rep("new_screen_ag", N_age_groups),1:N_age_groups), 1, paste, collapse=""), apply(cbind(rep("new_vac_ag", N_age_groups),1:N_age_groups), 1, paste, collapse="")) 
        rownames(new_vac)<-c(year_start:year_end)
        
        #Apply ageing & vaccination 
        post_age_vac<-AgeingVaccination(parms_mod=parms_mod, parms_gen=parms_gen, year=0, yH=Hini_vect, new_vac=new_vac, ageing_flag=TRUE)
        new_vac=post_age_vac$new_vac
        if(via_interface){GetStackInfo()}
        
        #Create vector for recording incidence
        inc_vect<-rep(0,N_age_groups*4*3*vac_levels)          #Incidence is recorded for each age group, each vaccination status, four serotype and three infection types (primary, secondary, post-secondary)
        
        #Aggregate the final vector with the initial state
        yini<-c(post_age_vac$yH, init_state$Vini, inc_vect)
        
        #Delete the initial state
        rm(init_state); gc()
        
        if(via_interface){GetStackInfo()}
        
        if(parms_gen$vac_switch){

          #------------------------------------------ 
          #--- Prepare the "Markov trace" -----------
          #------------------------------------------ 
          #If vaccination & efficacy boosting are included, this object is used to keep track of how long ago the individuals had their last boosting (vaccination or natural infection). 
          #With each boosting, the efficacy against symptoms (with or without hospitalization) is changed (goes to the next efficacy curve)
          #Without boosting, the efficacy values will be drawn from the same curve, for the year that corresponds to how long ago the individual was vaccinated
          #The object "trace" has 5 dimensions:
          #  1/ Nb of the vac cohort (cohorts of different age will be vaccinated and will experience breakthrough infections at different periods in time) --> nb levels = nb cohorts that are vaccinated (see below)
          #  2/ Nb infections before vaccination (0 to 4) --> nb levels = 5
          #  3/ Simulation year --> nb levels = nb years in the timeframe 
          #  4/ Current nb of infections (0 to 4) --> nb levels = 5 
          #     (distinction between the current nb of infections and nb of infections before vaccination is require to differentiate episodes post-vaccination; 
          #      for example the 2nd natural infection will be the second episode post-vaccination in someone who had no infections before being vaccinated; in someone who had one infection before being vaccinated, it will be the first episode post-vaccination. Etc.)
          #  5/ Nb of years since the last infection --> nb levels = nb years in the timeframe (equivalent of a tunnel state in a Markov trace)
          #The values in the object represent the nb of people in each state.
          #The values in the array represent the state at the beginning of each year. Where they are located in the tunnel states indicate which year of efficacy is applicable to them
          #For example, right after vaccination everyone is in the tunnel state #1 (for the corresponding nb of infections), they have efficacy of year 1
          #Next year they can either go to year 1 of the next infection (and have efficacy of year 1 of the next efficacy curve) or go to the year 2 of the same infection (and have efficacy of year 2 of the same efficacy curve)
          
          #NOTE: The structure of the object "trace" is illustrated in the Excel file "Structure of the object trace" located in the folder "4_docs"
          
          #--- Vaccinated cohorts -------------------
          #The 1st dimension of the objects described above is the nb of the vaccinated cohort
          #Multiple cohorts will be vaccinated during the simulations - some of them in the 1st year (routine and catch-up, if included) and some in other years (new routine cohorts)
          #Therefore, creating several objects:
          
          #TOTAL NB OF VACCINATED COHORTS
          vac_cohorts_nb<-NA                                                #Total nb of cohorts that will be vaccinated during the simulation
          vac_cohorts_nb<-parms_gen$timeframe                               #At a minimum, the number of vaccinated cohorts will be equal to the nb of simulated years (as one cohort will be routinely vaccinated every year)
          if(parms_mod$parms_vac_strat$catchup_switch){                     #If catch-up vaccination in the 1st year is included, estimating the additional nb of vaccinated cohorts
            age_catchup_high<-parms_mod$parms_vac_strat$age_catchup_high
            age_catchup_low <-parms_mod$parms_vac_strat$age_catchup_low
            nb_catchup_cohorts<-age_catchup_high-age_catchup_low+1
            vac_cohorts_nb<-vac_cohorts_nb+nb_catchup_cohorts               #... adding them to the previoulsy estimated nb
          }
          
          #AGE AT VACCINATION & YEAR OF VACCINATION
          #These two vectors have as many elements as vac_cohorts_nb and list, for each cohort, the age at vaccination and the year of vaccination
          #In the vectors the cohorts are ordered as follows: routine vaccination at year 1, catch-up vaccination at year 1 (if included), routine vaccination in the following years
          vac_cohorts_age_at_vac <-rep(NA, vac_cohorts_nb)
          vac_cohorts_year_of_vac<-rep(NA, vac_cohorts_nb)
          
          if(parms_mod$parms_vac_strat$catchup_switch){                                                           #If catch-up is included
            vac_cohorts_age_at_vac[1]<-parms_mod$parms_vac_strat$age_routine                                      #Cohort vaccinated as part of routine vaccination at year 1 (1st element of the vectors); age = routine
            vac_cohorts_year_of_vac[1]<-1                                                                         #...year = 1
            
            vac_cohorts_age_at_vac[2:(nb_catchup_cohorts+1)]<-c(age_catchup_low:age_catchup_high)                 #Cohorts vaccinated as part of catch-up vaccination at year 1 (starting with element 2 of the vector); listing all ages covered by catch-up
            vac_cohorts_year_of_vac[2:(nb_catchup_cohorts+1)]<-1                                                  #...year = 1
            
            vac_cohorts_age_at_vac[(nb_catchup_cohorts+2):vac_cohorts_nb]<-parms_mod$parms_vac_strat$age_routine  #Cohorts vaccinated as part of routine vaccination in the following years; age = routine
            vac_cohorts_year_of_vac[(nb_catchup_cohorts+2):vac_cohorts_nb]<-c(2:parms_gen$timeframe)              #... years = 2 and until the end of the simulated timeframe
          }else{   
            vac_cohorts_age_at_vac<-rep(parms_mod$parms_vac_strat$age_routine, parms_gen$timeframe)               #If catch-up is not included, all cohorts are vaccinated at the same age (age of routine vaccination)
            vac_cohorts_year_of_vac<-c(1:parms_gen$timeframe)                                                     #... years of vaccination are 1 to the end of the simulated timeframe
          }
          
          #--- Trace object ------------------------
          trace<-array(0, dim=c(vac_cohorts_nb,              #Nb vaccinated cohorts
                                5,                           #Nb infections before vaccination (0 to 4)
                                parms_gen$timeframe,         #Simulation year 
                                5,                           #Current nb of infections (0 to 4)
                                parms_gen$timeframe))        #Nb of years since the last infection ("tunnel" state)
          
          #--- New infections counts ---------------
          new_inf<-array(NA, dim=c(vac_cohorts_nb,           #Nb of new infections, for each vaccinated cohort, by simulation year and infection number (1st to 4th)
                                   parms_gen$timeframe,
                                   4))  
        }
        
        #------------------------------------------
        #--- Create objects to record parms -------
        #------------------------------------------
        #Some of the parameters related to vaccination (efficacy, proportion of symptomatic cases in vaccinated etc.) are estimated at each year, as simulation progresses
        #Some of them are required in later stages of the code, so creating larger objects to save these values generated at each year
        if(parms_gen$vac_switch){
          
          #Objects for efficacies
          eff_inf           <-array(NA, dim=c(101,4,2,parms_gen$timeframe))        #4D arrays with efficacy/effectiveness values. Dimensions are age (single year), serotype, serostatus at vaccination, simulation year
          eff_asympt        <-array(NA, dim=c(101,4,2,parms_gen$timeframe)) 
          eff_sympt_non_hosp<-array(NA, dim=c(101,4,2,parms_gen$timeframe)) 
          eff_sympt_hosp    <-array(NA, dim=c(101,4,2,parms_gen$timeframe))
          
          #Objects for severity distributions in vaccinated subjects
          prop_asympt_vac        <-array(NA, dim=c(N_age_groups,4,2,3,parms_gen$timeframe))   #% asymptomatic infections, by age, serotype, serostatus, infection type and simulation year
          prop_sympt_non_hosp_vac<-array(NA, dim=c(N_age_groups,4,2,3,parms_gen$timeframe))   #% non-hospitalized symptomatic infections
          prop_sympt_hosp_vac    <-array(NA, dim=c(N_age_groups,4,2,3,parms_gen$timeframe))   #% hospitalized symptomatic infections
          prop_sympt_vac         <-array(NA, dim=c(N_age_groups,4,2,3,parms_gen$timeframe))   #% total symptomatic infections (required for the calculation of average infectiousness)
          prop_severe_vac        <-array(NA, dim=c(N_age_groups,4,2,3,parms_gen$timeframe))   #% severe infections
        }
        
        #------------------------------------------
        #--- Simulate each year -------------------
        #------------------------------------------
        for(year in 1:parms_gen$timeframe){
          
          #message("Running simulation for year ", year)
          
          if(via_interface){progress_matrix<-ModelRun_UpdateProgress(queue=queue, progress_matrix=progress_matrix, current_sim=sim_nb, current_status=paste("Running transmission module (year ", year, "/", parms_gen$timeframe, ")", sep=""))}
          
          #--- Determine time points -------------------------
          #366 days are simulated at each loop - from Dec 31 of the previous year and until Dec 31 of the current year
          #The system size on Dec 31 is recorded and used in discrete processes (ageing & vaccination)
          #The resulting system size is then used as the initial state for the next year
          times_year<-times_sim[((year-1)*365+1):(year*365+1)]
          
          #--- Prepare model parameters ----------------------
          #List parms_C includes all the parameters that are passed to C codes 
          #The order of these parameters is crucial for C codes to run. If the order is changed, then C codes need to be adapted & re-compiled
          
          #Population & transmission parameters (always used in C codes)
          parms_C<-c(parms_mod$parms_epi_dyn_mod)
          if(via_interface){GetStackInfo()}
          
          #Vector control parameters
          if(parms_gen$vec_switch){
            parms_C$b<-parms_C$b*parms_mod$parms_vec_eff$vec_b_change
            parms_C$ratio_VH<-parms_C$ratio_VH*parms_mod$parms_vec_eff$vec_ratio_VH_change
            if(via_interface){GetStackInfo()}
          }
          
          #If vaccination is included
          if(parms_gen$vac_switch){
            
            #Convert the initial state for this year (for hosts only) into an array
            yini_6D<-array(yini[1:(vac_levels*NstatesH)], dim=c(N_age_groups,           #Nb age groups
                                                                5,5,5,5,                #Nb statuses for each serotype
                                                                vac_levels))            #Nb vac statuses

            #Get all the parameters for this year (theta, psy, proportion of symptomatic cases in vaccinated subjects etc.)
            parms_vac_eff_year<-ParmsVacEff_GetThetaAndPsyYear(year=year,                                        #Current simulation year
                                                               parms_epi=parms_mod$parms_epi_other,              #Population & transmission inputs
                                                               parms_vac_eff=parms_mod$parms_vac_eff,            #Vaccine efficacy inputs (previously aggregated and processed)
                                                               parms_vac_strat=parms_mod$parms_vac_strat,        #Vaccination strategy inputs
                                                               yini_6D=yini_6D,                                  #Array with the initial state for the hosts for the current year (post-ageing, post-vaccination)
                                                               trace=trace,                                      #Markov trace (required to estimate the efficacy if efficacy boosting is included),
                                                               vac_cohorts_nb=vac_cohorts_nb,                    #Nb of vaccinated cohorts,
                                                               vac_cohorts_age_at_vac=vac_cohorts_age_at_vac,    #Vector indicating age at vaccination, for each vaccinated cohort
                                                               vac_cohorts_year_of_vac=vac_cohorts_year_of_vac)  #Vector indicating calendar year of vaccination, for each vaccinated cohort  
            
            #Add objects for this year to the larger objects that will be used later
            eff_inf           [,,,year]<-parms_vac_eff_year$efficacies$eff_inf
            eff_asympt        [,,,year]<-parms_vac_eff_year$efficacies$eff_asympt
            eff_sympt_non_hosp[,,,year]<-parms_vac_eff_year$efficacies$eff_sympt_non_hosp
            eff_sympt_hosp    [,,,year]<-parms_vac_eff_year$efficacies$eff_sympt_hosp
            
            prop_asympt_vac        [,,,,year]<-parms_vac_eff_year$severities_vac$prop_asympt_vac
            prop_sympt_non_hosp_vac[,,,,year]<-parms_vac_eff_year$severities_vac$prop_sympt_non_hosp_vac
            prop_sympt_hosp_vac    [,,,,year]<-parms_vac_eff_year$severities_vac$prop_sympt_hosp_vac
            prop_sympt_vac         [,,,,year]<-parms_vac_eff_year$severities_vac$prop_sympt_vac
            prop_severe_vac        [,,,,year]<-parms_vac_eff_year$severities_vac$prop_severe_vac
            
            parms_vac_severities<-list(prop_asympt_vac=prop_asympt_vac,
                                       prop_sympt_non_hosp_vac=prop_sympt_non_hosp_vac,
                                       prop_sympt_hosp_vac=prop_sympt_hosp_vac,
                                       prop_sympt_vac=prop_sympt_vac,
                                       prop_severe_vac=prop_severe_vac)
            
            #If this is the last year, also add severity parameters to parms_mod so that they can be passed to further code
            if(year==parms_gen$timeframe){
              parms_mod<-c(parms_mod, parms_vac_severities)
            }
            
            #Add objects theta and psy to the list of parameters passed to the C code
            parms_C<-c(parms_C,
                       list(theta_vac_neg=parms_vac_eff_year$theta_vac_neg,           #theta for the hosts vaccinated as seronegative,
                            theta_vac_pos=parms_vac_eff_year$theta_vac_pos,           #theta for the hosts vaccinated as seropositive, for the current simulation year
                            psy_DENV1_vac_neg=parms_vac_eff_year$psy_DENV1_vac_neg,   #psy for the hosts vaccinated as seronegative, for the current simulation year
                            psy_DENV2_vac_neg=parms_vac_eff_year$psy_DENV2_vac_neg,
                            psy_DENV3_vac_neg=parms_vac_eff_year$psy_DENV3_vac_neg,
                            psy_DENV4_vac_neg=parms_vac_eff_year$psy_DENV4_vac_neg,
                            psy_DENV1_vac_pos=parms_vac_eff_year$psy_DENV1_vac_pos,   #psy for the hosts vaccinated as seropositive, for the current simulation year
                            psy_DENV2_vac_pos=parms_vac_eff_year$psy_DENV2_vac_pos,
                            psy_DENV3_vac_pos=parms_vac_eff_year$psy_DENV3_vac_pos,
                            psy_DENV4_vac_pos=parms_vac_eff_year$psy_DENV4_vac_pos))
            if(via_interface){GetStackInfo()}
          }
          
          #------ Run the ODE solver -------------------------
          out<-try(rk(y=yini,                    #Initial state vector (host compartments, vector compartments, incidence compartments)
                      times=times_year,          #Vector of time points to be estimated (day numbers)
                      parms=parms_C,             #List of parameters to be passed to C code
                      func="model",              #Name of the model function (as specified in the C code)
                      dllname=script_name,       #Name of the dll file (w/o extension)
                      initfunc="initmod",        #Name of the function initializing the model parameters (as specified in the C code)  
                      nout=nout,                 #Number of additional outputs to be recorded (length of the corresponding vector)  
                      method=ode_method,         #Integration method
                      hini=ode_hini,             #Initial step size to be attempted
                      atol=ode_atol),            #Absolute error tolerance
                   silent=TRUE)             
          
          #------ Get new initial state ----------------------
          #The object out has as many rows as estimation time points (elements in the vector times_year; normally 366 in a year).
          #It has the following columns:
          #-Times (nb columns = 1)
          #-Host compartments (nb columns = vac_levels*NstatesH)
          #-Vector compartments (nb columns = 9)
          #-Vector with recorded incidence (nb columns = N_age_groups*4*3*vac_levels)
          #-Vector yout with additional output (nb columns = nout)
          
          #The last row of the object out is taken to be used in discrete processes (and then as the initial state for the next year)
          #Only the columns with host compartments, vector compartments and incidence records are taken
          yini<-(out[length(out[,1]), 2:(1+vac_levels*NstatesH+9+N_age_groups*4*3*vac_levels)])  
          
          #Also, if the flag save_init_state is on, then the compartments for hosts and vectors are saved as of Dec 31 (last day of the current simulation year and before applying ageing)
          #These initial states are then used to run the model with interventions, which is run from the moment of intervention introduction using the appropriate initial state
          #This allows shortening runtime for the model with interventions (which are heavy due to a great number of compartments)
          #Note: the last line is saved from year 1 and until the end of the period used to test the intervention introduction 
          #comes from the calibration process (part "Calculation of initial state"). The length of the period used to test intervention introduction is set in the Excel file
          #list_model_timeframes stored in the folder "Data/00_lists"
          
          if(parms_gen$save_init_state && (year>=1) && (year<=(list_model_timeframes$timeframe_intervention_introduction-1))){
            save_yH_1D<-yini[1:(vac_levels*NstatesH)]                                        #Host compartments
            save_yV_1D<-yini[(vac_levels*NstatesH+1):(vac_levels*NstatesH+9)]                #Vector compartmenst
            
            save_yH_5D<-array(save_yH_1D, dim=c(N_age_groups,5,5,5,5))                       #Host compartments in vector form
            init_state<-list(Hini=save_yH_5D, Vini=save_yV_1D)
            save(init_state, file=file.path(path_inputs_epi, paste(parms_gen$set_epi, "_init_state_y", (year+1), ".Rdata", sep="")))    #Save as initial state for the NEXT simulation year
            rm(save_yH_1D,save_yV_1D,save_yH_5D,init_state); gc()
          }
          
          #------ Age & vaccinate ----------------------------
          post_age_vac<-AgeingVaccination(parms_mod=parms_mod, parms_gen=parms_gen, year=year, yH=yini[1:(vac_levels*NstatesH)], new_vac=new_vac, ageing_flag=TRUE)
          yini[1:(vac_levels*NstatesH)]<-post_age_vac$yH        #Inject new system for hosts into the initial state for the next year
          new_vac=post_age_vac$new_vac                          #Save the updated matrix with the nb of vaccinated/screened
          rm(post_age_vac);gc()
          
          #------ Fill the traces ----------------------------  #Doing this only if both vaccination and efficacy boosting are included
          #Based on the output from the solver, filling the "Markov trace" to obtain the distribution of each vaccinated cohort by their boosting status
          if(parms_gen$vac_switch){
            
            #Fill the trace
            res_fill_trace<-Model_FillTrace(trace=trace,                                     #Object with the "trace" that need to be filled
                                            new_inf=new_inf,                                 #Objects used to record the nb of new infections by number (1st-4th) and by year
                                            year=year,                                       #Nb of the current year
                                            vac_cohorts_nb=vac_cohorts_nb,                   #Nb of vaccinated cohorts,
                                            vac_cohorts_age_at_vac=vac_cohorts_age_at_vac,   #Vector indicating age at vaccination, for each vaccinated cohort
                                            vac_cohorts_year_of_vac=vac_cohorts_year_of_vac, #Vector indicating calendar year of vaccination, for each vaccinated cohort  
                                            timeframe=parms_gen$timeframe,                   #Simulation timeframe
                                            N_age_groups=N_age_groups,                       #Nb of age groups in the model,
                                            vac_levels=vac_levels,                           #Nb of vaccination statuses,
                                            parms_mod=parms_mod,                             #List of model parameters
                                            out_hosts=out[,2:(1+vac_levels*NstatesH)])       #Size of host compartments for each simulation day within the current year loop
            
            #Replace the objects for the next loop
            trace<-res_fill_trace$trace
            new_inf<-res_fill_trace$new_inf        
          }
          
          #------ If the last year, save results -------------
          #If vaccination & efficacy boosting are included, saving the generated effectiveness curves for each vaccinated cohort for later analyses
          if(year==parms_gen$timeframe && parms_gen$vac_switch){
            
            #Object listing all the vaccinated cohorts with the age at and year of vaccination
            vac_cohorts_list<-cbind((1:vac_cohorts_nb),                                  
                                    vac_cohorts_age_at_vac,
                                    vac_cohorts_year_of_vac)
            colnames(vac_cohorts_list)<-c("cohort_nb", "age_at_vac", "year_of_vac")
            vac_cohorts_list_df<-data.frame(vac_cohorts_list)
            #write.csv(vac_cohorts_list_df, file=file.path(path_results, paste(parms_gen$scen_ref, "_VacCohortsList.csv", sep="")), row.names=FALSE)
            
            #Efficacy & proportion symptomatic for each vaccinated cohort
            for(vac_cohort in 1:vac_cohorts_nb){                          #Loop through each vaccinated cohort
              
              age_at_vac_yrs<-vac_cohorts_age_at_vac[vac_cohort]          #Determine age at vaccination 
              age_at_vac_ind<-age_at_vac_yrs+1                            #Age index (position in the list of ages)
              year_of_vac<-vac_cohorts_year_of_vac[vac_cohort]            #Year of vaccination
              
              eff_inf_cohort           <-matrix(NA, nrow=parms_gen$timeframe, ncol=8)   #Matrix where each row is a year since vaccination; columns are efficacy by serotype and serostatus (some rows will remain NA if the cohort was not vaccinated at year 1)
              eff_asympt_cohort        <-matrix(NA, nrow=parms_gen$timeframe, ncol=8)
              eff_sympt_non_hosp_cohort<-matrix(NA, nrow=parms_gen$timeframe, ncol=8)
              eff_sympt_hosp_cohort    <-matrix(NA, nrow=parms_gen$timeframe, ncol=8)
              
              colnames(eff_inf_cohort)<-c("eff_inf_DENV1_neg","eff_inf_DENV2_neg","eff_inf_DENV3_neg","eff_inf_DENV4_neg","eff_inf_DENV1_pos","eff_inf_DENV2_pos","eff_inf_DENV3_pos","eff_inf_DENV4_pos")
              colnames(eff_asympt_cohort)<-c("eff_asympt_DENV1_neg","eff_asympt_DENV2_neg","eff_asympt_DENV3_neg","eff_asympt_DENV4_neg","eff_asympt_DENV1_pos","eff_asympt_DENV2_pos","eff_asympt_DENV3_pos","eff_asympt_DENV4_pos")
              colnames(eff_sympt_non_hosp_cohort)<-c("eff_sympt_non_hosp_DENV1_neg","eff_sympt_non_hosp_DENV2_neg","eff_sympt_non_hosp_DENV3_neg","eff_sympt_non_hosp_DENV4_neg","eff_sympt_non_hosp_DENV1_pos","eff_sympt_non_hosp_DENV2_pos","eff_sympt_non_hosp_DENV3_pos","eff_sympt_non_hosp_DENV4_pos")
              colnames(eff_sympt_hosp_cohort)<-c("eff_sympt_hosp_DENV1_neg","eff_sympt_hosp_DENV2_neg","eff_sympt_hosp_DENV3_neg","eff_sympt_hosp_DENV4_neg","eff_sympt_hosp_DENV1_pos","eff_sympt_hosp_DENV2_pos","eff_sympt_hosp_DENV3_pos","eff_sympt_hosp_DENV4_pos")
              
              prop_sympt_cohort<-matrix(NA, nrow=parms_gen$timeframe, ncol=4*2*3)   #Columns are % symptomatic infections by serotype, serostatus and type (primary, secondary, post-secondary)
              
              year_row<-1
              for(y in year_of_vac:parms_gen$timeframe){                  #Loop through all the years that cohort was vaccinated
                age_at_y_ind<-age_at_vac_ind+y-year_of_vac                #Index of this age group at year y
                
                
                #Efficacy
                eff_inf_cohort[year_row,(1:4)]<-eff_inf[age_at_y_ind,,1,y]                         #Efficacy against inf, by serotype, for seronegatives corresponding to the current year since vaccination and current age
                eff_inf_cohort[year_row,(5:8)]<-eff_inf[age_at_y_ind,,2,y]                         #... for seropositives
                
                eff_asympt_cohort[year_row,(1:4)]<-eff_asympt[age_at_y_ind,,1,y]                   #Similarly for the efficacy against asymptomatic dengue 
                eff_asympt_cohort[year_row,(5:8)]<-eff_asympt[age_at_y_ind,,2,y]   
                
                eff_sympt_non_hosp_cohort[year_row,(1:4)]<-eff_sympt_non_hosp[age_at_y_ind,,1,y]   #Similarly for the efficacy against symptomatic non-hospitalized dengue
                eff_sympt_non_hosp_cohort[year_row,(5:8)]<-eff_sympt_non_hosp[age_at_y_ind,,2,y]   
                
                eff_sympt_hosp_cohort[year_row,(1:4)]<-eff_sympt_hosp[age_at_y_ind,,1,y]           #Similarly for the efficacy against symptomatic hospitalized dengue
                eff_sympt_hosp_cohort[year_row,(5:8)]<-eff_sympt_hosp[age_at_y_ind,,2,y]   
                
                #Proportion symptomatic
                prop_sympt_cohort[year_row,(1:4)]  <-prop_sympt_vac[age_at_y_ind,,1,1,y]           #Seronegative, primary, by serotype
                prop_sympt_cohort[year_row,(5:8)]  <-prop_sympt_vac[age_at_y_ind,,1,2,y]           #Seronegative, secondary, by serotype
                prop_sympt_cohort[year_row,(9:12)] <-prop_sympt_vac[age_at_y_ind,,1,3,y]           #Seronegative, post-secondary, by serotype
                
                prop_sympt_cohort[year_row,(13:16)]<-prop_sympt_vac[age_at_y_ind,,2,1,y]           #Seropositive, primary, by serotype
                prop_sympt_cohort[year_row,(17:20)]<-prop_sympt_vac[age_at_y_ind,,2,2,y]           #Seropositive, secondary, by serotype
                prop_sympt_cohort[year_row,(21:24)]<-prop_sympt_vac[age_at_y_ind,,2,3,y]           #Seropositive, post-secondary, by serotype
                
                
                year_row<-year_row+1     #Update the row, to which to write the results, before the next loop
              }
              
              eff_cohort<-cbind(eff_inf_cohort, eff_asympt_cohort, eff_sympt_non_hosp_cohort, eff_sympt_hosp_cohort)
              eff_cohort_df<-data.frame(eff_cohort)
              
              prop_sympt_cohort_df<-data.frame(prop_sympt_cohort)
              
              #write.csv(eff_cohort_df,        file=file.path(path_results, paste(parms_gen$scen_ref, "_y", year_start, "_Efficacy_Cohort",  vac_cohort, ".csv", sep="")), row.names=FALSE)
              #write.csv(prop_sympt_cohort_df, file=file.path(path_results, paste(parms_gen$scen_ref, "_y", year_start, "_PropSympt_Cohort", vac_cohort, ".csv", sep="")), row.names=FALSE)
              
            }
            
            #Also saving the objects with traces and nb of infections
            #save(trace,   file=file.path(path_results, paste(parms_gen$scen_ref, "_y", year_start, "_Trace.Rdata", sep="")))
            #save(new_inf, file=file.path(path_results, paste(parms_gen$scen_ref, "_y", year_start, "_NewInf.Rdata", sep="")))
            
          }
          
          #------ Output results -----------------------------
          #In the object out_final the columns with the host compartments are discarded. 
          #For the first simulation year all the rows are taken. For the following years the first row is discarded, as it corresponds to the last day of the previous year
          #Note: the first row of the object out_final is the initial state, which should not be taken into account in the aggregation of the results
          
          if(year==1){
            out_final=out[,c(1,(vac_levels*NstatesH+2):(length(out[1,])))]                          #If first year --> create object out_final
          }else{
            out_final=rbind(out_final, out[-1,c(1,(vac_levels*NstatesH+2):(length(out[1,])))])      #Otherwise, add new results to the previously saved results
          }
          rm(out); gc()
          
        }
        
        model_out=list(res_daily=out_final,               #Simulation results, by day
                       new_vac=new_vac)                   #Number of newly vaccinated and screened hosts, by year
        
        #------------------------------------------
        #--- Aggregate & save epi results ---------
        #------------------------------------------
        #Split simulation results into categories (incidence, pop size, seroprevalence etc.)
        
        results_epi<-Simulation_SplitResults(res_daily=model_out$res_daily,
                                             parms_mod=parms_mod,
                                             parms_gen=parms_gen,
                                             keep_init_state=FALSE)                  #Flag to keep or not the first day of the simulation (represents the initial state)
        
        #Add information on the number of newly vaccinated/screened individuals
        results_epi<-c(results_epi,list(new_vac=model_out$new_vac))
        rm(model_out);gc()
        
        #Save results as R object & delete 
        save(results_epi, file=file.path(path_results, paste(parms_gen$scen_ref, "_y", year_start, "_Epi.Rdata", sep="")))
        rm(results_epi); gc()
        message("Results saved")
        
        #------------------------------------------
        #--- Unload C model -----------------------
        #------------------------------------------
        dyn.unload(file.path(path_scripts_src, dll_name))
        
      #------------------------------------------------------------------ #If the trasmission module does not need to run (parms_gen$run_epi=FALSE), this means that the results already exist for this combination of population & transmsission,
      #--- If the transmission module is not run, load old results ------ #vaccine effects, vaccination strategy and vector control effects. The reference of the scenario with this combination is identified when launching the model from the interface.
      #------------------------------------------------------------------ #It is saved as parms_gen$scen_ref scen_ref_res_epi. The epi results for this scenario are then loaded and saved under a new scenario reference    
      }else{
        
        #Load the already existing results (loaded as object results_epi)
        load(file=file.path(path_results, paste(parms_gen$scen_ref_res_epi, "_y", year_start, "_Epi.Rdata", sep="")))  
        
        #"Cut out" the required timeframe for each object
        nb_inf<-results_epi$nb_inf[1:parms_gen$timeframe,,,,,, drop=F]
        pop_host<-results_epi$pop_host[1:parms_gen$timeframe,,, drop=F]
        
        pop_vector<-data.frame(results_epi$pop_vector$pop_vector[1:parms_gen$timeframe])
        colnames(pop_vector)<-"pop_vector"
        
        seroprev<-results_epi$seroprev[1:parms_gen$timeframe,,, drop=F]
        new_vac<-results_epi$new_vac[1:parms_gen$timeframe,]
        rm(results_epi); gc()
        
        #Aggregate into a list & delete individual objects
        results_epi<-list(nb_inf=nb_inf,
                          pop_host=pop_host,
                          pop_vector=pop_vector,
                          seroprev=seroprev,
                          new_vac=new_vac)
        rm(nb_inf, pop_host, pop_vector, seroprev, new_vac); gc()
        
        save(results_epi, file=file.path(path_results, paste(parms_gen$scen_ref, "_y", year_start, "_Epi.Rdata", sep="")))
        rm(results_epi); gc()
      }
      
      #------------------------------------------------------------------ #Unlike the transmission module that is run only when needed, the economic module is run every time
      #--- Run economic module (always run) ----------------------------- 
      #------------------------------------------------------------------
      if(via_interface){progress_matrix<-ModelRun_UpdateProgress(queue=queue, progress_matrix=progress_matrix, current_sim=sim_nb, current_status="Running economic module")}
      
      #--- Load cost parameters ----------   #The required cost parameters are loaded and added to the lsit parms_mod
      #Cost of treatment (always loaded)
      parms_cost_ttt<-ParmsCostTtt_Load(path_inputs_cost_ttt, parms_gen$set_cost_ttt)                              
      parms_mod<-c(parms_mod, parms_cost_ttt)
      rm(parms_cost_ttt); gc()
      
      #Vaccination cost (if vaccination is included)
      if(parms_gen$vac_switch){
        parms_cost_vac<-ParmsCostVac_Load(path_inputs_cost_vac, parms_gen$set_cost_vac)
        parms_mod<-c(parms_mod, parms_cost_vac)
        rm(parms_cost_vac); gc()
      } 
      
      #Vector control cost (if vector control is included)   
      if(parms_gen$vec_switch){
        parms_cost_vec<-ParmsCostVec_Load(path_inputs_cost_vec, parms_gen$set_cost_vec)
        parms_mod<-c(parms_mod, parms_cost_vec)
        rm(parms_cost_vec); gc()
      }
      
      #--- Load QoL parameters -----------
      parms_QoL<-ParmsQoL_Load(path_inputs_QoL, parms_gen$set_QoL)
      parms_mod<-c(parms_mod, parms_QoL)
      rm(parms_QoL); gc()
      
      #--- Load & prepare epi results ----
      results_epi_for_eco<-ResultsEpi_PrepForEco(parms_mod, parms_gen, year_start)
      
      #--- Calculate costs & nb cases ----
      results_eco_temp<-ResultsEco_CalculateAll(parms_mod, parms_gen, results_epi_for_eco)
      
      #--- Prepare & save final output ---
      results_eco<-ResultsEco_PrepareOutput(parms_mod, parms_gen, results_epi_for_eco, results_eco_temp)
      rm(results_epi_for_eco, results_eco_temp); gc()
      save(results_eco, file=file.path(path_results, paste(parms_gen$scen_ref, "_y", year_start, "_Eco.Rdata", sep="")))
      rm(results_eco); gc()
      
      #------------------------------------------------------------------
      #--- Set the run success flag to 1 --------------------------------
      #------------------------------------------------------------------
      if(via_interface){progress_matrix<-ModelRun_UpdateProgress(queue=queue, progress_matrix=progress_matrix, current_sim=sim_nb, current_status="Successfully completed")}
      run_success[sim_nb]<-1
    }
    #------------------------------------------------------------------
    #--- Or to 0 if the model failed ----------------------------------
    #------------------------------------------------------------------
    ,error=function(e){
      
      #Print a user message to the console
      message("The simulation could not be completed due to the following error:")
      message(e)
      
      #Set success/failure flag to zero
      run_success[sim_nb]<-0
      return(F)
    }  
    #------------------------------------------------------------------
    #--- Force next simulation if the model failed --------------------
    #------------------------------------------------------------------   
    ,finally={
      #Update the status
      if(via_interface){
        if(run_success[sim_nb]==0){progress_matrix<-ModelRun_UpdateProgress(queue=queue, progress_matrix=progress_matrix, current_sim=sim_nb, current_status="Simulation failed")}
      }
      #Go onto the next simulation
      next
    })  
  }
  errmsg<-geterrmessage()
  message("Last error: ",errmsg)
  if(via_interface){
    return(list(progress_matrix=progress_matrix,
                run_success=run_success))
  }
}

############################################################
#This function fills the "Markov trace" which is described above (in the function Model_Run) and is required for the calculation of effectiveness values (i.e. taking into account the boosting of efficacy with each breakthrough infection)
Model_FillTrace<-function(trace,                    #5D object that needs to be filled
                          new_inf,                  #3D object with the nb of new infections (1st to 4th), for each vaccinated cohort and each simulation year
                          year,                     #Nb of the currently simualted year (the trace is filled at the end of the year for the next simulation year)
                          vac_cohorts_nb,           #Nb of vaccinated cohorts,
                          vac_cohorts_age_at_vac,   #Vector indicating age at vaccination, for each vaccinated cohort
                          vac_cohorts_year_of_vac,  #Vector indicating year of vaccination (as simulation year), for each vaccinated cohort  
                          timeframe,                #Simulation timeframe,
                          N_age_groups,             #Nb of age groups in the model,
                          vac_levels,               #Nb of vaccination statuses
                          parms_mod,                #List of model parameters,
                          out_hosts                 #Part of the object out (output by the solver) containing size of host compartments by day, that was simulated for the current simulation year
){
  
  #==================================================================
  #--- Convert solver output into an array --------------------------
  #==================================================================
  out_hosts_7D<-array(out_hosts, dim=c(length(out_hosts[,1]),  #Nb of simulation days (nb rows = length of the 1st column)
                                       N_age_groups,           #Nb age groups
                                       5,5,5,5,                #Nb statuses for each serotype
                                       vac_levels))            #Nb vac statuses
  
  
  #==================================================================
  #--- Prepare flags for negative infeciton counts ------------------
  #==================================================================
  #Some of the infection counts are going to be very small (especially for 3rd/4th infections in young cohorts)
  #Because of the precision loss during the integration, this may result in negative infection counts
  #In this case replacing them by zero and displaying a warning message in the console
  neg_inf1_flag<-FALSE
  neg_inf2_flag<-FALSE
  neg_inf3_flag<-FALSE
  neg_inf4_flag<-FALSE
  
  neg_inf1_min<-0
  neg_inf2_min<-0
  neg_inf3_min<-0
  neg_inf4_min<-0
  
  #==================================================================
  #--- Loop through each vaccinated cohort --------------------------
  #==================================================================
  for(vac_cohort in 1:vac_cohorts_nb){
    
    year_of_vac<-vac_cohorts_year_of_vac[vac_cohort]         #Determine the year when this cohort was (or will be) vaccinated
    years_since_vac<-year-year_of_vac+1                      #Determine which year since vaccination it is (e.g. someone vaccinated at year 3, will be at year 4 since vaccination at year 6; 6-3+1=4)
    
    #================================================
    #--- If cohort is vaccinated, continue ----------
    #================================================
    if(year_of_vac<=year){                                   #Year of vaccination is smaller than or equal to the current year --> cohort already vaccinated --> proceed (otherwise this cohort is ignored for this year and all values for it remain zero)
      
     #Determine current age of the cohort
      age_at_vac_yrs<-vac_cohorts_age_at_vac[vac_cohort]     #Age this cohort was at the moment of vaccination (years) --> taking this value from the previously constructed vector
      age_current_yrs<-age_at_vac_yrs+years_since_vac-1      #Calculating current age (in years) from the age at vaccination and years elapsed since vaccination --> differentiating age at vaccination and current age, as we want to track each vaccinated cohort over time (i.e. need to keep track of how old they are now to correctly position them within the compartments)
      
      age_at_vac_ind<-age_at_vac_yrs+1                       #Index of the age at vaccination (+1 as age starts with zero and indexing starts with 1)
      age_current_ind<-age_current_yrs+1                     #...same for the current age
      
      #Determine the probability of survival after this year (will be required further in the function)
      survival<-1-parms_mod$parms_epi_dyn_model$muHi[age_current_ind]
        
      if(parms_mod$parms_vac_strat$vac_strat_scope==1){       #If all individuals are vaccinated without testing then estimate object to extra check the calculations
          
        #Determine the size of the cohort (for checking the counts in the trace)
        coverage<-parms_mod$parms_vac_strat$coverage_by_age_year[age_at_vac_ind, year_of_vac]                #Coverage at which this cohort was vaccinated
        
        cohort_size_before_ageing<-parms_mod$parms_epi_other$NHi[age_current_ind]*coverage                   #Nb of vaccinated individuals in this cohort now (before ageing is applied)
        cohort_size_after_ageing <-parms_mod$parms_epi_other$NHi[age_current_ind+1]*coverage                 #... after ageing is applied
      }

      
      #================================================
      #-- Get the status at day 1 ---------------------
      #================================================
      #Estimate nb of hosts with each status at the first day of the year. To do so, from the solver output for the current year (object "out_hosts_7D") we take:
      #- First simulated day (initial state of the year 1 corresponds to situation right after the vaccination; DIMENSION 1)
      #- Age group that corresponds to the current vaccinated cohort (DIMENSION 2)
      #- Combination(s) of indices for the four serotypes that correnspond(s) to each number of infections before vaccination (DIMENSIONS 3-6; here -1 indicates "not equal to 1", i.e. excluding the state "Susceptible to a serotype")
      #- Vaccinated individiuals only (DIMENSION 7), as the probability of an infection in vaccinees can be different if the vaccine protects against infection
      
      x0_was<- sum(out_hosts_7D[1, age_current_ind, 1, 1, 1, 1,(2:vac_levels)])#History of no infections (index for each serotype is 1)
      
      x1_was<-(sum(out_hosts_7D[1, age_current_ind,-1, 1, 1, 1,(2:vac_levels)])+ #Hosts with the history of DENV-1 (index for DENV-1 is not 1),
               sum(out_hosts_7D[1, age_current_ind, 1,-1, 1, 1,(2:vac_levels)])+ #Hosts with the history of DENV-2 (index for DENV-2 is not 1),
               sum(out_hosts_7D[1, age_current_ind, 1, 1,-1, 1,(2:vac_levels)])+ #Hosts with the history of DENV-3 (index for DENV-3 is not 1),
               sum(out_hosts_7D[1, age_current_ind, 1, 1, 1,-1,(2:vac_levels)])) #Hosts with the history of DENV-4 (index for DENV-4 is not 1)
      
      x2_was<-(sum(out_hosts_7D[1, age_current_ind,-1,-1, 1, 1,(2:vac_levels)])+ #All possible combinations of serotypes that result in two previous infections; here DENV-1 & DENV-2
               sum(out_hosts_7D[1, age_current_ind,-1, 1,-1, 1,(2:vac_levels)])+ #... DENV-1 & DENV-3
               sum(out_hosts_7D[1, age_current_ind,-1, 1, 1,-1,(2:vac_levels)])+ #... DENV-1 & DENV-4
               sum(out_hosts_7D[1, age_current_ind, 1,-1,-1, 1,(2:vac_levels)])+ #... DENV-2 & DENV-3
               sum(out_hosts_7D[1, age_current_ind, 1,-1, 1,-1,(2:vac_levels)])+ #... DENV-2 & DENV-4
               sum(out_hosts_7D[1, age_current_ind, 1, 1,-1,-1,(2:vac_levels)])) #... DENV-3 & DENV-4
      
      x3_was<-(sum(out_hosts_7D[1, age_current_ind,-1,-1,-1, 1,(2:vac_levels)])+ #All combinations of serotypes that result in three previous infections (same logic as above)
               sum(out_hosts_7D[1, age_current_ind,-1,-1, 1,-1,(2:vac_levels)])+ 
               sum(out_hosts_7D[1, age_current_ind,-1, 1,-1,-1,(2:vac_levels)])+ 
               sum(out_hosts_7D[1, age_current_ind, 1,-1,-1,-1,(2:vac_levels)])) 
      
      x4_was<- sum(out_hosts_7D[1, age_current_ind,-1,-1,-1,-1,(2:vac_levels)]) #History of four infections (index for each serotype is not 1)
      
      #================================================
      #-- If just vaccinated, fill the trace top row --
      #================================================
      #If the cohort was vaccinated at the start of the current simulation year, filling the top row for them using the counts from day 1 (state at the beginning of the first year post-vaccination)
      
      if(year_of_vac==year){
        
        #On the left-hand side of the expressions below. Filling the "cell" for:
        #- The current vaccinated  cohort (DIMENSION 1 of the object "trace")
        #- Each nb of infections before vaccination (0 to 4 infections corresponding to levels 1 to 5 in the DIMENSION 2 of the object "trace")
        #- Start of the current simulation year (which is the first year of vaccination for the cohort in question; DIMENSION 3 of the object "trace")
        #  Note: the values for previous years for this cohort will remain zeros (they are not present in the trace until they are vaccinated)
        #- Current nb of the infections, which, right after the vaccination, is the same as the nb of infections before vaccination (i.e. DIMENSION 4 is equal to dimension 2)
        #- 1st ("tunnel") year since the last boosting (DIMENSION 5; since all these individuals were just vaccinated and vaccination serves as a boosting)
        
        trace[vac_cohort,1,year,1,1]<-x0_was 
        trace[vac_cohort,2,year,2,1]<-x1_was
        trace[vac_cohort,3,year,3,1]<-x2_was
        trace[vac_cohort,4,year,4,1]<-x3_was
        trace[vac_cohort,5,year,5,1]<-x4_was 
      }
      
      #================================================
      #-- Save the trace for this year ----------------
      #================================================
      trace_this_year<-trace[vac_cohort,,year,,]    #This "slice" of the trace is about to be modified, but later it needs to be reset, so saving it  
      
      #================================================
      #-- Day by day, fill the trace ------------------
      #================================================
      #Initialize (or reset) the counts for the cumulative counts of new infections in a given year
      new_inf1_total<-0    
      new_inf2_total<-0
      new_inf3_total<-0
      new_inf4_total<-0
      
      #Loop through each row in the solver output for the current year 
      #(366 rows where the 1st one is the last day of the previous year with ageing & vaccination applied to it)
      for(day in 2:366){  
        
        #---------------------------------------
        #-- Get status for the day ------------- 
        #---------------------------------------
        #Estimating the nb of hosts with each nb of infections at the current day (same logic as above)
        #Note: these variables are labelled "is" because they represent "today". At the end of the loop they will be re-labelled as "was" to be used as "yesterday" in the next loop
        
        x0_is<- sum(out_hosts_7D[day,age_current_ind, 1, 1, 1, 1,(2:vac_levels)])  
        
        x1_is<-(sum(out_hosts_7D[day,age_current_ind,-1, 1, 1, 1,(2:vac_levels)])+ 
                sum(out_hosts_7D[day,age_current_ind, 1,-1, 1, 1,(2:vac_levels)])+ 
                sum(out_hosts_7D[day,age_current_ind, 1, 1,-1, 1,(2:vac_levels)])+ 
                sum(out_hosts_7D[day,age_current_ind, 1, 1, 1,-1,(2:vac_levels)])) 
        
        x2_is<-(sum(out_hosts_7D[day,age_current_ind,-1,-1, 1, 1,(2:vac_levels)])+ 
                sum(out_hosts_7D[day,age_current_ind,-1, 1,-1, 1,(2:vac_levels)])+ 
                sum(out_hosts_7D[day,age_current_ind,-1, 1, 1,-1,(2:vac_levels)])+ 
                sum(out_hosts_7D[day,age_current_ind, 1,-1,-1, 1,(2:vac_levels)])+ 
                sum(out_hosts_7D[day,age_current_ind, 1,-1, 1,-1,(2:vac_levels)])+ 
                sum(out_hosts_7D[day,age_current_ind, 1, 1,-1,-1,(2:vac_levels)])) 
        
        x3_is<-(sum(out_hosts_7D[day,age_current_ind,-1,-1,-1, 1,(2:vac_levels)])+ 
                sum(out_hosts_7D[day,age_current_ind,-1,-1, 1,-1,(2:vac_levels)])+ 
                sum(out_hosts_7D[day,age_current_ind,-1, 1,-1,-1,(2:vac_levels)])+ 
                sum(out_hosts_7D[day,age_current_ind, 1,-1,-1,-1,(2:vac_levels)])) 
        
        x4_is<- sum(out_hosts_7D[day,age_current_ind,-1,-1,-1,-1,(2:vac_levels)])  
        
        #---------------------------------------
        #-- Estimate nb of infection this day --
        #---------------------------------------
        #Initialize or reset the counts
        new_inf1<-0
        new_inf2<-0
        new_inf3<-0
        new_inf4<-0
        
        #The nb of new infections that took place during this day can be estimated by comparing the size of compartments, one after another
        #(in the absence of mortality during the year, all exits from index 1 ("Susceptible") are due to infection acquisition)
        
        #The number of new first infections can be estimated by comparing the nb of people with zero infections today and yesterday
        #x0_is=x0_was-new_inf1, where new_inf1 are the new first infections that took place in this day
        new_inf1<-x0_was-x0_is
        
        #The number of new second infections can be estimated in a similar way
        #x1_is=x1_was-new_inf2+new_inf1 --> the nb of hosts with the history of one infection is determined by how many left (by acquiring a second infection) and how many came in (by acquiring a first infection; estimated above)
        new_inf2<-x1_was-x1_is+new_inf1
        
        #Similarly for the new third infections
        #x2_is=x2_was-new_inf3+new_inf2 --> the nb of hosts with the history of two infections is determined by how many left (by acquiring a third infection) and how many came in (by acquiring a second infection; estimated above)
        new_inf3<-x2_was-x2_is+new_inf2
        
        #Similarly for the new fourth infection
        #x3_is=x3_was-new_inf4+new_inf3
        new_inf4<-x3_was-x3_is+new_inf3
        
        #---------------------------------------
        #-- Manage negative counts -------------
        #---------------------------------------
        if(new_inf1<0){                                          #If the nb of new infections is negative
          neg_inf1_flag<-TRUE                                    #... set the flag to TRUE
          neg_inf1_min<-min(neg_inf1_min, new_inf1)              #... save the new min
          new_inf1<-0                                            #... replace the negative count by zero
        }
        if(new_inf2<0){
          neg_inf2_flag<-TRUE
          neg_inf2_min<-min(neg_inf2_min, new_inf2)
          new_inf2<-0
        }
        if(new_inf3<0){
          neg_inf3_flag<-TRUE
          neg_inf3_min<-min(neg_inf3_min, new_inf3)
          new_inf3<-0
        }
        if(new_inf4<0){
          neg_inf4_flag<-TRUE
          neg_inf4_min<-min(neg_inf4_min, new_inf4)
          new_inf4<-0
        }
        
        #---------------------------------------
        #-- Increase total counts --------------
        #---------------------------------------
        new_inf1_total<-new_inf1_total+new_inf1
        new_inf2_total<-new_inf2_total+new_inf2
        new_inf3_total<-new_inf3_total+new_inf3
        new_inf4_total<-new_inf4_total+new_inf4
        
        #---------------------------------------
        #-- Move from "is" to "was" ------------
        #---------------------------------------
        #Counts from today will become the counts from yesterday at the next loop
        x0_was<-x0_is
        x1_was<-x1_is
        x2_was<-x2_is
        x3_was<-x3_is
        x4_was<-x4_is
        
        #---------------------------------------
        #-- Fill the trace ---------------------
        #---------------------------------------
        
        if(year<timeframe){        #At each year the row of the trace is filled that corresponds to the NEXT year. Therefore, if the last simulation year, not filling the trace (as already in the last row of the object)
          
          #-- Move those with new infections ---------
          #As new infections occur, some individuals will move to the next "stratum" (tunnel year 1). The others will stay where they are.
          #More explanations of the process are provided in the code below
          #Note: As the simulation progresses, more and more "cells" will be filled (e.g. there will be multiple people susceptible to a secondary infection with different time since the primary infection, i.e. in different "tunnel" years)
          #      It is assumed that they all have an equal chance of being infected, however long it's been since their last infection (cross-protection here is, thus, ignored; it is assumed that it's implicitely taken into account in the estimation of the nb of new infections)
          #      Therefore, all the new infections will be distributed proportionally between the eligible "cells" based on how many people are in them
          #Note: For each cohort the values will only be filled once they are vaccinated. For example, if a cohort is vaccinated at year 3, rows for years 1 and 2 (dimension 3) will remain zero
          
          for(inf_nb in 1:4){         #Loop through all new infections (1st to 4th)
            
            #Determine the position in the array of the "stratum" that is at risk for an infection being considered in the current loop
            #Note: this value coincides with the current infection nb; for example if the new 1st infections are being considered in the loop, the stratum at risk are those with the history of 0 infections, i.e. level 1 of the 4th dimension of the object "trace"
            stratum_at_risk<-inf_nb       
            
            #Determine the position in the array of the next stratum
            next_stratum<-stratum_at_risk+1   
            
            #Get the nb of the new infections for the current loop
            if(inf_nb==1){nb_new_inf<-new_inf1} 
            if(inf_nb==2){nb_new_inf<-new_inf2}
            if(inf_nb==3){nb_new_inf<-new_inf3}
            if(inf_nb==4){nb_new_inf<-new_inf4}
            
            #Determine the total number of hosts who are at risk of such an infection (based on the current number of infections, whatever there status before vaccination was (dimension 2) and however long they were at risk (dimension 5))
            #Note: here the estimation of those at risk is based on the counts for the current year and the next year.
            #Those in the row for the current year (year) have been "sitting" there since the last year. 
            #Those in the row for the next year (year+1) already moved this year and are at risk again (accounting for the fact that some small proportion of individuals will have more than one infection in a year)
            pop_at_risk<-sum(trace[vac_cohort,,year:(year+1),stratum_at_risk,])
            
            #Display a warning message if any of the counts are illogical
            #Note: if the vaccine has a large impact on transmission, the size of the compartment will become very small, which may result in illogical counts (e.g. negative nb of infections) 
            #      or nb of infections that is greater than the population at risk. In this case, displaying warning messages and forcing a correction to allow for the model to run correctly.
            if(pop_at_risk<0){
              message("WARNING! Year ", year, "; day ", day, "; cohort ", vac_cohort, "; population at risk for infection # ", inf_nb, " is negative")
            }else if(pop_at_risk<nb_new_inf){
              message("WARNING! Year ", year, "; day ", day, "; cohort ", vac_cohort, "; population at risk for infection # ", inf_nb, " is smaller than the nb of infections")
              nb_new_inf<-pop_at_risk
            }
            
            #Loop through the dimensions that are not fixed (if pop_at_risk >0)
            if(pop_at_risk>0){
              weight_total<-0                            
              for(nb_inf_before_vac in 1:4){               #Nb infections before vaccination (0 to 3; ignoring those with 4 infections before vaccination because they are no longer at risk for a new one); this parameter does not change, it's fixed since the vaccination
                for(tunnel_year in 1:years_since_vac){     #"Tunnel" years (only looping in the years since vaccination as other cells will be empty; for example, someone who is currently in his third year since vaccination, cannot be in tunnel year 4)                                 
                  pop_cell_this_year<-trace[vac_cohort,nb_inf_before_vac, year   ,stratum_at_risk,tunnel_year]     #Nb of individuals in the "cell", who haven't moved this year yet (so they still sit in the row for this year)
                  pop_cell_next_year<-trace[vac_cohort,nb_inf_before_vac,(year+1),stratum_at_risk,tunnel_year]     #Nb of individuals in the same cell, who already moved this year (so they went to the row for the next year). Note: For those looping through each tunnel year as well, even if they can only be in tunnel year 1 (since they just experienced an infection this year)
                  
                  weight_cell_this_year<-pop_cell_this_year/pop_at_risk                                            #Weight of the "cell" in the total population at risk (to determine the proportion of the infections that will happen in this cell)
                  weight_cell_next_year<-pop_cell_next_year/pop_at_risk
                  
                  nb_inf_cell_this_year<-weight_cell_this_year*nb_new_inf                                          #Nb of infections that will happen from the current "cell" (i.e. among those who have not experienced a new infection yet)
                  nb_inf_cell_next_year<-weight_cell_next_year*nb_new_inf                                          #Nb of infections that will happen from the same "cell", but among those who already experience a new infeciton this year (and are at risk again)
                  
                  
                  #Move those with a new infection to the 1st tunnel year of the next stratum and the next simulation year (adding to those who are already there)
                  #Note: these new infections come from both who have and who have not experienced an infection this year already
                  trace[vac_cohort, nb_inf_before_vac, (year+1), next_stratum, 1]<-trace[vac_cohort, nb_inf_before_vac, (year+1), next_stratum, 1]+nb_inf_cell_this_year+nb_inf_cell_next_year
                  
                  #Remove them from the cells where they were before
                  trace[vac_cohort,nb_inf_before_vac, year   ,stratum_at_risk,tunnel_year]<-trace[vac_cohort,nb_inf_before_vac, year   ,stratum_at_risk,tunnel_year]-nb_inf_cell_this_year
                  trace[vac_cohort,nb_inf_before_vac,(year+1),stratum_at_risk,tunnel_year]<-trace[vac_cohort,nb_inf_before_vac,(year+1),stratum_at_risk,tunnel_year]-nb_inf_cell_next_year
                  
                  #Increase the count for total weight (for checks)
                  weight_total<-weight_total+weight_cell_this_year+weight_cell_next_year
                }  
              }
              
              #Check the total weight (i.e. that all the "cells" at risk were accounted for)
              if(round(weight_total,10)!=1){
                message("Total weight of pop at risk is not equal to 1")  
              }
            }
          }
          
          if(parms_mod$parms_vac_strat$vac_strat_scope==1){       #If all individuals are vaccinated without testing then estimate object to extra check the calculations
            #-- Check trace pop counts before ageing ---
            cohort_pop_trace<-sum(trace[vac_cohort,,year:(year+1),,])   
            if(round(cohort_pop_trace,0)!=round(cohort_size_before_ageing,0)){
              message("WARNING! Year ", year, "; day ", day, "; cohort ", vac_cohort, "; total cohort population is ", cohort_pop_trace, ", while it should be before", cohort_size_before_ageing)
            }
          }
          
          #When the size of the compartments is very small, sometimes values in the trace become very small but negative numbers
          #In this case forcing a correction by resetting those values to zero
          if(min(trace[vac_cohort,,year:(year+1),,])<0){
            message("WARNING! Year ", year, "; day ", day, "; cohort ", vac_cohort, "; there are NEGATIVE VALUES IN THE TRACE. Smallest values is equal to ", min(trace[vac_cohort,,year:(year+1),,])<0)
            trace[which(trace<0)]<-0
          }
          
          
          #-------------------------------------------
          #-- If last day ----------------------------
          #-------------------------------------------
          if(day==366){
            
            #-- Age everyone -----------
            #Loop through each cell in the row "year" and carry over everyone left there to the same strata (as they did not get infected this year), but next simulation year (year+1) and next tunnel year (tunnel_year+1)
            #Adding to those who are already there (although likely empty, because whoever had an infection will be in tunnel year 1 only)
            for(nb_inf_before_vac in 1:5){                      
              for(nb_inf_current in 1:5){
                for(tunnel_year in 1:years_since_vac){
                  trace[vac_cohort, nb_inf_before_vac, (year+1), nb_inf_current, (tunnel_year+1)]<-sum(trace[vac_cohort, nb_inf_before_vac, (year+1), nb_inf_current, (tunnel_year+1)]+   
                                                                                                         trace[vac_cohort, nb_inf_before_vac,  year   , nb_inf_current,  tunnel_year])
                }
              }
            }
            
            #-- Apply mortality --------
            trace[vac_cohort,,(year+1),,]<-trace[vac_cohort,,(year+1),,]*survival        
            
            #-- Check the counts in the trace ----------
            if(parms_mod$parms_vac_strat$vac_strat_scope==1){       #If all individuals are vaccinated without testing then estimate object to extra check the calculations
              cohort_pop_trace<-sum(trace[vac_cohort,,(year+1),,])   
              if(round(cohort_pop_trace,0)!=round(cohort_size_after_ageing,0)){
                message("WARNING! Year ", year, "; day ", day, "; cohort ", vac_cohort, "; total cohort population is ", cohort_pop_trace, ", while it should be after", cohort_size_after_ageing)
              }
            }
            
            if(min(trace[vac_cohort,,(year+1),,])<0){
              message("WARNING! Year ", year, "; day ", day, "; cohort ", vac_cohort, "; there are NEGATIVE VALUES IN THE TRACE")
            }
            
            #-- Restore the values for current year ----
            trace[vac_cohort,,year,,]<-trace_this_year
          }
        }
      }
      
      #---------------------------------------------------------------------------
      #--- Record the nb of new inf at this year ---------------------------------
      #---------------------------------------------------------------------------
      new_inf[vac_cohort, year, 1]<-new_inf1_total
      new_inf[vac_cohort, year, 2]<-new_inf2_total
      new_inf[vac_cohort, year, 3]<-new_inf3_total
      new_inf[vac_cohort, year, 4]<-new_inf4_total
    }
  }
  
  #-- Display a message about negative counts ----
  if(neg_inf1_flag){message("WARNING! Year ", year, ". Some of the 1st infection counts were negative and replaced by zero; lowest count equal to ", neg_inf1_min)}
  if(neg_inf2_flag){message("WARNING! Year ", year, ". Some of the 2nd infection counts were negative and replaced by zero; lowest count equal to ", neg_inf2_min)}
  if(neg_inf3_flag){message("WARNING! Year ", year, ". Some of the 3rd infection counts were negative and replaced by zero; lowest count equal to ", neg_inf3_min)}
  if(neg_inf4_flag){message("WARNING! Year ", year, ". Some of the 4th infection counts were negative and replaced by zero; lowest count equal to ", neg_inf4_min)}
  
  return(list(trace=trace,
              new_inf=new_inf))
}

############################################################
#This function displays and updates a progress window while the model is running

ModelRun_UpdateProgress<-function(queue,
                                  progress_matrix,
                                  current_sim,
                                  current_status){
  
  nb_sim<-length(queue[,"scen_ref"])                
  
  #--- Generate/update the progress matrix -----
  if(current_sim==0){                                   #If no simulations were run yet --> generate the matrix
    progress_matrix<-matrix(NA, nrow=nb_sim, ncol=3)
    for(sim in 1:nb_sim){
      progress_matrix[sim,1]<-paste("Scenario ", queue[sim, "scen_ref"], sep="")
      progress_matrix[sim,2]<-paste("(start year ", queue[sim, "year_start"], ")", sep="")
      progress_matrix[sim,3]<-": "
    }
  }else{                                                #Otherwise - update status for the current simulation
    progress_matrix<-progress_matrix
    progress_matrix[current_sim,3]<-paste(":", current_status, sep="")
  }
  
  #--- Convert in a tagList --------------------
  progress_list<-tagList()
  for(sim in 1:nb_sim){
    progress_list<-c(progress_list, tagList(paste(progress_matrix[sim, 1], progress_matrix[sim, 2], progress_matrix[sim, 3]), br()))
  }
  
  #--- Estimate runtime ------------------------
  timeframe<-as.numeric(queue[, "timeframe"])         #Timeframe for each simulation, years
  runtime_mins<-runtime_one_year*sum(timeframe)       #RUntime, minutes (based on the nb of simulation years to run and estimated runtime per one simulation year)
  runtime_hours<-runtime_mins/60                      #Runtime, hours
  if(runtime_mins<60){                                 
    runtime<-paste(format(round(runtime_mins, 0)), " min", sep="")
  }else{
    runtime<-paste(format(round(runtime_hours, 1)), " hours", sep="")
  }
  
  
  showModal(
    modalDialog(
      tags$b(paste("Total runtime (approx.): ", runtime, sep="" )), br(), br(),
      progress_list,
      title="Run progress",
      easyClose=FALSE,
      footer=NULL,
      fade=FALSE
    )
  )
  
  return(progress_matrix)
}

GetStackInfo<-function(currentStep=""){
	info<-Cstack_info()
	info["usage"]<-100*info["current"]/info["size"]
	if(currentStep!=""){message("Running: ", currentStep)}
	message("C Stack: ", info["current"], " used on ", info["size"], ". Usage: ", info["usage"],"%")
	rm(info); gc()
}

