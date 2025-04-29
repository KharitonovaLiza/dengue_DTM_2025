##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contains functions that import, select, initalise and load all the population & transmission parameters
###########################################################################################################################################################################################################

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function selects population & transmission inputs that are required to run the dynamic model or calibration

ParmsEpi_Select<-function(parms_epi_imported,     #Inputs imported from an Excel file Epi_XXXXX 
                          calibration){           #Calibration flag (if TRUE, the parameters required for the calibration process are output as well). 

  with(parms_epi_imported, {
  
    #----------------------------------------------------
    #--- LIST OF POPUALTION & TRANSMISSION PARAMETERS ---
    #----------------------------------------------------
    
    parms_epi=list(
    
      #Host population
      N_age_groups=as.integer(N_age_groups_model),                      #Number of age groups used in the model --> integer
      age_group_lower_bound=as.integer(age_group_lower_bound_model),    #Lower bounds of the age groups used in the model --> vector of integers, nb. elements = N_age_groups
      age_group_upper_bound=as.integer(age_group_upper_bound_model),    #Upper bounds of the age groups used in the model --> vector of integers, nb. elements = N_age_groups
      pop_size=as.numeric(pop_size_model),                 		          #Population size used in the model for each age group (adjusted to maintain constant population size) --> numeric vector, nb. elements = N_age_groups
      mortality=as.numeric(mortality_model),			                      #Mortality used in the model in each age group --> numeric vector, nb. elements = N_age_groups	
	    life_expectancy=as.numeric(life_expectancy),			                #Life expectancy --> numeric vector, nb. elements = N_age_groups
        
      #Vector population
      ratio_VH=as.numeric(ratio_VH),                                    #Number of adult female vectors per host (average without seasonal fluctuations) --> numeric
      LE_V=as.numeric(LE_V),                                            #Vector life expectancy (days) --> numeric
      
      #Seasonality
      seasonality=as.logical(seasonality),	                            #Inclusion of seasonality (TRUE/FALSE) --> boolean
      season_p1=as.numeric(season_p1),                                  #Amplitude of the sine function applied to the vector population size at each point in time --> numeric
      season_p2=as.numeric(season_p2),                                  #Horizontal shift of the sine function applied to the vector population size at each point in time --> numeric
      
      #Transmission in vectors
      betaHV=as.numeric(betaHV),                                        #Probability of virus transmission from an infectious host to a susceptible vector (given a bite)  --> numeric 
      b=as.numeric(b),                                                  #Number of bites an adult female vector makes per day --> numeric
      dur_latency_V=as.numeric(dur_latency_V),	                        #Duration of latency in vectors (days) --> numeric
      
      #Transmission in hosts
      betaVH=as.numeric(betaVH),                                                                 #Base of virus transmission from an infectious vector to a susceptible host (given a bite) --> numeric vector, nb. elements = 4 (one for each serotype)
      
      N_age_groups_betaVH_age_coef=as.integer(N_age_groups_betaVH_age_coef),                     #Number of age groups used to differentiate age-specific coefficient applied to betaVH (relative risk of the virus transmission from an infectious vector to a susceptible host, given a bite) --> integer
      betaVH_age_coeff_age_group_lower_bound=as.integer(betaVH_age_coeff_age_group_lower_bound), #Lower bounds of the age groups used to differentiate age-specific coefficient --> vector of integers, nb. elements = N_age_groups_betaVH_age_coef
      betaVH_age_coeff_age_group_upper_bound=as.integer(betaVH_age_coeff_age_group_upper_bound), #Upper bounds of the age groups used to differentiate age-specific coefficient --> vector of integers, nb. elements = N_age_groups_betaVH_age_coef
      betaVH_age_coef=as.numeric(betaVH_age_coef),                                               #Age-specific coefficient applied to betaVH --> numeric vector, nb. elements = N_age_groups_betaVH_age_coef
      
      level_infect_sympt_unvac=as.numeric(level_infect_sympt_unvac),    #Relative infectiousness of symptomatic infections in unvaccinated hosts (reference, always equal to 1) --> numeric
      level_infect_asympt_unvac=as.numeric(level_infect_asympt_unvac),	#Relative infectiousness of asymptomatic infections in unvaccinated hosts (relative to symptomatic) --> numeric
      dur_latency_H=as.numeric(dur_latency_H),                          #Duration of latency in hosts (days) --> numeric
      dur_vir=as.numeric(dur_vir_H),                                    #Duration of viremia in hosts (days) --> numeric     
      
      #Susceptibility enhancement characteristics
      dzetaCE=as.integer(dzetaCE),					                            #Inclusion of susceptibility enhancement and its frequency(0 = not considered, 1 = applied after the 1st infection only, 2 = applied after each infection) --> Integer 
      gammaCE=as.numeric(gammaCE),                                      #Susceptibility enhancement factor --> numeric
      
      #Cross-protection characteristics
      dzetaCP=as.logical(dzetaCP)[1],                                   #Inclusion of cross-protection (TRUE/FALSE) --> boolean (Note: for some reason R always sees it as a vector of two elements, so to avoid warning messages, only the first element is taken)
      gammaCP=as.numeric(gammaCP),                                      #Relative risk of infection during cross-protection period --> numeric
      dur_CP=as.numeric(dur_CP_days),                                   #Duration of cross-protection (days) --> numeric
      nb_infect_max=as.integer(nb_infect_max),	                        #Maximum allowed number of infections (up to 4) --> integer
      
      
      #Proportion of symptomatic infections (out of all infections)
      N_age_groups_prop_sympt_unvac=as.integer(N_age_groups_prop_sympt_unvac),                       #Number of age groups used to differentiate the probability of being symptomatic --> integer
      prop_sympt_unvac_age_group_lower_bound=as.integer(prop_sympt_unvac_age_group_lower_bound),     #... their lower bounds --> vector of integers, nb. elements = N_age_groups_prop_sympt_unvac
      prop_sympt_unvac_age_group_upper_bound=as.integer(prop_sympt_unvac_age_group_upper_bound),     #... their upper bounds --> vector of integers, nb. elements = N_age_groups_prop_sympt_unvac
      prop_sympt_1inf_DENV1_unvac =as.numeric(prop_sympt_1inf_DENV1_unvac),                          #Prop. of symptomatic infections out of all PRIMARY infections in unvaccinated hosts caused by DENV-1 --> numeric vector, nb. elements = N_age_groups_prop_sympt_unvac
      prop_sympt_2inf_DENV1_unvac =as.numeric(prop_sympt_2inf_DENV1_unvac),                          #Prop. of symptomatic infections out of all SECONDARY infections in unvaccinated hosts caused by DENV-1 --> numeric vector, nb. elements = N_age_groups_prop_sympt_unvac
      prop_sympt_34inf_DENV1_unvac=as.numeric(prop_sympt_34inf_DENV1_unvac),                         #Prop. of symptomatic infections out of all POST-SECONDARY infections in unvaccinated hosts caused by DENV-1 --> numeric vector, nb. elements = N_age_groups_prop_sympt_unvac
      prop_sympt_1inf_DENV2_unvac=as.numeric(prop_sympt_1inf_DENV2_unvac),                           #... same for DENV-2
      prop_sympt_2inf_DENV2_unvac=as.numeric(prop_sympt_2inf_DENV2_unvac),
      prop_sympt_34inf_DENV2_unvac=as.numeric(prop_sympt_34inf_DENV2_unvac),
      prop_sympt_1inf_DENV3_unvac=as.numeric(prop_sympt_1inf_DENV3_unvac),                           #... same for DENV-3
      prop_sympt_2inf_DENV3_unvac=as.numeric(prop_sympt_2inf_DENV3_unvac),
      prop_sympt_34inf_DENV3_unvac=as.numeric(prop_sympt_34inf_DENV3_unvac),
      prop_sympt_1inf_DENV4_unvac=as.numeric(prop_sympt_1inf_DENV4_unvac),                           #... same for DENV-4       
      prop_sympt_2inf_DENV4_unvac=as.numeric(prop_sympt_2inf_DENV4_unvac),
      prop_sympt_34inf_DENV4_unvac=as.numeric(prop_sympt_34inf_DENV4_unvac),
      
      #Proportion of severe infections (out of symptomatic infections)
      N_age_groups_prop_severe_sympt_unvac=as.integer(N_age_groups_prop_severe_sympt_unvac),
      prop_severe_sympt_unvac_age_group_lower_bound=as.integer(prop_severe_sympt_unvac_age_group_lower_bound),
      prop_severe_sympt_unvac_age_group_upper_bound=as.integer(prop_severe_sympt_unvac_age_group_upper_bound),
      prop_severe_sympt_1inf_DENV1_unvac =as.numeric(prop_severe_sympt_1inf_DENV1_unvac),
      prop_severe_sympt_2inf_DENV1_unvac =as.numeric(prop_severe_sympt_2inf_DENV1_unvac),
      prop_severe_sympt_34inf_DENV1_unvac=as.numeric(prop_severe_sympt_34inf_DENV1_unvac),
      prop_severe_sympt_1inf_DENV2_unvac=as.numeric(prop_severe_sympt_1inf_DENV2_unvac),
      prop_severe_sympt_2inf_DENV2_unvac=as.numeric(prop_severe_sympt_2inf_DENV2_unvac),
      prop_severe_sympt_34inf_DENV2_unvac=as.numeric(prop_severe_sympt_34inf_DENV2_unvac),
      prop_severe_sympt_1inf_DENV3_unvac=as.numeric(prop_severe_sympt_1inf_DENV3_unvac),
      prop_severe_sympt_2inf_DENV3_unvac=as.numeric(prop_severe_sympt_2inf_DENV3_unvac),
      prop_severe_sympt_34inf_DENV3_unvac=as.numeric(prop_severe_sympt_34inf_DENV3_unvac),
      prop_severe_sympt_1inf_DENV4_unvac=as.numeric(prop_severe_sympt_1inf_DENV4_unvac),
      prop_severe_sympt_2inf_DENV4_unvac=as.numeric(prop_severe_sympt_2inf_DENV4_unvac),
      prop_severe_sympt_34inf_DENV4_unvac=as.numeric(prop_severe_sympt_34inf_DENV4_unvac),
      
      #Proportion of hospitalized infections (out of symptomatic infections)
      N_age_groups_prop_hospit_sympt_unvac=as.integer(N_age_groups_prop_hospit_sympt_unvac),
      prop_hospit_sympt_unvac_age_group_lower_bound=as.integer(prop_hospit_sympt_unvac_age_group_lower_bound),
      prop_hospit_sympt_unvac_age_group_upper_bound=as.integer(prop_hospit_sympt_unvac_age_group_upper_bound),
      prop_hospit_sympt_1inf_DENV1_unvac =as.numeric(prop_hospit_sympt_1inf_DENV1_unvac),
      prop_hospit_sympt_2inf_DENV1_unvac =as.numeric(prop_hospit_sympt_2inf_DENV1_unvac),
      prop_hospit_sympt_34inf_DENV1_unvac=as.numeric(prop_hospit_sympt_34inf_DENV1_unvac),
      prop_hospit_sympt_1inf_DENV2_unvac=as.numeric(prop_hospit_sympt_1inf_DENV2_unvac),
      prop_hospit_sympt_2inf_DENV2_unvac=as.numeric(prop_hospit_sympt_2inf_DENV2_unvac),
      prop_hospit_sympt_34inf_DENV2_unvac=as.numeric(prop_hospit_sympt_34inf_DENV2_unvac),
      prop_hospit_sympt_1inf_DENV3_unvac=as.numeric(prop_hospit_sympt_1inf_DENV3_unvac),
      prop_hospit_sympt_2inf_DENV3_unvac=as.numeric(prop_hospit_sympt_2inf_DENV3_unvac),
      prop_hospit_sympt_34inf_DENV3_unvac=as.numeric(prop_hospit_sympt_34inf_DENV3_unvac),
      prop_hospit_sympt_1inf_DENV4_unvac=as.numeric(prop_hospit_sympt_1inf_DENV4_unvac),
      prop_hospit_sympt_2inf_DENV4_unvac=as.numeric(prop_hospit_sympt_2inf_DENV4_unvac),
      prop_hospit_sympt_34inf_DENV4_unvac=as.numeric(prop_hospit_sympt_34inf_DENV4_unvac),
      
      #Proportion of fatal infections (out of hospitalized severe infections)
      N_age_groups_prop_death_hospit_sympt_unvac=as.integer(N_age_groups_prop_death_hospit_sympt_unvac),
      prop_death_hospit_sympt_unvac_age_group_lower_bound=as.integer(prop_death_hospit_sympt_unvac_age_group_lower_bound),
      prop_death_hospit_sympt_unvac_age_group_upper_bound=as.integer(prop_death_hospit_sympt_unvac_age_group_upper_bound),
      prop_death_hospit_sympt_1inf_DENV1_unvac =as.numeric(prop_death_hospit_sympt_1inf_DENV1_unvac),
      prop_death_hospit_sympt_2inf_DENV1_unvac =as.numeric(prop_death_hospit_sympt_2inf_DENV1_unvac),
      prop_death_hospit_sympt_34inf_DENV1_unvac=as.numeric(prop_death_hospit_sympt_34inf_DENV1_unvac),
      prop_death_hospit_sympt_1inf_DENV2_unvac=as.numeric(prop_death_hospit_sympt_1inf_DENV2_unvac),
      prop_death_hospit_sympt_2inf_DENV2_unvac=as.numeric(prop_death_hospit_sympt_2inf_DENV2_unvac),
      prop_death_hospit_sympt_34inf_DENV2_unvac=as.numeric(prop_death_hospit_sympt_34inf_DENV2_unvac),
      prop_death_hospit_sympt_1inf_DENV3_unvac=as.numeric(prop_death_hospit_sympt_1inf_DENV3_unvac),
      prop_death_hospit_sympt_2inf_DENV3_unvac=as.numeric(prop_death_hospit_sympt_2inf_DENV3_unvac),
      prop_death_hospit_sympt_34inf_DENV3_unvac=as.numeric(prop_death_hospit_sympt_34inf_DENV3_unvac),
      prop_death_hospit_sympt_1inf_DENV4_unvac=as.numeric(prop_death_hospit_sympt_1inf_DENV4_unvac),
      prop_death_hospit_sympt_2inf_DENV4_unvac=as.numeric(prop_death_hospit_sympt_2inf_DENV4_unvac),
      prop_death_hospit_sympt_34inf_DENV4_unvac=as.numeric(prop_death_hospit_sympt_34inf_DENV4_unvac)

    )
    
    #----------------------------------------------------
    #--- LIST OF PARAMETERS FOR CALIBRATION -------------
    #----------------------------------------------------
    if(calibration){
      parms_calibr=list(
        
        #Host population  
        
        N_age_groups=as.integer(N_age_groups_calibr),                      #Number of age groups used in the calibration --> integer
        age_group_lower_bound=as.integer(age_group_lower_bound_calibr),    #Lower bounds of the age groups used in the calibration --> vector of integers, nb. elements = N_age_groups (in calibration)
        age_group_upper_bound=as.integer(age_group_upper_bound_calibr),    #Upper bounds of the age groups used in the calibration --> vector of integers, nb. elements = N_age_groups (in calibration)
        pop_size=as.numeric(pop_size_calibr),                            	 #Population size used in the calibration for each age group (adjusted to maintain constant population size) --> numeric vector, nb. elements = N_age_groups (in calibration)
        mortality=as.numeric(mortality_calibr), 									         #Mortality used in the calibration for the last year of each age group (calculated to maintant constant population size) --> numeric vector, nb. elements = N_age_groups (in calibration)

        #Incidence rate of symptomatic dengue, by age (observed)
        N_age_groups_incidence=as.integer(N_age_groups_incidence),                      #Number of age groups used to define data on observed incidence --> integer
        incidence_age_group_lower_bound=as.integer(incidence_age_group_lower_bound),    #Lower bounds of the age groups used to define data on observed incidence --> vector of integers, nb elements = N_age_groups_incidence
        incidence_age_group_upper_bound=as.integer(incidence_age_group_upper_bound),    #Upper bounds of the age groups used to define data on observed incidence --> vector of integers, nb elements = N_age_groups_incidence,
        incidence=as.numeric(incidence),                                                #Incidence rate of symptomatic dengue, by age (reported, not corrected for under-reporting)
        incidence_per=as.integer(incidence_per),                                        #Denominator for the incidence rate (e.g. per 100,000 people)
        
        #Incidence rate of symptomatic dengue, by month (observed),             #Incidence rate of symptomatic dengue observed, by calendar month --> numeric vector with 12 elements
        incidence_by_month=as.numeric(incidence_by_month),
        
        #Expansion factors
        expansion_factor_hospit=as.numeric(expansion_factor_hospit),            # Expansion (under-reporting) factor for hospitalized cases (real hosp cases=reported hosp cases * expansion_factor_hosp)
        expansion_factor_non_hospit=as.numeric(expansion_factor_non_hospit),    # Expansion (under-reporting) factor for non-hospitalized cases (real non-hosp cases=reported non-hosp cases * expansion_factor_non_hosp)
        
        #Seroprevalence data (used to calibrate betaVH)
        N_age_groups_seroprevalence=as.integer(N_age_groups_seroprevalence),                     #Number of age groups used to define seroprevalence data --> integer
        seroprevalence_age_group_lower_bound=as.integer(seroprevalence_age_group_lower_bound),   #Lower bounds of the age groups used to define seroprevalence data --> vector of integers, nb. elements = N_age_groups_seroprevalence
        seroprevalence_age_group_upper_bound=as.integer(seroprevalence_age_group_upper_bound),   #Upper bounds of the age groups used to define seroprevalence data --> vector of integers, nb. elements = N_age_groups_seroprevalence
        nx=as.integer(nx),                                                                       #Number of people susceptible to all four serotypes --> vector of integers, nb. elements = N_age_groups_seroprevalence
        nmz=as.integer(nmz))                                                                     #Number of people with history of DENV infection --> vector of integers, nb. elements = N_age_groups_seroprevalence
      
    }else{
      parms_calibr<-list()                         #If calibration=FALSE (i.e. this function is used to run the main model), the list parms_calibr remains empty
    }

    return(list(parms_epi=parms_epi,               #A list with population & transmission inputs
                parms_calibr=parms_calibr))        #A list with the parameters required for calibration
  })
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function formats the imported and selected population & transmission inputs, as well as calculates the additional model parameters based on the imported inputs
#NB! The order in which the parameters are listed in the list parms_epi_dyn_model is crucial for the infection process coded in C. 
#    If this order is changed, the C code needs to be adapted and recompiled.

ParmsEpi_Initialise<- function(parms_epi_selected,        #Selected and pre-formatted population & transmission inputs and calibration inputs
                               calibration){              #Calibration flag (if TRUE, the parameters required for the calibration process are output as well)
  
  
  #----------------------------------------------------
  #--- PARAMETERS FOR THE DYNAMIC MODEL ---------------
  #----------------------------------------------------
  
  #---------- Hosts -------------------
  if(calibration){											  																   #Whether the code is run for the main model or for the calibration determines which population inputs are loaded in the list parms_epi_dyn_model 
    N_age_groups<-parms_epi_selected$parms_calibr$N_age_groups
    age_group_lower_bound<-parms_epi_selected$parms_calibr$age_group_lower_bound
    age_group_upper_bound<-parms_epi_selected$parms_calibr$age_group_upper_bound
    NHi<-parms_epi_selected$parms_calibr$pop_size		                         #Population size in each age group used in the calibration (adjusted to maintain constant population size) --> numeric vector, nb. elements = N_age_groups (in calibration)
    NH <-sum(parms_epi_selected$parms_calibr$pop_size)                       #Population size, total --> numeric
    ai<-(1/(parms_epi_selected$parms_calibr$age_group_upper_bound-parms_epi_selected$parms_calibr$age_group_lower_bound+1))   #Proportion of the age group that ages every year (weight of the last year in the age group duration) --> numeric vector, nb. elements = N_age_groups (in calibration)
    muHi<-parms_epi_selected$parms_calibr$mortality               				   #Mortality for the last year of each age group used in the calibration (calculated to maintain constant population size) --> numeric vector, nb. elements = N_age_groups (in calibration)
  }else{
    N_age_groups<-parms_epi_selected$parms_epi$N_age_groups
    age_group_lower_bound<-parms_epi_selected$parms_epi$age_group_lower_bound
    age_group_upper_bound<-parms_epi_selected$parms_epi$age_group_upper_bound
    NHi<-parms_epi_selected$parms_epi$pop_size		          								 #Population size in each age group used in the model (adjusted to maintain constant population size) --> numeric vector, nb. elements = N_age_groups (in model)
    NH <-sum(parms_epi_selected$parms_epi$pop_size)                          #Population size, total --> numeric
    ai<-(1/(parms_epi_selected$parms_epi$age_group_upper_bound-parms_epi_selected$parms_epi$age_group_lower_bound+1))         #Proportion of the age group that ages every year (weight of the last year in the age group duration) --> numeric vector, nb. elements = N_age_groups (in model)
    muHi<-parms_epi_selected$parms_epi$mortality                             #Mortality for the last year of each age group used in the model (calculated to maintain constant population size) --> numeric vector, nb. elements = N_age_groups (in model)
  }
  
  xiH<-1/parms_epi_selected$parms_epi$dur_latency_H          #Daily transition rate from "Exposed" to "Infectious" (inverse of latency duration) --> numeric
  ro<-1/parms_epi_selected$parms_epi$dur_vir                 #Daily transition rate from "Infectious" to "Cross-protected" (or "Immune" if no cross-protection; inverse of viremia duration) --> numeric
  
  if(parms_epi_selected$parms_epi$dzetaCP){                  #Daily transition rate from "Cross-protected" to "Immune"... 
    phiCP<-1/parms_epi_selected$parms_epi$dur_CP             #... inverse of CP duration in days (if CP is included) --> numeric     
  }else{
    phiCP<-0                                                 #... zero if CP is not included
  }
  
  #Age-specific coefficient applied to betaVH
  #NOTE: This code only works if the model (or calibration) is run with 101 single-year age groups
  betaVH_age_coef<-rep(NA,N_age_groups)                                                           #Vector with one element per age group
  for(age_group in 1:parms_epi_selected$parms_epi$N_age_groups_betaVH_age_coef){                  #Populating each element of the vector with the value for the corresponding composite age group
    age_low<-parms_epi_selected$parms_epi$betaVH_age_coeff_age_group_lower_bound[age_group]+1     #Start & end age for the composite age group (+1 as age starts with zero)
    age_high<-parms_epi_selected$parms_epi$betaVH_age_coeff_age_group_upper_bound[age_group]+1
    betaVH_age_coef[age_low:age_high]<-parms_epi_selected$parms_epi$betaVH_age_coef[age_group]    #Plugging the value into the final vector
  }

  #Infection severity by single year of age
  sev_vector<-c('sympt','severe_sympt','hospit_sympt','death_hospit_sympt')                #Vector with the labels for severity
  inf_vector<-c('1','2','34')                                                              #Vector with the labels for infection type
  
  for(sev in 1:4){			             #Loop through severity labels/levels
    
    sev_id<-sev_vector[sev]          #Set the current severity label/level
    
    do.call('<-',list(paste0('prop_',sev_id,'_unvac'),array(NA, dim=c(parms_epi_selected$parms_epi$N_age_groups,4,3))))  #Create a  3D array with dimensions 101x3x4 and named by severity label (e.g. "prop_sympt_unvac"), which contains proportions of symptomatic infections by age, serotype and infection type
    
    N_age_groups_sev<-eval(parse(text=paste0('parms_epi_selected$parms_epi$N_age_groups_prop_',sev_id, '_unvac')))                       #Load the nb of age group used to define parameters for the selected severity level
    age_group_lower_bound_temp<-eval(parse(text=paste0('parms_epi_selected$parms_epi$prop_',sev_id, '_unvac','_age_group_lower_bound'))) #Load their lower age bounds
    age_group_upper_bound_temp<-eval(parse(text=paste0('parms_epi_selected$parms_epi$prop_',sev_id, '_unvac','_age_group_upper_bound'))) #Load their upper age bounds
    
    prop_to_get<-eval(parse(text=paste0('prop_',sev_id,'_unvac')))	                                                        #Create a temporary object
    for(infection in 1:3){                                                                                                  #Loop through all infection types: primary, secondary, post-secondary
      inf_id<-inf_vector[infection]                                                                                         #Set the current infection type/label
      for(serotype in 1:4){                                                                                                 #Loop through four serotypes
        prop_to_set<-eval(parse(text=paste0('parms_epi_selected$parms_epi$prop_',sev_id,'_',inf_id,'inf_DENV',serotype, '_unvac')))
        for(age_group in 1:N_age_groups_sev){
          age_low <-age_group_lower_bound_temp[age_group]+1     #Start & end age for the composite age group (+1 as age starts with zero)
          age_high<-age_group_upper_bound_temp[age_group]+1 
          prop_to_get[age_low:age_high,serotype,infection]<- prop_to_set[age_group]
        }
      }
    }
    do.call('<-',list(paste0('prop_',sev_id,'_unvac'),prop_to_get))
  }
  rm(sev_vector,inf_vector,sev_id,N_age_groups_sev,age_group_lower_bound_temp,age_group_upper_bound_temp,prop_to_get); gc()
  
  #Number of past infections for each combination of k,l and m --> 3D array (5x5x5) with integer values
  past_infections<-Parms_GetPastInfections()                                             
  
  #Theta, modifying factors for the force of infection --> 3D array (5x5x5) with numeric values
  theta_unvac<-Parms_GetThetaUnvac(gammaCP=parms_epi_selected$parms_epi$gammaCP,       
                                    dzetaCE=parms_epi_selected$parms_epi$dzetaCE, 
                                    gammaCE=parms_epi_selected$parms_epi$gammaCE, 
                                    nb_infect_max=parms_epi_selected$parms_epi$nb_infect_max)             
  
  #---- Get the object psy for unvaccinated individuals  ------------
  #The object psy captures relative infectiousness of vaccinated individuals, based on the infection type & presence of symptoms (which may depend on age)
  #More details can be found in the script "13_parms_psy.R"
  
  #Initialize the arrays psy for unvaccinated individuals (separate object for each serotype)
  #4D array (101x5x5x5) with numeric values; dimensions are age, status for three serotypes (other than the one being considered in a current infection process)
  psy_DENV1_unvac<-array(NA, dim=c(parms_epi_selected$parms_epi$N_age_groups,5,5,5))      
  psy_DENV2_unvac<-array(NA, dim=c(parms_epi_selected$parms_epi$N_age_groups,5,5,5))
  psy_DENV3_unvac<-array(NA, dim=c(parms_epi_selected$parms_epi$N_age_groups,5,5,5))
  psy_DENV4_unvac<-array(NA, dim=c(parms_epi_selected$parms_epi$N_age_groups,5,5,5))
  
  #Loop through the age dimensiom
  #Get the proportion of symptomatic infections for this age from the object "prop_sympt_unvac" (dimensions are age, serotype, infection type)
  for(age in 1:parms_epi_selected$parms_epi$N_age_groups){
    psy_DENV1_unvac[age,,,]<-Parms_GetPsy(prop_sympt=c(prop_sympt_unvac[age,1,1],prop_sympt_unvac[age,1,2],prop_sympt_unvac[age,1,3]),   #DENV-1 - middle index of prop_sympt_unvac = 1
                                          level_infect_sympt=parms_epi_selected$parms_epi$level_infect_sympt_unvac,
                                          level_infect_asympt=parms_epi_selected$parms_epi$level_infect_asympt_unvac)
    psy_DENV2_unvac[age,,,]<-Parms_GetPsy(prop_sympt=c(prop_sympt_unvac[age,2,1],prop_sympt_unvac[age,2,2],prop_sympt_unvac[age,2,3]),   #DENV-2 - middle index of prop_sympt_unvac = 2 etc.
                                          level_infect_sympt=parms_epi_selected$parms_epi$level_infect_sympt_unvac,
                                          level_infect_asympt=parms_epi_selected$parms_epi$level_infect_asympt_unvac)
    psy_DENV3_unvac[age,,,]<-Parms_GetPsy(prop_sympt=c(prop_sympt_unvac[age,3,1],prop_sympt_unvac[age,3,2],prop_sympt_unvac[age,3,3]),   
                                          level_infect_sympt=parms_epi_selected$parms_epi$level_infect_sympt_unvac,
                                          level_infect_asympt=parms_epi_selected$parms_epi$level_infect_asympt_unvac)
    psy_DENV4_unvac[age,,,]<-Parms_GetPsy(prop_sympt=c(prop_sympt_unvac[age,4,1],prop_sympt_unvac[age,4,2],prop_sympt_unvac[age,4,3]),
                                          level_infect_sympt=parms_epi_selected$parms_epi$level_infect_sympt_unvac,
                                          level_infect_asympt=parms_epi_selected$parms_epi$level_infect_asympt_unvac)
}
  #---------- Vectors -----------------
  NV<-NH*parms_epi_selected$parms_epi$ratio_VH               #Vector population size (calculated from the host population size and vector-to-host ratio) --> numeric
  muV<-1/parms_epi_selected$parms_epi$LE_V                   #Vector daily mortality rate (inverse of the life expectancy) --> numeric
  
  if (!parms_epi_selected$parms_epi$seasonality){            #Shape of the seasonality sine function. If seasonality if not included, both parameters are set to zero
    season_p1<-0                                             
    season_p2<-0
  }else{
    season_p1<-parms_epi_selected$parms_epi$season_p1
    season_p2<-parms_epi_selected$parms_epi$season_p2
  }
  
  xiV<-1/parms_epi_selected$parms_epi$dur_latency_V          #Daily transition rate from "Vector, exposed" to "Vector, infectious" (inverse of latency duration) --> numeric
  
  
  
  #---------- Parameter list ----------
  #Parameters used in the C part of the model
  parms_epi_dyn_model=list(ratio_VH=parms_epi_selected$parms_epi$ratio_VH,               #Numeric values
                           muV=muV,                         
                           b=parms_epi_selected$parms_epi$b,                             
                           ro=ro,                          
                           xiH=xiH, 
                           xiV=xiV, 
                           phiCP=phiCP, 
                           season_p1=season_p1,
                           season_p2=season_p2,  
                           
                           dzetaCP=parms_epi_selected$parms_epi$dzetaCP,                  #Booleans
                           seasonality=parms_epi_selected$parms_epi$seasonality, 
                           
                           betaVH=parms_epi_selected$parms_epi$betaVH,                    #Vectors
                           betaVH_age_coef=betaVH_age_coef,
                           betaHV=parms_epi_selected$parms_epi$betaHV, 
                           muHi=muHi, 
                           
                           theta_unvac=theta_unvac,                                       #Arrays
                           psy_DENV1_unvac=psy_DENV1_unvac, 
                           psy_DENV2_unvac=psy_DENV2_unvac, 
                           psy_DENV3_unvac=psy_DENV3_unvac, 
                           psy_DENV4_unvac=psy_DENV4_unvac, 
                           past_infections=past_infections)       
  
  #Parameters not used in the C part of the model but used throughout R codes
  parms_epi_other=list(N_age_groups=N_age_groups,
                       age_group_lower_bound=age_group_lower_bound,
                       age_group_upper_bound=age_group_upper_bound,
                       ai=ai,
                       NHi=NHi,
                       NH=NH,
                       NV=NV,
                       nb_infect_max=parms_epi_selected$parms_epi$nb_infect_max,
                       prop_sympt_unvac=prop_sympt_unvac,
                       prop_severe_sympt_unvac=prop_severe_sympt_unvac,
                       prop_hospit_sympt_unvac=prop_hospit_sympt_unvac,
                       prop_death_hospit_sympt_unvac=prop_death_hospit_sympt_unvac,
                       level_infect_sympt_unvac=parms_epi_selected$parms_epi$level_infect_sympt_unvac,
                       level_infect_asympt_unvac=parms_epi_selected$parms_epi$level_infect_asympt_unvac,
                       gammaCP=parms_epi_selected$parms_epi$gammaCP,
                       dzetaCE=parms_epi_selected$parms_epi$dzetaCE,
                       gammaCE=parms_epi_selected$parms_epi$gammaCE,
                       life_expectancy=parms_epi_selected$parms_epi$life_expectancy)
  
  #----------------------------------------------------
  #--- PARAMETERS FOR THE CALIBRATION -----------------
  #----------------------------------------------------
  
  if(calibration){
    
    with(parms_epi_selected$parms_calibr,{ 
      
      #--------- Incidence of symptomatic dengue ----------
      incidence_by_age<-cbind(incidence_age_group_lower_bound,
                              incidence_age_group_upper_bound,
                              incidence)
      
      #--------- Seasonality data -------------------------
      
      #--------- Seroprevalence data ----------------------
      #Identify missing values
      seroprev_temp=cbind(nx, nmz)                                   #Merge seroprevalence data into a single matrix
      for(i in 1:N_age_groups_seroprevalence){
        if(sum(seroprev_temp[i,])==0){                               #If the value for each serostatus is zero in a given age group --> no data, set to NA not to be considered in the calibration
          seroprev_temp[i,]<-NA 
        }
      }
      
      #Split into single-year age groups 
      seroprev_101<-matrix(NA, ncol=2, nrow=101)
      for(i in 1:N_age_groups_seroprevalence){
        if(length((seroprevalence_age_group_lower_bound[i]+1):(seroprevalence_age_group_upper_bound[i]+1))==1){  #if seroprevalence data are given for a single-year age group...
          seroprev_101[seroprevalence_age_group_lower_bound[i]+1,]<-seroprev_temp[i,]                            #... then copy it into a new matrix
        }else{
          length_age_group<-seroprevalence_age_group_upper_bound[i]-seroprevalence_age_group_lower_bound[i]+1    #Otherwise, calculate the duration of the used age group...
          for(j in (seroprevalence_age_group_lower_bound[i]+1):(seroprevalence_age_group_upper_bound[i]+1)){     #... and distribute uniformaly between all years in the age group
            seroprev_101[j,]<-seroprev_temp[i,]/length_age_group 
          }
        }
      }
      
      #Aggregate for the age groups used in the calibration
      seroprev_observed_2states<-matrix(NA, nrow=N_age_groups, ncol=2)                                                    
      for(i in 1:N_age_groups){
        if(length((age_group_lower_bound[i]+1):(age_group_upper_bound[i]+1)) == 1){                              #If a single-year age group
          seroprev_observed_2states[i,]<-seroprev_101[age_group_lower_bound[i]+1,]                               #... then copy the seroprevalence
        }else{                                                                                                   #Otherwise, sum for the composite age groups
          seroprev_observed_2states[i,]<-apply(seroprev_101[(age_group_lower_bound[i]+1):(age_group_upper_bound[i]+1),], 2, sum)
        }
      }
      
      #--------- Main model population --------------------                              #Population structure and size from the main model are required to re-calculate the initial state obtained in the calibration process
      pop_main_model<-list(NHi=parms_epi_selected$parms_epi$pop_size,	
                           NH=sum(parms_epi_selected$parms_epi$pop_size), 
                           ai=(1/(parms_epi_selected$parms_epi$age_group_upper_bound-parms_epi_selected$parms_epi$age_group_lower_bound+1)),
                           age_group_lower_bound_model=parms_epi_selected$parms_epi$age_group_lower_bound,
                           age_group_upper_bound_model=parms_epi_selected$parms_epi$age_group_upper_bound,
                           muHi=parms_epi_selected$parms_epi$mortality)                                      	   
      
      #--------- Parameter list ---------------------------
      parms_epi_other<-c(parms_epi_other, list(incidence_by_age_observed=incidence_by_age,
                                               incidence_observed_per=incidence_per, 
                                               incidence_by_month_observed=incidence_by_month,
                                               seroprev_observed_2states=seroprev_observed_2states,     #Calibration parameters are added to the list parms_epi_other
                                               pop_main_model=pop_main_model,
                                               expansion_factor_hospit=expansion_factor_hospit,
                                               expansion_factor_non_hospit=expansion_factor_non_hospit))
      
      return(list(parms_epi_dyn_model=parms_epi_dyn_model, 
                  parms_epi_other=parms_epi_other))
    })
  }else{
    return(list(parms_epi_dyn_model=parms_epi_dyn_model, 
                parms_epi_other=parms_epi_other))
    
  }
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Wrapper functions that imports, formats, calculates and outputs population & transmissions inputs

ParmsEpi_Load<-function(path_inputs_epi,            #Path to the folder containing the Excel file with the population & transmission inputs
                        set_epi,                    #Name of the Excel file containing the inputs (without extension)
                        calibration){               #Calibration flag (if TRUE, the parameters required for the calibration process will be output as well). 
  parms_epi_imported<-Parms_Import(path=path_inputs_epi, filename=set_epi)
  parms_epi_selected<-ParmsEpi_Select(parms_epi_imported=parms_epi_imported, calibration=calibration)
  parms_epi_final<-ParmsEpi_Initialise(parms_epi_selected=parms_epi_selected, calibration=calibration)
  
  return(parms_epi_final)
}