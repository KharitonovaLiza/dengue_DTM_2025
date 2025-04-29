##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contains functions that import, select, initalise and load all the vaccine effect parameters
###########################################################################################################################################################################################################

#This script contains the following functions:
# - ParmsVacEff_Select - This function generates a formatted list of all vaccine effects inputs from the inputs imported from an Excel file
# - ParmsVacEff_AggregateEff - This function aggregates the vaccine efficacy inputs to be used in the transmission module and outside it
# - ParmsVacEff_Initialise - This functions initializes all the efficacy-related inputs that will be used when running the model
# - ParmsVacEff_Load - Wrapper functions that imports, formats, calculates and outputs vaccination parameters
# - ParmsVacEff_GetThetaAndPsyYear - This function generates objects theta and psy for vaccinated individuals for the current simulation year
# - ParmsVacEff_GetEffectiveness - This function estimates effectiveness for a specific vaccinated cohort in a specific simulation year by taking into account the boosting of efficacy by each natural infection
# - ParmsVacEff_GetSevDistr - This function calculates the proportion of different severity of infections in vaccinated individuals for the year that is about to be simulated
# - ParmsVacEff_GetEffInf - This function estimates, for a given cohort, the efficacy/effectiveness against overall dengue infection based on the efficacy against each dengue outcome (asymptomatic, symptomatic non-hospitalized, symptomatic hospitalized) 
#                           and their relative probabilities in the given cohort for the current simulation year
# - ParmsVacEff_GetCohortDistr - This function estimates average weighted proportion (probability) of each dengue outcome (asymptomatic, symptomatic non-hospitalized, symptomatic hospitalized) 
#                                in a given cohort and for the current simulation year. It, thus, takes into account, what parts of the cohorts are at risk of 1st, 2nd, 3rd or 4th infection

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function generates a formatted list of all vaccine effects inputs from the inputs imported from an Excel file

ParmsVacEff_Select<-function(parms_vac_eff_imported){   #Inputs imported from the Excel file VacEff_XXXXX (output of the function Parms_Import)
  
  with(parms_vac_eff_imported,{ 
    
    parms_vac_eff<-list(
      
      #--- General settings --------------------------------------------------      
      vac_switch=as.logical(vac_switch),                                      #Inclusion of vaccination --> Boolean
      
      boosting_asymptomatic_switch=as.logical(boosting_asymptomatic_switch),  #Inclusion of natural boosting for the efficacy against asymptomatic dengue --> Boolean
      boosting_symptomatic_switch=as.logical(boosting_symptomatic_switch),    #Inclusion of natural boosting for the efficacy against  symptomatic dengue --> Boolean
      vac_as_silent_inf=as.logical(vac_as_silent_inf),                        #Vaccination acts as a silent infection --> Boolean
      
      sympt_becomes_asympt=as.logical(sympt_becomes_asympt),                  #Flag whether symptomatic infections avoided in vaccinated subjects become asymptomatic instead (TRUE) or are eliminated completely (FALSE) --> Boolean
      
      #--- Eff against asympt in seronegative --------------------------------
      eff_asympt_denv1_neg_ep1=as.numeric(eff_asympt_denv1_neg_ep1),  #Vaccine efficacy against asymptomatic infections with DENV-1 in seronegative at vaccination individuals (for 1st episode post-vaccination if boosting is included or overall if boosting is not included) --> Vector of numerics (nb of elements equal to max. allowed timeframe)
      eff_asympt_denv1_neg_ep2=as.numeric(eff_asympt_denv1_neg_ep2),  #... for 2nd episode post-vaccination if boosting is included; not used otherwise
      eff_asympt_denv1_neg_ep3=as.numeric(eff_asympt_denv1_neg_ep3),  #... for 3rd episode post-vaccination if boosting is included; not used otherwise
      eff_asympt_denv1_neg_ep4=as.numeric(eff_asympt_denv1_neg_ep4),  #... for 4th episode post-vaccination if boosting is included; not used otherwise
      
      eff_asympt_denv2_neg_ep1=as.numeric(eff_asympt_denv2_neg_ep1),  #Same for DENV-2 in seronegative
      eff_asympt_denv2_neg_ep2=as.numeric(eff_asympt_denv2_neg_ep2),  
      eff_asympt_denv2_neg_ep3=as.numeric(eff_asympt_denv2_neg_ep3),  
      eff_asympt_denv2_neg_ep4=as.numeric(eff_asympt_denv2_neg_ep4),  
      
      eff_asympt_denv3_neg_ep1=as.numeric(eff_asympt_denv3_neg_ep1),  #Same for DENV-3 in seronegative
      eff_asympt_denv3_neg_ep2=as.numeric(eff_asympt_denv3_neg_ep2),  
      eff_asympt_denv3_neg_ep3=as.numeric(eff_asympt_denv3_neg_ep3),  
      eff_asympt_denv3_neg_ep4=as.numeric(eff_asympt_denv3_neg_ep4),  
      
      eff_asympt_denv4_neg_ep1=as.numeric(eff_asympt_denv4_neg_ep1),  #Same for DENV-4 in seronegative
      eff_asympt_denv4_neg_ep2=as.numeric(eff_asympt_denv4_neg_ep2),  
      eff_asympt_denv4_neg_ep3=as.numeric(eff_asympt_denv4_neg_ep3),  
      eff_asympt_denv4_neg_ep4=as.numeric(eff_asympt_denv4_neg_ep4),  
      
      #--- Eff against asympt in seropositive --------------------------------      
      eff_asympt_denv1_pos_ep1=as.numeric(eff_asympt_denv1_pos_ep1),  #Same for DENV-1 in seropositive
      eff_asympt_denv1_pos_ep2=as.numeric(eff_asympt_denv1_pos_ep2),  
      eff_asympt_denv1_pos_ep3=as.numeric(eff_asympt_denv1_pos_ep3),  
      eff_asympt_denv1_pos_ep4=as.numeric(eff_asympt_denv1_pos_ep4),  
      
      eff_asympt_denv2_pos_ep1=as.numeric(eff_asympt_denv2_pos_ep1),  #Same for DENV-2 in seropositive
      eff_asympt_denv2_pos_ep2=as.numeric(eff_asympt_denv2_pos_ep2),  
      eff_asympt_denv2_pos_ep3=as.numeric(eff_asympt_denv2_pos_ep3),  
      eff_asympt_denv2_pos_ep4=as.numeric(eff_asympt_denv2_pos_ep4),  
      
      eff_asympt_denv3_pos_ep1=as.numeric(eff_asympt_denv3_pos_ep1),  #Same for DENV-3 in seropositive
      eff_asympt_denv3_pos_ep2=as.numeric(eff_asympt_denv3_pos_ep2),  
      eff_asympt_denv3_pos_ep3=as.numeric(eff_asympt_denv3_pos_ep3),  
      eff_asympt_denv3_pos_ep4=as.numeric(eff_asympt_denv3_pos_ep4),  
      
      eff_asympt_denv4_pos_ep1=as.numeric(eff_asympt_denv4_pos_ep1),  #Same for DENV-4 in seropositive
      eff_asympt_denv4_pos_ep2=as.numeric(eff_asympt_denv4_pos_ep2),  
      eff_asympt_denv4_pos_ep3=as.numeric(eff_asympt_denv4_pos_ep3),  
      eff_asympt_denv4_pos_ep4=as.numeric(eff_asympt_denv4_pos_ep4),  
      
      #--- Eff against sympt non-hosp in seronegative ------------------------
      eff_sympt_non_hosp_denv1_neg_ep1=as.numeric(eff_sympt_non_hosp_denv1_neg_ep1),  #Vaccine efficacy against symptomatic non-hospitalized infections with DENV-1 in seronegative at vaccination individuals (for 1st episode post-vaccination if boosting is included or overall if boosting is not included) --> Vector of numerics (nb of elements equal to max. allowed timeframe)
      eff_sympt_non_hosp_denv1_neg_ep2=as.numeric(eff_sympt_non_hosp_denv1_neg_ep2),  #... for 2nd episode post-vaccination if boosting is included; not used otherwise
      eff_sympt_non_hosp_denv1_neg_ep3=as.numeric(eff_sympt_non_hosp_denv1_neg_ep3),  #... for 3rd episode post-vaccination if boosting is included; not used otherwise
      eff_sympt_non_hosp_denv1_neg_ep4=as.numeric(eff_sympt_non_hosp_denv1_neg_ep4),  #... for 4th episode post-vaccination if boosting is included; not used otherwise
      
      eff_sympt_non_hosp_denv2_neg_ep1=as.numeric(eff_sympt_non_hosp_denv2_neg_ep1),  #Same for DENV-2 in seronegative
      eff_sympt_non_hosp_denv2_neg_ep2=as.numeric(eff_sympt_non_hosp_denv2_neg_ep2),  
      eff_sympt_non_hosp_denv2_neg_ep3=as.numeric(eff_sympt_non_hosp_denv2_neg_ep3),  
      eff_sympt_non_hosp_denv2_neg_ep4=as.numeric(eff_sympt_non_hosp_denv2_neg_ep4),  
      
      eff_sympt_non_hosp_denv3_neg_ep1=as.numeric(eff_sympt_non_hosp_denv3_neg_ep1),  #Same for DENV-3 in seronegative
      eff_sympt_non_hosp_denv3_neg_ep2=as.numeric(eff_sympt_non_hosp_denv3_neg_ep2),  
      eff_sympt_non_hosp_denv3_neg_ep3=as.numeric(eff_sympt_non_hosp_denv3_neg_ep3),  
      eff_sympt_non_hosp_denv3_neg_ep4=as.numeric(eff_sympt_non_hosp_denv3_neg_ep4),  
      
      eff_sympt_non_hosp_denv4_neg_ep1=as.numeric(eff_sympt_non_hosp_denv4_neg_ep1),  #Same for DENV-4 in seronegative
      eff_sympt_non_hosp_denv4_neg_ep2=as.numeric(eff_sympt_non_hosp_denv4_neg_ep2),  
      eff_sympt_non_hosp_denv4_neg_ep3=as.numeric(eff_sympt_non_hosp_denv4_neg_ep3),  
      eff_sympt_non_hosp_denv4_neg_ep4=as.numeric(eff_sympt_non_hosp_denv4_neg_ep4),  
      
      #--- Eff against sympt non-hosp in seropositive ------------------------
      eff_sympt_non_hosp_denv1_pos_ep1=as.numeric(eff_sympt_non_hosp_denv1_pos_ep1),  #Same for DENV-1 in seropositive
      eff_sympt_non_hosp_denv1_pos_ep2=as.numeric(eff_sympt_non_hosp_denv1_pos_ep2),  
      eff_sympt_non_hosp_denv1_pos_ep3=as.numeric(eff_sympt_non_hosp_denv1_pos_ep3),  
      eff_sympt_non_hosp_denv1_pos_ep4=as.numeric(eff_sympt_non_hosp_denv1_pos_ep4),  
      
      eff_sympt_non_hosp_denv2_pos_ep1=as.numeric(eff_sympt_non_hosp_denv2_pos_ep1),  #Same for DENV-2 in seropositive
      eff_sympt_non_hosp_denv2_pos_ep2=as.numeric(eff_sympt_non_hosp_denv2_pos_ep2),  
      eff_sympt_non_hosp_denv2_pos_ep3=as.numeric(eff_sympt_non_hosp_denv2_pos_ep3),  
      eff_sympt_non_hosp_denv2_pos_ep4=as.numeric(eff_sympt_non_hosp_denv2_pos_ep4),  
      
      eff_sympt_non_hosp_denv3_pos_ep1=as.numeric(eff_sympt_non_hosp_denv3_pos_ep1),  #Same for DENV-3 in seropositive
      eff_sympt_non_hosp_denv3_pos_ep2=as.numeric(eff_sympt_non_hosp_denv3_pos_ep2),  
      eff_sympt_non_hosp_denv3_pos_ep3=as.numeric(eff_sympt_non_hosp_denv3_pos_ep3),  
      eff_sympt_non_hosp_denv3_pos_ep4=as.numeric(eff_sympt_non_hosp_denv3_pos_ep4),  
      
      eff_sympt_non_hosp_denv4_pos_ep1=as.numeric(eff_sympt_non_hosp_denv4_pos_ep1),  #Same for DENV-4 in seropositive
      eff_sympt_non_hosp_denv4_pos_ep2=as.numeric(eff_sympt_non_hosp_denv4_pos_ep2),  
      eff_sympt_non_hosp_denv4_pos_ep3=as.numeric(eff_sympt_non_hosp_denv4_pos_ep3),  
      eff_sympt_non_hosp_denv4_pos_ep4=as.numeric(eff_sympt_non_hosp_denv4_pos_ep4),  
      
      #--- Eff against sympt hosp in seronegative ----------------------------
      eff_sympt_hosp_denv1_neg_ep1=as.numeric(eff_sympt_hosp_denv1_neg_ep1),          #Vaccine efficacy against symptomatic hospitalized infections with DENV-1 in seronegative at vaccination individuals (for 1st episode post-vaccination if boosting is included or overall if boosting is not included) --> Vector of numerics (nb of elements equal to max. allowed timeframe)
      eff_sympt_hosp_denv1_neg_ep2=as.numeric(eff_sympt_hosp_denv1_neg_ep2),          #... for 2nd episode post-vaccination if boosting is included; not used otherwise
      eff_sympt_hosp_denv1_neg_ep3=as.numeric(eff_sympt_hosp_denv1_neg_ep3),          #... for 3rd episode post-vaccination if boosting is included; not used otherwise
      eff_sympt_hosp_denv1_neg_ep4=as.numeric(eff_sympt_hosp_denv1_neg_ep4),          #... for 4th episode post-vaccination if boosting is included; not used otherwise
      
      eff_sympt_hosp_denv2_neg_ep1=as.numeric(eff_sympt_hosp_denv2_neg_ep1),          #Same for DENV-2 in seronegative
      eff_sympt_hosp_denv2_neg_ep2=as.numeric(eff_sympt_hosp_denv2_neg_ep2),  
      eff_sympt_hosp_denv2_neg_ep3=as.numeric(eff_sympt_hosp_denv2_neg_ep3),  
      eff_sympt_hosp_denv2_neg_ep4=as.numeric(eff_sympt_hosp_denv2_neg_ep4),  
      
      eff_sympt_hosp_denv3_neg_ep1=as.numeric(eff_sympt_hosp_denv3_neg_ep1),          #Same for DENV-3 in seronegative
      eff_sympt_hosp_denv3_neg_ep2=as.numeric(eff_sympt_hosp_denv3_neg_ep2),  
      eff_sympt_hosp_denv3_neg_ep3=as.numeric(eff_sympt_hosp_denv3_neg_ep3),  
      eff_sympt_hosp_denv3_neg_ep4=as.numeric(eff_sympt_hosp_denv3_neg_ep4),  
      
      eff_sympt_hosp_denv4_neg_ep1=as.numeric(eff_sympt_hosp_denv4_neg_ep1),          #Same for DENV-4 in seronegative
      eff_sympt_hosp_denv4_neg_ep2=as.numeric(eff_sympt_hosp_denv4_neg_ep2),  
      eff_sympt_hosp_denv4_neg_ep3=as.numeric(eff_sympt_hosp_denv4_neg_ep3),  
      eff_sympt_hosp_denv4_neg_ep4=as.numeric(eff_sympt_hosp_denv4_neg_ep4),  
      
      #--- Eff against sympt hosp in seropositive ----------------------------
      eff_sympt_hosp_denv1_pos_ep1=as.numeric(eff_sympt_hosp_denv1_pos_ep1),          #Same for DENV-1 in seropositive
      eff_sympt_hosp_denv1_pos_ep2=as.numeric(eff_sympt_hosp_denv1_pos_ep2),  
      eff_sympt_hosp_denv1_pos_ep3=as.numeric(eff_sympt_hosp_denv1_pos_ep3),  
      eff_sympt_hosp_denv1_pos_ep4=as.numeric(eff_sympt_hosp_denv1_pos_ep4),  
      
      eff_sympt_hosp_denv2_pos_ep1=as.numeric(eff_sympt_hosp_denv2_pos_ep1),          #Same for DENV-2 in seropositive
      eff_sympt_hosp_denv2_pos_ep2=as.numeric(eff_sympt_hosp_denv2_pos_ep2),  
      eff_sympt_hosp_denv2_pos_ep3=as.numeric(eff_sympt_hosp_denv2_pos_ep3),  
      eff_sympt_hosp_denv2_pos_ep4=as.numeric(eff_sympt_hosp_denv2_pos_ep4),  
      
      eff_sympt_hosp_denv3_pos_ep1=as.numeric(eff_sympt_hosp_denv3_pos_ep1),          #Same for DENV-3 in seropositive
      eff_sympt_hosp_denv3_pos_ep2=as.numeric(eff_sympt_hosp_denv3_pos_ep2),  
      eff_sympt_hosp_denv3_pos_ep3=as.numeric(eff_sympt_hosp_denv3_pos_ep3),  
      eff_sympt_hosp_denv3_pos_ep4=as.numeric(eff_sympt_hosp_denv3_pos_ep4),  
      
      eff_sympt_hosp_denv4_pos_ep1=as.numeric(eff_sympt_hosp_denv4_pos_ep1),          #Same for DENV-4 in seropositive
      eff_sympt_hosp_denv4_pos_ep2=as.numeric(eff_sympt_hosp_denv4_pos_ep2),  
      eff_sympt_hosp_denv4_pos_ep3=as.numeric(eff_sympt_hosp_denv4_pos_ep3),  
      eff_sympt_hosp_denv4_pos_ep4=as.numeric(eff_sympt_hosp_denv4_pos_ep4),  
      
      #--- Relative infectiousness of vaccinees ------------------------------
      level_infect_sympt_vac=as.numeric(level_infect_sympt_vac),   				    #Infectiousness (contribution into the force of infection) of a vaccinated individual with symptomatic infection (relative to an unvaccinated individual with asymptomatic infection) --> Numeric
      level_infect_asympt_vac=as.numeric(level_infect_asympt_vac)   			    #Infectiousness (contribution into the force of infection) of a vaccinated individual with asymptomatic infection (relative to an unvaccinated individual with asymptomatic infection) --> Numeric
    )   																		 			 			  
    
    return(parms_vac_eff)
  })
}


#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function aggregates the vaccine efficacy inputs into larger arrays to be used in the other functions

ParmsVacEff_AggregateEff<-function(parms_vac_eff_selected,         #Pre-selected vaccine effects inputs
                                   timeframe){                     #Simulation timeframe
  
  with(parms_vac_eff_selected, {
    
    #-------------------------------------------------
    #--- Eff. against asymptomatic dengue ------------
    #-------------------------------------------------
    if(!(boosting_asymptomatic_switch)){                   #If natural boosting of efficacy against asymptomatic dengue is not included, only using the efficacy against 1st episode (and applying it based on the time elapsed since vaccination)
      eff_asympt<-array(NA, dim=c(4,2,timeframe))          #Efficacy by serotype (4), serostatus (2), year post-vaccination     
      for(serotype in 1:4){                                             
        for(status in 1:2){
          if(status==1){status_suffix="neg"}
          if(status==2){status_suffix="pos"}
          eff_asympt[serotype,status,]<-get(paste("eff_asympt_denv",serotype,"_",status_suffix,"_ep1",sep=""))[1:timeframe]
        }
      }
    }else{                                                 #...otherwise adding another dimension to the array to track efficacy for each episode post-vaccination
      eff_asympt<-array(NA, dim=c(4,2,4,timeframe))        #Efficacy by serotype (4), serostatus (2), episode post-vaccination (4), year post-boosting
      for(serotype in 1:4){                                             
        for(status in 1:2){
          if(status==1){status_suffix="neg"}
          if(status==2){status_suffix="pos"}
          for(episode in 1:4){
            eff_asympt[serotype,status,episode,]<-get(paste("eff_asympt_denv",serotype,"_",status_suffix,"_ep",episode,sep=""))[1:timeframe]
          }
        }
      }
    }
    
    #-------------------------------------------------             #Same general logic as for the efficacy against asymptomatic infections above
    #--- Eff. against sympt non-hospitalized dengue --
    #-------------------------------------------------
    if(!(boosting_symptomatic_switch)){                            #Note the change in the variable (switches for boosting are separate for asymptomatic and symptomatic infections)
      eff_sympt_non_hosp<-array(NA, dim=c(4,2,timeframe))
      for(serotype in 1:4){                                             
        for(status in 1:2){
          if(status==1){status_suffix="neg"}
          if(status==2){status_suffix="pos"}
          eff_sympt_non_hosp[serotype,status,]<-get(paste("eff_sympt_non_hosp_denv",serotype,"_",status_suffix,"_ep1",sep=""))[1:timeframe]
        }
      }
    }else{
      eff_sympt_non_hosp<-array(NA, dim=c(4,2,4,timeframe))
      for(serotype in 1:4){                                             
        for(status in 1:2){
          if(status==1){status_suffix="neg"}
          if(status==2){status_suffix="pos"}
          for(episode in 1:4){
            eff_sympt_non_hosp[serotype,status,episode,]<-get(paste("eff_sympt_non_hosp_denv",serotype,"_",status_suffix,"_ep",episode,sep=""))[1:timeframe]
          }
        }
      }
    }
    
    #-------------------------------------------------             #Same logic as above
    #--- Eff. against sympt hospitalized dengue ------
    #-------------------------------------------------
    if(!(boosting_symptomatic_switch)){
      eff_sympt_hosp<-array(NA, dim=c(4,2,timeframe))          
      for(serotype in 1:4){                                             
        for(status in 1:2){
          if(status==1){status_suffix="neg"}
          if(status==2){status_suffix="pos"}
          eff_sympt_hosp[serotype,status,]<-get(paste("eff_sympt_hosp_denv",serotype,"_",status_suffix,"_ep1",sep=""))[1:timeframe]
        }
      }
    }else{                                                          
      eff_sympt_hosp<-array(NA, dim=c(4,2,4,timeframe))        
      for(serotype in 1:4){                                             
        for(status in 1:2){
          if(status==1){status_suffix="neg"}
          if(status==2){status_suffix="pos"}
          for(episode in 1:4){
            eff_sympt_hosp[serotype,status,episode,]<-get(paste("eff_sympt_hosp_denv",serotype,"_",status_suffix,"_ep",episode,sep=""))[1:timeframe]
          }
        }
      }
    }
    
    return(list(eff_asympt=eff_asympt,
                eff_sympt_non_hosp=eff_sympt_non_hosp,
                eff_sympt_hosp=eff_sympt_hosp))
  })
}



#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This functions initializes all the efficacy-related inputs that will be used when running the model
#Note: additional parameters will be calculated at each year of the simulation, using this inputs

ParmsVacEff_Initialise<- function(parms_vac_eff_selected,     #Selected and pre-formatted vaccine effects inputs
                                  parms_vac_strat,            #Vaccination strategy parameters
                                  timeframe){
  
  #Arrays with efficacy inputs
  eff_inputs<-ParmsVacEff_AggregateEff(parms_vac_eff_selected=parms_vac_eff_selected, 
                                       timeframe=timeframe)
  
  parms_vac_eff<-list(eff_asympt=eff_inputs$eff_asympt,
                      eff_sympt_non_hosp=eff_inputs$eff_sympt_non_hosp,
                      eff_sympt_hosp=eff_inputs$eff_sympt_hosp,
                      vac_cohorts=parms_vac_strat$parms_vac_strat$vac_cohorts,                           #Matrix indicating whether a cohort is vaccinated or not 
                      vac_cohorts_cov=parms_vac_strat$parms_vac_strat$vac_cohorts_cov,
                      level_infect_sympt_vac=parms_vac_eff_selected$level_infect_sympt_vac,
                      level_infect_asympt_vac=parms_vac_eff_selected$level_infect_asympt_vac,
                      sympt_becomes_asympt=parms_vac_eff_selected$sympt_becomes_asympt,
                      boosting_asymptomatic_switch=parms_vac_eff_selected$boosting_asymptomatic_switch,
                      boosting_symptomatic_switch=parms_vac_eff_selected$boosting_symptomatic_switch,
                      vac_as_silent_inf=parms_vac_eff_selected$vac_as_silent_inf)
  
  return(list(parms_vac_eff=parms_vac_eff))
  
}



#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Wrapper functions that imports, formats, calculates and outputs vaccination parameters

ParmsVacEff_Load<-function(path_inputs_vac_eff,        #Path to the folder containing the Excel file with the vaccination effects inputs
                           set_vac_eff,                #Name of the Excel file containing the inputs (without extension) 
                           parms_vac_strat,            #Vaccination strategy parameters
                           timeframe){
  
  parms_vac_eff_imported<-Parms_Import(path=path_inputs_vac_eff, filename=set_vac_eff)        
  parms_vac_eff_selected<-ParmsVacEff_Select(parms_vac_eff_imported=parms_vac_eff_imported)
  parms_vac_eff_final<-ParmsVacEff_Initialise(parms_vac_eff_selected=parms_vac_eff_selected, parms_vac_strat=parms_vac_strat, timeframe=timeframe)
  return(parms_vac_eff_final)
}


#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function generates objects theta and psy for vaccinated individuals for the current simulation year

ParmsVacEff_GetThetaAndPsyYear<-function(year,                     #Current simulation year
                                         parms_epi,                #Population & transmission inputs
                                         parms_vac_eff,            #Vaccine efficacy inputs (previously aggregated and processed)
                                         parms_vac_strat,          #Vaccination strategy inputs
                                         yini_6D,                  #Array with the initial state for the hosts for the current year (post-ageing, post-vaccination)
                                         trace){                   #Markov trace (required to estimate the efficacy if efficacy boosting is included),
  
  #-----------------------------------------------------------------------------------------------  
  #---- Get efficacy values for the current simulation year --------------------------------------
  #-----------------------------------------------------------------------------------------------  
  #Previously we aggregated efficacy inputs into larger arrays, where the values are efficacy by serotype, serostatus at vaccination and year post-vaccination (and by episode post-vaccination if boosting is included)
  #However, all cohorts are vaccinated at different years, so year post-vaccination is not equal to calendar (simulation) year
  #We need to estimate, for the currently simulated year, the level of protection in each simulated cohort
  #For that we use the previously created matrix vac_cohorts_ind
  
  #-----------------------------------
  #--- Initialize the objects --------
  #-----------------------------------
  eff_asympt        <-array(0, dim=c(101,4,2)) #3D arrays with efficacy/effectiveness values. Dimensions are age (single year), serotype, serostatus at vaccination; for the currently simulated year
  eff_sympt_non_hosp<-array(0, dim=c(101,4,2)) #The values represent efficacy (if boosting is not included) or effectiveness (if boosting is included)
  eff_sympt_hosp    <-array(0, dim=c(101,4,2))
  eff_inf           <-array(0, dim=c(101,4,2)) #Additional object with the efficacy/effectiveness against infection, which will be calculated from the efficacy against different types of dengue (asympt, sympt non-hosp, sympt hosp)
  
  #-----------------------------------
  #--- Get values for the year -------
  #-----------------------------------
  
  #Loop through every year of age (age_ind). Note: here age_ind corresponds to the index of an age group (e.g. age_ind=1 indicates 0-year-olds)
  for(age_ind in 1:101){                                      
    
    #If the cohort aged age_ind is not vaccinated in the current simulation year (flag value in vac_cohorts_ind is 0), then zero efficacy for all serotypes and serostatuses
    if(parms_vac_strat$vac_cohorts_ind[age_ind,year]==0){       
      eff_asympt        [age_ind,,]<-0    
      eff_sympt_non_hosp[age_ind,,]<-0           
      eff_sympt_hosp    [age_ind,,]<-0
      eff_inf           [age_ind,,]<-0
      
    #If cohort is vaccinated, continue
    }else{
      
      #---------------
      #---------------
      #Get age at vaccination, year since vaccination etc. (may be required to define the current level of efficacy)
      
      #If vaccinated  (flag value in vac_cohorts_ind is 1)
      if(parms_vac_strat$vac_cohorts_ind[age_ind,year]==1){
        age_at_vac<-parms_vac_strat$vac_cohorts_age_at_vac[age_ind,year]   #...age at vaccination
        year_of_vac<-parms_vac_strat$vac_cohorts_year_of_vac[age_ind,year] #...year of vaccination
        year_since_vac<-year-year_of_vac+1                                 #...year since vaccination is the difference between current year and year of vaccination
      }   
      
      #---------------
      #---------------
      #Get EFFICACY/EFFECTIVENESS
      #If this is THE FIRST YEAR SINCE THE COHORT'S VACCINATION, then their efficacy is not (yet) determined by boosting
      #--> applying the efficacy of the first year (and against the 1st episode if boosting is included; as the cohort was just vaccinated and vaccination serves as the first boosting)
      #Note: Splitting the definition into sub-cases because the dimensions of the object parms_vac_eff$eff_XXX depend on whether or not the boosting is not included (3D if not, 4D if yes)
      #      In both cases first two dimensions are serotype and serostatus
      if(year==year_of_vac){
        
        #Get efficacy for ASYMPTOMATIC
        if(!parms_vac_eff$boosting_asymptomatic_switch){                                 #If boosting is not included
          eff_asympt[age_ind,,]<-parms_vac_eff$eff_asympt[,,1,   drop=FALSE]             #...taking values for each serotype & serostatus, for the 1st year  
        }else{                                                                           #If boosting is included
          eff_asympt[age_ind,,]<-parms_vac_eff$eff_asympt[,,1,1, drop=FALSE]             #...taking values for each serotype & serostatus, for the 1st year of the efficacy against the 1st episode  
        }
        
        #Get efficacy for SYMPTOMATIC (similar logic as for asymptomatic)
        if(!parms_vac_eff$boosting_symptomatic_switch){
          eff_sympt_non_hosp[age_ind,,]<-parms_vac_eff$eff_sympt_non_hosp[,,1,   drop=FALSE]    
          eff_sympt_hosp    [age_ind,,]<-parms_vac_eff$eff_sympt_hosp    [,,1,   drop=FALSE] 
        }else{
          eff_sympt_non_hosp[age_ind,,]<-parms_vac_eff$eff_sympt_non_hosp[,,1,1, drop=FALSE]              
          eff_sympt_hosp    [age_ind,,]<-parms_vac_eff$eff_sympt_hosp    [,,1,1, drop=FALSE]
        }

        
      #If this is NOT THE FIRST YEAR SINCE THE COHORT'S VACCINATION
      }else{
        #If boosting of efficacy (symptomatic or asymptomatic) is included, get effectiveness values for each type of dengue + the distribution of cohort by the type of dengue they are at risk at
        #Calling the function ParmsVacEff_GetEffectiveness for this
        if(parms_vac_eff$boosting_asymptomatic_switch || parms_vac_eff$boosting_symptomatic_switch){
          
          effectiveness<-ParmsVacEff_GetEffectiveness(year=year,
                                                      trace=trace,
                                                      vac_cohorts_nb=parms_vac_strat$trac_vac_cohorts_nb,
                                                      vac_cohorts_age_at_vac=parms_vac_strat$trac_vac_cohorts_age_at_vac,
                                                      vac_cohorts_year_of_vac=parms_vac_strat$trac_vac_cohorts_year_of_vac,
                                                      age_at_vac_loop=age_at_vac,
                                                      year_of_vac_loop=year_of_vac,
                                                      eff_asympt=parms_vac_eff$eff_asympt,
                                                      eff_sympt_non_hosp=parms_vac_eff$eff_sympt_non_hosp,
                                                      eff_sympt_hosp=parms_vac_eff$eff_sympt_hosp,
                                                      boosting_asymptomatic_switch=parms_vac_eff$boosting_asymptomatic_switch,
                                                      boosting_symptomatic_switch=parms_vac_eff$boosting_symptomatic_switch)
        }
        
        #Get efficacy/effectiveness against ASYMPTOMATIC dengue
        if(!parms_vac_eff$boosting_asymptomatic_switch){       #If boosting is not included, the efficacy depends on the time since vaccination
          for(year_since_vac_id in 1:(year_since_vac)){
            eff_asympt[age_ind,,]<-eff_asympt[age_ind,,]+parms_vac_eff$eff_asympt[,,year_since_vac_id]*parms_vac_strat$vac_cohorts_cov[age_ind,year,year_since_vac_id]/sum(parms_vac_strat$vac_cohorts_cov[age_ind,year,])
          }
          
        }else{                                                 #If boosting is included, the effectiveness was calculated already by the function ParmsVacEff_GetEffectiveness and is stored in the object effectiveness
          eff_asympt[age_ind,,]<-effectiveness$effect_asympt         
        }
        
        #Get efficacy/effectiveness against SYMPTOMATIC dengue (similar logic as for asymptomatic above)
        if(!parms_vac_eff$boosting_symptomatic_switch){      
          for(year_since_vac_id in 1:(year_since_vac)){
            eff_sympt_non_hosp[age_ind,,]<-eff_sympt_non_hosp[age_ind,,]+parms_vac_eff$eff_sympt_non_hosp[,,year_since_vac_id]*parms_vac_strat$vac_cohorts_cov[age_ind,year,year_since_vac_id]/sum(parms_vac_strat$vac_cohorts_cov[age_ind,year,])
            eff_sympt_hosp    [age_ind,,]<-eff_sympt_hosp    [age_ind,,]+parms_vac_eff$eff_sympt_hosp    [,,year_since_vac_id]*parms_vac_strat$vac_cohorts_cov[age_ind,year,year_since_vac_id]/sum(parms_vac_strat$vac_cohorts_cov[age_ind,year,])  
          }
        }else{                                                 
          eff_sympt_non_hosp[age_ind,,]<-effectiveness$effect_sympt_non_hosp
          eff_sympt_hosp    [age_ind,,]<-effectiveness$effect_sympt_hosp
        }
      }
      
      for(serostatus in 1:2){
        #Estimate the population at risk of a new infection, for a given serotype and a given serostatus
        for(serotype in 1:4){
          if     (serotype==1){pop_at_risk_total<-sum(yini_6D[age_ind,1, , , ,(serostatus+1)])}
          else if(serotype==2){pop_at_risk_total<-sum(yini_6D[age_ind, ,1, , ,(serostatus+1)])}
          else if(serotype==3){pop_at_risk_total<-sum(yini_6D[age_ind, , ,1, ,(serostatus+1)])}
          else                {pop_at_risk_total<-sum(yini_6D[age_ind, , , ,1,(serostatus+1)])}
          if(pop_at_risk_total==0){
            eff_asympt        [age_ind,serotype,serostatus]<-0
            eff_sympt_non_hosp[age_ind,serotype,serostatus]<-0 
            eff_sympt_hosp    [age_ind,serotype,serostatus]<-0 
          }
        }
      }
      
      #---------------
      #---------------
      #Get COHORT DISTRIBUTION by each type of dengue outcome (asympt, sympt non-hosp, sympt hosp) to be able to calculated weighted average efficacy/effectiveness
      cohort_distr<-ParmsVacEff_GetCohortDistr(age_ind=age_ind,     
                                               yini_6D=yini_6D,
                                               parms_epi=parms_epi)
      
      
      #---------------
      #---------------
      #Calculate efficacy/effectiveness against DENGUE INFECTION, input it in the final object
      eff_inf[age_ind,,]<-ParmsVacEff_GetEffInf(sympt_becomes_asympt=parms_vac_eff$sympt_becomes_asympt,
                                                eff_asympt=eff_asympt[age_ind,,],
                                                eff_sympt_non_hosp=eff_sympt_non_hosp[age_ind,,],
                                                eff_sympt_hosp=eff_sympt_hosp[age_ind,,],
                                                prop_asympt=cohort_distr$prop_asympt,
                                                prop_sympt_non_hosp=cohort_distr$prop_sympt_non_hosp,
                                                prop_sympt_hosp=cohort_distr$prop_sympt_hosp)
    }
  }
  
  #---------------
  #Aggregate all new objects in a single list
  efficacies<-list(eff_inf=eff_inf,
                   eff_asympt=eff_asympt,
                   eff_sympt_non_hosp=eff_sympt_non_hosp,
                   eff_sympt_hosp=eff_sympt_hosp)
  
  
  #-----------------------------------------------------------------------------------------------  
  #---- Get object theta for the current simulation year -----------------------------------------
  #-----------------------------------------------------------------------------------------------  
  #The object theta depends (among other things, see more details in the script "12_parms_theta.R") on the effectiveness against dengue infection
  #Preparing two such objects (separately for seronegative and seropositive vaccinees) for the year that is about to be simulated
  
  #--- Format for the transmission module ------
  #In the transmission module the values of efficacy against infection are stored in two arrays - separately for the individuals who are seronegative and seropositive at vaccination.
  
  #Initiate the arrays
  eta_vac_neg<-array(0, dim=c(parms_epi$N_age_groups,4))        
  eta_vac_pos<-array(0, dim=c(parms_epi$N_age_groups,4))        
  
  #Fill the arrays
  eta_vac_neg<-eff_inf[,,1]   #Take a "slice" of eff_inf for seronegative individuals (third index = 1)
  eta_vac_pos<-eff_inf[,,2]   #... for seropositive individuals (third index = 2)
  
  #--- Calculate object theta ------------------
  theta_vac_neg<-Parms_GetThetaVac(N_age_groups=parms_epi$N_age_groups,
                                   gammaCP=parms_epi$gammaCP,                         #Theta for the hosts vaccinated as seronegative
                                   dzetaCE=parms_epi$dzetaCE, 
                                   gammaCE=parms_epi$gammaCE, 
                                   nb_infect_max=parms_epi$nb_infect_max, 
                                   eta=eta_vac_neg) 
  
  theta_vac_pos<-Parms_GetThetaVac(N_age_groups=parms_epi$N_age_groups,
                                   gammaCP=parms_epi$gammaCP,                         #Theta for the hosts vaccinated as seropositive
                                   dzetaCE=parms_epi$dzetaCE, 
                                   gammaCE=parms_epi$gammaCE, 
                                   nb_infect_max=parms_epi$nb_infect_max, 
                                   eta=eta_vac_pos)
  
  #-----------------------------------------------------------------------------------------------  
  #---- Estimate the severity distribution in the vaccinated subjects ----------------------------
  #-----------------------------------------------------------------------------------------------  
  #The severity distribution (proportion of infections that are symptomatic, hospitalized etc) depends on the efficacy parameters
  #Calculating such parameters for two reasons:
  #- To determine the object psy (relative infectiousness which depends on the chance of being symptomatic, as symptomatic infections can be more transmissible)
  #- To store for later use (once the number of infections is generated by the transmission module, these proportions are used to estimated the number of cases of different severity)
  #These parameters are calculated for the current simulation year
  
  severities_vac<-ParmsVacEff_GetSevDistr(parms_epi=parms_epi, 
                                          parms_vac_eff=parms_vac_eff,
                                          efficacies=efficacies)
  
  #-----------------------------------------------------------------------------------------------  
  #---- Get the object psy for vaccinated individuals for the current simulation year ------------
  #-----------------------------------------------------------------------------------------------  
  #The object psy captures relative infectiousness of vaccinated individuals
  #More details can be found in the script "13_parms_psy.R"
  
  #Initialize the arrays psy for vaccinated individuals (separate object for each serotype and serostatus at vaccination)
  #Dimensions are age, status for three serotypes (other than the one being considered in a current infection process)
  psy_DENV1_vac_neg<-array(NA, dim=c(parms_epi$N_age_groups,5,5,5))      
  psy_DENV2_vac_neg<-array(NA, dim=c(parms_epi$N_age_groups,5,5,5))
  psy_DENV3_vac_neg<-array(NA, dim=c(parms_epi$N_age_groups,5,5,5))
  psy_DENV4_vac_neg<-array(NA, dim=c(parms_epi$N_age_groups,5,5,5))
  
  psy_DENV1_vac_pos<-array(NA, dim=c(parms_epi$N_age_groups,5,5,5))      
  psy_DENV2_vac_pos<-array(NA, dim=c(parms_epi$N_age_groups,5,5,5))
  psy_DENV3_vac_pos<-array(NA, dim=c(parms_epi$N_age_groups,5,5,5))
  psy_DENV4_vac_pos<-array(NA, dim=c(parms_epi$N_age_groups,5,5,5))
  
  #Create the objects
  for(age in 1:parms_epi$N_age_groups){
    
    #Reminder: in the object prop_sympt_vac[,,,] 2nd index indicates serotype; 3rd index indicates serostatus at vaccination (1=seronegative; 2=seropositive) 
    
    #For those seronegative at vaccination
    psy_DENV1_vac_neg[age,,,]<-Parms_GetPsy(prop_sympt=severities_vac$prop_sympt_vac[age,1,1,],level_infect_sympt=parms_vac_eff$level_infect_sympt_vac,level_infect_asympt=parms_vac_eff$level_infect_asympt_vac)
    psy_DENV2_vac_neg[age,,,]<-Parms_GetPsy(prop_sympt=severities_vac$prop_sympt_vac[age,2,1,],level_infect_sympt=parms_vac_eff$level_infect_sympt_vac,level_infect_asympt=parms_vac_eff$level_infect_asympt_vac)
    psy_DENV3_vac_neg[age,,,]<-Parms_GetPsy(prop_sympt=severities_vac$prop_sympt_vac[age,3,1,],level_infect_sympt=parms_vac_eff$level_infect_sympt_vac,level_infect_asympt=parms_vac_eff$level_infect_asympt_vac)
    psy_DENV4_vac_neg[age,,,]<-Parms_GetPsy(prop_sympt=severities_vac$prop_sympt_vac[age,4,1,],level_infect_sympt=parms_vac_eff$level_infect_sympt_vac,level_infect_asympt=parms_vac_eff$level_infect_asympt_vac)
    
    #For those seropositive at vaccination
    psy_DENV1_vac_pos[age,,,]<-Parms_GetPsy(prop_sympt=severities_vac$prop_sympt_vac[age,1,2,],level_infect_sympt=parms_vac_eff$level_infect_sympt_vac,level_infect_asympt=parms_vac_eff$level_infect_asympt_vac)
    psy_DENV2_vac_pos[age,,,]<-Parms_GetPsy(prop_sympt=severities_vac$prop_sympt_vac[age,2,2,],level_infect_sympt=parms_vac_eff$level_infect_sympt_vac,level_infect_asympt=parms_vac_eff$level_infect_asympt_vac)
    psy_DENV3_vac_pos[age,,,]<-Parms_GetPsy(prop_sympt=severities_vac$prop_sympt_vac[age,3,2,],level_infect_sympt=parms_vac_eff$level_infect_sympt_vac,level_infect_asympt=parms_vac_eff$level_infect_asympt_vac)
    psy_DENV4_vac_pos[age,,,]<-Parms_GetPsy(prop_sympt=severities_vac$prop_sympt_vac[age,4,2,],level_infect_sympt=parms_vac_eff$level_infect_sympt_vac,level_infect_asympt=parms_vac_eff$level_infect_asympt_vac)
  }
  
  #-----------------------------------------------------------------------------------------------  
  #---- Return all the objects required for running the simulation -------------------------------
  #-----------------------------------------------------------------------------------------------
  return(list(efficacies=efficacies,
              theta_vac_neg=theta_vac_neg,
              theta_vac_pos=theta_vac_pos,
              severities_vac=severities_vac,
              psy_DENV1_vac_neg=psy_DENV1_vac_neg,
              psy_DENV2_vac_neg=psy_DENV2_vac_neg,
              psy_DENV3_vac_neg=psy_DENV3_vac_neg,
              psy_DENV4_vac_neg=psy_DENV4_vac_neg,
              psy_DENV1_vac_pos=psy_DENV1_vac_pos,
              psy_DENV2_vac_pos=psy_DENV2_vac_pos,
              psy_DENV3_vac_pos=psy_DENV3_vac_pos,
              psy_DENV4_vac_pos=psy_DENV4_vac_pos))
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function estimates effectiveness for a specific vaccinated cohort in a specific simulation year by taking into account the boosting of efficacy by each natural infection
#Reminder to facilitate reading the code: the object "trace" has 5 dimensions:
#  1/ Nb of the vac cohort (cohorts of different age will be vaccinated and will experience breakthrough infections at different periods in time) --> nb levels = nb cohorts that are vaccinated (see below)
#  2/ Nb infections before vaccination (0 to 4) --> nb levels = 5
#  3/ Simulation year --> nb levels = nb years in the timeframe 
#  4/ Current nb of infections (0 to 4) --> nb levels = 5 
#     (distinction between the current nb of infections and nb of infections before vaccination is require to differentiate episodes post-vaccination; 
#      for example the 2nd natural infection will be the second episode post-vaccination in someone who had no infections before being vaccinated; in someone who had one infection before being vaccinated, it will be the first episode post-vaccination. Etc.)
#  5/ Nb of years since the last infection --> nb levels = nb years in the timeframe (equivalent of a tunnel state in a Markov trace)

ParmsVacEff_GetEffectiveness<-function(year,
                                       trace,
                                       vac_cohorts_nb,
                                       vac_cohorts_age_at_vac,     #Matrix indicating the age of vaccination (in years) for each vaccinated cohort
                                       vac_cohorts_year_of_vac,    #Matrix indicating the year of vaccination for each vaccinated cohort
                                       age_at_vac_loop,            #Age at vaccination (in years) for the cohort considered in this loop
                                       year_of_vac_loop,           #Year of vaccination for the cohort considered in this loop
                                       eff_asympt,
                                       eff_sympt_non_hosp,
                                       eff_sympt_hosp,
                                       boosting_asymptomatic_switch,
                                       boosting_symptomatic_switch){
  
  for(vac_cohort in 1:vac_cohorts_nb){                #Loop through all the cohorts included in the object "trace"
    age_at_vac<-vac_cohorts_age_at_vac[vac_cohort]    #Age this cohort was at the moment of vaccination (years)
    year_of_vac<-vac_cohorts_year_of_vac[vac_cohort]  #Determine the year when this cohort was (or will be) vaccinated
    
    if(age_at_vac==age_at_vac_loop && year_of_vac==year_of_vac_loop){     #If this cohort corresponds to the one of interest
      #-------------------------------------------------------------------
      #--- Initialize the objects ----------------------------------------
      #-------------------------------------------------------------------
      #Initialize the objects to record effectiveness, by serotype & serostatus
      #Here naming the objects "effect" to indicate that this is no longer efficacy, but real-life (field) effectiveness (as it depends on the boosting, i.e. transmission intensity)
      effect_asympt        <-matrix(0, nrow=4, ncol=2)
      effect_sympt_non_hosp<-matrix(0, nrow=4, ncol=2)
      effect_sympt_hosp    <-matrix(0, nrow=4, ncol=2)
      #-------------------------------------------------------------------
      #--- Estimate effectiveness for those seroneg. at vaccination ------
      #-------------------------------------------------------------------
      #Isolate the trace for seronegatives for the current cohort and currently simulated year
      nb_inf_before_vac<-0
      neg_trace_year<-trace[vac_cohort,(nb_inf_before_vac+1),year,(1:4),] #2D object, only those who had zero infections before vaccination (i.e. seronegative at vaccination; 1st level of 2nd dimension) and has 0 to 3 infections now (ignoring those with 4 infections because they are not at risk anymore; +1 to all indices because infection counts start with 0 and indexing starts with 1)
      if(sum(neg_trace_year) >0){
        neg_trace_year_prop<-neg_trace_year/sum(neg_trace_year)            #Express values in the trace as proportions of the entire population who were seronegative at vaccination; proportions are required for weighing the efficacy
      }else{
        neg_trace_year_prop<-neg_trace_year
      }
      
      #Loop through serotypes
      for(serotype in 1:4){                                               #Loop through all serotypes
        for(nb_inf_current in nb_inf_before_vac:3){                       #Loop through the number of current infections starting with the nb of infections before vac and until a maximum of 3 (people with 4 infections do not need any efficacy as they are no longer at risk of a new infection)
          
          nb_episode<-nb_inf_current-nb_inf_before_vac+1                  #Calculate the episode for which the individual is at risk; this defines efficacy
          
          #Estimate effectiveness for the current loop
          #Estimating effectiveness by multiplying the proportion of population with a specific nb of current infections (by tunnel year, i.e. by how long they've been at this state) by the efficacy curve for the episode for which they are at risk
          
          #For asymptomatic infections (if boosting is included, otherwise the object effect_asympt will remain filled with zeros, but it will not be used any further)
          if(boosting_asymptomatic_switch){
            #Estimate efficacy for the current "fraction" of the population (loop)
            effect_asympt_loop<-sum(neg_trace_year_prop[(nb_inf_current+1),]*eff_asympt[serotype,1,nb_episode,])    
            
            #Check that the values of effectiveness are comprised between 0 and 1
            if((effect_asympt_loop<0)||(effect_asympt_loop>1)){message("WARNING: Effectiveness not comprised between 0 and 1")}
            
            #Add it to the effectiveness calculated in previous loops (i.e. parts of the cohort with different nb of current infections), for this serotype
            effect_asympt[serotype,1]<-effect_asympt[serotype,1]+effect_asympt_loop
            
            #Zero out the counts for the next loop
            effect_asympt_loop<-0
          }
          
          #Similarly for symptomatic infections (if boosting is included)
          if(boosting_symptomatic_switch){
            effect_sympt_non_hosp_loop<-sum(neg_trace_year_prop[(nb_inf_current+1),]*eff_sympt_non_hosp[serotype,1,nb_episode,])       
            effect_sympt_hosp_loop    <-sum(neg_trace_year_prop[(nb_inf_current+1),]*eff_sympt_hosp    [serotype,1,nb_episode,])   
            
            if((effect_sympt_non_hosp_loop<0)||(effect_sympt_non_hosp_loop>1)){message("WARNING: Effectiveness not comprised between 0 and 1")}
            if((effect_sympt_hosp_loop    <0)||(effect_sympt_hosp_loop    >1)){message("WARNING: Effectiveness not comprised between 0 and 1")}
            
            effect_sympt_non_hosp[serotype,1]<-effect_sympt_non_hosp[serotype,1]+effect_sympt_non_hosp_loop
            effect_sympt_hosp    [serotype,1]<-effect_sympt_hosp    [serotype,1]+effect_sympt_hosp_loop
            
            effect_sympt_non_hosp_loop<-0
            effect_sympt_hosp_loop    <-0
          }
        }
      }
      #-------------------------------------------------------------------
      #--- Estimate effectiveness for those seropos. at vaccination ------
      #-------------------------------------------------------------------
      #Subjects who were seropositive at vaccination include several "strata" - those who had 1 infection before vaccination, 2 infections etc.
      #These strata have to be treated separately as the nb of of infections before vaccination (when compared to the current number of infection) will determine the nb of the episode post-vaccination, for which they are at risk and, hence, the efficacy to be applied
      
      #Isolate the trace for seropositives for the current cohort and currently simulated year
      pos_trace_year<-trace[vac_cohort,(2:4),year,(2:4),]                 #3D object for those who had 1 to 3 infections before vaccination and has 1 to 3 infections now (converting to 2:4 as indexing starts with 1 and not zero); ignoring those with 4 infections (before vaccination or now) as they are no longer at risk of a new infection
      if(sum(pos_trace_year)>0){
        pos_trace_year_prop<-pos_trace_year/sum(pos_trace_year)             #Express values in the trace as proportions of the entire population who were seropositive at vaccination (but at risk of a new infection, i.e. ignoring those who had 4 infections); proportions are required for weighing the efficacy
      }else{
        pos_trace_year_prop<-pos_trace_year
      }
      
      #Loop through serotypes
      for(serotype in 1:4){                                                 #Loop through all serotypes
        for(nb_inf_before_vac in 1:3){                                      #Loop through the "strata" (i.e. categories of vaccinees with different nb of infections before vaccination)
          for(nb_inf_current in nb_inf_before_vac:3){                       #Loop through the number of current infections starting with the nb of infections before vac and until a maximum of 3 (people with 4 infections do not need any efficacy as they are no longer at risk of a new infection)
            
            nb_episode<-nb_inf_current-nb_inf_before_vac+1                  #Calculate the episode for which the individual is at risk; this defines efficacy
            
            #Estimate effectiveness for the current loop
            #Estimating effectiveness by multiplying the proportion of population in a specific "stratum" with a specific nb of current infections (by tunnel year, i.e. by how long they've been at this state) by the efficacy curve for the episode for which they are at risk
            #Note: here in the objects eff_XXX taking index 2 for 2nd dimension (seropositives)
            
            #For asymptomatic infections (if boosting is included)
            if(boosting_asymptomatic_switch){
              
              #Estimate efficacy for the current "fraction" of the population (loop)
              effect_asympt_loop<-sum(pos_trace_year_prop[nb_inf_before_vac,nb_inf_current,]*eff_asympt[serotype,2,nb_episode,])   #Note: unlike in the section for seronegatives above, not adding 1 to the nb_inf_current, because the object pos_trace_year_prop only contains those who were seropositive (i.e. first level corresponds to history of 1 infection already)
              
              #Check that the values of effectiveness are comprised between 0 and 1 
              if((effect_asympt_loop<0)||(effect_asympt_loop>1)){message("WARNING: Effectiveness not comprised between 0 and 1")}
              
              #Add it to the effectiveness calculated in previous loops
              effect_asympt[serotype,2]<-effect_asympt[serotype,2]+effect_asympt_loop
              
              #Zero out the counts for the next loop
              effect_asympt_loop<-0
            }
            
            #For symptomatic infections (if boosting is included; similarly to the above)
            if(boosting_symptomatic_switch){
              effect_sympt_non_hosp_loop<-sum(pos_trace_year_prop[nb_inf_before_vac,nb_inf_current,]*eff_sympt_non_hosp[serotype,2,nb_episode,])
              effect_sympt_hosp_loop    <-sum(pos_trace_year_prop[nb_inf_before_vac,nb_inf_current,]*eff_sympt_hosp    [serotype,2,nb_episode,]) 
              
              if((effect_sympt_non_hosp_loop<0)||(effect_sympt_non_hosp_loop>1)){message("WARNING: Effectiveness not comprised between 0 and 1")}
              if((effect_sympt_hosp_loop    <0)||(effect_sympt_hosp_loop    >1)){message("WARNING: Effectiveness not comprised between 0 and 1")}
              
              effect_sympt_non_hosp[serotype,2]<-effect_sympt_non_hosp[serotype,2]+effect_sympt_non_hosp_loop
              effect_sympt_hosp    [serotype,2]<-effect_sympt_hosp    [serotype,2]+effect_sympt_hosp_loop
              
              effect_sympt_non_hosp_loop<-0
              effect_sympt_hosp_loop    <-0
            }
          }
        }
      }
      
      break                                            #Exit the loop to avoid looping through all the remaining cohorts (the cohort of interest is already found)
    }
    
  }
  
  return(list(effect_asympt=effect_asympt,
              effect_sympt_non_hosp=effect_sympt_non_hosp,
              effect_sympt_hosp=effect_sympt_hosp))
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function estimates average weighted proportion (probability) of each dengue outcome (asymptomatic, symptomatic non-hospitalized, symptomatic hospitalized)
#in a given cohort and for the current simulation year. It, thus, takes into account, what parts of the cohorts are at risk of 1st, 2nd, 3rd or 4th infection

ParmsVacEff_GetCohortDistr<-function(age_ind,       #Age (index) of the cohort in question
                                     yini_6D,       #Array with the initial state for the current simulation year (post-ageing, post-vaccination)
                                     parms_epi      #List with severity parameters (sub-list parms_epi_other of the list parms_mod)
){
  
  #-------------------------------------------------------------------
  #--- Initialize the objects ----------------------------------------
  #-------------------------------------------------------------------
  #Initialize the objects to record average weighted proportion of each outcome in this cohort, by serotype (as prevalence over time will differ and severity with each serotype may differ too),
  #and serostatus (as efficacy/effectiveness may differ)
  prop_asympt        <-matrix(0, nrow=4, ncol=2)         
  prop_sympt_non_hosp<-matrix(0, nrow=4, ncol=2)
  prop_sympt_hosp    <-matrix(0, nrow=4, ncol=2)
  
  #-------------------------------------------------------------------
  #--- Estimate the proportions for the cohort -----------------------
  #-------------------------------------------------------------------
  #Loop through each serotype
  for(serotype in 1:4){
    
    #Loop through each serostatus at vaccination
    for(serostatus in 1:2){
      #Estimate the population at risk of a new infection, for a given serotype and a given serostatus
      #To be considered as at-risk, the index for the serotype in question should be 1 (susceptible); the values for other serotypes do not matter 
      #(i.e. the subjects considered in the loop may or may not have history of other serotypes; it is of no importance as long as they are susceptible to the serotype in question)
      
      if     (serotype==1){pop_at_risk_total<-sum(yini_6D[age_ind,1, , , ,(serostatus+1)])}
      else if(serotype==2){pop_at_risk_total<-sum(yini_6D[age_ind, ,1, , ,(serostatus+1)])}
      else if(serotype==3){pop_at_risk_total<-sum(yini_6D[age_ind, , ,1, ,(serostatus+1)])}
      else                {pop_at_risk_total<-sum(yini_6D[age_ind, , , ,1,(serostatus+1)])}
      
      if(pop_at_risk_total>0){
        #Initialize the variable "weight" used to make sure that the entire population considered in the current loop is accounted for
        weight_total<-0
        
        #Loop through the possible combinations of statuses for the serotypes other than the one considered in the current loop
        for(k in 1:5){
          for(l in 1:5){
            for(m in 1:5){
              
              #Determine the nb of past infections
              values_klm<-c(k,l,m)                                      #A vector with the current values of indices k, l and m
              nb_past_inf<-sum((values_klm >= 2) & (values_klm <= 5))   #Nb of already experienced infections (0 to 3; by definition cannot have 4 past infection as then this fraction of the cohort is no longer at risk)
              
              #Determine the size of the compartment
              #Keeping the index 1 (susceptible) for the serotype in question, and using the current combination of indices k, l and m for the remaining three serotypes
              if     (serotype==1){pop_at_risk_loop<-sum(yini_6D[age_ind,1,k,l,m,(serostatus+1)])}
              else if(serotype==2){pop_at_risk_loop<-sum(yini_6D[age_ind,k,1,l,m,(serostatus+1)])}
              else if(serotype==3){pop_at_risk_loop<-sum(yini_6D[age_ind,k,l,1,m,(serostatus+1)])}
              else                {pop_at_risk_loop<-sum(yini_6D[age_ind,k,l,m,1,(serostatus+1)])}
              
              #Determing the compartment weight in the total population at risk
              weight_loop<-pop_at_risk_loop/pop_at_risk_total
              
              #Add the weight of the current compartment to the total weight (to be able to check later that the whole population at risk, for a given serotype & serostatus, was taken into account)
              weight_total<-weight_total+weight_loop
              #Get the input values for the probability of each outcome for the nb of infections and serotype considered in this loop
              prop_sympt_loop<-eval(parse(text=paste0('parms_epi$prop_sympt_unvac[',age_ind,',',serotype,',',min((nb_past_inf+1),3),']')))           #Symptomatic, given infection
              prop_asympt_loop<-1-prop_sympt_loop                                                                                                    #Asymptomatic, given infection
              prop_hosp_loop<-eval(parse(text=paste0('parms_epi$prop_hospit_sympt_unvac[',age_ind,',',serotype,',',min((nb_past_inf+1),3),']')))     #Hospitalized, given symptomatic infection
              prop_sympt_non_hosp_loop<-prop_sympt_loop*(1-prop_hosp_loop)                                                                           #Non-hospitalized, given infection    
              prop_sympt_hosp_loop<-prop_sympt_loop*prop_hosp_loop                                                                                   #Hospitalized, given infection
              
              #Check that the proportions add up to 1
              if(round((prop_asympt_loop+prop_sympt_non_hosp_loop+prop_sympt_hosp_loop),5)!=1){message("WARNING: Input severity proportions do not add up to 1")}
              
              #Add weighted values from the current loop to the values calculated in the previous loops
              prop_asympt        [serotype,serostatus]<-prop_asympt        [serotype,serostatus]+weight_loop*prop_asympt_loop
              prop_sympt_non_hosp[serotype,serostatus]<-prop_sympt_non_hosp[serotype,serostatus]+weight_loop*prop_sympt_non_hosp_loop
              prop_sympt_hosp    [serotype,serostatus]<-prop_sympt_hosp    [serotype,serostatus]+weight_loop*prop_sympt_hosp_loop
              
              #Reset the values for the next loop
              prop_sympt_loop         <-0
              prop_asympt_loop        <-0
              prop_hosp_loop          <-0
              prop_sympt_non_hosp_loop<-0
              prop_sympt_hosp_loop    <-0
            }
          }
        }
        
        #Check that all the little parts of the cohorts were accounted for
        if(round(weight_total,5)!=1){message("WARNING: Sum of weights does not add up to 1")}
        
        #Check that the proportions for the current serotype & serostatus add up to 1
        if(serostatus==1){
          if(round((prop_asympt[serotype,serostatus]+prop_sympt_non_hosp[serotype,serostatus]+prop_sympt_hosp[serotype,serostatus]),5)!=1){message("WARNING: Proportions for DENV-", serotype, " in seronegatives do not add up to 1")}
        }else{
          if(round((prop_asympt[serotype,serostatus]+prop_sympt_non_hosp[serotype,serostatus]+prop_sympt_hosp[serotype,serostatus]),5)!=1){message("WARNING: Proportions for DENV-", serotype, " in seropositives do not add up to 1")}
        }
      }
    }
  }
  
  return(list(prop_asympt=prop_asympt,
              prop_sympt_non_hosp=prop_sympt_non_hosp,
              prop_sympt_hosp=prop_sympt_hosp))
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function estimates, for a given cohort, the efficacy/effectiveness against overall dengue infection based on the efficacy against each dengue outcome 
#(asymptomatic, symptomatic non-hospitalized, symptomatic hospitalized) and their relative probabilities in the given cohort, for the current simulation year

ParmsVacEff_GetEffInf<-function(sympt_becomes_asympt,    #Flag whether symptomatic infections avoided with vaccination become asymptomatic (TRUE) or are completely eliminated (FALSE)
                                eff_asympt,              #Efficacy/effectiveness against asymptomatic dengue, for a given cohort (estimated in the function ParmsVacEff_GetThetaAndPsyYear); matrix 4x2 (serotype x serostatus at vac)
                                eff_sympt_non_hosp,      #...same for symptomatic non-hospitalized dengue
                                eff_sympt_hosp,          #... same for hospitalized dengue  
                                prop_asympt,             #Weighted average proportion (probability) of asymptomatic dengue in a given cohort; matrix 4x2 as well
                                prop_sympt_non_hosp,     #... same for symptomatic non-hospitalized dengue
                                prop_sympt_hosp          #... same for hospitalized dengue
){
  
  #------------
  #------------
  #Initialize the object for efficacy against infection --> matrix 4x2 (serotype x serostatus), for a given age
  eff_inf<-matrix(NA, nrow=4, ncol=2)
  
  #------------
  #------------
  #If it is assumed that the symptomatic infections avoided with vaccination become asymptomatic,
  #then first calculating the new share of asymptomatic and only then applying the efficacy against asymptomatic onto it
  
  if(sympt_becomes_asympt){
    
    #Calculate new share of asymptomatic
    prop_asympt_new<-(prop_asympt+                                #Proportion of infections that were asymptomatic to begin with,
                        prop_sympt_non_hosp*eff_sympt_non_hosp+     #... proportion of symptomatic non-hospitalized that will become asymptomatic with vaccination
                        prop_sympt_hosp*eff_sympt_hosp)             #... same with symptomatic hospitalized
    
    #Calculate efficacy against infection (equal to avoided asymptomatic)
    eff_inf<-prop_asympt_new*eff_asympt
    
  }
  
  #------------
  #------------
  #If it is assumed that the symptomatic infections avoided with vaccination are eliminated (disappear),
  #then the efficacy against infection is simply a weighted average of all type of efficacy
  
  if(!sympt_becomes_asympt){
    eff_inf<-(prop_asympt        *eff_asympt+                   #Asymptomatic infections that will disappear,
                prop_sympt_non_hosp*eff_sympt_non_hosp+           #... symptomatic non-hospitalized that will disappear, and
                prop_sympt_hosp    *eff_sympt_hosp)               #... symptomatic hospitalized that will disappear
  }
  
  return(eff_inf)
}


#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function calculates the proportion of different severity of infections in vaccinated individuals for the year that is about to be simulated
#These calculations take into account the vaccine efficacy/effectiveness by:
#- Type (asymptomatic, symptomatic non-hospitalized or symptomatic hospitalized dengue)
#- Serotype
#- Serostatus at vaccination
#- Individual's age
#The calculated proportions are then used to estimate the number of different outcomes (once the overall number of infections in protected individuals is simulated in the transmission module)
#In addition to that, the calculated proportion of symptomatic infections is used to calculate weighted average infectiousness to be used in the transmission module 
#(as it depends on whether the infection is symptomatic or asymptomatic; see function ParmsVacEff_Initialise, calculation of objects psy)

#The calculations use the objects with efficacy/effectiveness calculated for the current year (in the function ParmsVacEff_GetThetaAndPsyYear)
#The calculation is based on the following formula (Note: The described process is an example for one type of outcome only. It also does not reflect the fact that
#the probability of being symptomatic in an unvaccinated individual (prop_sympt_unvac) is different for primary, secondary and post-secondary infections)

#nb_inf_vac = nb_inf_unvac x (1-eff_inf)
#nb_sympt_vac = nb_sympt_unvac x (1-eff_sympt)       
#             = nb_inf_unvac x prop_sympt_unvac x (1-eff_sympt)

#                 nb_sympt_vac   nb_inf_unvac x prop_sympt_unvac x (1-eff_sympt)         1
#prop_sympt_vac = ------------ = ------------------------------------------------ = ----------- x prop_sympt_unvac x (1-eff_sympt) [here the terms nb_inf_unvac in numerator and denominator are cancelled out]
#                  nb_inf_vac            nb_inf_unvac x (1-eff_inf)                  1-eff_inf 


#This logic can be reproduced separately for symptomatic cases that are hospitalized or not hospitalized

ParmsVacEff_GetSevDistr<-function(parms_epi,       #Population & transmission parameters 
                                  parms_vac_eff,
                                  efficacies){     #List of objects with efficacy values
  
  with(c(parms_epi,
         efficacies),{
           
           #-----------
           #Get proportion of symptomatic infections (hospitalised and non-hospitalised) and severe infections in unvaccinated individuals, for each serotype
           #In the original inputs the proportion of hospitalized dengue infections is expressed as % of symptomatic (same for severe infections)
           #They will all now be expressed as % of infections (except for the probability of severe, see below)
           prop_sympt_hosp_unvac    <-array(NA, dim=c(N_age_groups,4,3))        #Dimensions are age, serotype & infection type
           prop_sympt_non_hosp_unvac<-array(NA, dim=c(N_age_groups,4,3))
           prop_severe_unvac        <-array(NA, dim=c(N_age_groups,4,3))
           
           prop_sympt_hosp_unvac<-prop_sympt_unvac*prop_hospit_sympt_unvac
           prop_sympt_non_hosp_unvac<-prop_sympt_unvac*(1-prop_hospit_sympt_unvac)
           prop_severe_unvac<-prop_severe_sympt_unvac                              #The final % of severe in the overall infections will depend on the new % of symptomatic (after the efficacy is applied)
           #                                                                       #Therefore, not multiplying by the proportion of symptomatic in unvaccinated. 
           #                                                                       #These values, thus, represent % severe out of symptomatic infection (not overall infections)
           
           #-----------
           #Initialize the objects with the proportions of infections of different severity
           #Dimensions of each object are age group, serotype, serostatus at vaccination, infection type (primary, secondary, post-secondary)
           #Note: Some of the elements of this arrrays are irrelevant. For example % of asymptomatic primary infections in someone who is seropositive at vaccination 
           #(being seropositive at vaccination, they will never have a primary infection as vaccinated). However, for the sake of simplicity these values are calculated, but they will not be used later on
           prop_asympt_vac        <-array(NA, dim=c(N_age_groups,4,2,3))   #% asymptomatic infections
           prop_sympt_non_hosp_vac<-array(NA, dim=c(N_age_groups,4,2,3))   #% non-hospitalized symptomatic infections
           prop_sympt_hosp_vac    <-array(NA, dim=c(N_age_groups,4,2,3))   #% hospitalized symptomatic infections
           prop_sympt_vac         <-array(NA, dim=c(N_age_groups,4,2,3))   #% total symptomatic infections (required for the calculation of average infectiousness)
           prop_severe_vac        <-array(NA, dim=c(N_age_groups,4,2,3))   #% severe infections
           
           #-----------
           #Calculate the proportion of each type of infection infection in effectively vaccinated individuals by looping through each element of the array
           #Note: This section of the code may use some optimization (loop is not the most efficient way), but it's only executed once for each simulation
           for(age in 1:N_age_groups){
             for(serotype in 1:4){
               for(serostatus in 1:2){
                 for(type in 1:3){
                   
                   #if vaccination acts as a silent infection then vaccination "shift" the underlying probability of being symptomatic/hospitalized/severe
                   shift=0
                   if(parms_vac_eff$vac_as_silent_inf){
                     if(type<3){shift=1}
                   }
                   #Calculate % of each outcome
                   #Symptomatic hospitalized & non-hospitalized
                   if((1-eff_inf[age,serotype,serostatus])>0){
                     prop_sympt_non_hosp_vac[age,serotype,serostatus,type]<-1/(1-eff_inf[age,serotype,serostatus])*prop_sympt_non_hosp_unvac[age,serotype,type+shift]*(1-eff_sympt_non_hosp[age,serotype,serostatus])
                     prop_sympt_hosp_vac    [age,serotype,serostatus,type]<-1/(1-eff_inf[age,serotype,serostatus])*prop_sympt_hosp_unvac    [age,serotype,type+shift]*(1-eff_sympt_hosp    [age,serotype,serostatus])
                     #Total symptomatic & asymptomatic
                     prop_sympt_vac [age,serotype,serostatus,type]<-prop_sympt_hosp_vac[age,serotype,serostatus,type]+prop_sympt_non_hosp_vac[age,serotype,serostatus,type]
                     prop_asympt_vac[age,serotype,serostatus,type]<-1-prop_sympt_vac [age,serotype,serostatus,type]
                     #Severe (calculating last as it depends on the new proportion of symptomatic)
                     prop_severe_vac[age,serotype,serostatus,type]<-prop_sympt_vac[age,serotype,serostatus,type]*prop_severe_unvac[age,serotype,type+shift]       #Note: Here prop_severe_unvac is the % of severe infections OUT OF SYMPTOMATIC infections
                   }else{
                     prop_sympt_non_hosp_vac[age,serotype,serostatus,type]<-0
                     prop_sympt_hosp_vac    [age,serotype,serostatus,type]<-0
                     prop_sympt_vac [age,serotype,serostatus,type]<-0
                     prop_asympt_vac[age,serotype,serostatus,type]<-0
                     prop_severe_vac[age,serotype,serostatus,type]<-0
                   }
                   
                 }
               }
             }
           }
           
           return(list(prop_sympt_non_hosp_vac=prop_sympt_non_hosp_vac,
                       prop_sympt_hosp_vac=prop_sympt_hosp_vac,
                       prop_severe_vac=prop_severe_vac,
                       prop_sympt_vac=prop_sympt_vac,
                       prop_asympt_vac=prop_asympt_vac))
         })
}

