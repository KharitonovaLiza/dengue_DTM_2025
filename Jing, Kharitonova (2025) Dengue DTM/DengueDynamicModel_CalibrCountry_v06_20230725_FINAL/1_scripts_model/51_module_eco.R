##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contains functions that calculate economic outcomes (costs, nb of cases)
###########################################################################################################################################################################################################

ResultsEpi_PrepForEco<-function(parms_mod,      #Epi and cost inputs
                                parms_gen,      #Scenario general inputs
                                year_start){    #Simulation start year
  
  #---------------------------------------------------------------------------------------
  #---- Import the results of the transmission module ------------------------------------
  #---------------------------------------------------------------------------------------
  load(file=file.path(path_results, paste(parms_gen$scen_ref, "_y", year_start, "_Epi.Rdata", sep="")))  #Load the object with epi results (loaded under name results_epi)
  epi_results_all<-results_epi
  
  #---------------------------------------------------------------------------------------
  #---- Aggregate results ----------------------------------------------------------------
  #---------------------------------------------------------------------------------------
  
  #Host population size, by age
  pop_host<-matrix(NA, nrow=parms_gen$timeframe, ncol=parms_gen$N_age_groups)           #Matrix with one row per year and one column per age group
  pop_host<-apply(epi_results_all$pop_host[,,], c(1,2), sum)                            #In the results epi pop_host is a 3D array with the population size by year, age and vacination status. Sum across all vaccination statuses

  #Nb. of infections
  incidence<-array(NA,dim=c(parms_gen$timeframe, parms_gen$N_age_groups, 6))            #3D array where each element is the number of infections by year, age and severity
  incidence[,,]<-apply(epi_results_all$nb_inf[,,,,,, drop=F], c(1,2,6), sum, drop=F)    #Nb_inf is a 6D array with the nb of infection by year, age, serotype, type, vaccination status & severity --> summing it across all serotypes, types and vac statuses for each year, age and severity
 
  #Newly vaccinated
  new_vac<-epi_results_all$new_vac[,(parms_gen$N_age_groups+1):(parms_gen$N_age_groups*2)]      #Number of newly vaccinated (second half of the columns in the object new_vac)
  
  #Newly screened for vaccination
  new_screen<-epi_results_all$new_vac[,1:(parms_gen$N_age_groups)]    
  
  return(list(pop_host=pop_host,
              incidence=incidence,
              new_vac=new_vac,
              new_screen=new_screen))
}

#===================================================================================================================================
#===================================================================================================================================
#Function that calculates the costs in each category 
ResultsEco_CalculateCosts<-function(parms_mod,
                                    parms_gen,
                                    results_epi){
  
 with(c(parms_mod, parms_gen, results_epi),{
  
   #-----------------------------------------------------------
   #-------- Cost of treatment --------------------------------
   #-----------------------------------------------------------
   
   #Initialise matrices for each type of costs
   labels_age<-apply(cbind(rep("ag", N_age_groups),1:N_age_groups), 1, paste, collapse="") 
   labels_year=c(1:timeframe)
   
   cost_med_non_hospit       <-array(NA,dim=c(timeframe,N_age_groups),dimnames=list(labels_year,labels_age))   #Direct medical cost (including the cost of persistent dengue)
   cost_med_hospit_mild      <-array(NA,dim=c(timeframe,N_age_groups),dimnames=list(labels_year,labels_age))
   cost_med_hospit_severe    <-array(NA,dim=c(timeframe,N_age_groups),dimnames=list(labels_year,labels_age))
   
   cost_non_med_non_hospit   <-array(NA,dim=c(timeframe,N_age_groups),dimnames=list(labels_year,labels_age))   #Direct non-medical cost
   cost_non_med_hospit_mild  <-array(NA,dim=c(timeframe,N_age_groups),dimnames=list(labels_year,labels_age))
   cost_non_med_hospit_severe<-array(NA,dim=c(timeframe,N_age_groups),dimnames=list(labels_year,labels_age))
   
   income_lost_non_hospit    <-array(NA,dim=c(timeframe,N_age_groups),dimnames=list(labels_year,labels_age))   #Income lost
   income_lost_hospit_mild   <-array(NA,dim=c(timeframe,N_age_groups),dimnames=list(labels_year,labels_age))
   income_lost_hospit_severe <- array(NA,dim=c(timeframe,N_age_groups),dimnames=list(labels_year,labels_age))
   
   absenteeism_non_hospit    <-array(NA,dim=c(timeframe,N_age_groups),dimnames=list(labels_year,labels_age))   #School absenteeism
   absenteeism_hospit_mild   <-array(NA,dim=c(timeframe,N_age_groups),dimnames=list(labels_year,labels_age))
   absenteeism_hospit_severe <- array(NA,dim=c(timeframe,N_age_groups),dimnames=list(labels_year,labels_age))
   
   #Estimate costs
   for(j in 1:parms_cost_ttt$N_age_groups_eco){                 #Loop through every age group used to define costs   
     
     age_low<-parms_cost_ttt$age_group_lower_bound_eco[j]+1     #Get the first and the last year included in the composite age group (+1 as the age starts with zero and indexing starts with 1)
     age_high<-parms_cost_ttt$age_group_upper_bound_eco[j]+1
     
     for(i in age_low:age_high){                                #Loop through every year included in the age group
       
       #-------- Direct medical costs (incl. persistent dengue) ---
       cost_med_non_hospit   [,i]<-(incidence[,i,2]+incidence[,i,4])*(parms_cost_ttt$costs_non_hosp_med   [j]+parms_cost_ttt$long_term_proportion[j]*parms_cost_ttt$long_term_cost[j]*parms_QoL$dur_long_term)
       cost_med_hospit_mild  [,i]<- incidence[,i,3]                 *(parms_cost_ttt$costs_hosp_mild_med  [j]+parms_cost_ttt$long_term_proportion[j]*parms_cost_ttt$long_term_cost[j]*parms_QoL$dur_long_term)
       cost_med_hospit_severe[,i]<-((incidence[,i,5]+incidence[,i,6])*parms_cost_ttt$costs_hosp_severe_med[j]+                                                                                                  #Direct medical cost of treatment (applicable for both fatal and non-fatal cases, i.e. indices 5 and 6)
                                     incidence[,i,5]                 *                                        parms_cost_ttt$long_term_proportion[j]*parms_cost_ttt$long_term_cost[j]*parms_QoL$dur_long_term)  #Cost of persistent dengue (applicable for non-fatal cases only, i.e. index 5)
       
       
       #-------- Direct non-medical costs -------------------------
       cost_non_med_non_hospit   [,i]<-(incidence[,i,2]+incidence[,i,4])*parms_cost_ttt$costs_non_hosp_non_med   [j]
       cost_non_med_hospit_mild  [,i]<- incidence[,i,3]                 *parms_cost_ttt$costs_hosp_mild_non_med  [j]
       cost_non_med_hospit_severe[,i]<-(incidence[,i,5]+incidence[,i,6])*parms_cost_ttt$costs_hosp_severe_non_med[j]
       
       #-------- Income losses ------------------------------------
       income_lost_non_hospit   [,i]<-(incidence[,i,2]+incidence[,i,4])*parms_cost_ttt$income_lost_non_hosp   [j]
       income_lost_hospit_mild  [,i]<- incidence[,i,3]                 *parms_cost_ttt$income_lost_hosp_mild  [j]
       income_lost_hospit_severe[,i]<-(incidence[,i,5]+incidence[,i,6])*parms_cost_ttt$income_lost_hosp_severe[j]
       
       #-------- School absenteeism cost --------------------------
       absenteeism_non_hospit   [,i]<-(incidence[,i,2]+incidence[,i,4])*parms_cost_ttt$sd_lost_non_hosp   [j]*parms_cost_ttt$absenteeism_cost
       absenteeism_hospit_mild  [,i]<- incidence[,i,3]                 *parms_cost_ttt$sd_lost_hosp_mild  [j]*parms_cost_ttt$absenteeism_cost
       absenteeism_hospit_severe[,i]<-(incidence[,i,5]+incidence[,i,6])*parms_cost_ttt$sd_lost_hosp_severe[j]*parms_cost_ttt$absenteeism_cost   
     }
   }
   
   #-----------------------------------------------------------
   #-------- Cost of vaccination & screening ------------------
   #-----------------------------------------------------------
   cost_vac_screening<-array(0,dim=c(timeframe,N_age_groups))           #Initialising objects for vaccination-related costs
   cost_vac_vaccine<-array(0,dim=c(timeframe,N_age_groups))             #Originally all values are equal to zero (and they remain equal to zero if there's no vaccination)
   cost_vac_other<-array(0,dim=c(timeframe,N_age_groups))
   
   if(vac_switch){
     #Cost of screening
     cost_vac_screening[,]<-new_screen[,]*parms_cost_vac$cost_screening_per_capita                                  
     
     #Cost of vaccination (vaccine and administration)
     cost_vac_vaccine[,]<-new_vac[,]*parms_cost_vac$nb_doses*(parms_cost_vac$cost_vaccine_per_dose+parms_cost_vac$cost_admin_per_dose)   
     
     #Misc. cost (per dispensed dose)
     cost_vac_other[,]<-new_vac[,]*parms_cost_vac$nb_doses*parms_cost_vac$cost_extra_per_dose
   }
   
   #-----------------------------------------------------------
   #-------- Vector control cost ------------------------------
   #-----------------------------------------------------------
   if(!vec_switch){                                                              #If vector control is not included, then the cost is zero for each year
     cost_vector_control<-rep(0, timeframe)
   }else{
     cost_vector_control<-rep(parms_cost_vec$vec_costs_annual, timeframe)        #If vector control is included, then annual cost of vector control (same at each year)
   }
   
   
   #-----------------------------------------------------------
   #-------- Return final objects with costs ------------------
   #-----------------------------------------------------------
   return(list(cost_med_non_hospit=cost_med_non_hospit,
               cost_med_hospit_mild=cost_med_hospit_mild,
               cost_med_hospit_severe=cost_med_hospit_severe,
               cost_non_med_non_hospit=cost_non_med_non_hospit,
               cost_non_med_hospit_mild=cost_non_med_hospit_mild,
               cost_non_med_hospit_severe=cost_non_med_hospit_severe,
               income_lost_non_hospit=income_lost_non_hospit,
               income_lost_hospit_mild=income_lost_hospit_mild,
               income_lost_hospit_severe=income_lost_hospit_severe,
               absenteeism_non_hospit=absenteeism_non_hospit,
               absenteeism_hospit_mild=absenteeism_hospit_mild,
               absenteeism_hospit_severe=absenteeism_hospit_severe,
               cost_vac_screening=cost_vac_screening,
               cost_vac_vaccine=cost_vac_vaccine,
               cost_vac_other=cost_vac_other,
               cost_vector_control=cost_vector_control))
 })
}


#===================================================================================================================================
#===================================================================================================================================
#Function that calculates the number of cases (infections) for each year and each severity

ResultsEco_CalculateNbCases<-function(results_epi){
  with(results_epi,{
    return(list(cases_asympt                        =incidence[,,1],
                cases_sympt                         =incidence[,,2]+incidence[,,3]+incidence[,,4]+incidence[,,5]+incidence[,,6],
                cases_sympt_mild                    =incidence[,,2]+incidence[,,3],
                cases_sympt_mild_non_hospit         =incidence[,,2],
                cases_sympt_mild_hospit             =incidence[,,3],
                cases_sympt_severe                  =incidence[,,4]+incidence[,,5]+incidence[,,6],
                cases_sympt_severe_non_hospit       =incidence[,,4],
                cases_sympt_severe_hospit           =incidence[,,5]+incidence[,,6],
                cases_hospit_severe_hospit_non_death=incidence[,,5],
                cases_hospit_severe_hospit_death    =incidence[,,6],
                cases_sympt_hospit                  =incidence[,,3]+incidence[,,5]+incidence[,,6], 
                cases_all                           =incidence[,,1]+incidence[,,2]+incidence[,,3]+incidence[,,4]+incidence[,,5]+incidence[,,6])) 
  })
}


#===================================================================================================================================
#===================================================================================================================================
#Wrapper function calculating all the economic outcomes
ResultsEco_CalculateAll<-function(parms_mod,
                                  parms_gen,
                                  results_epi){
  
  costs<-ResultsEco_CalculateCosts(parms_mod=parms_mod, parms_gen=parms_gen, results_epi=results_epi)  #Calculate costs
  nb_cases<-ResultsEco_CalculateNbCases(results_epi=results_epi)                                       #Calculate nb of cases (infections)
  
  return(list(costs=costs, cases=nb_cases))
}


##########################################################################################################################################################################################################
# Function that prepare the output of economic module
ResultsEco_PrepareOutput<-function(parms_mod,
                                   parms_gen,
                                   results_epi,
                                   results_eco){
  
  with(c(parms_gen, parms_mod$parms_epi_other, results_epi, results_eco),{
    
    labels_year=c(1:timeframe)
    
    #--- Cost of treatment & vaccination and nb cases ----------------
    labels_age<-apply(cbind(rep("ag", N_age_groups),1:N_age_groups), 1, paste, collapse="") 
    labels_eco<-c("pop_host",
                  "cost_med_non_hospit","cost_med_hospit_mild","cost_med_hospit_severe",
                  "cost_non_med_non_hospit","cost_non_med_hospit_mild","cost_non_med_hospit_severe",
                  "income_lost_non_hospit","income_lost_hospit_mild","income_lost_hospit_severe",
                  "absenteeism_non_hospit","absenteeism_hospit_mild","absenteeism_hospit_severe",
                  "cost_vac_screening","cost_vac_vaccine","cost_vac_other",
                  "cases_asympt","cases_sympt","cases_sympt_mild","cases_sympt_mild_non_hospit",
                  "cases_sympt_mild_hospit","cases_sympt_severe","cases_sympt_severe_non_hospit",
                  "cases_sympt_severe_hospit","cases_hospit_severe_hospit_non_death", 
                  "cases_hospit_severe_hospit_death","cases_sympt_hospit","cases_all")
    costs_cases<-array(NA, dim=c(timeframe, N_age_groups, 28),dimnames = list(labels_year,labels_age,labels_eco))
    costs_cases[,,1 ]<-pop_host
    costs_cases[,,2 ]<-costs$cost_med_non_hospit
    costs_cases[,,3 ]<-costs$cost_med_hospit_mild
    costs_cases[,,4 ]<-costs$cost_med_hospit_severe
    costs_cases[,,5 ]<-costs$cost_non_med_non_hospit
    costs_cases[,,6 ]<-costs$cost_non_med_hospit_mild
    costs_cases[,,7 ]<-costs$cost_non_med_hospit_severe
    costs_cases[,,8 ]<-costs$income_lost_non_hospit
    costs_cases[,,9 ]<-costs$income_lost_hospit_mild
    costs_cases[,,10]<-costs$income_lost_hospit_severe
    costs_cases[,,11]<-costs$absenteeism_non_hospit
    costs_cases[,,12]<-costs$absenteeism_hospit_mild
    costs_cases[,,13]<-costs$absenteeism_hospit_severe
    costs_cases[,,14]<-costs$cost_vac_screening
    costs_cases[,,15]<-costs$cost_vac_vaccine
    costs_cases[,,16]<-costs$cost_vac_other
    costs_cases[,,17]<-cases$cases_asympt
    costs_cases[,,18]<-cases$cases_sympt
    costs_cases[,,19]<-cases$cases_sympt_mild
    costs_cases[,,20]<-cases$cases_sympt_mild_non_hospit
    costs_cases[,,21]<-cases$cases_sympt_mild_hospit
    costs_cases[,,22]<-cases$cases_sympt_severe
    costs_cases[,,23]<-cases$cases_sympt_severe_non_hospit
    costs_cases[,,24]<-cases$cases_sympt_severe_hospit
    costs_cases[,,25]<-cases$cases_hospit_severe_hospit_non_death
    costs_cases[,,26]<-cases$cases_hospit_severe_hospit_death
    costs_cases[,,27]<-cases$cases_sympt_hospit
    costs_cases[,,28]<-cases$cases_all

    return(list(costs_cases=costs_cases,                     #Costs of treatment & vaccination; nb of cases of each severity
                vector_control=costs$cost_vector_control))   #Cost of vector control (separated because it does not have age dimension)  
    
  })
  
}    



