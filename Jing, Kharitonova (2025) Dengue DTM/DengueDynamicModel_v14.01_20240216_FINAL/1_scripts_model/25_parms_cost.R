##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contains functions that import and load all the cost inputs
###########################################################################################################################################################################################################


####################################################
#Vaccination cost
ParmsCostVac_Load<-function(path_inputs_cost_vac,
                            set_cost_vac){   
  
  parms_cost_vac_imported<-Parms_Import(path=path_inputs_cost_vac, filename=set_cost_vac)
  parms_cost_vac_final<-list(  
    nb_doses=as.numeric(parms_cost_vac_imported$nb_doses),                                          #Number of vaccine doses required to vaccinated one person --> Numeric  
    cost_vaccine_per_dose=as.numeric(parms_cost_vac_imported$cost_vaccine_per_dose),   							#Cost of one vaccine dose --> Numeric
    cost_admin_per_dose=as.numeric(parms_cost_vac_imported$cost_admin_per_dose),   					        #Administration cost of one vaccine dose  --> Numeric
    cost_extra_per_dose=as.numeric(parms_cost_vac_imported$cost_extra_per_dose),                    #Extra cost per dispensed dose (vaccine storage, awareness campaigns, logistics etc) --> Numeric
    cost_screening_per_capita=as.numeric(parms_cost_vac_imported$cost_screening_per_capita)   			#Cost of screening per screened person (only if vaccination by serostatus is enabled) --> Numeric
  )
  return(list(parms_cost_vac=parms_cost_vac_final))
}

####################################################
#Cost of vector control
ParmsCostVec_Load<-function(path_inputs_cost_vec,         #Path to the folder containing the Excel file with the vector control inputs
                            set_cost_vec){                #Name of the Excel file containing the inputs (without extension) 
  
  parms_cost_vec_imported<-Parms_Import(path=path_inputs_cost_vec, filename=set_cost_vec) 
  parms_cost_vec_final<-list(vec_costs_annual=as.numeric(parms_cost_vec_imported$vec_costs_annual))  			#Annual cost of all vector control policies --> Numeric
  return(list(parms_cost_vec=parms_cost_vec_final))
}

####################################################
#Cost of treatment
ParmsCostTtt_Load<-function(path_inputs_cost_ttt,
                            set_cost_ttt){
  
  #Import the inputs
  parms_cost_ttt_imported<-Parms_Import(path=path_inputs_cost_ttt, filename=set_cost_ttt) 
  
  #Extend the proportion of cases with persistent (long-term) dengue to single-year age groups
  long_term_proportion_101y<-rep(NA,101)
  for(j in 1:as.integer(parms_cost_ttt_imported$N_age_groups_eco)){                 #Loop through every age group used to define costs   
    age_low <-as.integer(parms_cost_ttt_imported$age_group_lower_bound_eco)[j]+1    #Get the first and the last year included in the composite age group (+1 as the age starts with zero and indexing starts with 1)
    age_high<-as.integer(parms_cost_ttt_imported$age_group_upper_bound_eco)[j]+1
    for(i in age_low:age_high){ 
      long_term_proportion_101y[i]<-as.numeric(parms_cost_ttt_imported$long_term_proportion)[j]
    }
  }
  
  #Aggregate the final list
  parms_cost_ttt_final = list(
    cost_by_age=as.logical(parms_cost_ttt_imported$cost_by_age),										          #True if infection cost are splitted by age --> Boolean
    N_age_groups_eco=as.integer(parms_cost_ttt_imported$N_age_groups_eco),                    #Number of age groups used to define costs --> integer
    age_group_lower_bound_eco=as.integer(parms_cost_ttt_imported$age_group_lower_bound_eco),  #Lower bounds of the age groups used to define costs --> vector of integeres with as many elements as age groups
    age_group_upper_bound_eco=as.integer(parms_cost_ttt_imported$age_group_upper_bound_eco),  #Upper bounds of the age groups used to define costs --> vector of integeres with as many elements as age groups
    
    costs_non_hosp_med=as.numeric(parms_cost_ttt_imported$costs_non_hosp_med),                #Direct medical costs of non-hospitalised cases --> numeric vector with as many elements as age groups
    costs_non_hosp_non_med=as.numeric(parms_cost_ttt_imported$costs_non_hosp_non_med),        #Direct non-medical costs of non-hospitalised cases --> numeric vector with as many elements as age groups
    income_lost_non_hosp=as.numeric(parms_cost_ttt_imported$income_lost_non_hosp),            #Income lost due to non-hospitalised cases --> numeric vector with as many elements as age groups
    wd_lost_non_hosp=as.numeric(parms_cost_ttt_imported$wd_lost_non_hosp),								    #Work days lost due to mild non-hospitalised cases --> numeric vector with as many elements as age groups
    sd_lost_non_hosp=as.numeric(parms_cost_ttt_imported$sd_lost_non_hosp),								    #School days lost due to mild non-hospitalised cases --> numeric vector with as many elements as age groups
    
    costs_hosp_mild_med=as.numeric(parms_cost_ttt_imported$costs_hosp_mild_med),              #Direct medical cost of mild hospitalised cases --> numeric vector with as many elements as age groups
    costs_hosp_mild_non_med=as.numeric(parms_cost_ttt_imported$costs_hosp_mild_non_med),      #Direct non-medical cost of mild hospitalised cases --> numeric vector with as many elements as age groups
    income_lost_hosp_mild=as.numeric(parms_cost_ttt_imported$income_lost_hosp_mild),          #Income lost due to mild hospitalised cases --> numeric vector with as many elements as age groups
    wd_lost_hosp_mild=as.numeric(parms_cost_ttt_imported$wd_lost_hosp_mild),							    #Work days lost due to mild hospitalised cases --> numeric vector with as many elements as age groups
    sd_lost_hosp_mild=as.numeric(parms_cost_ttt_imported$sd_lost_hosp_mild),							    #School days lost due to mild hospitalised cases --> numeric vector with as many elements as age groups
    
    costs_hosp_severe_med=as.numeric(parms_cost_ttt_imported$costs_hosp_severe_med),          #Direct medical cost of severe hospitalised cases --> numeric vector with as many elements as age groups
    costs_hosp_severe_non_med=as.numeric(parms_cost_ttt_imported$costs_hosp_severe_non_med),  #Direct non-medical cost of severe hospitalised cases --> numeric vector with as many elements as age groups
    income_lost_hosp_severe=as.numeric(parms_cost_ttt_imported$income_lost_hosp_severe),      #Income lost due to severe hospitalised cases --> numeric vector with as many elements as age groups
    
    wd_lost_hosp_severe=as.numeric(parms_cost_ttt_imported$wd_lost_hosp_severe),						  #Work days lost due to severe hospitalised cases --> numeric vector with as many elements as age groups
    sd_lost_hosp_severe=as.numeric(parms_cost_ttt_imported$sd_lost_hosp_severe),						  #School days lost due to severe hospitalised cases --> numeric vector with as many elements as age groups
    
    long_term_proportion=as.numeric(parms_cost_ttt_imported$long_term_proportion),
    long_term_proportion_101y=long_term_proportion_101y,                                      #% cases resulting in persistent dengue --> numeric vector with 101 elements
    long_term_cost=as.numeric(parms_cost_ttt_imported$long_term_cost),                        #Cost of persistent dengue, per month --> numeric vector with as many elements as age groups
    
    daily_wage=as.numeric(parms_cost_ttt_imported$daily_wage),											                         #Average wage --> scalar
    mortality_related_productivity=as.logical(parms_cost_ttt_imported$mortality_related_productivity),       #Inclusion of productivity loss cost due to mortality  --> Boolean
    workdays_annual=as.numeric(parms_cost_ttt_imported$workdays_annual),									               #Nb of workdays in year --> scalar
    age_work_start=as.numeric(parms_cost_ttt_imported$age_work_start),											                 #Age at work start --> integer
    age_retirement=as.numeric(parms_cost_ttt_imported$age_retirement),											                 #Age of retirement --> integer
    emp_rate=as.numeric(parms_cost_ttt_imported$emp_rate),										                               #Employment rate --> scalar
    absenteeism_cost=as.numeric(parms_cost_ttt_imported$absenteeism_cost))										               #Average wage --> scalar
  
  return(list(parms_cost_ttt=parms_cost_ttt_final))
}


