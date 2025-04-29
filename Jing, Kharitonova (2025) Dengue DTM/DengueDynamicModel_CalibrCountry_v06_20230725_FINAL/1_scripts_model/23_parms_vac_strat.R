##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contains functions that import, select, initalise and load all the vaccination strategy parameters
###########################################################################################################################################################################################################

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function selects the vaccination strategy inputs that are required to run the model out of all inputs imported from an Excel file
ParmsVacStrat_Select<-function(parms_vac_strat_imported){   #Inputs imported from the Excel file VacStrat_XXXXX (output of the function Parms_Import)
  
  with(parms_vac_strat_imported,{ 
    
    parms_vac_strat_selected<-list(
      test_before_vac=as.logical(test_before_vac),   							     #Flag for serostatus-based vaccination (only seropositive individuals are vaccinated). The coverage rates will be considered as screening rate --> Boolean
      test_sensitivity=as.numeric(test_sensitivity),   	 	             #Probability to identify a true positive --> Numeric
      test_specificity=as.numeric(test_specificity),	                 #Probability to identify a true negative --> Numeric
      
      age_routine=as.integer(age_routine),												     #Age for routine vaccination, years --> Integer
      coverage_routine=as.numeric(coverage_routine),   								 #Coverage rate for the routine vaccination (or screening if enabled) during the year 1, year 2, year 3 and years 4+ --> Numeric vector (4 elements)
      
      catchup_switch=as.logical(catchup_switch),										   #Flag for the inclusion of catch-up campaign --> Boolean
      age_catchup_low=as.integer(age_catchup_low),										 #Lower bound for the age included in the catch-up campaign, years --> Integer
      age_catchup_high=as.integer(age_catchup_high),									 #Upper bound for the age included in the catch-up campaign, years --> Integer
      coverage_catchup=as.numeric(coverage_catchup)   								 #Coverage rate for catch-up vaccination  --> Numeric
    )
    
    return(parms_vac_strat_selected)
  })
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function calculates the vaccination coverage rate for each age group and at each year of the model
ParmsVacStrat_GetCov<-function(parms_epi,                             #Population & transmission parameters,
                               parms_vac_strat_selected){              #Pre-selected vaccination parameters
                               
  
 with(c(parms_epi$parms_epi_other, 
         parms_vac_strat_selected),{
           
           #------ Get coverage rates by single year of age for each year ------------------------
           
           #Initialise objects
           coverage_by_age_year_101<-matrix(0, nrow=101, ncol=100)                               #Matrix with 101 rows (one per a single year of age) and one column per simulation year
           
           #Routine vaccination
           coverage_by_age_year_101[age_routine+1,1]=coverage_routine[1]                         #Year 1 of the routine vaccination (+1 to row number, because the age starts with 0 and indexing starts with 1)
           coverage_by_age_year_101[age_routine+1,2]=coverage_routine[2]                         #Year 2 of the routine vaccination
           coverage_by_age_year_101[age_routine+1,3]=coverage_routine[3]                         #Year 3 of the routine vaccination
           for(i in 4:100){
             coverage_by_age_year_101[age_routine+1,i]=coverage_routine[4]                       #Years 4+ of the routine vaccination
           }
           
           #Catch-up vaccination (if included)
           if(catchup_switch){
             coverage_by_age_year_101[(age_catchup_low+1):(age_catchup_high+1),1]=coverage_catchup      #Catch-up vaccination coverage (only at the year of vaccination introduction)
           }
           
           #------ Re-calculate for model age groups (if not single-year groups) -----------------
           coverage_by_age_year<-matrix(0,nrow=N_age_groups, ncol=100)                            #Matrix with as many rows as age groups used in the model
           
           if(N_age_groups!=101){
             #Loop through the model age groups
             for(age_gr in 1:N_age_groups){
               if(ai[age_gr]==1){                                                                                                           #If the age group only includes one year...
                 coverage_sum<-coverage_by_age_year_101[(age_group_lower_bound[age_gr]+1),]                                                 #... then copy the corresponding value from the matrix with single years
               }else{
                 coverage_sum<-apply(coverage_by_age_year_101[(age_group_lower_bound[age_gr]+1):(age_group_upper_bound[age_gr]+1),],2,sum)  #Otherwise, calculate the sum of the coverage rates for each single year included in the age group
               }
               coverage_by_age_year[age_gr,]<-coverage_sum*ai[age_gr]                                                                       #Weigh the sum of the coverage rates by the weight of each year in the age group 
             }             
           }else{
             coverage_by_age_year<-coverage_by_age_year_101
           }
           

           
           return(list(coverage_by_age_year_101=coverage_by_age_year_101,
                       coverage_by_age_year=coverage_by_age_year)) 
         })
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function initializes all the parameters for vaccination strategy
ParmsVacStrat_Initialise<- function(parms_epi,                	#Population & transmission parameters
                                    parms_vac_strat_selected){  #Selected and pre-formatted vaccination inputs (output of the function ParmsVacStrat_Select)
  
  with(c(parms_epi, parms_vac_strat_selected),{
    
    #--- Determine vaccination scope -------------------------------------
    #Determine which individuals are to be vaccinated (in the appropriate age groups):
    # vac_strat_scope=1: Everyone 
    # vac_strat_scope=2: Only the seropositive individuals (all of them)
           
    if(test_before_vac){          #If vaccination of seropositive only...
      vac_strat_scope=2           #... --> vac_strat_scope = 2
    }else{
      vac_strat_scope=1           #Otherwise, vaccination of everyone (seropositive or seronegative)     
    }
           
    #--- Determine coverage for each year --------------------------------       
    coverage<-ParmsVacStrat_GetCov(parms_epi=parms_epi, 
                                   parms_vac_strat_selected=parms_vac_strat_selected)
    
    #--- Aggregate vaccinatino strategy parameters -----------------------   
    parms_vac_strat=list(vac_strat_scope=vac_strat_scope,                         			     #Vaccination scope --> integer
                         test_sensitivity=test_sensitivity,   	 	                           #Probability to identify a true positive
                         test_specificity=test_specificity,	                                 #Probability to identify a true negative
                         age_routine=parms_vac_strat_selected$age_routine,                   #Age of routine vaccination --> integer
                         coverage_routine=parms_vac_strat_selected$coverage_routine,         #Coverage of routine vaccination for years 1, 2, 3 and 4+ --> numeric vector with four elements
                         catchup_switch=parms_vac_strat_selected$catchup_switch,             #Presence of catch-up campaign --> boolean
                         age_catchup_low=parms_vac_strat_selected$age_catchup_low,           #Lower bound of the age group vaccinated during catch-up --> integer
                         age_catchup_high=parms_vac_strat_selected$age_catchup_high,         #Upper bound of the age group vaccinated during catch-up --> integer
                         coverage_catchup=parms_vac_strat_selected$coverage_catchup,         #Coverage of catch-up vaccination --> numeric
                         coverage_by_age_year=coverage$coverage_by_age_year)                 #Coverage rate for each age group and each year --> numeric matrix N_age_groups x nb years
    
    #------ Return lists with parameters -------------------------------------------------
    return(list(parms_vac_strat=parms_vac_strat))
    
  })
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Wrapper functions that imports, formats, calculates and outputs vaccination strategy parameters
ParmsVacStrat_Load<-function(path_inputs_vac_strat,        #Path to the folder containing the Excel file with the vaccination inputs
                             set_vac_strat,                #Name of the Excel file containing the inputs (without extension) 
                             parms_epi){                    #Population & transmission parameters
  
  
  parms_vac_strat_imported<-Parms_Import(path=path_inputs_vac_strat, filename=set_vac_strat) 
  parms_vac_strat_selected<-ParmsVacStrat_Select(parms_vac_strat_imported=parms_vac_strat_imported)
  parms_vac_strat_final<-ParmsVacStrat_Initialise(parms_epi=parms_epi, parms_vac_strat_selected=parms_vac_strat_selected)
  return(parms_vac_strat_final)
}

