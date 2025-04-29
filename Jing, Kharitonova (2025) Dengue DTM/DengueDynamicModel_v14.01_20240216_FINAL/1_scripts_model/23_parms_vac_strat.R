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
    
    cov_matrix<-matrix(NA, ncol=101, nrow=30)                     #Initialize a matrix with coverage rates for each year of age (columns) and each simulation year (rows)
    for(age in 0:100){                                            #Loop through all ages
      cov<-as.numeric(eval(parse(text=paste0("cov_age",age))))    #Load values with coverage for specific age (numeric)
      cov<-c(cov,rep(0,max(30-length(cov),0)))                    #Create a vector with coverage for specific age (appendix zeroes as needed)
      cov_matrix[,age+1]<-cov                                     #Input into the final coverage matrix (+1 as age starts with 0 and column counts start with 1)
    }
    
    parms_vac_strat_selected<-list(
      test_before_vac=as.logical(test_before_vac),   							     #Flag for serotesting before vaccination (only those who test positive are vaccinated). The coverage rates will be considered as screening rate --> Boolean
      test_sensitivity=as.numeric(test_sensitivity),   	 	             #Sensitivity (probability to identify a true positive) --> Numeric
      test_specificity=as.numeric(test_specificity),	                 #Specificify (probability to identify a true negative) --> Numeric
      cov_matrix=cov_matrix                                            #Coverage rate for each age and timeframe year --> Numeric
    )
    
    return(parms_vac_strat_selected)
  })
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function aggregates the vaccination coverage rate for each age and each year of the model
#It also creates a matrix indicating, for each year of age and each simulation year, whether this cohort is vaccinated or not

ParmsVacStrat_GetCov<-function(parms_epi,                       #Population & transmission parameters,
                               parms_vac_strat_selected,        #Pre-selected vaccination parameters
                               timeframe){                      #Simulation timeframe
                               
  
 with(c(parms_epi$parms_epi_other, 
         parms_vac_strat_selected),{
           
           #------ Get coverage rates by year of age & simulation year ----------------------------------------------------------------
           coverage_by_age_year<-matrix(0, nrow=101, ncol=timeframe)           #Matrix with 101 rows (one per a single year of age) and one column per simulation year
           for(year in 1:timeframe){
             coverage_by_age_year[,year]<-t(cov_matrix[year,])
           }

           #------ Prepare additional objects describing proportions of vaccinated ----------------------------------------------------
           vac_cohorts_cov        <-array(0,dim=c(101,timeframe,timeframe))    #Array with the proportion of a cohort of specific age that is vaccinated at a specific simulation year and year since vaccination
           vac_cohorts            <-array(0,dim=c(101,timeframe))              #Matrix with the proportion of a cohort of specific age that is vaccinated at a specific simulation year
           unvac_cohorts          <-array(1,dim=c(101,timeframe))              #Matrix with the proportion of a cohort of specific age that is unvaccinated at a specific simulation year
           vac_cohorts_age_at_vac <-array(0,dim=c(101,timeframe))              #Age at first vaccination for a specific cohort
           vac_cohorts_year_of_vac<-array(0,dim=c(101,timeframe))              #Year of first vaccination for specific cohort
           
           for(year in 1:timeframe){
             
             #First save the information for newly vaccinated at this year
             if(year==1){
               vac_cohorts    [,year  ]<-coverage_by_age_year  [,year]
               vac_cohorts_cov[,year,1]<-coverage_by_age_year  [,year]
               unvac_cohorts  [,year  ]<-1-coverage_by_age_year[,year]
             }else{
               vac_cohorts    [,year  ]<-vac_cohorts  [,year]+unvac_cohorts[,year]*coverage_by_age_year[,year]
               vac_cohorts_cov[,year,1]<-unvac_cohorts[,year]*coverage_by_age_year[,year]
               unvac_cohorts  [,year  ]<-unvac_cohorts[,year]*(1-coverage_by_age_year[,year])
             }
             
             #Then update the cells for the following years based on this year's vaccination
             for(age_ind in 1:100){
               if(year<timeframe){
                 vac_cohorts[age_ind+1,year+1]<-vac_cohorts[age_ind+1,year+1]+vac_cohorts[age_ind,year]
                 for(year1 in 1:(timeframe-year)){
                   if(age_ind+year1<=101){vac_cohorts_cov[age_ind+year1,year+year1,year1+1]<-vac_cohorts_cov[age_ind,year,1]}
                 }
                 unvac_cohorts[age_ind+1,year+1]<-1-vac_cohorts[age_ind+1,year+1]
               }
             }
           }
           
           #Create matrices with information at which age and year a specific cohort was first vaccinated
           trac_vac_cohorts_nb<-0               #Total nb of cohorts that will be vaccinated during the simulation
           trac_vac_cohorts_age_at_vac<-c()     #Age at first vaccination, for each vaccinated cohort
           trac_vac_cohorts_year_of_vac<-c()    #Simulation year of first vaccination, for each vaccinated cohort
           
           for(ind in 1:101){                                              #Loop through ages (0 to 100, i.e. age indices 1 to 101)
             for(year in 1:timeframe){                                     #Loop through simulation years
               if(vac_cohorts[ind,year]>0){                                #If the cohort was vaccinated (proportion of vaccinated is greater than zero)
                 if(ind==1 || year==1){                                    #If it is the first year of timeframe or the first year of the cohort's life then the cohort could not be vaccinated previously
                   trac_vac_cohorts_nb<-trac_vac_cohorts_nb+1                            #Increase the counter of vaccinate cohorts
                   trac_vac_cohorts_age_at_vac <-c(trac_vac_cohorts_age_at_vac,ind-1)    #Append the age of vaccination (0 years as ind==1) to the vector
                   trac_vac_cohorts_year_of_vac<-c(trac_vac_cohorts_year_of_vac,year)    #Append the year of vaccination (1st year) to the vector
                   
                   for(i in 1:(min(timeframe-year+1,101-ind+1))){                  #Save information of age at vac and year of vac for all cells in the matrix dedicated for this cohort
                     vac_cohorts_age_at_vac [ind+i-1,year+i-1]<-ind-1
                     vac_cohorts_year_of_vac[ind+i-1,year+i-1]<-year
                   }
                   
                 }else{
                   if(vac_cohorts[ind-1,year-1]==0){                                      #If the cohort was not vaccinated the year before then this is the first year of vaccination (we know that the cohort WAS vaccinated as the value of vac_cohorts is >0)
                     trac_vac_cohorts_nb<-trac_vac_cohorts_nb+1                           #Increase the counter
                     trac_vac_cohorts_age_at_vac <-c(trac_vac_cohorts_age_at_vac,ind-1)   #Append the age of vaccination
                     trac_vac_cohorts_year_of_vac<-c(trac_vac_cohorts_year_of_vac,year)   #Append the year of vaccination
                     
                     for(i in 1:(min(timeframe-year+1,101-ind+1))){                #Set information for all cells in matrices from this year to the end of timeframe
                       vac_cohorts_age_at_vac [ind+i-1,year+i-1]<-ind-1
                       vac_cohorts_year_of_vac[ind+i-1,year+i-1]<-year
                     }
                     
                     for(i in 1:(min(ind,year))){                                  #Set information for all cells in matrices before this year
                       vac_cohorts_age_at_vac [ind-i+1,year-i+1]<-ind-1
                       vac_cohorts_year_of_vac[ind-i+1,year-i+1]<-year
                     }
                   }
                 }
                }
              }
            }

           #Identify all vaccinated cohorts at each year of vaccination campaign
           #A flag for each cohort (total of 101 cohorts aged 0 to 100 years) at each year of vaccination campaign:
           #=0: Unvaccinated
           #=1: Vaccinated 
           vac_cohorts_ind<-ifelse(vac_cohorts>0,1,0) 

           
           return(list(coverage_by_age_year=coverage_by_age_year,
                       vac_cohorts=vac_cohorts,
                       vac_cohorts_age_at_vac=vac_cohorts_age_at_vac,
                       vac_cohorts_year_of_vac=vac_cohorts_year_of_vac,
                       vac_cohorts_ind=vac_cohorts_ind,
                       vac_cohorts_cov=vac_cohorts_cov,
                       trac_vac_cohorts_nb=trac_vac_cohorts_nb,
                       trac_vac_cohorts_age_at_vac=trac_vac_cohorts_age_at_vac,
                       trac_vac_cohorts_year_of_vac=trac_vac_cohorts_year_of_vac)) 
         })
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function initializes all the parameters for vaccination strategy
ParmsVacStrat_Initialise<- function(parms_epi,                 #Population & transmission parameters
                                    parms_vac_strat_selected,  #Selected and pre-formatted vaccination inputs (output of the function ParmsVacStrat_Select)
                                    timeframe){                #Simulation timeframe  
  
  with(c(parms_epi, parms_vac_strat_selected),{
    
    #--- Determine vaccination scope -------------------------------------
    #Determine which individuals are to be vaccinated (in the appropriate age groups):
    # vac_strat_scope=1: Everyone 
    # vac_strat_scope=2: Only those who tested seropositive (true or false)
           
    if(test_before_vac){          #If vaccination of seropositive only...
      vac_strat_scope=2           #... --> vac_strat_scope = 2
    }else{
      vac_strat_scope=1           #Otherwise, vaccination of everyone (seropositive or seronegative)     
    }
           
    #--- Determine coverage for each year --------------------------------       
    coverage<-ParmsVacStrat_GetCov(parms_epi=parms_epi, 
                                   parms_vac_strat_selected=parms_vac_strat_selected,
                                   timeframe=timeframe)
    
    #--- Aggregate vaccination strategy parameters -----------------------   
    parms_vac_strat=list(vac_strat_scope=vac_strat_scope,                         			     #Vaccination scope --> integer
                         test_sensitivity=test_sensitivity,   	 	                           #Probability to identify a true positive
                         test_specificity=test_specificity,	                                 #Probability to identify a true negative
                         coverage_by_age_year=coverage$coverage_by_age_year,                 #Coverage rate for each age group and each year --> numeric matrix N_age_groups x nb years
                         vac_cohorts_ind=coverage$vac_cohorts_ind,
                         vac_cohorts=coverage$vac_cohorts,
                         vac_cohorts_cov=coverage$vac_cohorts_cov,
                         vac_cohorts_age_at_vac=coverage$vac_cohorts_age_at_vac,
                         vac_cohorts_year_of_vac=coverage$vac_cohorts_year_of_vac,
                         trac_vac_cohorts_nb=coverage$trac_vac_cohorts_nb,
                         trac_vac_cohorts_age_at_vac=coverage$trac_vac_cohorts_age_at_vac,
                         trac_vac_cohorts_year_of_vac=coverage$trac_vac_cohorts_year_of_vac)
    
    #------ Return lists with parameters -------------------------------------------------
    return(list(parms_vac_strat=parms_vac_strat))
    
  })
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Wrapper functions that imports, formats, calculates and outputs vaccination strategy parameters
ParmsVacStrat_Load<-function(path_inputs_vac_strat,        #Path to the folder containing the Excel file with the vaccination inputs
                             set_vac_strat,                #Name of the Excel file containing the inputs (without extension) 
                             parms_epi,                    #Population & transmission parameters
                             timeframe){
  
  parms_vac_strat_imported<-Parms_Import(path=path_inputs_vac_strat, filename=set_vac_strat) 
  parms_vac_strat_selected<-ParmsVacStrat_Select(parms_vac_strat_imported=parms_vac_strat_imported)
  parms_vac_strat_final<-ParmsVacStrat_Initialise(parms_epi=parms_epi, parms_vac_strat_selected=parms_vac_strat_selected,timeframe=timeframe)
  return(parms_vac_strat_final)
}

