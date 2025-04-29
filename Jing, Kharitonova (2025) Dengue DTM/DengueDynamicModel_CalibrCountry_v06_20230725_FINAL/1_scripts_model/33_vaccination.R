##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contain a function that re-calculates the size of the system for hosts after discrete vaccination
###########################################################################################################################################################################################################

Vaccination<-function(list_stateH,             #List of vectors with the size of host compartments (one vector per vaccination status)
                      parms_mod,               #List of all model parameters
                      parms_gen,               #List of simulation characteristics
                      year){                   #The number of the last completed simulation year
  
  with(c(list_stateH, parms_mod$parms_vac_strat, parms_mod$parms_epi_other, parms_mod$parms_vac_other, parms_gen),{
    
    #----------------------------------------------------------------------------------------------------------------------------
    #------------ Initialise the vectors ----------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------  
    #The vectors below will contain the differential nb. of hosts in each compartment (for each vaccination status)
    #This number can be positive (if more hosts enter the compartment than leave the compartment) or negative (if more hosts leave the compartment than enter the compartment)
    #One element of each vector is calculated at each loop (see below)
    
    dH_unvac<-rep(0,NstatesH)   		#Unvaccinated --> always used
    dH_vac_neg<-rep(0,NstatesH)    	#Vaccinated as seronegative --> may remain unused if a model w/o vaccination
    dH_vac_pos<-rep(0,NstatesH)   	#Vaccinated as seropositive --> may remain unused if a model w/o vaccination
    
    #Nb.of hosts who are newly vaccinated (or screened for vaccination)
    new_vac<-rep(0,(2*N_age_groups))                #Vector with the cumulative numbers for the year (number of screened, by age; number of vaccinated, by age)
    
    #----------------------------------------------------------------------------------------------------------------------------
    #------------ Calculate population movements due to vaccination -------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------  
    next_year=year+1                                                #Simulation year that is about to start (determines the coverage rates to be applied)
    
    #------ Loop through all age groups --------------------------------------
    for(i in 1:N_age_groups){        
      
      #------ Get the coverage rates -------------------------------------------
      cov_rate<-coverage_by_age_year[i,next_year]     #The coverage rate represents either the vaccination rate (if everybody is vaccinated) or the screening rate (if only seropositive individuals are vaccinated)
                                                      #Taking the coverage rate for age group i and the simulation year that is about to start
      
      if(cov_rate!=0){                                #If the coverage rate for the current age group is not zero, then continue
        
        #------ Loop through all compartments ------------------------------------
        for(j in 1:5){
          for(k in 1:5){
            for(l in 1:5){
              for(m in 1:5){
                
                pos<-((i-1)+N_age_groups*((j-1)+5*((k-1)+5*((l-1)+5*(m-1)))))+1    #Determine the position of the current compartment in the vector with all host compartments
                
                #------ Determine whether the compartment is to be vaccinated ---------------------

                #Determine the nb. of previous infections in this compartment
                values_jklm<-c(j,k,l,m)                                            #With the serostatus-based vaccination, the number of previous infections determine whether or not the compartment is to be vaccinated
                nb_prev_inf<-sum(values_jklm >= 2)                                 #Counting the nb of indices >=2 (i.e. not susceptible)
                
                #If universal vaccination (both seronegative & seropositive are vaccinated)
                if(vac_strat_scope==1){                                       
                  
                  #If serostatus-based vaccination
                }else{                                                           
                  
                  #Record the new cumulative number of screened hosts (among unvaccinated)
                  new_vac[i]<-new_vac[i]+cov_rate*stateH_unvac[pos]                #Nb. of hosts screened in this compartment is added to the nb. of hosts screened in other compartments (in the age group i)             

                }
                
                #------ Vaccinate (if needed) -----------------------------------------------------

                if(nb_prev_inf==0){                                           # if seronegative at vaccination
                  if(vac_strat_scope==1){                                     # when individuals are not tested to be seropositive..
                    multiply_par<-1                                           # ...then vaccinate all
                  }else{                                                      # otherwise...
                    multiply_par<-(1-test_specificity)                        # ...vaccine only false positives (identified as positive although are negative)
                  }
                  #If they are seronegative (zero previous infections), they move to an identical compartment for hosts vaccinated as seronegative 
                  #Record the new cumulative number of vaccinated hosts (among unvaccinated)
                  new_vac[i+N_age_groups]<-new_vac[i+N_age_groups]+cov_rate*stateH_unvac[pos]*multiply_par     #Nb. of hosts vaccinated in this compartment is added to the nb. of hosts vaccinated in other compartments
                  #The calculation is the same as for the number of screened hosts 
                  #Calculate change in the compartment size
                  dH_unvac[pos]<-dH_unvac[pos]-cov_rate*stateH_unvac[pos]*multiply_par                         #Newly vaccinated hosts leave the compartment jlkm for unvaccinated hosts
                  dH_vac_neg[pos]<-dH_vac_neg[pos]+cov_rate*stateH_unvac[pos]*multiply_par 
                  
                }else{  #If they are seropositive, they move to an identical compartment for hosts vaccinated as seropositive 
                  if(vac_strat_scope==1){                                     # when individuals are not tested to be seropositive..
                    multiply_par<-1                                           # ...then vaccinate all
                  }else{                                                      # otherwise...
                    multiply_par<-test_sensitivity                            # ...vaccine only true positives (identified as positive when really positive)
                  }
                  #Record the new cumulative number of vaccinated hosts (among unvaccinated)
                  new_vac[i+N_age_groups]<-new_vac[i+N_age_groups]+cov_rate*stateH_unvac[pos]*multiply_par     #Nb. of hosts vaccinated in this compartment is added to the nb. of hosts vaccinated in other compartments
                  #The calculation is the same as for the number of screened hosts 
                  #Calculate change in the compartment size
                  dH_unvac[pos]<-dH_unvac[pos]-cov_rate*stateH_unvac[pos]*multiply_par                         #Newly vaccinated hosts leave the compartment jlkm for unvaccinated hosts
                  dH_vac_pos[pos]<-dH_vac_pos[pos]+cov_rate*stateH_unvac[pos]*multiply_par

                }
              }
            }
          }      
        }        
      }          
    }
    
    #----------------------------------------------------------------------------------------------------------------------------
    #------------ Apply population movements due to vaccination to the previous system for the hosts ----------------------------
    #---------------------------------------------------------------------------------------------------------------------------- 
    stateH_unvac <- stateH_unvac + dH_unvac
    stateH_vac_neg <- stateH_vac_neg + dH_vac_neg
    stateH_vac_pos <- stateH_vac_pos + dH_vac_pos
    
    rm(dH_unvac, dH_vac_neg, dH_vac_pos); gc()
    
    #----------------------------------------------------------------------------------------------------------------------------
    #------------ Aggregate & return the new system for the hosts ---------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------  
    if(!vac_switch){yH_vac<-c(stateH_unvac)}
    if(vac_switch) {yH_vac<-c(stateH_unvac, 
                              stateH_vac_neg,
                              stateH_vac_pos)}
    
    rm(stateH_unvac, stateH_vac_neg, stateH_vac_pos); gc()
    
    return(list(yH_vac=yH_vac,
                new_vac=new_vac))
  })         
}








