##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contain a function that re-calculates the size of the system for hosts after discrete ageing
###########################################################################################################################################################################################################

Ageing<-function(list_stateH,              #List of vectors with the size of host compartments (one vector per vaccination status)
                 parms_mod,                #List of model parameters
                 parms_gen){               #List of simulation characteristics
  
 with(c(list_stateH, parms_mod, parms_mod$parms_epi_other, parms_mod$parms_epi_dyn_model, parms_gen),{
  
    #----------------------------------------------------------------------------------------------------------------------------
    #------------ Initialise the vectors ----------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------  
    #The vectors below will contain the differential nb. of hosts in each compartment (for each vaccination status)
    #This number can be positive (if more hosts enter the compartment than leave the compartment) or negative (if more hosts leave the compartment than enter the compartment)
    #One element of each vector is calculated at each loop (see below)
    
    dH_unvac<-rep(0,NstatesH)   	#Unvaccinated --> always used
    dH_vac_neg<-rep(0,NstatesH)   #Vaccinated as seronegative --> may remain unused if a model w/o vaccination
    dH_vac_pos<-rep(0,NstatesH)   #Vaccinated as seropositive --> may remain unused if a model w/o vaccination
    
    #----------------------------------------------------------------------------------------------------------------------------
    #------------ Calculate population movements due to ageing ------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------  
    
    #Loop through all the compartments
    for(i in 1:N_age_groups){
      for(j in 1:5){
        for(k in 1:5){
          for(l in 1:5){
            for(m in 1:5){
              
              pos<-((i-1)+N_age_groups*((j-1)+5*((k-1)+5*((l-1)+5*(m-1)))))+1        #Determine the position of the current compartment in the vector (Note: the name of the variable pos might be misleading; not to confuse with "pos" for "seropositive")
              
              #--------- First age group -----------------------------------------   #First age group is affected by the births and ageing out (no ageing in)
              if(i==1){
                
                #Unvaccinated individuals (only the fully susceptible compartment) --> births 
                if((j*k*l*m)==1){dH_unvac[pos]<-dH_unvac[pos]+ai[i]*NHi[i]}          #Fully susceptible compartment --> all indices are equal to 1 (j*k*l*m=1) 
                                                                                     #The number of births is equal to the number of hosts aged 0 years (fraction ai[1] of the age group size NHi[1])              
                #Unvaccinated individuals (all) --> ageing out
                dH_unvac[pos]<-dH_unvac[pos]-ai[i]*stateH_unvac[pos]                        
                
                #Vaccinated individuals --> ageing out
                if(vac_switch){														
                  dH_vac_neg[pos]<- dH_vac_neg[pos] - ai[i]*stateH_vac_neg[pos]
                  dH_vac_pos[pos]<- dH_vac_pos[pos] - ai[i]*stateH_vac_pos[pos]
                }                                                      
                
                #--------- Other age groups ----------------------------------------  #Other age groups are affected by the ageing in and ageing out
                                                                                      #Ageing in includes hosts who age AND survive in the process
              }else{
                
                pos_prev_age<-((i-2)+N_age_groups*((j-1)+5*((k-1)+5*((l-1)+5*(m-1)))))+1   #Determine the position of the current compartment for the previous age group (to calculate ageing in) 
                
                #Unvaccinated individuals --> ageing in & out
                dH_unvac[pos]<-dH_unvac[pos]+(ai[i-1]*stateH_unvac[pos_prev_age]*(1-muHi[i-1]))-(ai[i]*stateH_unvac[pos])   
                
                #Vaccinated individuals --> ageing in & out                                                                  
                if(vac_switch){
                  dH_vac_neg[pos]<-dH_vac_neg[pos] + (ai[i-1]*stateH_vac_neg[pos_prev_age]*(1-muHi[i-1])) - (ai[i]*stateH_vac_neg[pos])
                  dH_vac_pos[pos]<-dH_vac_pos[pos] + (ai[i-1]*stateH_vac_pos[pos_prev_age]*(1-muHi[i-1])) - (ai[i]*stateH_vac_pos[pos])
                }              
              }  
            }
          }
        }
      }
    }  
    
    #----------------------------------------------------------------------------------------------------------------------------
    #------------ Apply population movements due to ageing to the previous system for the hosts ---------------------------------
    #----------------------------------------------------------------------------------------------------------------------------  
    stateH_unvac <- stateH_unvac + dH_unvac
    stateH_vac_neg <- stateH_vac_neg + dH_vac_neg
    stateH_vac_pos <- stateH_vac_pos + dH_vac_pos
    
    rm(dH_unvac, dH_vac_neg, dH_vac_pos); gc()
    
    #----------------------------------------------------------------------------------------------------------------------------
    #------------ Aggregate & return the new system for the hosts ---------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------  
    if(!vac_switch){yH_ageing<-c(stateH_unvac)}
    if(vac_switch) {yH_ageing<-c(stateH_unvac, 
                                 stateH_vac_neg,
                                 stateH_vac_pos)}
    
    rm(stateH_unvac, stateH_vac_neg, stateH_vac_pos); gc()
    
    return(yH_ageing)
  })
}
