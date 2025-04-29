##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  SCRIPT NAME: 3_ageing.R (for calibration!)
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA,   Olivier CRISTEAU , Aurélien JAMOTTE
###########################################################################################################################################################################################################

Process_Age<-function(stateH_NV,                #Vector with the sizes of host compartments (unvaccinated hosts only)
                      parms_epi){               #Population & transmission inputs
  
  #----------------------------------------------------------------------------------------------------------------------------
  #------------ Get required parameters ---------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------------------  
  muHi=parms_epi$parms_epi_dyn_model$muHi
  ai=parms_epi$parms_epi_other$ai
  NHi=parms_epi$parms_epi_other$NHi
  
  #----------------------------------------------------------------------------------------------------------------------------
  #------------ Initialise the vector with the incremental nb. of individuals -------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------------------  
  #This number can be positive (if more hosts enter the compartment than leave the compartment) or negative (if more hosts leave the compartment than enter the compartment)
  #One element in each vector is calculated at each loop (see below)
  
  dH_NV <-rep(0,NstatesH)   #Incremental nb. of individuals in each compartment --> one element of this vector is calculated at each loop (see below)

  
  #----------------------------------------------------------------------------------------------------------------------------
  #------------ Calculate population movements due to ageing ------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------------------------------------- 
  
  #Loop through all the compartments
  for(i in 1:N_age_groups){
    for(j in 1:5){
      for(k in 1:5){
        for(l in 1:5){
          for(m in 1:5){
            
            vpos<-((i-1)+N_age_groups*((j-1)+5*((k-1)+5*((l-1)+5*(m-1)))))+1        #Determine the position of the current compartment in the vector with all host compartments
            
            #--------- First age group -----------------------------------------    #First age group is concerned by the births and ageing out (no ageing in)
            if(i==1){
              if((j*k*l*m)==1){dH_NV[vpos] <- dH_NV[vpos] + ai[i]*NHi[i]}           #Births (fully susceptible individuals)
              dH_NV[vpos] <- dH_NV[vpos] - ai[i]*stateH_NV[vpos]                    #Ageing out
              
              #--------- Other age groups ----------------------------------------  #Other age groups are concerned by the ageing in and ageing out
                                                                                    #Ageing in includes individuals who age and survive in the process
            } else {
              pos_prev_agegp<-((i-2)+N_age_groups*((j-1)+5*((k-1)+5*((l-1)+5*(m-1)))))+1                                #Determine the position of the current compartment for the previous age group (to calculate ageing in) 
              dH_NV[vpos] <- dH_NV[vpos] + (ai[i-1]*stateH_NV[pos_prev_agegp]*(1-muHi[i-1]))- (ai[i]*stateH_NV[vpos])   #Ageing in (those who survive) and ageing out
            }  
          }
        }
      }
    }
  }  
  
  #Apply population movements due to waning to the previous system for the hosts
  stateH_NV<-stateH_NV + dH_NV
  return(yH_ageing=stateH_NV)
}