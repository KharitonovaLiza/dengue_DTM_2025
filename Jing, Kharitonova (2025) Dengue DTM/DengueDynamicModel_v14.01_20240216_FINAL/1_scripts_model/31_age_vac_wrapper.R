##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script launches the discrete vaccination & ageing
###########################################################################################################################################################################################################

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This wrapper function applies discrete ageing and vaccination to the results of the dynamic simulation. Ageing and vaccination are assumed to occur between the simulation years.
#They are, thus, applied to the system estimated on the last day of the current simulation year.The resulting system is then used as the initial state for the next simulation year.
#Ageing is applied before vaccination (to ensure that the hosts change the age groups before being vaccinated)


AgeingVaccination <-function(parms_mod,       #List of all model parameters 
                             parms_gen,       #List of simulation characteristics
                             year,            #The number of the last completed simulation year
                             yH,              #A vector with the size of host compartments as of the last day of the last completed simulation year
                             new_vac,         #A matrix with the nb. of newly vaccinated/screened individuals for each simulation year
                             ageing_flag){    #A boolean indicating whether ageing has to be applied

 with(c(parms_mod, parms_gen),{
  
  yH<-yH
  new_vac<-new_vac
  
  #------- If model without vaccination -----------------------       #Ageing is applied at the end of each simulation year & before the simulation starts (as the initial states represents Dec 31). . NOTE: In the RCT calibration ageing is now always applied (see comments in the script launch.R and simulation.R)
  if(!vac_switch){
    if(ageing_flag){
      list_stateH<-SplitStateH(yH=yH, parms_gen=parms_gen)            #Split the vector with the host compartments into sub-vectors (one per vaccination status)
      rm(yH); gc()
      
      yH<-Ageing(list_stateH=list_stateH,                             #Apply ageing 
                 parms_mod=parms_mod,
                 parms_gen=parms_gen)
      rm(list_stateH); gc()
    }
    
  #------- If model with vaccination --------------------------       #Vaccination is applied before each simulation year (including the first year). NOTE: In the RCT calibration ageing is now always applied (see comments in the script launch.R and simulation.R)
  }else if(vac_switch){
    #Apply ageing (if necessary)
    if(ageing_flag){
      list_stateH<-SplitStateH(yH=yH, parms_gen=parms_gen)            #Split the vector with the host compartments         
      rm(yH); gc()
      
      yH<-Ageing(list_stateH=list_stateH,                             #Apply ageing
                 parms_mod=parms_mod,
                 parms_gen=parms_gen)
      rm(list_stateH); gc()      
    }
    
    #Apply vaccination (if at least one simulation year left)
    if((year+1)<=timeframe){
      
      list_stateH<-SplitStateH(yH=yH, parms_gen=parms_gen)           #Split the vector with the host compartments      
      rm(yH); gc()
      
      post_vac<-Vaccination(list_stateH=list_stateH,                 #Apply vaccination
                            parms_mod=parms_mod, 
                            parms_gen=parms_gen,
                            year)
      
      yH=post_vac$yH_vac                                             #Get new system for the hosts
      new_vac[(year+1),]<-post_vac$new_vac                           #Update the matrix with the nb of newly vaccinated/screened
    }
  }
   
  #------- Return results -------------------------------------
  return(list(yH=yH,                                 
              new_vac=new_vac))                      
  })
}

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function splits the system for the hosts into sub-vectors (by vaccination status)

SplitStateH<-function(yH,                   #System for the hosts (a vector with the size of each compartment)
                      parms_gen){           #List of simulation characteristics
  
  with(parms_gen, {
    #---------- Initialise all vectors -------------- 
    stateH_unvac<-rep(0,NstatesH)                  #Unvaccinated --> always used
    stateH_vac_neg<-rep(0,NstatesH)                #Vaccinated as seronegative --> may remain unused if a model w/o vaccination
    stateH_vac_pos<-rep(0,NstatesH)                #Vaccinated as seropositive --> may remain unused if a model w/o vaccination
    
    #---------- Invaccinated individuals ------------
    stateH_unvac<-yH[1:NstatesH]                         
    
    #---------- Vaccinated individuals --------------                         #If a model with vaccination, define the vectors for vaccinated hosts
    if(vac_switch){              		                                          #Otherwise they only contain zeros
      stateH_vac_neg <- yH[(NstatesH*1+1):(NstatesH*2)]
      stateH_vac_pos <- yH[(NstatesH*2+1):(NstatesH*3)]
    }         
    
    #---------- Return results ----------------------
    return(list(stateH_unvac=stateH_unvac,                                       
                stateH_vac_neg=stateH_vac_neg,   
                stateH_vac_pos=stateH_vac_pos))    
  })
}
