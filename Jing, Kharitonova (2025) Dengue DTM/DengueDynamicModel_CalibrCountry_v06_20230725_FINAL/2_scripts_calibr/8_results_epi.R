##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  SCRIPT NAME: 8_results_epi.R (for calibration!)
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA,   Olivier CRISTEAU , Aurélien JAMOTTE
###########################################################################################################################################################################################################

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Function that splits the ODE solver output into categories (incidence, population size, seroprevalence etc.)

Simulation_SplitResults<-function(out,                   #ODE solver output with or without information on the size of host compartments (returned by the function Simulation_Run as the object res_daily)
                                  parms,                 #List with all the population, transmission & intervention inputs
                                  hosts_included,        #Boolean indicating whether the object out (ODE solver output) contains information on the size of host compartments
                                  keep_first){           #Boolean indicating whether the first row of the object out should be kept (the first row represents the initial state, i.e. the day before the simulation started)
  
  with(c(parms, parms$parms_dyn_model,parms$other_parms), {
    
    times<-out[,1]           
    length_times<-length(times)
    
    #-------------------------------------------------------------------------------------------------------------
    #------------------- Identify the first column (index) for each category of output ---------------------------
    #-------------------------------------------------------------------------------------------------------------
    #The object out is organized as follows:
    # - Host compartments (if included) --> nb. columns = (size of the system) x (nb. vaccination levels)
    # - Vector compartments --> nb. columns = 9
    # - Incidence --> nb. columns = (nb. age groups) x (4 serotypes) x (3 infection types) x (nb. vaccination levels)
    # - Host population size --> nb. columns = (nb. age groups) 
    # - Vector population size --> nb. columns = 1
    # - Seroprevalence --> all remaining columns
    
    #Vector compartments
    if(hosts_included){ 
      ind_StateH<-2
      ind_StateV<-ind_StateH+NstatesH*vac_levels
    }else{
      ind_StateV<-2
    }
    
    #Incidence
    ind_incidence<-ind_StateV+9
    
    #Host population size
    ind_NHit<-ind_incidence+N_age_groups*4*3*vac_levels
    
    #Vector population size
    ind_NV<-ind_NHit+N_age_groups*(nb_vac+1)       #Population size by vaccination status is tracked for each vaccine (+ 1 level for unvaccinated population)
    
    #Seroprevalence
    ind_seroprevalence<-ind_NV+1
    
    
    #-------------------------------------------------------------------------------------------------------------
    #------------------- Aggregate the results -------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    
    #---- Vector compartments -------
    resultV<-out[,ind_StateV:(ind_incidence-1)]    
    
    #---- Incidence -----------------
    #Incidence --> nb. columns = (nb. age groups) x (4 serotypes) x (3 infection types) 
    #Determine dimensions for the incidence array
    dim_incid_array<-c(length_times,N_age_groups,4,3)
    
    
    #Extract incidence from the object out    
    incidence<-out[,ind_incidence:(ind_NHit-1)]
    
    #Transform into an array
    #NB! During the integration, the incidence is recorded as cumulative values --> difference between the adjacent lines (days) represents the nb. of new cases
    incidence=array(rbind(incidence[1,],incidence[-1,]-incidence[-length_times,]), dim=dim_incid_array)
    
    #---- Host population size ------
    NHit<-out[,ind_NHit:(ind_NV-1)]
    
    #---- Vector population size ----
    NV<-out[,ind_NV:(ind_seroprevalence-1)]
    
    #---- Seroprevalence ------------
    seroprevalence<-array(out[,ind_seroprevalence:length(out[1,])], dim=c(length_times,N_age_groups,8))
    
    #-------------------------------------------------------------------------------------------------------------
    #------------------- Remove the first line (if needed) -------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    if(!keep_first){
      times=times[-1]
      resultV=resultV[-1,]
      incidence=incidence[-1,,,]
      NHit=NHit[-1,]
      NV=NV[-1]
      seroprevalence=seroprevalence[-1,,]
    }
    
    #-------------------------------------------------------------------------------------------------------------
    #------------------- Return results --------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------
    
    return(list(times=times,
                resultV=resultV,
                incidence=incidence,
                NHit=NHit,
                NV=NV,
                seroprevalence=seroprevalence))
  })
}

