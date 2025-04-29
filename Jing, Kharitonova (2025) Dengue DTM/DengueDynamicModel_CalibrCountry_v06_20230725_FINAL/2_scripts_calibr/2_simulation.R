##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  SCRIPT NAME: 2_simulation.R (for calibration!)
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA,   Olivier CRISTEAU , Aurélien JAMOTTE
###########################################################################################################################################################################################################

#Simulation_Run_TwoPeriods - This function runs the simulation separately for the burn-in and analysis periods. 
# In both cases the simulation is run year by year and the last estimation in each loop is used as the initial state for the next loop. 
# The results from the burn-in period are discarded (with the exception of the estimations from the last time point, which are used to initialise the simulation for the analysis period).
# The results from the analysis period are used to estimate the model fit.
# The initial state of the model is simulated

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function runs the simulation separately for the burn-in and analysis periods. In both cases the simulation is run year by year and the last estimation in each loop is used as the initial state for
#the next loop. The results from the burn-in period are discarded (with the exception of the estimations from the last time point, which are used to initialise the simulation for the analysis period).
#The results from the analysis period are used to estimate the model fit. The initial state is created.

Simulation_Run_TwoPeriods<-function(yini,                   #Vector with the initial state
                                    parms_epi){                 #List of population & transmission inputs
  #----------------------------------------------------------------------
  #------- PREPARE THE MODEL PARAMETERS ---------------------------------
  #----------------------------------------------------------------------
  parms_C<-c(parms_epi$parms_epi_dyn_model)

  #----------------------------------------------------------------------
  #------- LOAD C MODEL -------------------------------------------------
  #----------------------------------------------------------------------
  script_name<-"M1_no_vac"                                                 #Name of the C script (w/o extension)
  dll_name<-paste(script_name,.Platform$dynlib.ext, sep="")                #Name of the dll file
  dyn.load(file.path(path_scripts_src, dll_name), local=TRUE, now=FALSE)   #Load the model
    
  #----------------------------------------------------------------------
  #------- SIMULATE BURN-IN PERIOD --------------------------------------
  #----------------------------------------------------------------------
  for(year in 1:dur_period1){
    
    #print(paste("Running year ", year, "/", dur_period1, sep=""))
    
    #------ Determine current time points --
    #For the burn-in period, 13 time points are simulated at each loop (12 for each month of the current year and 1 for the first day of the next year)
    #The exception is the last loop (last year in the burn-in period), which has 12 points for each month of the last year and 1 additional point for the last day of the last year (initial state for the analysis period)
    
    times_current<-times_period1_year[[year]]                                   #List of time points to be simulated at the current loop (13 elements)
    
    out<-rk(y=yini,                    #Initial state vector (host compartments, vector compartments, incidence compartments)
            times=times_current,       #Vector of time points to be estimated (day numbers)
            parms=parms_C,             #List of parameters to be passed to C code
            func="model",              #Name of the model function (as specified in the C code)
            dllname=script_name,       #Name of the dll file (w/o extension)
            initfunc="initmod",        #Name of the function initialising the model parameters (as specified in the C code)  
            nout=(9*N_age_groups+1),   #Number of additional outputs to be recorded (length of the corresponding vector)  
            method=ode_method,         #Integration method
            hini=ode_hini,             #Initial step size to be attempted
            atol=ode_atol)             #Absolute error tolerance
    
    #------ Get new initial state ----------
    yini<- out[length(out[,1]), 2:(1+NstatesH+9+N_age_groups*4*3)]             #Last estimation point in the current loop represents the the initial state for the next year
    
    #------ Apply ageing -------------------
    yH_ageing<-Process_Age(stateH_NV=yini[1:NstatesH], parms_epi=parms_epi)    #Calculate the system for hosts after ageing
    yini[1:NstatesH]<-yH_ageing                                                #Inject new system for hosts into the initial state for the next year
    
    #------ Output compartments ------------
    #The state of the system is output at the end of a specified year. It is then used as the initial state when running the core model with these parameters
    
     if(year==dur_period1){
      
      row_nb<-length(out[,1])
      # print(paste("Year is ", year, ". Outputing the ",row_nb, "th line of this year to be used as the initial state"))
  
      #For hosts
      out_hosts<-out[row_nb,2:(NstatesH+1)]                              #Taking last line of the solver output, as it represents the situation at the end of the year
      out_hosts_array<-array(out_hosts, dim=c(N_age_groups,5,5,5,5))     #Transforming into an array
      
      #For vectors
      out_vectors<-out[row_nb,(NstatesH+2):(NstatesH+2+9-1)]             #Also taking last line
    }
    
  }
  
  #----------------------------------------------------------------------
  #------- SIMULATE MAIN PERIOD -----------------------------------------
  #----------------------------------------------------------------------
  
  for (year in 1:dur_period2){
    
    #print(paste("Running year ", year, "/", dur_period2, sep=""))
    
    #------ Determine current time points ---
    #For the analysis period, 366 time points are simulated at each loop (365 for each day of the current year and 1 for the last day of the previous year)
    
    times_current=times_period2_year[[year]]     #List of time points to be simulated at the current loop (366 elements, 365 elements for the last simulation year)
    
    #------ Run the ODE solver -------------
    out<-rk(y=yini,                    #Initial state vector (host compartments, vector compartments, incidence compartments)
            times=times_current,       #Vector of time points to be estimated (day numbers)
            parms=parms_C,             #List of parameters to be passed to C code
            func="model",              #Name of the model function (as specified in the C code)
            dllname=script_name,       #Name of the dll file (w/o extension)
            initfunc="initmod",        #Name of the function initialising the model parameters (as specified in the C code)  
            nout=(9*N_age_groups+1),   #Number of additional outputs to be recorded (length of the corresponding vector)  
            method=ode_method,         #Integration method
            hini=ode_hini,             #Initial step size to be attempted
            atol=ode_atol)             #Absolute error tolerance
    
    #------ Get new initial state ----------
    yini<-out[length(out[,1]), 2:(1+NstatesH+9+N_age_groups*4*3)]              #Last estimation point in the current loop represents the initial state for the next year
    
    #------ Apply ageing -------------------
    yH_ageing<-Process_Age(stateH_NV=yini[1:NstatesH], parms_epi=parms_epi)    #Calculate new state of the system after ageing
    yini[1:NstatesH]<-yH_ageing                                                #Inject new system for hosts into the initial state for the next year
    
    #------ Output results -----------------
    if(year==1){                                                                     #If the first simulation year, the object with the results does not yet exist...
      out_wo_hosts=out[,c(1,(NstatesH+2):(length(out[1,])))]                         #...create it by saving solver output (without host compartments, which are not required for the calibration)
    }else{                                                                           #If not the first year...
      out_wo_hosts=rbind(out_wo_hosts, out[-1,c(1,(NstatesH+2):(length(out[1,])))])  #... add new results to the existing object (without the first estimation point from the current loop, which was already included in the results for the previous loop)
    }

    rm(out)
  }
  
  #----------------------------------------------------------------------
  #------- UNLOAD C MODEL -----------------------------------------------
  #----------------------------------------------------------------------  
  dyn.unload(file.path(path_scripts_src, dll_name))
  
  return(list(out_wo_hosts=out_wo_hosts,
              out_hosts_array=out_hosts_array,
              out_vectors=out_vectors))  
}



