##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  TASK: Calibrating the model for a county based on observed incidence
#  SCRIPT NAME: 0_calibration_launch.R
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA
###########################################################################################################################################################################################################

#The code below allows fitting the model to age-specific incidence reported in the country of interest
#The parameters that are being calibrated:
#- betaVH (probability of virus transmission from an infectious vector to a susceptible host given a bite)
#- seasonality parameters (amplitude and horizontal shift of the sine function applied to the size of vector population)

#The parameters are being calibrated one after another (i.e. first betaVH is calibrated given some plausible values of seasonality; then seasonality is calibrated given the best-fitting value of beta)
#If needed, this process can be repeated (i.e. betaVH can be re-calibrated with the best-fitting values of seasonality)

#To calibrate, the model is run over a long period of time to let it equilibrate.
#The last interval of the simulation is then used to assess the average predicted incidence and to compare it to empirical (observed) data
#The goodness of fit is expressed as log-likelihood.
#Based on the goodness of fit, the algorithm adjusts the value of betaVH and runs again until the fit can no longer be significantly improved.

#NOTE ON LIMITATIONS:
#- This code does not work if composite age groups are used for the calibration (i.e. less than 101)
#- This code does not work if betaVH is different for different serotypes

rm(list=ls())   #Clearing the workspace before running a new calibration

#=======================================================================================================================
#--- PROVIDE ALL PARAMETERS & SETTINGS (CAN BE MODIFIED) ---------------------------------------------------------------
#=======================================================================================================================

#-- FILE WITH EPI INPUTS --------------
epi_file_name<-"calibr_inputs"         #Name of the Excel file with population & transmission inputs (without extension .xlsx)

#-- SIMULATION DURATIONS --------------
dur_period1<-650                       #Duration of the period used to equilibrate the model (burn-in period)
dur_period2<-100                       #Duration of the period used to estimate the average predicted seroprevalence (and to search for a period with the target serotype distribution if this is required, see below)

#-- ODE SOLVER SETTINGS ---------------
ode_method<-"ode23"                    #"ode23" - an integration method for ODE systems using 2nd and 3rd order Runge-Kutta-Fehlberg formulas with automatic step size
#                                      #"euler" - an integration method with fixed step size
ode_hini<-0                            #Initial step size to be attempted. If 0, the initial step size is determined automatically by solvers with flexible time step
ode_atol<-1e-06                        #Absolute error tolerance. Value of 1e-06 means an absolute error of 0.000001 is tolerated for each compartment during the integration.
#                                      #As the population consist of thousands (or millions) of people, this level of accuracy was deemed sufficient

#-- OPTIM FUNCTION SETTINGS -----------
optim_reltol=0.001                     #Relative convergence tolerance. The algorithm will stop if it is unable to reduce the value by a factor of reltol*(abs(val)+reltol) at a step.

#-- NB AGE GROUPS ---------------------
N_age_groups<-101
NstatesH<-N_age_groups*5*5*5*5

#-- SEROTYPE WEIGHTS ------------------
prop_serotype<-c(0.1, 0.2, 0.3, 0.4)   #Serotype weights are used to destabilize the number of infectious vectors in the beginning of the calibration 
#                                      #(if betaVH is the same for each serotype this is required to avoid identical incidence curves for each serotype and to reproduce serotype replacement)

#-- PARMS TO CALIBRATE ----------------
calibr_betaVH<-TRUE                   #Flag whether BetaVH is calibrated
calibr_seasonality<-FALSE               #Flag indicating whether seasonality parameters are calibrated

#--- IMITATE BETA CALIBRATION TO SAVE BEST_FIT FILES
save_SP_data_only<-FALSE
                                     
#=======================================================================================================================
#--- PREPARE FOR CALIBRATION (THIS CODE SHOULD NOT BE MODIFIED) --------------------------------------------------------
#=======================================================================================================================

#--------------------------------------
#-- VACCINATION STATUSES --------------
#--------------------------------------
vac_levels<-1                           #The calibration code re-uses some of the codes from the core model 
nb_vac<-0                               #These parameters indicate that the model is being run without vaccination

#--------------------------------------
#-- ESTIMATION TIME POINTS ------------
#--------------------------------------
#Total period duration
dur_total=dur_period1+dur_period2

#Time points for each period
times_period1<-round(((0:11)+rep(0:(dur_period1-1), each=12)*12)*(365/12)+1,digits=0)   #Burn-in period, by month (in order to save memory and reduce runtime)
times_period2<-(((dur_period1)*365):((dur_period1+dur_period2)*365))                    #Analysis period, by day
times_total<-c(times_period1, times_period2)                      						          #Total period, by month in burn-in period / by day in analysis period(in order to run the model with the best-fitting parameters)
times_period1=c(times_period1, times_period2[1])                                        #Adding the first day of the analysis period to be simulated in the burn-in period (to obtain the initial state for the analysis period)

#Time points for each year within each period
times_period1_year<-list()
times_period2_year<-list()
for (i in 1:dur_period1){times_period1_year[[i]]<-times_period1[((i-1)*12+1):(i*12+1)]}
for (i in 1:dur_period2){times_period2_year[[i]]<-times_period2[((i-1)*365+1):(i*365+1)]}
times_total_year<-times_period1_year
times_total_year<-c(times_total_year, times_period2_year)

times_total_year[[dur_total]]<-times_total_year[[dur_total]][-366]                       #For the last year to be simulated, remove the last estimation point (no need to simulate the initial state for the next year)


#--------------------------------------
#-- PATHS TO SUB-FOLDERS --------------
#--------------------------------------
path_wd<-getwd()                                                  #Working directory
path_scripts_model<<-file.path(path_wd, "1_scripts_model")        #Scripts for transmission model (R) 
path_scripts_src<<-file.path(path_scripts_model, "src")           #Script for infection process (C)
path_scripts_calibr<-file.path(path_wd, "2_scripts_calibr")       #Scripts for calibration
path_data_calibr<-file.path(path_wd, "0_data_calibr")             #Data for calibration

#--------------------------------------
#-- REQUIRED FUNCTIONS ----------------
#--------------------------------------
#Scripts for the dynamic model
debugSource(file.path(path_scripts_model, "11_parms_common.R"))             
debugSource(file.path(path_scripts_model, "12_parms_theta.R"))
debugSource(file.path(path_scripts_model, "13_parms_psy.R"))
debugSource(file.path(path_scripts_model, "21_parms_epi.R"))

#Scripts for calibration
debugSource(file.path(path_scripts_calibr, "2_simulation.R"))
debugSource(file.path(path_scripts_calibr, "3_ageing.R"))
debugSource(file.path(path_scripts_calibr, "5_calibration_dynamic.R"))
debugSource(file.path(path_scripts_calibr, "8_results_epi.R"))

#--------------------------------------
#-- REQUIRED PACKAGES -----------------
#--------------------------------------
options(repos="http://cran.us.r-project.org")             #set the repository to download the new packages from
if (!("openxls" %in% rownames(installed.packages()))){install.packages("openxls")}; require(openxls)
if (!("deSolve" %in% rownames(installed.packages()))){install.packages("deSolve")}; require(deSolve)
if (!("ggplot2" %in% rownames(installed.packages()))){install.packages("ggplot2")}; require(ggplot2)
if (!("reshape2" %in% rownames(installed.packages()))){install.packages("reshape2")}; require(reshape2)
if (!("dplyr" %in% rownames(installed.packages()))){install.packages("dplyr")}; require(dplyr)

#--------------------------------------
#-- EPI PARAMETERS --------------------
#--------------------------------------
parms_epi_imported<-Parms_Import(path=path_data_calibr, filename=epi_file_name)                  #Import
parms_epi_selected<-ParmsEpi_Select(parms_epi_imported=parms_epi_imported,calibration=TRUE)      #Pre-select
rm(parms_epi_imported)                                                                           #Remove the object that is no longer required

# print(paste("Dur CP: ", parms_epi_selected$parms_epi$dur_CP))
# print(paste("Infectiousness sympt: ", parms_epi_selected$parms_epi$level_infect_sympt_unvac))
# print(paste("Infectiousness asympt: ", parms_epi_selected$parms_epi$level_infect_asympt_unvac))
# print(paste("Beta_VH: ", parms_epi_selected$parms_epi$betaVH[1]))
#--------------------------------------
#-- RUN THE CALIBRATION ---------------
#--------------------------------------

CalibrDyn_Run()

