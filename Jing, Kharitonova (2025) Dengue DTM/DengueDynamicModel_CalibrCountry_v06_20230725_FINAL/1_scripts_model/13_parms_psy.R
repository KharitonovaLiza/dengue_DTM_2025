##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contains functions that calculate the object psy (infectiousness)
###########################################################################################################################################################################################################

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Psy is a modifying factor used in the calculation of the force of infection for vectors. It represents the relative infectiousness of different categories of hosts (i.e. their relative contribution into 
#the force of infection). The infectiousness depends on:
#- Whether the infection is symptomatic or asymptomatic, and
#- Whether the host is vaccinated or not (the model allows exploring a change in the transmission potential of vaccinated individuals, including sterile immunity).

#Since the model does not include explicit compartments for symptomatic and asymptomatic infections, the value of psy represents the average infectiousness weighted by the proportions of symptomatic and
#asymptomatic infections. These proportions, in their turn, depend on the infection type (primary, secondary or post-secondary), which can be determined by observing the number of previous infections.
#The function below outputs a 3D array where each elements is the relative infectiousness for each combination of the indices k, l and m. Similarly to the calculation of other array objects,
#here the indices k, l and m DO NOT always represent the status for DENV-2, DENV-3 and DENV-4 (as it is the case in the general notation of the model equations). Instead, they represent the three serotypes
#other than the one being considered in the current infection process. For instance, in the infection process for DENV-1, the indices k, l and m will indeed represent the status for DENV-2, DENV-3 and DENV-4, respectively. 
#However, in the infection process for DENV-2, they will represent the status for DENV-1, DENV-3 and DENV-4, respectively. Same for the other serotypes.

#The proportions of symptomatic and asymptomatic infections may as well depend on the host's vaccination status. The present function is, thus, called multiple times and a different array psy is calculated
#for each population category, including:
#-Unvaccinated hosts (called by the function ParmsEpi_Initialise)
#-Vaccinated hosts (seronegative at vaccination; called by the function ParmsVac_Initialise)
#-Vaccinated hosts (seropositive at vaccination; called by the function ParmsVac_Initialise)

Parms_GetPsy<-function(prop_sympt,              #Proportion of symptomatic infections with primary, secondary and post-secondary infections (numeric vector with 3 elements)
                       level_infect_sympt,      #Infectiousness of a symptomatic infection (relative to an asymptomatic infection in an unvaccinated individual)
                       level_infect_asympt){    #Infectiousness of an asymptomatic infection (relative to an asymptomatic infection in an unvaccinated individual) 
                                                #Note: the infectiousness of an asymptomatic infection in an unvaccinated individual is used as a reference and is always equal to 1
  
  #Calculate average weighted infectiousness --> numeric vector with 3 elements
  infectiousness<-prop_sympt[]*level_infect_sympt+(1-prop_sympt[])*level_infect_asympt    
  
  
  #Calculate the array psy
  psy<-Parms_GetInfSpecValues(value_primary=infectiousness[1], 
                              value_secondary=infectiousness[2], 
                              value_post_secondary=infectiousness[3])
  return(psy)
}


#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This generic function calculates, for each combination of indices k,l and m, the values, which depend on the number of past infections 

Parms_GetInfSpecValues<-function(value_primary,               #Value of the parameter in question for the primary infection
                                 value_secondary,             #Value of the parameter in question for the secondary infection
                                 value_post_secondary){       #Value of the parameter in question for the post-secondary infection
  
  past_infections<-Parms_GetPastInfections()                  #3D array with the number of ongoing or completed infections for each combination of indices k,l, and m
  
  inf_spec_values<-array(NA, dim=c(5,5,5))          #3D array with 5 levels for each dimension
  
  inf_spec_values<-((past_infections==0)*value_primary +             #If zero past infections --> apply parameter value for the primary infection (as the infection process that is being considered may result in a primary infection)
                    (past_infections==1)*value_secondary +           #If one past infection --> apply parameter value for the secondary infection
                    (past_infections==2)*value_post_secondary +      #If two or three past infections --> apply parameter value for the post-secondary infection
                    (past_infections==3)*value_post_secondary)
  
  return(inf_spec_values)
}
