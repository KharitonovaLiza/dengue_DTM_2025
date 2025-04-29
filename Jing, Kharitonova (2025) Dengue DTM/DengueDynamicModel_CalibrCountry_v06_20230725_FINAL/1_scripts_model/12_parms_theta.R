##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contains functions that calculate the object theta
###########################################################################################################################################################################################################

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Theta is a modifying factor applied to the force of infection for hosts in order to take into account:
# - The impossibility of co-infection,
# - The presence of cross-protection & relative risk of being infected during the cross-protection period,
# - The maximum allowed number of infections, and
# - The presence of susceptibility enhancement & increase in susceptibility after the first (or each) infection

#This function calculates the array theta for UNVACCINATED HOSTS.
#It outputs a 3D array containing the value of theta for each combination of k, l and m. Here, the indices k, l and m DO NOT always represent the status for DENV-2, DENV-3 and DENV-4 
#(as it is the case in the general notation of the model equations). Instead, they represent the three serotypes that are not being considered in the current infection process. 
#For instance, in the infection process for DENV-1, the indices k, l and m will indeed represent the status for DENV-2, DENV-3 and DENV-4, respectively. 
#However, in the infection process for DENV-2, they will represent the status for DENV-1, DENV-3 and DENV-4, respectively. Same for the other serotypes.

#The IMPOSSIBILITY OF CO-INFECTION is taken into account by setting the value of theta to zero if at least one of the indices k,l or m is equal to 2 or 3 (which indicates that the
#individual is currently exposed to or infected with another serotype). The theta equal to zero is then applied to the force of infection for the serotype being considered, 
#thus making an infection with this serotype impossible. This is done in the infection process for each serotype, thus guaranteeing that no new infection can occur while another infection is ongoing.

#Similarly, the IMPACT OF CROSS-PROTECTION is taken into account by setting the value of theta equal to the parameter gammaCP (relative risk of being infected during the cross-protection period). 
#This is done if k,l or m are equal to 4, which indicates that the individual is currently cross-protected following an earlier infection with another serotype.

#If the indices k,l and m indicate that the host has already had the MAXIMUM ALLOWED NUMBER OF INFECTIONS, theta is set to zero to make a new infection impossible.

#Finally, the value of theta can be greater than 1 in the presence of SUSCEPTIBILITY ENHANCEMENT (after the first or after each infection).

Parms_GetThetaUnvac<-function(gammaCP,             	  #Relative risk of being infected during the cross-protection period
                              dzetaCE,                #Presence & frequency of susceptibility enhancement (type of cross-enhancement, CE)
                                                      # =0 --> not included
                                                      # =1 --> applied after the 1st infection only
                                                      # =2 --> applied after each infection
                              gammaCE,                #Susceptibility enhancement factor
                              nb_infect_max){         #Maximum allowed number of infections
  
  #------------------------------
  #Initialize array theta
  #------------------------------
  theta=array(NA, dim=c(5,5,5))                           #3D array with 5 levels for each dimension (as k, l and m can have values from 1 to 5)
  
  #------------------------------
  #Account for cross-protection & impossibility of co-infection
  #------------------------------
  sigma=c(1,             #k/l/m = 1: Susceptible to a serotype --> no modification of the force of infection for another serotype (force of infection for this other serotype is multiplied by 1)
          0,             #k/l/m = 2: Exposed to a serotype --> infection with another serotype is impossible (force of infection for this other serotype is multiplied by 0)
          0,             #k/l/m = 3: Infected with a serotype --> infection with another serotype is impossible (force of infection for this other serotype is multiplied by 0)
          gammaCP,       #k/l/m = 4: Cross-protected following an infection with a serotype --> the probability of infection with another serotype is reduced or zero 
                         #(force of infection for this other serotype is multiplied by the relative risk gammaCP comprised between 0 and 1)
          1)             #k/l/m = 5: Immune to a serotype --> no modification of the infection process for another serotype (force of infection for this other serotype is multiplied by 1)
  
  
  for(k in 1:5){                    #Loop through all elements of the array theta
    for(l in 1:5){
      for(m in 1:5){
        
        values_klm<-c(k,l,m)        #A vector with the current values of indices k, l and m

        #------------------------------
        #Account for the maximum number of infections
        #------------------------------
        if(sum((values_klm >= 2) & (values_klm <= 5))>=nb_infect_max){      #Count the number of previous infections (ongoing or completed; indices equal to 2, 3, 4 or 5) and compare to the maximum allowed number of infections
          theta[k,l,m]<-0                                                   #If the number of previous infections is already equal to the maximum, a new infection is impossible (theta equal to 0)
        }else{                                                              
          theta[k,l,m]<-min(sigma[k], sigma[l], sigma[m])                   #Otherwise, theta is equal to the minimum of the three values sigma (before susceptibility enhancement is taken into account)
          
          #------------------------------
          #Account for the susceptibility enhancement                     
          #------------------------------
          #If susceptibility enhancement is applied (dzetaCE!=0) --> apply after the 1st infection
          if(dzetaCE!=0){                                                     #Here, susceptibility enhancement is only applied if values_klm = c(5,1,1) (in any order). This is the only combination of indices that indicates:
            if((sum(values_klm==1)==2) & (sum(values_klm==5)==1)){            # - History of exactly one infection,
              theta[k,l,m]<-theta[k,l,m]*gammaCE                              # - No ongoing infection with another serotype (in which case theta should be equal to 0), and
            }                                                                 # - No ongoing cross-protection from a previous infection (in which case theta should be equal to gammaCP)
          }                                                                   #The value of theta is then multiplied by the susceptbility enhancement factor 
          
          #If susceptibility enhancement is applied after each infction --> also apply after second & third infections
          if(dzetaCE==2){                                                     #Here, susceptibility enhancement is only applied if values_klm = c(5,5,1) (in any order) or values_klm = c(5,5,5)
            if(((sum(values_klm==1)==1) & (sum(values_klm==5)==2))||(sum(values_klm==5)==3)){                                     
              theta[k,l,m]<-theta[k,l,m]*gammaCE                                           
            }                                                                   
          }        
        }
      }
    }
  }
  
  return(theta)      
}


#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function calculates the arrays theta for VACCINATED INDIVIDUALS.
#Two separate arrays are calculated - for the individuals who were seronegative and seropositive at the moment of vaccination. The calculation process is identical for the 
#both arrays, but the function is called twice with the different sets of vaccine efficacy values (object eta). 

#In addition to the effects considered for unvaccinated population (impossibility of co-infection, presence of cross-protection etc.), the efficacy of the vaccine is taken into account.
#To do so, the array theta now has two additional dimensions for:
#-Infecting serotype (as the vaccine efficacy may be serotype-specific),
#-Age group (determining, for each cohort, how long ago vaccination took place and, thus, remaining efficacy).
#These additional effects are captured in the array eta (for the calculation details see function ParmsVac_GetEta in the script 13_parms_vac.R)

#The function below, thus, returns a 5D array with the values of theta for each age group, serotype and combination of indices k, l and m.

Parms_GetThetaVac<-function(N_age_groups,   #Number of age groups used in the model
                            gammaCP,        #Relative risk of being infected during the cross-protection period
                            dzetaCE,        #Presence & frequency of susceptibility enhancement (type of cross-enhancement, CE)
                                            # =0 --> not included
                                            # =1 --> applied after the 1st infection only
                                            # =2 --> applied after each infection
                            gammaCE,        #Susceptibility enhancement factor
                            nb_infect_max,	#Maximum allowed number of infections,
                            eta){           #Vaccine efficacy (output of the function ParmsVac_GetEta)
  
  theta<-array(NA, dim=c(N_age_groups,4,5,5,5))  #5D array with the values of theta for each age group (N_age_groups), serotype (4) and combination of indices k, l and m (5x5x5)

  for(i in 1:N_age_groups){                      #Loop through all elements of the array theta
    for(s in 1:4){
      for(k in 1:5){
        for(l in 1:5){
          for(m in 1:5){
            #------------------
            #Account for cross-protection, impossibility of co-infection & vaccine efficacy
            #------------------
            sigma=c(1-eta[i,s],       #k/l/m = 1: Susceptible to a serotype --> the force of infection for another serotype is only modified by the vaccine efficacy eta
                    0,             		#k/l/m = 2: Exposed to a serotype --> infection with another serotype is impossible (force of infection for this other serotype is multiplied by 0)
                    0,             		#k/l/m = 3: Infected with a serotype --> infection with another serotype is impossible (force of infection for this other serotype is multiplied by 0)
                    gammaCP,       		#k/l/m = 4: Cross-protected following an infection with a serotype --> the probability of infection with another serotype is reduced or zero 
                    #                             (force of infection for this other serotype is multiplied by the relative risk gammaCP comprised between 0 and 1)
                    1-eta[i,s])       #k/l/m = 5: Immune to a serotype --> the force of infection for another serotype is only modified by the vaccine efficacy eta
            
            values_klm<-c(k,l,m)            #A vector with the current values of indices k, l and m
            
            #------------------
            #Account for the maximum number of infections
            #------------------          
            if(sum((values_klm >= 2) & (values_klm <= 5))>=nb_infect_max){  #Count the number of previous infections (ongoing or completed; indices equal to 2, 3, 4 or 5) and compare to the maximum allowed number of infections
              theta[i,s,k,l,m]<-0                                         #If the number of previous infections is already equal to the maximum, a new infection is impossible (theta equal to 0)
            }else{
              theta[i,s,k,l,m]<-min(sigma[k], sigma[l], sigma[m])         #Otherwise, theta is equal to the minimum of the three values sigma (before susceptibility enhancement is taken into account)
              
              #------------------
              #Account for the susceptibility enhancement   
              #------------------
              #If susceptibility enhancement is applied (dzetaCE!=0) --> apply after the 1st infection
              if(dzetaCE!=0){                                                #Here, susceptibility enhancement is only applied if values_klm = c(5,1,1) (in any order). This is the only combination of indices that indicates:
                if((sum(values_klm==1)==2) & (sum(values_klm==5)==1)){       # - History of exactly one infection,
                  theta[i,s,k,l,m]<-theta[i,s,k,l,m]*gammaCE                 # - No ongoing infection with another serotype (in which case theta should be equal to 0), and
                }                                                            # - No ongoing cross-protection from a previous infection (in which case theta should be equal to gammaCP)
              }                                                              #The value of theta is then multiplied by the susceptibility enhancement factor
              
              #If susceptibility enhancement is applied after each infection --> also apply after second & third infections                                          
              if(dzetaCE==2){                                                #Here, susceptibility enhancement is only applied if values_klm = c(5,5,1) (in any order) or values_klm = c(5,5,5)
                if(((sum(values_klm==1)==1) & (sum(values_klm==5)==2))||(sum(values_klm==5)==3)){ 
                  theta[i,s,k,l,m]<-theta[i,s,k,l,m]*gammaCE
                }                                                                   
              } 
            }              
          }
        }
      }
    }
  }
  
  return(theta)
}

