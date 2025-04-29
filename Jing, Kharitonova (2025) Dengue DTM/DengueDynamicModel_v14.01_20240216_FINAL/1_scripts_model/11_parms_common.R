##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contains functions that are used throughout the model
###########################################################################################################################################################################################################

#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#This function imports all inputs from a specified CSV file

Parms_Import<-function(path,        #The path to the folder containing the CSV file to be imported
                       filename){   #The name of the CSV file to be imported (without the folder path and without the file extension .csv)
  
  #Import parameters
  parms<-read.csv(file=file.path(path, paste(filename, ".csv", sep="")), stringsAsFactors=FALSE)
  parms_list<-as.list(parms)
  
  #Remove N/A values that appear because of the different number of elements in each column
  nb_parms<-length(parms_list)
  for(i in 1:nb_parms){
    parms_list[[i]]<-parms_list[[i]][!is.na(parms_list[[i]])]
  }
  return(parms_list) #A list of inputs as list (one column in the CSV file = one element of the list)
}


#=========================================================================================================================================================================================================
#=========================================================================================================================================================================================================
#Some of the model parameters depend on the infection type (primary, secondary or post-secondary). This function counts the number of ongoing or already completed infections to determine which value of 
#the input should be used in the current infection process. Here, the indices k, l and m DO NOT represent the status DENV-2, DENV-3 and DENV-4 (as it is the case in the general notation of the model).
#Instead, they represent the three serotypes that are not considered in the current infection process. For instance, in the infection process for DENV-1, the indices k, l and m will indeed represent the status 
#for DENV-2, DENV-3 and DENV-4, respectively. In the infection process for DENV-2, they will represent the status for DENV-1, DENV-3 and DENV-4, respectively. Same for the other serotypes.

Parms_GetPastInfections<-function(){
  
  #Create a 3D array with 5 levels for each dimension
  past_infections<-array(NA, dim=c(5,5,5))
  
  #Loop through each combination of k,l and m
  for(k in 1:5){
    for(l in 1:5){
      for(m in 1:5){
        values_klm<-c(k,l,m)
        past_infections[k,l,m]<-sum((values_klm >= 2) & (values_klm <= 5))      #Count the number of ongoing or completed infections (indices equal to 2,3,4 or 5)                       
      }
    }
  }
  return(past_infections)      #A 3D array where each element represents the number of previous infections for each combination of indices k, l amd m
}


