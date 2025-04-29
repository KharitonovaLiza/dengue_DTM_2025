##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contains function that import and load the vector control effects inputs
###########################################################################################################################################################################################################

#This function selects & formats the imported inputs for vector control effects
ParmsVecEff_Load<-function(path_inputs_vec_eff,         #Path to the folder containing the Excel file with the vector control inputs
                           set_vec_eff){                #Name of the Excel file containing the inputs (without extension) 
  
  parms_vec_eff_imported<-Parms_Import(path=path_inputs_vec_eff, filename=set_vec_eff) 
  
  parms_vec_eff<-list(
    vec_switch=as.logical(parms_vec_eff_imported$vec_switch),                     #Inclusion of vector control (required in the model interface) --> Boolean
    vec_ratio_VH_change=as.numeric(parms_vec_eff_imported$vec_ratio_VH_change),  	#Number of adult female vectors per host after the introduction of vector control (in % of original value) --> Numeric
    vec_b_change=as.numeric(parms_vec_eff_imported$vec_b_change)  					      #Daily biting rate of an adult female vector after the introduction of vector control (in % of original value) --> Numeric
  )
  
  return(list(parms_vec_eff=parms_vec_eff))
}