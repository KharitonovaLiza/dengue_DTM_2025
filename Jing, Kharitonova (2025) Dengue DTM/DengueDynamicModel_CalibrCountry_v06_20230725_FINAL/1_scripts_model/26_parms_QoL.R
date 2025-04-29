##########################################################################################################################################################################################################
#  PROJECT: Dengue dynamic model
#  AUTHORS: Elizaveta KHARITONOVA, Anna TYTULA, Olivier CRISTEAU, Aurélien JAMOTTE
#  DESCRIPTION: This script contains a function that imports and loads the quality-of-life inputs
###########################################################################################################################################################################################################

#This function loads the quality-of-life inputs
ParmsQoL_Load<-function(path_inputs_QoL,
                        set_QoL){
  
  parms_QoL_imported<-Parms_Import(path=path_inputs_QoL, filename=set_QoL) 
  
  parms_QoL_final = list(
    disab_weight_mild=as.numeric(parms_QoL_imported$disab_weight_mild),                #Disability weight for a mild infection --> numeric
    disab_weight_severe=as.numeric(parms_QoL_imported$disab_weight_severe),            #Disability weight for a severe infection --> numeric
    disab_weight_long_term=as.numeric(parms_QoL_imported$disab_weight_long_term),      #Disability weight for long-term (persistent) dengue --> numeric
    
    disab_dur_mild=as.numeric(parms_QoL_imported$disab_dur_mild),                      #Duration of a mild infection (days) --> numeric
    disab_dur_severe=as.numeric(parms_QoL_imported$disab_dur_severe),                  #Duration of a severe infection (days) --> numeric
    dur_long_term=as.numeric(parms_QoL_imported$dur_long_term),                        #Duration of long-term (persistent) dengue (months) --> numeric
    
    DALY_weigh_by_age=as.logical(parms_QoL_imported$DALY_weigh_by_age),                #Weighing of DALYs by age --> boolean
    DALYs_param_C=as.numeric(parms_QoL_imported$DALYs_param_C),                        #C parameter for the DALY function --> numeric
    DALYs_param_b=as.numeric(parms_QoL_imported$DALYs_param_b))                        #b parameter for the DALY function --> numeric
  
  return(list(parms_QoL=parms_QoL_final))
}
