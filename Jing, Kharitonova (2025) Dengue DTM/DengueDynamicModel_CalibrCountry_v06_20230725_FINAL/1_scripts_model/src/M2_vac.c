/*
 * Project: Dengue dynamic model
 * Author(s): Elizaveta KHARITONOVA, Anna TYTULA, Aurelien JAMOTTE, Olivier CRISTEAU
 * Description: This file contains the source code for the infection process without vaccination used in the calibration (composite age groups)
 */

/*This file contains the source code of the functions defining the model*/


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>
#include <stdio.h>
#include "M2_vac.h"  /*header of the file*/

/********* DEFINE FUNCTIONS USED IN THE MODEL ********************************************/

/*Position of an element in the vector StateH*/
int get_elt_StateH(int i,int j, int k, int l, int m){
  return (i+N_AGE_GROUPS*(j+NN*(k+NN*(l+NN*m))));
}

/*Position of an element in a 3D array with dimensions NN x NN x NN*/
int get_elt_3D_serotype(int j, int k, int l){ 
  return (j+NN*(k+NN*l));
}

/*Position of an element in the incidence vector*/
int get_elt_incidence(int i, int s, int status){ 
  return (i+N_AGE_GROUPS*(s+NS*status));
}

/*Position of an element in the seroprevalence vector*/
int get_elt_seroprevalence(int i, int status){ 
  return (i+N_AGE_GROUPS*status);
}

/*Position of an element in an object theta for vaccinated individuals */
int get_elt_theta_vac(int i, int s, int j, int k, int l){
  return (i+N_AGE_GROUPS*(s+NS*(j+NN*(k+NN*l))));
}

/*Position of an element in an object psy for vaccinated individual */
int get_elt_psy_vac(int i, int j, int k, int l){
  return (i+N_AGE_GROUPS*(j+NN*(k+NN*l)));
}
    

/*Sum of the elements of a vector*/
double Vsum(double vector1[], int length_vector){
  int i;
  double sum=0;
  for(i=0; i<length_vector; i++){
    sum+=vector1[i];
  }
  return(sum);
}

/*Product of the corresponding elements of two vectors of same length*/
void Vprod(double vector1[], double vector2[], int length_vector, double result[]){
  int i;
  for(i=0; i<length_vector; i++){
    result[i]=vector1[i]*vector2[i];
  }
}

/*Sum of the corresponding elements of two vectors of same length*/
void Vadd(double vector1[], double vector2[], int length_vector, double result[]){
  int i;
  for(i=0; i<length_vector; i++){
    result[i]=vector1[i]+vector2[i];
  }
}

/*Difference of the corresponding elements of two vectors of same length*/
void Vdiff(double vector1[], double vector2[], int length_vector, double result[]){
  int i;
  for(i=0; i<length_vector; i++){
    result[i]=vector1[i]-vector2[i];
  }
}

/*Sum of the products of the corresponding elements of two vectors of same length*/
double Vsum_prod(double vector1[], double vector2[], int length_vector){
  double result = 0;
  int i;
  FOR_I result += vector1[i]*vector2[i];
  return result;
}

/*Sum of all compartments for each age group*/
void sum_overi_StateH(double *StateH, double *NHi){
  int i,j,k,l,m;
  double sum;
  FOR_I {
    sum=0;
    FOR_J FOR_K FOR_L FOR_M sum+=StateH[G_ELT1(i,j,k,l,m)];
    NHi[i]=sum;
  }
}

/*Sum of the 4D array (in a vector form), for each age group */
void sum_overi_psy_StateH(double psy_stateH[][NN][NN][NN], double result[]){
  int i, j, k, l;
  FOR_I {
    result[i]=0;
    FOR_J FOR_K FOR_L result[i]+=psy_stateH[i][j][k][l];
  }
}

/*Sum over age of several 4D arrays*/
void sum_overi_psy_StatesH(double psy_stateH_1[][NN][NN][NN],double psy_stateH_2[][NN][NN][NN],double psy_stateH_3[][NN][NN][NN], double result[]){
  int i, j, k, l;
  FOR_I {
    result[i]=0;
    FOR_J FOR_K FOR_L result[i]+=psy_stateH_1[i][j][k][l]+psy_stateH_2[i][j][k][l]+psy_stateH_3[i][j][k][l];
  }
}
/*Multiplication of two objects - 4D array psy (infectiousness) and a vector StateH (size of all compartments for a given vaccination status)*/
/*pos_ind - serotype in question; val_ind - status for the serotype in question */
void Vprod_psy_StateH(double psy[N_AGE_GROUPS][NN][NN][NN], double StateH[], int pos_ind, int val_ind, double result[][NN][NN][NN]){
  int i,j,k,l,m;
  switch(pos_ind){
  case 1: 
    FOR_I FOR_K FOR_L FOR_M result[i][k][l][m]=psy[i][k][l][m]*StateH[G_ELT1(i,val_ind,k,l,m)];
    break;
  case 2: 
    FOR_I FOR_J FOR_L FOR_M result[i][j][l][m]=psy[i][j][l][m]*StateH[G_ELT1(i,j,val_ind,l,m)];
    break;
  case 3:
    FOR_I FOR_J FOR_K FOR_M result[i][j][k][m]=psy[i][j][k][m]*StateH[G_ELT1(i,j,k,val_ind,m)];
    break;
  case 4:
    FOR_I FOR_J FOR_K FOR_L result[i][j][k][l]=psy[i][j][k][l]*StateH[G_ELT1(i,j,k,l,val_ind)];
    break;
  default:
    printf("Error of value for pos_ind in function Vprod_psy: must be 1,2,3 or 4 while the value was %d", pos_ind);
  exit(EXIT_FAILURE);
  break;
  }
}

/********* FUNCTION INITIALIZING THE MODEL PARAMETERS ************************************/
/*This functions allows to initialize the model parameters passed to C code as C objects to be used in the model*/
/*This parameters are included in the list parms_C, that is passed to the R solver procedure as the argument "parms" */

void initmod(){
  int nparms;
  DL_FUNC get_deSolve_gparms;
  SEXP gparms;
  get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");
  gparms = get_deSolve_gparms(); 
  nparms = LENGTH(gparms);
  
  /*Check the number of parameters passed to the model */
  if (nparms!=N_PARMS){  
    printf("The number of parameters in the list parms_C should be %d while it is %d", N_PARMS, nparms);
    PROBLEM "The number of elements in the list parms_C is not as expected"
    ERROR;
    
  }else{
    /*If the number of parameters is correct, initialize each model parameter*/
    /*Note: In C the elements are numbered from 0*/
    /*Note: prefix p_ indicates vectors or arrays*/
    
    ratio_VH=REAL(VECTOR_ELT(gparms, 0))[0];                                //Scalar of type double
    muV=REAL(VECTOR_ELT(gparms, 1))[0];                                     //Scalar of type double
    b=REAL(VECTOR_ELT(gparms, 2))[0];                                       //Scalar of type double
    ro=REAL(VECTOR_ELT(gparms, 3))[0];                                      //Scalar of type double
    xiH=REAL(VECTOR_ELT(gparms, 4))[0];                                     //Scalar of type double
    xiV=REAL(VECTOR_ELT(gparms, 5))[0];                                     //Scalar of type double
    phiCP=REAL(VECTOR_ELT(gparms, 6))[0];                                   //Scalar of type double
    season_p1=REAL(VECTOR_ELT(gparms, 7))[0];                               //Scalar of type double
    season_p2=REAL(VECTOR_ELT(gparms, 8))[0];                               //Scalar of type double
    dzetaCP=LOGICAL(VECTOR_ELT(gparms, 9))[0];                              //Boolean (logical), as integer
    seasonality=LOGICAL(VECTOR_ELT(gparms, 10))[0];                         //Boolean (logical), as integer
    double *p_betaVH=REAL(VECTOR_ELT(gparms, 11));                          //Vector of doubles
    double *p_betaVH_age_coef=REAL(VECTOR_ELT(gparms, 12));                 //Vector of doubles
    betaHV=REAL(VECTOR_ELT(gparms, 13))[0];                                 //Scalar of type double
    double *p_muHi=REAL(VECTOR_ELT(gparms, 14));                            //Vector of doubles
    
    double *p_theta_unvac=REAL(VECTOR_ELT(gparms, 15));                     //3D array of doubles; array theta for unvaccinated individuals
    double *p_psy_DENV1_unvac=REAL(VECTOR_ELT(gparms, 16));                 //3D array of doubles; array psy for unvaccinated individuals
    double *p_psy_DENV2_unvac=REAL(VECTOR_ELT(gparms, 17));                 //3D array of doubles
    double *p_psy_DENV3_unvac=REAL(VECTOR_ELT(gparms, 18));                 //3D array of doubles
    double *p_psy_DENV4_unvac=REAL(VECTOR_ELT(gparms, 19));                 //3D array of doubles
    
    int *p_past_infections=INTEGER(VECTOR_ELT(gparms, 20));                 //3D array of integers
    
    double *p_theta_vac_neg=REAL(VECTOR_ELT(gparms,21));                   //5D array of doubles; array theta for hosts vaccinated as seronegative
    double *p_theta_vac_pos=REAL(VECTOR_ELT(gparms,22));                   //5D array of doubles; array theta for hosts vaccinated as seropositive
    
    double *p_psy_DENV1_vac_neg=REAL(VECTOR_ELT(gparms,23));               //4D array of doubles; array psy for serostatus & vaccination, by serotype
    double *p_psy_DENV2_vac_neg=REAL(VECTOR_ELT(gparms,24));
    double *p_psy_DENV3_vac_neg=REAL(VECTOR_ELT(gparms,25));
    double *p_psy_DENV4_vac_neg=REAL(VECTOR_ELT(gparms,26));
    
    double *p_psy_DENV1_vac_pos=REAL(VECTOR_ELT(gparms,27)); 
    double *p_psy_DENV2_vac_pos=REAL(VECTOR_ELT(gparms,28));
    double *p_psy_DENV3_vac_pos=REAL(VECTOR_ELT(gparms,29));
    double *p_psy_DENV4_vac_pos=REAL(VECTOR_ELT(gparms,30));
    
    
    /* Get all elements of vectors */
    int i,s;
    FOR_I muHi[i]=p_muHi[i];
    FOR_S betaVH[s]=p_betaVH[s];
    FOR_I betaVH_age_coef[i]=p_betaVH_age_coef[i];
    
    /* Get all elements of arrays */
    int j,k,l;
    FOR_J FOR_K FOR_L past_infections[j][k][l]=p_past_infections[G_ELT2(j,k,l)];                   //3D array past_infections; dimensions [NN][NN][NN]
    
    FOR_I FOR_J FOR_K FOR_L psy_DENV1_unvac[i][j][k][l]=p_psy_DENV1_unvac[G_ELT6(i,j,k,l)];        //4D array psy_DENV1_unvac for unvaccinated individuals 
    FOR_I FOR_J FOR_K FOR_L psy_DENV2_unvac[i][j][k][l]=p_psy_DENV2_unvac[G_ELT6(i,j,k,l)];        
    FOR_I FOR_J FOR_K FOR_L psy_DENV3_unvac[i][j][k][l]=p_psy_DENV3_unvac[G_ELT6(i,j,k,l)];
    FOR_I FOR_J FOR_K FOR_L psy_DENV4_unvac[i][j][k][l]=p_psy_DENV4_unvac[G_ELT6(i,j,k,l)];
    
    FOR_I FOR_J FOR_K FOR_L psy_DENV1_vac_neg[i][j][k][l]=p_psy_DENV1_vac_neg[G_ELT6(i,j,k,l)];    //4D arrays psy for vaccinated individuals
    FOR_I FOR_J FOR_K FOR_L psy_DENV2_vac_neg[i][j][k][l]=p_psy_DENV2_vac_neg[G_ELT6(i,j,k,l)];    //Dimensions [N_AGE_GROUPS][NN][NN][NN]
    FOR_I FOR_J FOR_K FOR_L psy_DENV3_vac_neg[i][j][k][l]=p_psy_DENV3_vac_neg[G_ELT6(i,j,k,l)];    //These arrays are calculated in R code and provided to C code
    FOR_I FOR_J FOR_K FOR_L psy_DENV4_vac_neg[i][j][k][l]=p_psy_DENV4_vac_neg[G_ELT6(i,j,k,l)]; 
    
    FOR_I FOR_J FOR_K FOR_L psy_DENV1_vac_pos[i][j][k][l]=p_psy_DENV1_vac_pos[G_ELT6(i,j,k,l)];
    FOR_I FOR_J FOR_K FOR_L psy_DENV2_vac_pos[i][j][k][l]=p_psy_DENV2_vac_pos[G_ELT6(i,j,k,l)];
    FOR_I FOR_J FOR_K FOR_L psy_DENV3_vac_pos[i][j][k][l]=p_psy_DENV3_vac_pos[G_ELT6(i,j,k,l)];
    FOR_I FOR_J FOR_K FOR_L psy_DENV4_vac_pos[i][j][k][l]=p_psy_DENV4_vac_pos[G_ELT6(i,j,k,l)];
    
    
    FOR_J FOR_K FOR_L theta_unvac[j][k][l]=p_theta_unvac[G_ELT2(j,k,l)];                             //3D array theta for unvaccinated individuals
    FOR_I FOR_S FOR_J FOR_K FOR_L theta_vac_neg[i][s][j][k][l]=p_theta_vac_neg[G_ELT5(i,s,j,k,l)];   //5D array theta for vaccinated individuals
    FOR_I FOR_S FOR_J FOR_K FOR_L theta_vac_pos[i][s][j][l][k]=p_theta_vac_pos[G_ELT5(i,s,j,k,l)];
  }
}

/********* MODEL FUNCTION ****************************************************************/
/*Function for the infection process, which contains all the core calculations including the ODEs. The signature of the function should not be modified. 
 Function arguments:
- neq represents the number of equations (compartments) in the model, it is informed by the ODE solver
- t is the current time being evaluated
- y is a vector with the current size of each compartment. 
- ydot is a vector with the values of derivatives for each compartment. This is the main model output. 
- yout is a vector with additional outputs 
- ip is an array informed by the solver. Depending on the algorithm it can contain different elements

Note: The function does not return anything since the results are stored in the parameters yout and ydot via pointers
Note: The vectors used in this function have the following dimensions and structure:
- Vectors y, ydot include three categories (in the given order)
  - Host compartments (nb elements = N_STATES_H), 
  - Vector compartments (nb elements = 9),
  - Additional compartments used to record incidence (nb elements = N_age_groups*4*3*vac_levels, i.e. one compartment for each age group, serotype, type (primary, secondry or post-secondary) and vaccination status).
- Vector yout includes the following categories:
  - Host population size (nb elements = 1 for each age group and each vaccination status)
  - Vector population size (nb elements = 1)
  - Seroprevalence (nb elements = 8 for each age group)
*/

void model(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
  if (ip[0]<1) error("nout should be at least 1");								                  //The first element of ip is the number of additional outputs (argument nout)
  
  /****** Initialise ydot and yout **********************************/
  int a;
  for(a=0; a<(*neq); a++) ydot[a]=0;
  for(a=0; a<(ip[0]); a++) yout[a]=0;
  
  /****** Define pointers to split each vector **********************/
  /*Pointers for the vector y (current size of each compartment)*/
  /*The vector has three sections - one per vaccination status - unvaccinated, vaccinated as seronegative, vaccinated as seropositive*/
  /*Each section of the vector has N_STATES_H elements (one per compartment)*/
  double *stateH_unvac=&y[0];                    //Pointer to the first element of a sub-vector with the size of compartments for unvaccinated hosts
  double *stateH_vac_neg=&y[N_STATES_H];         //...for hosts vaccinated as seronegative (goes after compartments for unvaccinated); not adding 1 as the indexing started with 0
  double *stateH_vac_pos=&y[(N_STATES_H*2)];     //...for hosts vaccinated as seropositive (goes after compartments for vaccinated as seronegative)

  double *StateV=&y[(N_STATES_H*3)];            //Pointer to the first element of the vector compartments (goes after all host compartments)

  /*Pointers for the vector ydot (compartment derivatives)*/
  /*The vector structure is the same as for vector y*/
  double *dH_unvac=&ydot[0];                    //Pointers to the first element of a sub-vector with derivatives for unvaccinated hosts
  double *dH_vac_neg=&ydot[N_STATES_H];         //...for hosts vaccinated as seronegative (goes after compartments for unvaccinated)
  double *dH_vac_pos=&ydot[(N_STATES_H*2)];     //...for hosts vaccinated as seropositive (goes after compartments for vaccinated as seronegative)
  
  double *dV=&ydot[(N_STATES_H*3)];             //Pointer to the first element of a sub-vector with derivatives for vector compartments (goes after all host compartments)
  
  /*Pointers for the vector user to record incidence (appended to the vector y and, hence, ydot; starts after all the host and vector compartments, i.e. at the element (N_STATES_H*3+9))*/
  /*The vector is split into five sections - one per vaccination status*/
  /*Each section contains elements to record the nb of cases by age, serotype and type, i.e. total of N_AGE_GROUPS*NS*3 elements*/
  double *d_incidence_unvac=&ydot[(N_STATES_H*3+9)];                              // Unvaccinated
  double *d_incidence_vac_neg=&ydot[(N_STATES_H*3+9+N_AGE_GROUPS*NS*3)];          // Vaccinated as seronegative
  double *d_incidence_vac_pos=&ydot[(N_STATES_H*3+9+N_AGE_GROUPS*NS*3*2)];        // Vaccinated as seropositive

  /*Pointers for the vector yout (additional output)*/
  /*The vector yout has several sections - to record the size of host population with each vac status, to record the total size of vector population and the seroprevalence (by age and serostatus x 8)*/
  double *NHit_unvac=&yout[0];                          //Pointer to the first element of a sub-vector with the population size for unvaccinated hosts (nb elements = N_AGE_GROUPS)
  double *NHit_vac_neg=&yout[N_AGE_GROUPS];             //...for hosts vaccinated as seronegative
  double *NHit_vac_pos=&yout[N_AGE_GROUPS*2];           //...for hosts vaccinated as seropositive

  double *NV_record=&yout[(N_AGE_GROUPS*3)];            //Pointer to the first element of a sub-vector with the population size for vector (nb elements = 1)
  
  double *seroprevalence=&yout[(N_AGE_GROUPS*3+1)];     //Pointer to the first element with seroprevalence (comes after host and vector population)
  
  
  /*Population size, by age, for each vaccination status (sum of the compartments in the corresponding stateH_...)*/
  /*Each vector has length N_age_groups */
  sum_overi_StateH(stateH_unvac,   NHit_unvac);       //Unvaccinated (sum of stateH_unvac)
  sum_overi_StateH(stateH_vac_neg, NHit_vac_neg);     //Vaccinated as seronegative (sum of stateH_vac_neg)
  sum_overi_StateH(stateH_vac_pos, NHit_vac_pos);     //Vaccinated as seropositive (sum of stateH_vac_pos)
  
  /*Population size, by age (sum of the vectors calculated above)*/
  /*Each vector has length N_age_groups*/
  /*Function Vadd allows summing the corresponding elements of two vectors of the same length*/
  
  double NHit_vac[N_AGE_GROUPS]={0};             //Total size of vaccinated host population, by age
  Vadd(NHit_vac_neg, NHit_vac_pos, N_AGE_GROUPS, NHit_vac);   

  double NHit[N_AGE_GROUPS]={0};                  //Total size of host population, by age
  Vadd(NHit_unvac,NHit_vac,N_AGE_GROUPS,NHit);
  
  /*Host population size, total*/
  double NH=Vsum(NHit, N_AGE_GROUPS);             
  
  /*Vector population size, total (sum of corresponding compartments)*/
  (*NV_record)=Vsum(StateV,9);                    

    
  /****** Frequently used variables *********************************/
  int i,j,k,l,m;                     //Counters for loops
  const int not_dzetaCP=(1-dzetaCP); //Inverse of dzetaCP
  
  
  /****** Force of infection for vectors ****************************/
  /*
  The computation of the force of infection for vectors, for each serotype, is comprised of several steps:
  1) Multiply the number of hosts in each compartment by value of infectiousness for this compartment (i.e. by the corresponding element of the array psy).
  Here, infectiousness represents the relative contribution of a host with given characteristics into the force of infection. This step is done separately for each vaccination status
  2) Sum, for each age group, the number of infectious hosts with each vaccination status (obtained at step 1)
  3.1) Sum across all age groups
  3.2) Divide the obtained number by the total host population size to obtain the overall proportion of infectious hosts (i.e. the probability that a bite falls on an infectious host)
  3.3) Multiply the obtained number by the probability of virus transmission from an infectious host to a susceptible vector given a bite (betaHV)
  3.4) Multiply the obtained number by the daily biting rate of an adult female vector to obtain the force of infection
  */ 
 
  /*Initialise objects*/
  double lambV[4];                                             //Final vector with four values of FoI for vectors (one for each serotype)
  double temp_result_step1_1[N_AGE_GROUPS][NN][NN][NN];        //Nb of infectious unvaccinated hosts, weighted by relative infectiousness (by age and status for the three serotypes that are not considered in the current infection process)
  double temp_result_step1_2[N_AGE_GROUPS][NN][NN][NN];        //Same for hosts vaccinated as seronegative
  double temp_result_step1_3[N_AGE_GROUPS][NN][NN][NN];        //Same for hosts vaccinated as seropositive
 
  double temp_result_step2[N_AGE_GROUPS];
  
  /* FoI for DENV-1 */
  //Here, pos_ind=1 as serotype in question is DENV-1; val_ind=2 as it is the value indicating the status "Infected" (Corresponds to 3 in R codes) 
  Vprod_psy_StateH(psy_DENV1_unvac,   stateH_unvac,  1,2, temp_result_step1_1); 
  Vprod_psy_StateH(psy_DENV1_vac_neg, stateH_vac_neg,1,2, temp_result_step1_2);   
  Vprod_psy_StateH(psy_DENV1_vac_pos, stateH_vac_pos,1,2, temp_result_step1_3);   

  sum_overi_psy_StatesH(temp_result_step1_1,temp_result_step1_2, temp_result_step1_3, temp_result_step2);     //Note: Aggregating for each age group (before aggregating the overall number) is left-over from older code, when each age group had a different probability of being bitten
  lambV[0]=(betaHV*b/NH)*Vsum(temp_result_step2, N_AGE_GROUPS); 
  
  /* FoI for DENV-2 */
  //Here, pos_ind=2 as serotype in question is DENV-1; val_ind=2 as it is the value indicating the status "Infected" (Corresponds to 3 in R codes) 
  Vprod_psy_StateH(psy_DENV2_unvac,   stateH_unvac,  2,2, temp_result_step1_1); 
  Vprod_psy_StateH(psy_DENV2_vac_neg, stateH_vac_neg,2,2, temp_result_step1_2);
  Vprod_psy_StateH(psy_DENV2_vac_pos, stateH_vac_pos,2,2, temp_result_step1_3);   

  sum_overi_psy_StatesH(temp_result_step1_1,temp_result_step1_2,temp_result_step1_3, temp_result_step2); 
  lambV[1]=(betaHV*b/NH)*Vsum(temp_result_step2, N_AGE_GROUPS); 
  
  /* FoI for DENV-3 */
  //Here, pos_ind=3 as serotype in question is DENV-1; val_ind=2 as it is the value indicating the status "Infected" (Corresponds to 3 in R codes) 
  Vprod_psy_StateH(psy_DENV3_unvac,   stateH_unvac,  3,2, temp_result_step1_1); 
  Vprod_psy_StateH(psy_DENV3_vac_neg, stateH_vac_neg,3,2, temp_result_step1_2);
  Vprod_psy_StateH(psy_DENV3_vac_pos, stateH_vac_pos,3,2, temp_result_step1_3);   
  
  sum_overi_psy_StatesH(temp_result_step1_1,temp_result_step1_2,temp_result_step1_3, temp_result_step2); 
  lambV[2]=(betaHV*b/NH)*Vsum(temp_result_step2, N_AGE_GROUPS); 
  
  /* FoI for DENV-4 */
  //Here, pos_ind=4 as serotype in question is DENV-1; val_ind=2 as it is the value indicating the status "Infected" (Corresponds to 3 in R codes) 
  Vprod_psy_StateH(psy_DENV4_unvac,   stateH_unvac,  4,2, temp_result_step1_1); 
  Vprod_psy_StateH(psy_DENV4_vac_neg, stateH_vac_neg,4,2, temp_result_step1_2);
  Vprod_psy_StateH(psy_DENV4_vac_pos, stateH_vac_pos,4,2, temp_result_step1_3);   
  
  sum_overi_psy_StatesH(temp_result_step1_1,temp_result_step1_2,temp_result_step1_3, temp_result_step2); 
  lambV[3]=(betaHV*b/NH)*Vsum(temp_result_step2, N_AGE_GROUPS); 
  
  
  /****** ODE for vectors *******************************************/
  
  /*Calculate seasonality factor for the current t*/
  /*This seasonality coefficient is applied to the vector recruitment rate */
  double season_coefficient=season_p1*sin(((*t)+1)*2*Pi/365.0+season_p2); 
  if(season_coefficient<-1) season_coefficient=-1;                             /*Check that the coefficient is not smaller than -1 to avoid negative recruitment rate */
  
  /*Calculate the number of new recruited vectors*/
  /*The number of new recruited vectors is calculated so that the total number of vectors is equal to the target number (given the seasonality coefficient). This includes compensating for the vector mortality*/
  double BV=ratio_VH*NH*(1+seasonality*season_coefficient)+((*NV_record)*muV)-(*NV_record);
  
  /*Calculate the total exit rate from the exposed compartment*/
  double sum_muV_xiV=muV+xiV; 
  
  /*Calculate derivatives for the vector compartments*/
  dV[0]=BV-(Vsum(lambV, 4)+muV)*StateV[0];        //Susceptible vectors
  dV[1]=lambV[0]*StateV[0]-sum_muV_xiV*StateV[1]; //Vectors exposed to DENV1
  dV[2]=lambV[1]*StateV[0]-sum_muV_xiV*StateV[2]; //Vectors exposed to DENV2
  dV[3]=lambV[2]*StateV[0]-sum_muV_xiV*StateV[3]; //Vectors exposed to DENV3
  dV[4]=lambV[3]*StateV[0]-sum_muV_xiV*StateV[4]; //Vectors exposed to DENV4
  dV[5]=xiV*StateV[1]-muV*StateV[5];              //Vectors infected/infectious with DENV1
  dV[6]=xiV*StateV[2]-muV*StateV[6];              //Vectors infected/infectious with DENV2
  dV[7]=xiV*StateV[3]-muV*StateV[7];              //Vectors infected/infectious with DENV3
  dV[8]=xiV*StateV[4]-muV*StateV[8];              //Vectors infected/infectious with DENV4
  
  /****** Force of infection for hosts ******************************/
  double lambH[N_AGE_GROUPS][NS];                    //Final vector of FOI for hosts (one value per serotype & per age)
  
  FOR_I{
    lambH[i][0]=(betaVH[0]*betaVH_age_coef[i]*b*StateV[5])/NH; //DENV-1 (for age group i; taking into account the age-specific relative risk applied to betaVH for this age group)
    lambH[i][1]=(betaVH[1]*betaVH_age_coef[i]*b*StateV[6])/NH; //DENV-2 
    lambH[i][2]=(betaVH[2]*betaVH_age_coef[i]*b*StateV[7])/NH; //DENV-3 
    lambH[i][3]=(betaVH[3]*betaVH_age_coef[i]*b*StateV[8])/NH; //DENV-4 
  }
  
  /****** Infection process for hosts *******************************/
  /*Initialise vectors containing derivatives for each compartment, by serotype & vaccination status*/
  /*Unvaccinated hosts*/
  double (*I1_unvac)[NN][NN][NN][NN];                  // DENV1 - Define pointer
  I1_unvac = calloc(N_AGE_GROUPS, sizeof * I1_unvac);  // Allocate and zero-init memory
  double (*I2_unvac)[NN][NN][NN][NN];                  // DENV2 - Define pointer
  I2_unvac = calloc(N_AGE_GROUPS, sizeof * I2_unvac);  // Allocate and zero-init memory
  double (*I3_unvac)[NN][NN][NN][NN];                  // DENV3 - Define pointer
  I3_unvac = calloc(N_AGE_GROUPS, sizeof * I3_unvac);  // Allocate and zero-init memory
  double (*I4_unvac)[NN][NN][NN][NN];                  // DENV4 - Define pointer
  I4_unvac = calloc(N_AGE_GROUPS, sizeof * I4_unvac);  // Allocate and zero-init memory  
  
  /*Hosts vaccinated as seronegative */
  double (*I1_vac_neg)[NN][NN][NN][NN];                   
  I1_vac_neg = calloc(N_AGE_GROUPS, sizeof * I1_vac_neg); 
  double (*I2_vac_neg)[NN][NN][NN][NN];                   
  I2_vac_neg = calloc(N_AGE_GROUPS, sizeof * I2_vac_neg); 
  double (*I3_vac_neg)[NN][NN][NN][NN];                   
  I3_vac_neg = calloc(N_AGE_GROUPS, sizeof * I3_vac_neg); 
  double (*I4_vac_neg)[NN][NN][NN][NN];                   
  I4_vac_neg = calloc(N_AGE_GROUPS, sizeof * I4_vac_neg); 
  	
  /*Hosts vaccinated as seropositive */
  double (*I1_vac_pos)[NN][NN][NN][NN];                   
  I1_vac_pos = calloc(N_AGE_GROUPS, sizeof * I1_vac_pos); 
  double (*I2_vac_pos)[NN][NN][NN][NN];                   
  I2_vac_pos = calloc(N_AGE_GROUPS, sizeof * I2_vac_pos); 
  double (*I3_vac_pos)[NN][NN][NN][NN];                   
  I3_vac_pos = calloc(N_AGE_GROUPS, sizeof * I3_vac_pos); 
  double (*I4_vac_pos)[NN][NN][NN][NN];                   
  I4_vac_pos = calloc(N_AGE_GROUPS, sizeof * I4_vac_pos); 
	
  /*Initialise variables used in the calculations*/
  int elt_i0klm, elt_i1klm, elt_i2klm, elt_i3klm; 
  int elt_ij0lm, elt_ij1lm, elt_ij2lm, elt_ij3lm;
  int elt_ijk0m, elt_ijk1m, elt_ijk2m, elt_ijk3m; 
  int elt_ijkl0, elt_ijkl1, elt_ijkl2, elt_ijkl3; 
  
  /* Initialise product of theta and FOI, for each vaccination status and each serotype */
  /*Unvaccinated hosts*/
  double theta_unvac_klm_lambH0;
  double theta_unvac_jlm_lambH1;
  double theta_unvac_jkm_lambH2;
  double theta_unvac_jkl_lambH3;
  
  /*Hosts vaccinated as seronegative */
  double theta_vac_neg_klm_lambH0;
  double theta_vac_neg_jlm_lambH1;
  double theta_vac_neg_jkm_lambH2;
  double theta_vac_neg_jkl_lambH3;
  
  /*Hosts vaccinated as seropositive */
  double theta_vac_pos_klm_lambH0;
  double theta_vac_pos_jlm_lambH1;
  double theta_vac_pos_jkm_lambH2;
  double theta_vac_pos_jkl_lambH3;
  
  FOR_I{
    
    /** Infection process for DENV-1 **/
    FOR_K FOR_L FOR_M {
      
      /*Helper functions*/
      elt_i0klm=G_ELT1(i,0,k,l,m);
      elt_i1klm=G_ELT1(i,1,k,l,m);
      elt_i2klm=G_ELT1(i,2,k,l,m);
      elt_i3klm=G_ELT1(i,3,k,l,m);
      
      /** Unvaccinated hosts **/
      /*Product of theta and FoI for the serotype in question*/
      theta_unvac_klm_lambH0=theta_unvac[k][l][m]*lambH[i][0]; 
      
      /*ODE*/
      I1_unvac[i][0][k][l][m]=-theta_unvac_klm_lambH0*stateH_unvac[elt_i0klm];
      I1_unvac[i][1][k][l][m]=theta_unvac_klm_lambH0*stateH_unvac[elt_i0klm]-xiH*stateH_unvac[elt_i1klm];
      I1_unvac[i][2][k][l][m]=xiH*stateH_unvac[elt_i1klm]-ro*stateH_unvac[elt_i2klm];
      I1_unvac[i][3][k][l][m]=dzetaCP*(ro*stateH_unvac[elt_i2klm]-phiCP*stateH_unvac[elt_i3klm]);
      I1_unvac[i][4][k][l][m]=dzetaCP*phiCP*stateH_unvac[elt_i3klm]+not_dzetaCP*ro*stateH_unvac[elt_i2klm];
	  
      /** Hosts vaccinated as seronegative **/
      /*Product of theta and FoI for the serotype in question*/
      theta_vac_neg_klm_lambH0=theta_vac_neg[i][0][k][l][m]*lambH[i][0]; 
      
      /*ODE*/
      I1_vac_neg[i][0][k][l][m]=-theta_vac_neg_klm_lambH0*stateH_vac_neg[elt_i0klm];
      I1_vac_neg[i][1][k][l][m]=theta_vac_neg_klm_lambH0*stateH_vac_neg[elt_i0klm]-xiH*stateH_vac_neg[elt_i1klm];
      I1_vac_neg[i][2][k][l][m]=xiH*stateH_vac_neg[elt_i1klm]-ro*stateH_vac_neg[elt_i2klm];
      I1_vac_neg[i][3][k][l][m]=dzetaCP*(ro*stateH_vac_neg[elt_i2klm]-phiCP*stateH_vac_neg[elt_i3klm]);
      I1_vac_neg[i][4][k][l][m]=dzetaCP*phiCP*stateH_vac_neg[elt_i3klm]+not_dzetaCP*ro*stateH_vac_neg[elt_i2klm];      

      /** Hosts vaccinated as seropositive **/
      /*Product of theta and FoI for the serotype in question*/
      theta_vac_pos_klm_lambH0=theta_vac_pos[i][0][k][l][m]*lambH[i][0]; 
      
      /*ODE*/
      I1_vac_pos[i][0][k][l][m]=-theta_vac_pos_klm_lambH0*stateH_vac_pos[elt_i0klm];
      I1_vac_pos[i][1][k][l][m]=theta_vac_pos_klm_lambH0*stateH_vac_pos[elt_i0klm]-xiH*stateH_vac_pos[elt_i1klm];
      I1_vac_pos[i][2][k][l][m]=xiH*stateH_vac_pos[elt_i1klm]-ro*stateH_vac_pos[elt_i2klm];
      I1_vac_pos[i][3][k][l][m]=dzetaCP*(ro*stateH_vac_pos[elt_i2klm]-phiCP*stateH_vac_pos[elt_i3klm]);
      I1_vac_pos[i][4][k][l][m]=dzetaCP*phiCP*stateH_vac_pos[elt_i3klm]+not_dzetaCP*ro*stateH_vac_pos[elt_i2klm]; 
      
    }
    
    /** Infection process for DENV-2 **/
    FOR_J FOR_L FOR_M {
      /*Helper functions*/
      elt_ij0lm=G_ELT1(i,j,0,l,m);
      elt_ij1lm=G_ELT1(i,j,1,l,m);
      elt_ij2lm=G_ELT1(i,j,2,l,m);
      elt_ij3lm=G_ELT1(i,j,3,l,m);
      
      /** Unvaccinated hosts **/
      /*Product of theta and FoI for the serotype in question*/
      theta_unvac_jlm_lambH1=theta_unvac[j][l][m]*lambH[i][1];
      
      /*ODE*/
      I2_unvac[i][j][0][l][m]=-theta_unvac_jlm_lambH1*stateH_unvac[elt_ij0lm];
      I2_unvac[i][j][1][l][m]=theta_unvac_jlm_lambH1*stateH_unvac[elt_ij0lm]-xiH*stateH_unvac[elt_ij1lm];
      I2_unvac[i][j][2][l][m]=xiH*stateH_unvac[elt_ij1lm]-ro*stateH_unvac[elt_ij2lm];
      I2_unvac[i][j][3][l][m]=dzetaCP*(ro*stateH_unvac[elt_ij2lm]-phiCP*stateH_unvac[elt_ij3lm]);
      I2_unvac[i][j][4][l][m]=dzetaCP*phiCP*stateH_unvac[elt_ij3lm]+not_dzetaCP*ro*stateH_unvac[elt_ij2lm];
      
   	  /** Hosts vaccinated as seronegative **/
      /*Product of theta and FoI for the serotype in question*/
      theta_vac_neg_jlm_lambH1=theta_vac_neg[i][1][j][l][m]*lambH[i][1];  
      
      /*ODE*/
      I2_vac_neg[i][j][0][l][m]=-theta_vac_neg_jlm_lambH1*stateH_vac_neg[elt_ij0lm];
      I2_vac_neg[i][j][1][l][m]=theta_vac_neg_jlm_lambH1*stateH_vac_neg[elt_ij0lm]-xiH*stateH_vac_neg[elt_ij1lm];
      I2_vac_neg[i][j][2][l][m]=xiH*stateH_vac_neg[elt_ij1lm]-ro*stateH_vac_neg[elt_ij2lm];
      I2_vac_neg[i][j][3][l][m]=dzetaCP*(ro*stateH_vac_neg[elt_ij2lm]-phiCP*stateH_vac_neg[elt_ij3lm]);
      I2_vac_neg[i][j][4][l][m]=dzetaCP*phiCP*stateH_vac_neg[elt_ij3lm]+not_dzetaCP*ro*stateH_vac_neg[elt_ij2lm];
	  
      /** Hosts vaccinated as seropositive **/
      /*Product of theta and FoI for the serotype in question*/
      theta_vac_pos_jlm_lambH1=theta_vac_pos[i][1][j][l][m]*lambH[i][1];
      
      /*ODE*/
      I2_vac_pos[i][j][0][l][m]=-theta_vac_pos_jlm_lambH1*stateH_vac_pos[elt_ij0lm];
      I2_vac_pos[i][j][1][l][m]=theta_vac_pos_jlm_lambH1*stateH_vac_pos[elt_ij0lm]-xiH*stateH_vac_pos[elt_ij1lm];
      I2_vac_pos[i][j][2][l][m]=xiH*stateH_vac_pos[elt_ij1lm]-ro*stateH_vac_pos[elt_ij2lm];
      I2_vac_pos[i][j][3][l][m]=dzetaCP*(ro*stateH_vac_pos[elt_ij2lm]-phiCP*stateH_vac_pos[elt_ij3lm]);
      I2_vac_pos[i][j][4][l][m]=dzetaCP*phiCP*stateH_vac_pos[elt_ij3lm]+not_dzetaCP*ro*stateH_vac_pos[elt_ij2lm];	  
      
    }
    
    /** Infection process for DENV-3 **/
    FOR_J FOR_K FOR_M {
      /*Helper functions*/
      elt_ijk0m=G_ELT1(i,j,k,0,m);
      elt_ijk1m=G_ELT1(i,j,k,1,m);
      elt_ijk2m=G_ELT1(i,j,k,2,m);
      elt_ijk3m=G_ELT1(i,j,k,3,m);
      
      /** Unvaccinated hosts **/
      /*Product of theta and FoI for the serotype in question*/
      theta_unvac_jkm_lambH2=theta_unvac[j][k][m]*lambH[i][2];
      
      /*ODE*/
      I3_unvac[i][j][k][0][m]=-theta_unvac_jkm_lambH2*stateH_unvac[elt_ijk0m];
      I3_unvac[i][j][k][1][m]=theta_unvac_jkm_lambH2*stateH_unvac[elt_ijk0m]-xiH*stateH_unvac[elt_ijk1m];
      I3_unvac[i][j][k][2][m]=xiH*stateH_unvac[elt_ijk1m]-ro*stateH_unvac[elt_ijk2m];
      I3_unvac[i][j][k][3][m]=dzetaCP*(ro*stateH_unvac[elt_ijk2m]-phiCP*stateH_unvac[elt_ijk3m]);
      I3_unvac[i][j][k][4][m]=dzetaCP*phiCP*stateH_unvac[elt_ijk3m]+not_dzetaCP*ro*stateH_unvac[elt_ijk2m];
      
      /** Hosts vaccinated as seronegative **/
      /*Product of theta and FoI for the serotype in question*/
      theta_vac_neg_jkm_lambH2=theta_vac_neg[i][2][j][k][m]*lambH[i][2];
      
      /*ODE*/
      I3_vac_neg[i][j][k][0][m]=-theta_vac_neg_jkm_lambH2*stateH_vac_neg[elt_ijk0m];
      I3_vac_neg[i][j][k][1][m]=theta_vac_neg_jkm_lambH2*stateH_vac_neg[elt_ijk0m]-xiH*stateH_vac_neg[elt_ijk1m];
      I3_vac_neg[i][j][k][2][m]=xiH*stateH_vac_neg[elt_ijk1m]-ro*stateH_vac_neg[elt_ijk2m];
      I3_vac_neg[i][j][k][3][m]=dzetaCP*(ro*stateH_vac_neg[elt_ijk2m]-phiCP*stateH_vac_neg[elt_ijk3m]);
      I3_vac_neg[i][j][k][4][m]=dzetaCP*phiCP*stateH_vac_neg[elt_ijk3m]+not_dzetaCP*ro*stateH_vac_neg[elt_ijk2m];

     /** Hosts vaccinated as seropositive **/
      /*Product of theta and FoI for the serotype in question*/
      theta_vac_pos_jkm_lambH2=theta_vac_pos[i][2][j][k][m]*lambH[i][2];
      
      /*ODE*/
      I3_vac_pos[i][j][k][0][m]=-theta_vac_pos_jkm_lambH2*stateH_vac_pos[elt_ijk0m];
      I3_vac_pos[i][j][k][1][m]=theta_vac_pos_jkm_lambH2*stateH_vac_pos[elt_ijk0m]-xiH*stateH_vac_pos[elt_ijk1m];
      I3_vac_pos[i][j][k][2][m]=xiH*stateH_vac_pos[elt_ijk1m]-ro*stateH_vac_pos[elt_ijk2m];
      I3_vac_pos[i][j][k][3][m]=dzetaCP*(ro*stateH_vac_pos[elt_ijk2m]-phiCP*stateH_vac_pos[elt_ijk3m]);
      I3_vac_pos[i][j][k][4][m]=dzetaCP*phiCP*stateH_vac_pos[elt_ijk3m]+not_dzetaCP*ro*stateH_vac_pos[elt_ijk2m];
	  
    }
    
    /** Infection process for DENV-4 **/
    FOR_J FOR_K FOR_L {
      /*Helper functions*/
      elt_ijkl0=G_ELT1(i,j,k,l,0);
      elt_ijkl1=G_ELT1(i,j,k,l,1);
      elt_ijkl2=G_ELT1(i,j,k,l,2);
      elt_ijkl3=G_ELT1(i,j,k,l,3);
      
      /** Unvaccinated hosts **/
      /*Product of theta and FoI for the serotype in question*/
      theta_unvac_jkl_lambH3=theta_unvac[j][k][l]*lambH[i][3];
      
      /*ODE*/
      I4_unvac[i][j][k][l][0]=-theta_unvac_jkl_lambH3*stateH_unvac[elt_ijkl0];
      I4_unvac[i][j][k][l][1]=theta_unvac_jkl_lambH3*stateH_unvac[elt_ijkl0]-xiH*stateH_unvac[elt_ijkl1];
      I4_unvac[i][j][k][l][2]=xiH*stateH_unvac[elt_ijkl1]-ro*stateH_unvac[elt_ijkl2];
      I4_unvac[i][j][k][l][3]=dzetaCP*(ro*stateH_unvac[elt_ijkl2]-phiCP*stateH_unvac[elt_ijkl3]);
      I4_unvac[i][j][k][l][4]=dzetaCP*phiCP*stateH_unvac[elt_ijkl3]+not_dzetaCP*ro*stateH_unvac[elt_ijkl2];
      
      /** Hosts vaccinated as seronegative **/   
      /*Product of theta and FoI for the serotype in question*/
      theta_vac_neg_jkl_lambH3=theta_vac_neg[i][3][j][k][l]*lambH[i][3];
      
      /*ODE*/
      I4_vac_neg[i][j][k][l][0]=-theta_vac_neg_jkl_lambH3*stateH_vac_neg[elt_ijkl0];
      I4_vac_neg[i][j][k][l][1]=theta_vac_neg_jkl_lambH3*stateH_vac_neg[elt_ijkl0]-xiH*stateH_vac_neg[elt_ijkl1];
      I4_vac_neg[i][j][k][l][2]=xiH*stateH_vac_neg[elt_ijkl1]-ro*stateH_vac_neg[elt_ijkl2];
      I4_vac_neg[i][j][k][l][3]=dzetaCP*(ro*stateH_vac_neg[elt_ijkl2]-phiCP*stateH_vac_neg[elt_ijkl3]);
      I4_vac_neg[i][j][k][l][4]=dzetaCP*phiCP*stateH_vac_neg[elt_ijkl3]+not_dzetaCP*ro*stateH_vac_neg[elt_ijkl2];
       
      /** Hosts vaccinated as seropositive **/ 
      /*Product of theta and FoI for the serotype in question*/
      theta_vac_pos_jkl_lambH3=theta_vac_pos[i][3][j][k][l]*lambH[i][3];
      
      /*ODE*/
      I4_vac_pos[i][j][k][l][0]=-theta_vac_pos_jkl_lambH3*stateH_vac_pos[elt_ijkl0];
      I4_vac_pos[i][j][k][l][1]=theta_vac_pos_jkl_lambH3*stateH_vac_pos[elt_ijkl0]-xiH*stateH_vac_pos[elt_ijkl1];
      I4_vac_pos[i][j][k][l][2]=xiH*stateH_vac_pos[elt_ijkl1]-ro*stateH_vac_pos[elt_ijkl2];
      I4_vac_pos[i][j][k][l][3]=dzetaCP*(ro*stateH_vac_pos[elt_ijkl2]-phiCP*stateH_vac_pos[elt_ijkl3]);
      I4_vac_pos[i][j][k][l][4]=dzetaCP*phiCP*stateH_vac_pos[elt_ijkl3]+not_dzetaCP*ro*stateH_vac_pos[elt_ijkl2];	  
       
    }
  }
  
  /** Sum of all infection processes for hosts **/
  FOR_I FOR_J FOR_K FOR_L FOR_M {
    dH_unvac[G_ELT1(i,j,k,l,m)]=I1_unvac[i][j][k][l][m] + I2_unvac[i][j][k][l][m] + I3_unvac[i][j][k][l][m] + I4_unvac[i][j][k][l][m] ;
  }
  free(I1_unvac);
  free(I2_unvac);
  free(I3_unvac);
  free(I4_unvac);
  
  FOR_I FOR_J FOR_K FOR_L FOR_M {
    dH_vac_neg[G_ELT1(i,j,k,l,m)]=I1_vac_neg[i][j][k][l][m] + I2_vac_neg[i][j][k][l][m] + I3_vac_neg[i][j][k][l][m] + I4_vac_neg[i][j][k][l][m] ;
  }
  free(I1_vac_neg);
  free(I2_vac_neg);
  free(I3_vac_neg);
  free(I4_vac_neg); 
  
  FOR_I FOR_J FOR_K FOR_L FOR_M {
    dH_vac_pos[G_ELT1(i,j,k,l,m)]=I1_vac_pos[i][j][k][l][m] + I2_vac_pos[i][j][k][l][m] + I3_vac_pos[i][j][k][l][m] + I4_vac_pos[i][j][k][l][m] ;
  }
  free(I1_vac_pos);
  free(I2_vac_pos);
  free(I3_vac_pos);
  free(I4_vac_pos);

  
    /****** Computation of incidence compartments *********************/
  /* Incidence is recorded for each age group, serotype and infection type (primary, secondary, post-secondary)*/
  /* Thus the length of the array is N_AGE_GROUPS*NS*3 */
  /* The recorded number of cases is cumulative and is obtained by adding the number of hosts moving from the state "Exposed" to the state "Infectious"*/
  FOR_J FOR_K FOR_L {                                         //Here j, k and l represent the status for the three serotypes other than the serotype in question
    switch(past_infections[j][k][l]){                         //Reminder: object past_infections contains the number of previous infections for the three serotypes
    case 0:  //First infection (primary)
      /* Unvaccinated */
      FOR_I d_incidence_unvac[G_ELT3(i,0,0)]+=xiH*stateH_unvac[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence_unvac[G_ELT3(i,1,0)]+=xiH*stateH_unvac[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence_unvac[G_ELT3(i,2,0)]+=xiH*stateH_unvac[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence_unvac[G_ELT3(i,3,0)]+=xiH*stateH_unvac[G_ELT1(i,j,k,l,1)];
      /* Vaccinated as seronegative*/    
      FOR_I d_incidence_vac_neg[G_ELT3(i,0,0)]+=xiH*stateH_vac_neg[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence_vac_neg[G_ELT3(i,1,0)]+=xiH*stateH_vac_neg[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence_vac_neg[G_ELT3(i,2,0)]+=xiH*stateH_vac_neg[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence_vac_neg[G_ELT3(i,3,0)]+=xiH*stateH_vac_neg[G_ELT1(i,j,k,l,1)]; 
      /* Vaccinated as seropositive*/    
      FOR_I d_incidence_vac_pos[G_ELT3(i,0,0)]+=xiH*stateH_vac_pos[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence_vac_pos[G_ELT3(i,1,0)]+=xiH*stateH_vac_pos[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence_vac_pos[G_ELT3(i,2,0)]+=xiH*stateH_vac_pos[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence_vac_pos[G_ELT3(i,3,0)]+=xiH*stateH_vac_pos[G_ELT1(i,j,k,l,1)]; 
      break;
    case 1:  //Second infection (secondary)
      /* Unvaccinated */
      FOR_I d_incidence_unvac[G_ELT3(i,0,1)]+=xiH*stateH_unvac[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence_unvac[G_ELT3(i,1,1)]+=xiH*stateH_unvac[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence_unvac[G_ELT3(i,2,1)]+=xiH*stateH_unvac[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence_unvac[G_ELT3(i,3,1)]+=xiH*stateH_unvac[G_ELT1(i,j,k,l,1)];
      /* Vaccinated as seronegative*/    
      FOR_I d_incidence_vac_neg[G_ELT3(i,0,1)]+=xiH*stateH_vac_neg[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence_vac_neg[G_ELT3(i,1,1)]+=xiH*stateH_vac_neg[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence_vac_neg[G_ELT3(i,2,1)]+=xiH*stateH_vac_neg[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence_vac_neg[G_ELT3(i,3,1)]+=xiH*stateH_vac_neg[G_ELT1(i,j,k,l,1)]; 
      /* Vaccinated as seropositive*/    
      FOR_I d_incidence_vac_pos[G_ELT3(i,0,1)]+=xiH*stateH_vac_pos[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence_vac_pos[G_ELT3(i,1,1)]+=xiH*stateH_vac_pos[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence_vac_pos[G_ELT3(i,2,1)]+=xiH*stateH_vac_pos[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence_vac_pos[G_ELT3(i,3,1)]+=xiH*stateH_vac_pos[G_ELT1(i,j,k,l,1)]; 
      break;
    case 2: //Third infection (post-secondary)
      /* Unvaccinated */
      FOR_I d_incidence_unvac[G_ELT3(i,0,2)]+=xiH*stateH_unvac[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence_unvac[G_ELT3(i,1,2)]+=xiH*stateH_unvac[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence_unvac[G_ELT3(i,2,2)]+=xiH*stateH_unvac[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence_unvac[G_ELT3(i,3,2)]+=xiH*stateH_unvac[G_ELT1(i,j,k,l,1)];
      /* Vaccinated as seronegative*/    
      FOR_I d_incidence_vac_neg[G_ELT3(i,0,2)]+=xiH*stateH_vac_neg[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence_vac_neg[G_ELT3(i,1,2)]+=xiH*stateH_vac_neg[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence_vac_neg[G_ELT3(i,2,2)]+=xiH*stateH_vac_neg[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence_vac_neg[G_ELT3(i,3,2)]+=xiH*stateH_vac_neg[G_ELT1(i,j,k,l,1)]; 
      /* Vaccinated as seropositive*/    
      FOR_I d_incidence_vac_pos[G_ELT3(i,0,2)]+=xiH*stateH_vac_pos[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence_vac_pos[G_ELT3(i,1,2)]+=xiH*stateH_vac_pos[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence_vac_pos[G_ELT3(i,2,2)]+=xiH*stateH_vac_pos[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence_vac_pos[G_ELT3(i,3,2)]+=xiH*stateH_vac_pos[G_ELT1(i,j,k,l,1)]; 
      break;
    case 3: //Fourth infection (post-secondary)
      /* Unvaccinated */
      FOR_I d_incidence_unvac[G_ELT3(i,0,2)]+=xiH*stateH_unvac[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence_unvac[G_ELT3(i,1,2)]+=xiH*stateH_unvac[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence_unvac[G_ELT3(i,2,2)]+=xiH*stateH_unvac[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence_unvac[G_ELT3(i,3,2)]+=xiH*stateH_unvac[G_ELT1(i,j,k,l,1)];
      /* Vaccinated as seronegative*/    
      FOR_I d_incidence_vac_neg[G_ELT3(i,0,2)]+=xiH*stateH_vac_neg[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence_vac_neg[G_ELT3(i,1,2)]+=xiH*stateH_vac_neg[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence_vac_neg[G_ELT3(i,2,2)]+=xiH*stateH_vac_neg[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence_vac_neg[G_ELT3(i,3,2)]+=xiH*stateH_vac_neg[G_ELT1(i,j,k,l,1)]; 
      /* Vaccinated as seropositive*/    
      FOR_I d_incidence_vac_pos[G_ELT3(i,0,2)]+=xiH*stateH_vac_pos[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence_vac_pos[G_ELT3(i,1,2)]+=xiH*stateH_vac_pos[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence_vac_pos[G_ELT3(i,2,2)]+=xiH*stateH_vac_pos[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence_vac_pos[G_ELT3(i,3,2)]+=xiH*stateH_vac_pos[G_ELT1(i,j,k,l,1)]; 
      break;
    default:
      break;
    }
  }
  
  /****** Computation of seroprevalence *****************************/
  /*Seroprevalence for hosts is recorded for each  age group and for each of the following statuses:
  0) Susceptible to all four serotypes
  1) History of one infection (DENV-1)
  2) History of one infection (DENV-2)
  3) History of one infection (DENV-3)
  4) History of one infection (DENV-4)
  5) History of two infections
  6) History of three infections
  7) History of four infections
  */
  
  int o1, o2, o3, o4;
  FOR_I{
    seroprevalence[G_ELT4(i,0)]+=stateH_unvac[G_ELT1(i,0,0,0,0)]+stateH_vac_neg[G_ELT1(i,0,0,0,0)]+stateH_vac_pos[G_ELT1(i,0,0,0,0)];
    
    /*History of one infection (by serotype) */
    for(o1=1; o1<5; o1++){
      seroprevalence[G_ELT4(i,1)]+=stateH_unvac[G_ELT1(i,o1,0,0,0)]+stateH_vac_neg[G_ELT1(i,o1,0,0,0)]+stateH_vac_pos[G_ELT1(i,o1,0,0,0)]; //DENV1
      seroprevalence[G_ELT4(i,2)]+=stateH_unvac[G_ELT1(i,0,o1,0,0)]+stateH_vac_neg[G_ELT1(i,0,o1,0,0)]+stateH_vac_pos[G_ELT1(i,0,o1,0,0)]; //DENV2
      seroprevalence[G_ELT4(i,3)]+=stateH_unvac[G_ELT1(i,0,0,o1,0)]+stateH_vac_neg[G_ELT1(i,0,0,o1,0)]+stateH_vac_pos[G_ELT1(i,0,0,o1,0)]; //DENV3
      seroprevalence[G_ELT4(i,4)]+=stateH_unvac[G_ELT1(i,0,0,0,o1)]+stateH_vac_neg[G_ELT1(i,0,0,0,o1)]+stateH_vac_pos[G_ELT1(i,0,0,0,o1)]; //DENV4
      
      /*History of two infections */
      for(o2=1; o2<5; o2++){ 
        seroprevalence[G_ELT4(i,5)]+=
          stateH_unvac[G_ELT1(i,o1,o2,0,0)]+
          stateH_unvac[G_ELT1(i,o1,0,o2,0)]+
          stateH_unvac[G_ELT1(i,o1,0,0,o2)]+
          stateH_unvac[G_ELT1(i,0,o1,o2,0)]+
          stateH_unvac[G_ELT1(i,0,o1,0,o2)]+
          stateH_unvac[G_ELT1(i,0,0,o1,o2)]+
          stateH_vac_neg[G_ELT1(i,o1,o2,0,0)]+stateH_vac_pos[G_ELT1(i,o1,o2,0,0)]+
          stateH_vac_neg[G_ELT1(i,o1,0,o2,0)]+stateH_vac_pos[G_ELT1(i,o1,0,o2,0)]+
          stateH_vac_neg[G_ELT1(i,o1,0,0,o2)]+stateH_vac_pos[G_ELT1(i,o1,0,0,o2)]+
          stateH_vac_neg[G_ELT1(i,0,o1,o2,0)]+stateH_vac_pos[G_ELT1(i,0,o1,o2,0)]+
          stateH_vac_neg[G_ELT1(i,0,o1,0,o2)]+stateH_vac_pos[G_ELT1(i,0,o1,0,o2)]+
          stateH_vac_neg[G_ELT1(i,0,0,o1,o2)]+stateH_vac_pos[G_ELT1(i,0,0,o1,o2)];
        for(o3=1; o3<5; o3++){ //3 infections
          seroprevalence[G_ELT4(i,6)]+=
            stateH_unvac[G_ELT1(i,o1,o2,o3,0)]+
            stateH_unvac[G_ELT1(i,o1,o2,0,o3)]+
            stateH_unvac[G_ELT1(i,o1,0,o2,o3)]+
            stateH_unvac[G_ELT1(i,0,o1,o2,o3)]+
            stateH_vac_neg[G_ELT1(i,o1,o2,o3,0)]+stateH_vac_pos[G_ELT1(i,o1,o2,o3,0)]+
            stateH_vac_neg[G_ELT1(i,o1,o2,0,o3)]+stateH_vac_pos[G_ELT1(i,o1,o2,0,o3)]+
            stateH_vac_neg[G_ELT1(i,o1,0,o2,o3)]+stateH_vac_pos[G_ELT1(i,o1,0,o2,o3)]+
            stateH_vac_neg[G_ELT1(i,0,o1,o2,o3)]+stateH_vac_pos[G_ELT1(i,0,o1,o2,o3)];
          for(o4=1; o4<5; o4++){ // 4 infections
            seroprevalence[G_ELT4(i,7)]+= stateH_unvac[G_ELT1(i,o1,o2,o3,o4)]+stateH_vac_neg[G_ELT1(i,o1,o2,o3,o4)]+stateH_vac_pos[G_ELT1(i,o1,o2,o3,o4)];
          }
        }
      }
    }
  }
  
}
