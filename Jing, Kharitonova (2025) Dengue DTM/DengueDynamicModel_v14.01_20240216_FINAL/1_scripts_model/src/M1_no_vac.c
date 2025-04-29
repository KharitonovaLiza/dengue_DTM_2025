/*
 * Project: Dengue dynamic model
 * Author(s): Elizaveta KHARITONOVA, Anna TYTULA, Aurelien JAMOTTE, Olivier CRISTEAU
 * Description: This file contains the source code for the infection process without vaccination used in the calibration (composite age groups)
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>
#include <stdio.h>
#include "M1_no_vac.h"   //Header file


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

/*Multiplication of two objects - 3D array and a vector */
/*pos_int - serotype in question; val_ind - status for the serotype in question */
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
    
    ratio_VH=REAL(VECTOR_ELT(gparms, 0))[0];                //Scalar of type double
    muV=REAL(VECTOR_ELT(gparms, 1))[0];                     //Scalar of type double
    b=REAL(VECTOR_ELT(gparms, 2))[0];                       //Scalar of type double
    ro=REAL(VECTOR_ELT(gparms, 3))[0];                      //Scalar of type double
    xiH=REAL(VECTOR_ELT(gparms, 4))[0];                     //Scalar of type double
    xiV=REAL(VECTOR_ELT(gparms, 5))[0];                     //Scalar of type double
    phiCP=REAL(VECTOR_ELT(gparms, 6))[0];                   //Scalar of type double
    season_p1=REAL(VECTOR_ELT(gparms, 7))[0];               //Scalar of type double
    season_p2=REAL(VECTOR_ELT(gparms, 8))[0];               //Scalar of type double
    dzetaCP=LOGICAL(VECTOR_ELT(gparms, 9))[0];              //Boolean (logical), as integer
    seasonality=LOGICAL(VECTOR_ELT(gparms, 10))[0];         //Boolean (logical), as integer
    double *p_betaVH=REAL(VECTOR_ELT(gparms, 11));          //Vector of doubles
    double *p_betaVH_age_coef=REAL(VECTOR_ELT(gparms, 12)); //Vector of doubles
    betaHV=REAL(VECTOR_ELT(gparms, 13))[0];                 //Scalar of type double
    double *p_muHi=REAL(VECTOR_ELT(gparms, 14));            //Vector of doubles
    double *p_theta=REAL(VECTOR_ELT(gparms, 15));           //3D array of doubles
    double *p_psy_DENV1=REAL(VECTOR_ELT(gparms, 16));             //3D array of doubles
    double *p_psy_DENV2=REAL(VECTOR_ELT(gparms, 17));             //3D array of doubles
    double *p_psy_DENV3=REAL(VECTOR_ELT(gparms, 18));             //3D array of doubles
    double *p_psy_DENV4=REAL(VECTOR_ELT(gparms, 19));             //3D array of doubles
    int *p_past_infections=INTEGER(VECTOR_ELT(gparms, 20)); //3D array of integers
    
    
    /* Get all elements of vectors */
    int i,s;
    FOR_I muHi[i]=p_muHi[i];
    FOR_S betaVH[s]=p_betaVH[s];
    FOR_I betaVH_age_coef[i]=p_betaVH_age_coef[i];
    
    /* Get all elements of arrays */
    int j,k,l;
    FOR_J FOR_K FOR_L theta[j][k][l]=p_theta[G_ELT2(j,k,l)];  
    FOR_I FOR_J FOR_K FOR_L psy_DENV1[i][j][k][l]=p_psy_DENV1[G_ELT5(i,j,k,l)];  
    FOR_I FOR_J FOR_K FOR_L psy_DENV2[i][j][k][l]=p_psy_DENV2[G_ELT5(i,j,k,l)]; 
    FOR_I FOR_J FOR_K FOR_L psy_DENV3[i][j][k][l]=p_psy_DENV3[G_ELT5(i,j,k,l)]; 
    FOR_I FOR_J FOR_K FOR_L psy_DENV4[i][j][k][l]=p_psy_DENV4[G_ELT5(i,j,k,l)]; 
    FOR_J FOR_K FOR_L past_infections[j][k][l]=p_past_infections[G_ELT2(j,k,l)];  
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
Note: The vectors used in this function has the following dimensions and structure:
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
  if (ip[0]<1) error("nout should be equal to at least 1");       //The first element of ip is the number of additional outputs (argument nout)
  
  
  /****** Initialise ydot and yout **********************************/
  int a;
  for(a=0; a<(*neq); a++) ydot[a]=0;         
  for(a=0; a<(ip[0]); a++) yout[a]=0;        
  
  
  /****** Define pointers to split each vector **********************/
  /*Pointers for the vector y (current size of each compartment)*/
  double *StateH=&y[0];                      //Pointer to the first element of a sub-vector with the size of host compartments (nb elements in the sub-vector = N_STATES_H)
  double *StateV=&y[N_STATES_H];             //pointer to the first element of a sub-vector with the size of vector compartments (nb elements in the sub-vector = 9)
  
  
  /*Pointers for the vector ydot (compartment derivatives)*/
  double *dH=&ydot[0];                       //Pointer to the first element of a sub-vector with the derivatives of the host compartments (nb elements in the sub-vector = N_STATES_H)
  double *dV=&ydot[N_STATES_H];              //Pointer to the first element of a sub-vector with the derivatives of the vector compartments (nb elements in the sub-vector = 9)
  double *d_incidence=&ydot[(N_STATES_H+9)]; //Pointer to the first element of a sub-vector with the derivatives of the incidence "compartments" (i.e. nb of new infections)
  
  /*Pointers for the vector yout (additional outputs)*/
  double *NHit=&yout[0];                           //Pointer to the first element of a sub-vector with the host population size (nb elements = N_AGE_GROUPS)
  double *NV_record=&yout[N_AGE_GROUPS];           //Pointer to the first element of a sub-vector with the vector population size (nb elements = 1)
  double *seroprevalence=&yout[(N_AGE_GROUPS+1)];  //Pointer to the first element of a sub-vector with the seroprevalence (nb elements = 8*N_AGE_GROUPS)
  
  
  /****** Sum of the host and vector compartments *******************/
  sum_overi_StateH(StateH, NHit);      //Total number of hosts, by age
  double NH=Vsum(NHit, N_AGE_GROUPS);  //Total number of hosts
  (*NV_record)=Vsum(StateV,9);         //Total number of vectors
  
  
  /****** Frequently used variables *********************************/
  int i,j,k,l,m;                     //Counters for loops
  const int not_dzetaCP=(1-dzetaCP); //Inverse of dzetaCP
  
  
  /****** Force of infection for vectors ****************************/
  
  /*
  The computation of the force of infection for vectors, for each serotype, is comprised of several steps:
  1) Multiply the number of hosts in each compartment by value of infectiousness for this compartment (i.e. by the corresponding element of the array psy).
     Here, infectiousness represents the relative contribution of a host with given characteristics into the force of infection
  2) Sum the obtained values for each age group to obtain the total number of infectious hosts in each age group (weighed by their relative contribution into the force of infection)
     Here, the sums are done across age groups to make use of the pre-developed functions
  3.1) Sum across all age groups
  3.2) Divide the obtained number by the total host population size to obtain the overall proportion of infectious hosts (i.e. the probability that a bite falls on an infectious host)
  3.3) Multiply the obtained number by the probability of virus transmission from an infectious host to a susceptible vector given a bite (betaHV)
  3.4) Multiply the obtained number by the daily biting rate of an adult female vector to obtain the force of infection
  */
  
  /*Initialise objects */
  double lambV[4];                                       //Final vector with four values of FoI for vectors (one for each serotype)
  double temp_result_step1[N_AGE_GROUPS][NN][NN][NN];    //Intermediate results
  double temp_result_step2[N_AGE_GROUPS];
  
  /* FoI for DENV-1 */
  Vprod_psy_StateH(psy_DENV1, StateH, 1, 2, temp_result_step1);           //Here, pos_ind=1 as serotype in question is DENV-1; val_ind=2 as it is the value indicating the status "Infected" (Corresponds to 3 in R codes)
  sum_overi_psy_StateH(temp_result_step1, temp_result_step2); 
  lambV[0]=(betaHV*b/NH)*Vsum(temp_result_step2, N_AGE_GROUPS);  
  
  /* FoI for DENV-2 */
  Vprod_psy_StateH(psy_DENV2, StateH, 2, 2, temp_result_step1); 
  sum_overi_psy_StateH(temp_result_step1, temp_result_step2); 
  lambV[1]=(betaHV*b/NH)*Vsum(temp_result_step2, N_AGE_GROUPS);  
  
  /* FoI for DENV-3 */
  Vprod_psy_StateH(psy_DENV3, StateH, 3, 2, temp_result_step1); 
  sum_overi_psy_StateH(temp_result_step1, temp_result_step2); 
  lambV[2]=(betaHV*b/NH)*Vsum(temp_result_step2, N_AGE_GROUPS);  
  
  /* FoI for DENV-4 */
  Vprod_psy_StateH(psy_DENV4, StateH, 4, 2, temp_result_step1); 
  sum_overi_psy_StateH(temp_result_step1, temp_result_step2); 
  lambV[3]=(betaHV*b/NH)*Vsum(temp_result_step2, N_AGE_GROUPS);                               
  
  
  /****** ODE for vectors *******************************************/
  
  /*Calculate seasonality factor for the current t*/
  /*This seasonality coefficient is applied to the vector recruitment rate */
  double season_coefficient=season_p1*sin(((*t)+1)*2*Pi/365.0+season_p2); 
  if(season_coefficient<-1) season_coefficient=-1;                              /*Check that the coefficient is not smaller than -1 to avoid negative recruitment rate */
  
  /*Calculate the number of new recruited vectors*/
  /*The number of new recruited vectors is calculated so that the total number of vectors is equal to the target number (given the seasonality coefficient). THis includes compensating for the vector mortality*/
  double BV=ratio_VH*NH*(1+seasonality*season_coefficient)+((*NV_record)*muV)-(*NV_record);
  
  /*Calculate the total exit rate from the exposed compartment*/
  double sum_muV_xiV=muV+xiV; 
  
  /*Calculate derivatives for the vector compartments*/
  dV[0]=BV-(Vsum(lambV, 4)+muV)*StateV[0];        //Susceptible
  dV[1]=lambV[0]*StateV[0]-sum_muV_xiV*StateV[1]; //Exposed to DENV-1
  dV[2]=lambV[1]*StateV[0]-sum_muV_xiV*StateV[2]; //Exposed to DENV-2
  dV[3]=lambV[2]*StateV[0]-sum_muV_xiV*StateV[3]; //Exposed to DENV-3
  dV[4]=lambV[3]*StateV[0]-sum_muV_xiV*StateV[4]; //Exposed to DENV-4
  dV[5]=xiV*StateV[1]-muV*StateV[5];              //Infected with DENV-1
  dV[6]=xiV*StateV[2]-muV*StateV[6];              //Infected with DENV-2
  dV[7]=xiV*StateV[3]-muV*StateV[7];              //Infected with DENV-3
  dV[8]=xiV*StateV[4]-muV*StateV[8];              //Infected with DENV-4
  
  
  
  /****** Force of infection for hosts ******************************/
  double lambH[N_AGE_GROUPS][NS];                    //Final vector of FOI for hosts (one value per serotype & per age)
  
  FOR_I{
    lambH[i][0]=(betaVH[0]*betaVH_age_coef[i]*b*StateV[5])/NH; //DENV-1 (for age group i; taking into account the age-specific relative risk applied to betaVH for this age group)
    lambH[i][1]=(betaVH[1]*betaVH_age_coef[i]*b*StateV[6])/NH; //DENV-2 
    lambH[i][2]=(betaVH[2]*betaVH_age_coef[i]*b*StateV[7])/NH; //DENV-3 
    lambH[i][3]=(betaVH[3]*betaVH_age_coef[i]*b*StateV[8])/NH; //DENV-4 
  }
  
  
  /****** Infection process for hosts *******************************/
  /*Initialise vectors containing derivatives for each compartment, by serotype*/
  double I1[N_AGE_GROUPS][NN][NN][NN][NN]={0}; //DENV-1
  double I2[N_AGE_GROUPS][NN][NN][NN][NN]={0}; //DENV-2
  double I3[N_AGE_GROUPS][NN][NN][NN][NN]={0}; //DENV-3
  double I4[N_AGE_GROUPS][NN][NN][NN][NN]={0}; //DENV-4
  
  
  /*Initialise variables used in the calculations*/
  int elt_i0klm, elt_i1klm, elt_i2klm, elt_i3klm; 
  int elt_ij0lm, elt_ij1lm, elt_ij2lm, elt_ij3lm;
  int elt_ijk0m, elt_ijk1m, elt_ijk2m, elt_ijk3m; 
  int elt_ijkl0, elt_ijkl1, elt_ijkl2, elt_ijkl3; 
  
  double thetaklm_lambH0;     //Initialising the object, which is a product of theta and FoI for DENV-1 (and the current age group)
  double thetajlm_lambH1;     //Product of theta and FoI for DENV-2
  double thetajkm_lambH2;     //Product of theta and FoI for DENV-3
  double thetajkl_lambH3;     //Product of theta and FoI for DENV-4
  
  FOR_I{
    
    /** Infection process for DENV-1 **/
    
    FOR_K FOR_L FOR_M {
      
      /*Helper functions*/
      elt_i0klm=G_ELT1(i,0,k,l,m);  
      elt_i1klm=G_ELT1(i,1,k,l,m);
      elt_i2klm=G_ELT1(i,2,k,l,m);
      elt_i3klm=G_ELT1(i,3,k,l,m);
      
      /*Product of theta and FoI for DENV-1*/
      thetaklm_lambH0=theta[k][l][m]*lambH[i][0]; 
      
      /*ODE*/
      I1[i][0][k][l][m]=-thetaklm_lambH0*StateH[elt_i0klm];
      I1[i][1][k][l][m]=thetaklm_lambH0*StateH[elt_i0klm]-xiH*StateH[elt_i1klm];
      I1[i][2][k][l][m]=xiH*StateH[elt_i1klm]-ro*StateH[elt_i2klm];
      I1[i][3][k][l][m]=dzetaCP*(ro*StateH[elt_i2klm]-phiCP*StateH[elt_i3klm]);
      I1[i][4][k][l][m]=dzetaCP*phiCP*StateH[elt_i3klm]+not_dzetaCP*ro*StateH[elt_i2klm];
    }
    
    /** Infection process for DENV-2 **/
    FOR_J FOR_L FOR_M {
      
      /*Helper functions*/
      elt_ij0lm=G_ELT1(i,j,0,l,m);
      elt_ij1lm=G_ELT1(i,j,1,l,m);
      elt_ij2lm=G_ELT1(i,j,2,l,m);
      elt_ij3lm=G_ELT1(i,j,3,l,m);
      
      /*Product of theta and FoI for DENV-2*/
      thetajlm_lambH1=theta[j][l][m]*lambH[i][1];
      
      /*ODE*/
      I2[i][j][0][l][m]=-thetajlm_lambH1*StateH[elt_ij0lm];
      I2[i][j][1][l][m]=thetajlm_lambH1*StateH[elt_ij0lm]-xiH*StateH[elt_ij1lm];
      I2[i][j][2][l][m]=xiH*StateH[elt_ij1lm]-ro*StateH[elt_ij2lm];
      I2[i][j][3][l][m]=dzetaCP*(ro*StateH[elt_ij2lm]-phiCP*StateH[elt_ij3lm]);
      I2[i][j][4][l][m]=dzetaCP*phiCP*StateH[elt_ij3lm]+not_dzetaCP*ro*StateH[elt_ij2lm];
    }
    
    /** Infection process for DENV-3 **/
    FOR_J FOR_K FOR_M {
      
      /*Helper functions*/
      elt_ijk0m=G_ELT1(i,j,k,0,m);
      elt_ijk1m=G_ELT1(i,j,k,1,m);
      elt_ijk2m=G_ELT1(i,j,k,2,m);
      elt_ijk3m=G_ELT1(i,j,k,3,m);
      
      /*Product of theta and FoI for DENV-3*/
      thetajkm_lambH2=theta[j][k][m]*lambH[i][2];
      
      /*ODE*/
      I3[i][j][k][0][m]=-thetajkm_lambH2*StateH[elt_ijk0m];
      I3[i][j][k][1][m]=thetajkm_lambH2*StateH[elt_ijk0m]-xiH*StateH[elt_ijk1m];
      I3[i][j][k][2][m]=xiH*StateH[elt_ijk1m]-ro*StateH[elt_ijk2m];
      I3[i][j][k][3][m]=dzetaCP*(ro*StateH[elt_ijk2m]-phiCP*StateH[elt_ijk3m]);
      I3[i][j][k][4][m]=dzetaCP*phiCP*StateH[elt_ijk3m]+not_dzetaCP*ro*StateH[elt_ijk2m];
    }
    
    /** Infection process for DENV-4 **/
    FOR_J FOR_K FOR_L {
      
      /*Helper functions*/
      elt_ijkl0=G_ELT1(i,j,k,l,0);
      elt_ijkl1=G_ELT1(i,j,k,l,1);
      elt_ijkl2=G_ELT1(i,j,k,l,2);
      elt_ijkl3=G_ELT1(i,j,k,l,3);
      
      /*Product of theta and FoI for DENV-4*/
      thetajkl_lambH3=theta[j][k][l]*lambH[i][3];
      
      /*ODE*/
      I4[i][j][k][l][0]=-thetajkl_lambH3*StateH[elt_ijkl0];
      I4[i][j][k][l][1]=thetajkl_lambH3*StateH[elt_ijkl0]-xiH*StateH[elt_ijkl1];
      I4[i][j][k][l][2]=xiH*StateH[elt_ijkl1]-ro*StateH[elt_ijkl2];
      I4[i][j][k][l][3]=dzetaCP*(ro*StateH[elt_ijkl2]-phiCP*StateH[elt_ijkl3]);
      I4[i][j][k][l][4]=dzetaCP*phiCP*StateH[elt_ijkl3]+not_dzetaCP*ro*StateH[elt_ijkl2];
    }
  }
  
  /** Sum of all infectious processes for hosts **/
  FOR_I FOR_J FOR_K FOR_L FOR_M {
    dH[G_ELT1(i,j,k,l,m)]=I1[i][j][k][l][m] + I2[i][j][k][l][m] + I3[i][j][k][l][m] + I4[i][j][k][l][m];
  }
  
  /****** Computation of incidence compartments *********************/
  /* Incidence is recorded for each age group, serotype and infection type (primary, secondary, post-secondary)*/
  /* Thus the length of the array is N_AGE_GROUPS*NS*3 */
  
  FOR_J FOR_K FOR_L {                                                         //Here j, k and l represent the status for the three serotypes other than the serotype in question
    switch(past_infections[j][k][l]){                                         //Reminder: object past_infections contains the number of previous infections for the three serotypes
    case 0:  //First infection (primary)
      FOR_I d_incidence[G_ELT3(i,0,0)]+=xiH*StateH[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence[G_ELT3(i,1,0)]+=xiH*StateH[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence[G_ELT3(i,2,0)]+=xiH*StateH[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence[G_ELT3(i,3,0)]+=xiH*StateH[G_ELT1(i,j,k,l,1)];
      break;
    case 1:  //Second infection (secondary)
      FOR_I d_incidence[G_ELT3(i,0,1)]+=xiH*StateH[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence[G_ELT3(i,1,1)]+=xiH*StateH[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence[G_ELT3(i,2,1)]+=xiH*StateH[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence[G_ELT3(i,3,1)]+=xiH*StateH[G_ELT1(i,j,k,l,1)];
      break;
    case 2: //Third infection (post-secondary)
      FOR_I d_incidence[G_ELT3(i,0,2)]+=xiH*StateH[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence[G_ELT3(i,1,2)]+=xiH*StateH[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence[G_ELT3(i,2,2)]+=xiH*StateH[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence[G_ELT3(i,3,2)]+=xiH*StateH[G_ELT1(i,j,k,l,1)];
      break;
    case 3: //Fourth infection (post-secondary)
      FOR_I d_incidence[G_ELT3(i,0,2)]+=xiH*StateH[G_ELT1(i,1,j,k,l)];
      FOR_I d_incidence[G_ELT3(i,1,2)]+=xiH*StateH[G_ELT1(i,j,1,k,l)];
      FOR_I d_incidence[G_ELT3(i,2,2)]+=xiH*StateH[G_ELT1(i,j,k,1,l)];
      FOR_I d_incidence[G_ELT3(i,3,2)]+=xiH*StateH[G_ELT1(i,j,k,l,1)];
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
  
  The length of the seroprevalence array is N_AGE_GROUPS*8 */
  
  int o1, o2, o3, o4 ;
  FOR_I{
    seroprevalence[G_ELT4(i,0)]+=StateH[G_ELT1(i,0,0,0,0)];
    
    /*History of one infection */
    for(o1=1; o1<5; o1++){
      seroprevalence[G_ELT4(i,1)]+=StateH[G_ELT1(i,o1,0,0,0)]; //DENV-1
      seroprevalence[G_ELT4(i,2)]+=StateH[G_ELT1(i,0,o1,0,0)]; //DENV-2
      seroprevalence[G_ELT4(i,3)]+=StateH[G_ELT1(i,0,0,o1,0)]; //DENV-3
      seroprevalence[G_ELT4(i,4)]+=StateH[G_ELT1(i,0,0,0,o1)]; //DENV-4
      
      /*History of two infections */
      for(o2=1; o2<5; o2++){ 
        seroprevalence[G_ELT4(i,5)]+=
          StateH[G_ELT1(i,o1,o2,0,0)]+
          StateH[G_ELT1(i,o1,0,o2,0)]+
          StateH[G_ELT1(i,o1,0,0,o2)]+
          StateH[G_ELT1(i,0,o1,o2,0)]+
          StateH[G_ELT1(i,0,o1,0,o2)]+
          StateH[G_ELT1(i,0,0,o1,o2)];
        
        /*History of three infections */
        for(o3=1; o3<5; o3++){ 
          seroprevalence[G_ELT4(i,6)]+=
            StateH[G_ELT1(i,o1,o2,o3,0)]+
            StateH[G_ELT1(i,o1,o2,0,o3)]+
            StateH[G_ELT1(i,o1,0,o2,o3)]+
            StateH[G_ELT1(i,0,o1,o2,o3)];
          
          /*History of four infections */
          for(o4=1; o4<5; o4++){ 
            seroprevalence[G_ELT4(i,7)]+= StateH[G_ELT1(i,o1,o2,o3,o4)] ;
            
          }
        }
      }
    }            
  }
}