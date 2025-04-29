/*
 * Project: Dengue dynamic model
 * Author(s): Elizaveta KHARITONOVA, Anna TYTULA, Aurelien JAMOTTE, Olivier CRISTEAU
 */

/*=====================================================================================================================*/
/*MODEL IDENTIFIERS*/
#ifndef M1_NO_VAC
#define M1_NO_VAC

/*MAIN CHARACTERISTICS*/
#define N_PARMS 21                     //Nb. of parameters passed to the C code (list parms_C in R code)
#define N_AGE_GROUPS 101               //Nb. of age groups used to run the dynamic model
#define NN 5                           //Nb. of possible immunological statuses for each serotype
#define NS 4                           //Nb. of serotypes
#define Pi 3.14159265358979323846      //Value of Pi (used in the definition of seasonality)

/*ALIASES FOR LOOPS*/
#define FOR_S for(s = 0; s < NS; s++) 
#define FOR_I for(i = 0; i < N_AGE_GROUPS; i++)
#define FOR_J for(j = 0; j < NN; j++) 
#define FOR_K for(k = 0; k < NN; k++) 
#define FOR_L for(l = 0; l < NN; l++) 
#define FOR_M for(m = 0; m < NN; m++) 

/*NB. OF HOST COMPARTMENTS*/
const int N_STATES_H=N_AGE_GROUPS*NN*NN*NN*NN;  


/*=====================================================================================================================*/
/*DECLARATION OF MODEL PARAMETERS*/
/*Here the model parameters are listed in the same order as in the list parms_C created in R code*/

static double ratio_VH;
static double muV;
static double b;
static double ro;
static double xiH;
static double xiV;
static double phiCP;
static double season_p1;
static double season_p2;
static int dzetaCP;
static int seasonality;
static double betaVH[NS];
static double betaVH_age_coef[N_AGE_GROUPS];
static double betaHV; 
static double muHi[N_AGE_GROUPS];
static double theta[NN][NN][NN];
static double psy_DENV1[N_AGE_GROUPS][NN][NN][NN];
static double psy_DENV2[N_AGE_GROUPS][NN][NN][NN];
static double psy_DENV3[N_AGE_GROUPS][NN][NN][NN];
static double psy_DENV4[N_AGE_GROUPS][NN][NN][NN];
static int past_infections[NN][NN][NN];


/*=====================================================================================================================*/
/*SIGNATURES OF THE FUNCTIONS USED IN THE CORRESPONDING .C FILE*/

/*Position of an element in the vector StateH*/
#define G_ELT1 get_elt_StateH                             
int get_elt_StateH(int i,int j, int k, int l, int m);     

/*Position of an element in a 3D array with dimensions NN x NN x NN*/
#define G_ELT2 get_elt_3D_serotype 
int get_elt_3D_serotype(int j, int k, int l);             

/*Position of an element in the incidence vector*/
#define G_ELT3 get_elt_incidence 
int get_elt_incidence(int i, int s, int status);          

/*Position of an element in the seroprevalence vector*/
#define G_ELT4 get_elt_seroprevalence 
int get_elt_seroprevalence(int i, int status);    

/*Position of an element in an array psy for unvaccinated hosts (dimensions N_AGE_GROUPS x NN x NN x NN) */
#define G_ELT5 get_elt_psy_vac
int get_elt_psy_vac(int i, int j, int k, int l);

/*Sum of the elements of a vector*/
double Vsum(double vector1[], int length_vector);         

/*Product of the corresponding elements of two vectors of same length*/
void Vprod(double vector1[], double vector2[], int length_vector, double result[]);

/*Sum of the corresponding elements of two vectors of same length*/
void Vadd(double vector1[], double vector2[], int length_vector, double result[]);

/*Difference of the corresponding elements of two vectors of same length*/
void Vdiff(double vector1[], double vector2[], int length_vector, double result[]); 

/*Sum of the products of the corresponding elements of two vectors of same length*/
double Vsum_prod(double vector1[], double vector2[], int length_vector); 

/*Sum of all compartments for each age group*/
void sum_overi_StateH(double StateH[], double NHi[]); 


/*Multiplication of two objects:
* -3D array (in a vector form) where each element is the value of a parameter for a given combination of statuses for three serotype (for example, array psy)
* -Vector where each element is the nb. of hosts of a given age i with the given statuses for each serotype j,k,l,m

*Here, pos_ind indicates the serotype in question (different in each serotype-specific infection process)
*      val_ind indicates the immunological status for the serotype in question

*The function returns a 4D array (in a vector form) where each element if the value of a parameter for a given age and three immunological statuses*/

void Vprod_psy_StateH(double psy[N_AGE_GROUPS][NN][NN][NN], double StateH[], int pos_ind, int val_ind, double result[][NN][NN][NN]);

/*Sum of the 4D array (in a vector form) described above, for each age group*/
void sum_overi_psy_StateH(double psy_stateH[][NN][NN][NN], double result[]);

#endif

