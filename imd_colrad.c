#include "imd.h"
#include <sys/time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

// ***************************************************
// * TODO:
// * mpi2 und mpi3 anpassen weil energies delE und kbTe in eV
// * SPEEDUP pragma omp   auch in ydot
// * floats wo möglich
// ****************************************************


#ifdef LOADBALANCE
#define node  l1[i]
#define node2 l2[i]
#else
#define node  l1[i][j][k]
#define node2 l2[i][j][k]
#endif


//#define USEFLOAT  // hauptsächlich in der funktion genexptint. Profiling zeigte, dass
                  // hier die meiste zeit verbraucht wird -> float verdoppelt performance

#ifdef USEFLOAT
typedef float realtype;
#define REALTYPE MPI_FLOAT
#else
typedef double realtype;
#define REALTYPE MPI_DOUBLE
#endif


#ifdef USEFLOAT
  #define EXPR     exp //exp zu floaten ist eine ganz mieeese idee
  #define SQRTR    sqrtf
  #define POWR     powf
  #define LOGR     logf
#else
  #define EXPR     exp
  #define SQRTR    sqrt
  #define POWR     pow
#define LOGR       log
#endif

#define OMP
#define LAPACK
//#define MULTIPHOTON
// #define SPONT  //<-- spontante emission, Kaum effekt 
//#define STARK  //<-- reabsorption via stark effect
#define DOIPD    //
#define RHOMIN 50

#ifdef OMP
#include <omp.h>
#endif

#define TIMING      //time-elapsed to stdout

#define MAXLEVEL 4 // bis zu welchem Ionisationsgrad?

// #ifdef USEFLOAT          //Dieser Versuch ist fehlgeschlagen...
// #define EXPR(x) expf(x)  //floats an wichtigen stellen erschweren 
// #else                    //konvergenz für solver
// #define EXPR(x) exp(x)
// #endif



const int MAX_LINE_LENGTH=3000; //wird von colrad_read benuthzt
                                //Hängt davon ab wieviele states betrachtet werden
                                //sollte eigentlich besser dynamisch berechnet werden


// **********************************************************************************************
// *          PAR INTEGRAL STUFF
// **********************************************************************************************
realtype  integral_simpson(realtype (*f)(realtype, void*), realtype a, realtype b,int n,void* p);

int simpson_error;
const realtype tolmax=1e-20;
#define INITIAL_STACK_SIZE 128   /* initial size of new stacks */



realtype  
integral_simpson_par(
     realtype (*f)(realtype,void*), /* function to integrate */
     realtype ah,           /* left interval boundary  */
     realtype bh,           /* right interval boundary */
     realtype tolh,
     realtype S,
     realtype fa,
     realtype fb,
     realtype fm,
     int rec,
     void* p); 


/* the stack structure */
struct stack_s{               
  int el_count;            /* count of elements on stack */    
  int el_size;             /* size of an element */
  int mem_reserve;         /* allocated memory for stack */
  void* elements;          /* pointer to begin of stack */
};

typedef struct _work_t{
  realtype a;
  realtype b;
  realtype tol;
  realtype S;
  realtype fa;
  realtype fb;
  realtype fm;
  realtype rec;
} work_t;


typedef struct stack_s* stack_t;


void create_stack(stack_t* stack, int element_size);
int  empty_stack(stack_t stack);
void push_stack(stack_t stack, void* element);
void pop_stack(stack_t stack, void* element);




// *********************************************************
// PHYSICAL CONSTANTS
// *********************************************************
// const double eV2J=1.6021766E-19;
const realtype eV2H=0.03674932; //eV to Hartree
const realtype colrad_reltol=1e-4;
const realtype colrad_abstol=10.0;

// const realtype J2eV=6.2415091E18;
const realtype planck=6.62607004E-34; // J/s
const realtype bohr_radius=0.52917721067E-10; // m
const realtype bohr_radius_sq=2.800285202924816e-21;
const realtype hbar_cub=1.172812163789953e-102; //hbar^3 
const realtype double_emass_pow_3_2 = 2.459112949719466e-45; // (2*emass)^3/2
const int MAXLINE = 255;
const realtype  pi=3.141592653589793;
const realtype  pi_sq=9.869604401089358;
const realtype E_ion_H=13.6; // eV
const realtype E_ion_H_J=2.178960176000000e-18; // J
const realtype E_ion_H_sq_J=4.747867448593952e-36;

const realtype colrad_tequi=1e-12;//TEST// 1e-12; //bei initial equi ohne Temperatur-variation erst einmal 
                                //die Saha-besetzungsdichten equilibrieren

//const double  LIGHTSPEED=2.997925458e8; // m/s
realtype  HBAR;
realtype  LASERFREQ;


int num_threads;
//const double  EMASS=9.10938356e-31;     // kg
//const double  ECONST=8.854187817e-12;   // As/Vm
//const double  BOLTZMAN=1.38064852e-23;  // J/K
//const double  ECHARGE=1.60217662e-19;  // C
//const double  AMU=1.66053904020e-27;   // atomic mass unit

// ******************************************************************************
// *                CROSS SECTION INTEGRATION STUFF
// ******************************************************************************

gsl_integration_workspace * winteg_inner=NULL;
gsl_integration_workspace * winteg_outer=NULL;
gsl_integration_workspace * winteg_fermi=NULL;
gsl_integration_workspace * winteg_exc=NULL; //excitation

// gsl_integration_romberg_workspace * winteg_rb_inner=NULL;
// gsl_integration_romberg_workspace * winteg_rb_outer=NULL;


struct my_f_params { realtype ne; realtype T;realtype mu; realtype E;realtype DeltaE; int allowed;};
struct my_f_params fparams_inner; //For inner integrand
struct my_f_params fparams_outer; //outer integrand
struct my_f_params fparams_fermi; 
struct my_f_params fparams_exc; 

realtype inner_integrand_ionization(realtype x, void *p); // integrate along E'
realtype outer_integrand_ionization(realtype x,void *p);  // integrate along E
realtype double_integral_ionization(realtype ne,realtype T, realtype mu, realtype DeltaE); //evaluates double integral


realtype inner_integrand_recombination(realtype x, void *p);
realtype outer_integrand_recombination(realtype x,void *p);
realtype double_integral_recombination(realtype ne,realtype T, realtype mu, realtype DeltaE);

realtype integrand_excitation(realtype x,void *p);
realtype eval_excitation_integral(realtype ne,realtype T,realtype mu, realtype DeltaE, int allowed); 

realtype eval_dexcitation_integral(realtype ne,realtype T,realtype mu, realtype DeltaE, int allowed);
realtype integrand_deexcitation(realtype x,void *p);


realtype fermi_integrand(realtype x, void *p);
realtype eval_fermi_integrand(realtype ne,realtype T, realtype mu);
const realtype integ_reltol = 1e-2;

const realtype integ_abstol = 1.0e-6; //1e-6; 
const realtype integ_abstol_ioniz= 1.0e-4; //-6; 
const realtype integ_abstol_recomb= 1.0e-4; //1e-6; 

const int    integ_meshdim =1500;

const realtype alpha_i =0.3; //0.05; //0.3;
const realtype beta_i  =0.9; //4.0 ; //0.9;
const realtype ioniz_const = 1.573949440579906e+71; // konstanten zus. gefasst und aus dem doppelintegral gezogen
const realtype recomb_const= 6.213703330335829e+72; // selbes für 3b-recomb.

#define muINF (50*eV2J) //ACHTUNG: Unebdingt was sinnvolleres finden (mu wird kleiner mit Te)

//wenn geschätzte max. ioniz.rate kleiner als das --> integrale garnicht erst berechnen 
//Schätzungen mit Hilfe von matlab-script DESCHAUD_TEST.m
#define MINRATE  1e16            // in 1/s/m^3 --> alles unter 10 Ionisationsereignissen pro fs pro m^3 wird ignoriert
                                 // habe ich also eine ionenkonz. von 10^26/m^3 und eine rate von  10/m^3/fs
                                 // so ändert sich pro md-step die konz. um den Bruchteil 10/1E26 = 1E-25
                                 // -->wayne

#define RATEIONIZMAX  1.5252e-13 // in m^3/s, das entspricht dem doppelintegral 
#define RATERECOMBMAX 6e-35 //6.1052e-42 // in m^6/s, das entspricht dem doppelintegral
#define MINCONC 1e-60   //im Saha-init sollen zu kleine konzentrationen ignoriert werden

realtype k_EE_MAX, k_EI_MAX, k_EE_REV_MAX, k_EI_REV_MAX; //DEBUG PURPOSE

// *********************************************************
//             CVODE-STRUCT FOR SOLVER PARAMS
// *********************************************************
typedef struct {
  realtype It; //Intesity
  realtype IPD0,IPD1,IPD2,IPD3,IPD4;
  realtype EF;
  realtype P_EE,P_EI,P_MPI2,P_MPI3,P_RR;
  realtype P_TOTAL; //komplette colrad-Leistungsdichte, für eng-filge
  realtype dens; //weil cv F(dens,temp) 
  realtype ni; //für Z=ne/ni
  bool initial_equi;
  realtype Tinit; //initial Temp. during equi must not change!
  realtype Heatrate; //um temp.sprünge zu vermeiden //TEST!
} *colrad_UserData;
colrad_UserData  cdata;

// *********************************************************
//                    MAIN
// *********************************************************
void do_colrad(realtype dt)
{

  k_EE_MAX= k_EI_MAX= k_EE_REV_MAX= k_EI_REV_MAX=0.0;  

  int flag;
  realtype t;
  realtype tout=dt;
  int i,j,k;
  N_Vector y;
  realtype Te0,Ti0,rho0,ni0,ne0;
  colrad_ptotal=0.0;

  if(myid==0 && cdata->initial_equi)
    printf("COLRAD performs pre-equilibration of Saha-distribution\nfor t=%.4e s...This may take some time.\n",colrad_tequi);

#ifdef TIMING
  //if(myid==0)
  // {
    struct timeval start, end;
    gettimeofday(&start, NULL);   
  // }
#endif


  for(i=1;i<local_fd_dim.x-1;i++)
  {

#ifndef LOADBALANCE
    for(j=1;j<local_fd_dim.y-1;j++)
    {
      for(k=1;k<local_fd_dim.z-1;k++)
      {        
#endif   

        //clear -> irgendwas stimmt nicht!
        // node.ne=node2.ne=0.0;
        // node.Z=node2.Z=0.0;
        // node.P_EE=node2.P_EE=0.0;
        // node.P_EI=node2.P_EI=0.0;
        // node.P_MPI2=node2.P_MPI2=0.0;
        // node.P_MPI3=node2.P_MPI3=0.0;
        // node.P_RR=node2.P_RR=0.0;

        if(node.natoms < fd_min_atoms || node.dens < RHOMIN) continue;
        y=node.y;

        Te0=node.temp*11604.5;
        Ti0=node.md_temp*11604.5;        
        rho0=node.dens;
        ni0=rho0/AMU/26.9185; //1e28; //1e26/m^3 entspricht etwa 1e-4/Angtrom^3   

        if(cdata->initial_equi==true)     
        {
          realtype Zmean=MeanCharge(Te0, rho0, atomic_charge, atomic_weight,i,j,k);          
          ne0= Zmean* node.dens / (atomic_weight * AMU);          
          node.ne=ne0; //Saha init greift darauf zu
          colrad_Saha_init(i, j, k);                 
        }
        else //NORMAL
        {
          ne0=node.ne;          
        }

        realtype Told=Ith(y,0); //DEBUG //letzter step, checke wie groß temp.sprung
        Ith(y,0)=Told; //TEST//Te0;

        cdata->Tinit=Told;

        Ith(y,1)=Ti0;
        Ith(y,2)=ne0;

        flag = CVodeReInit(cvode_mem, 0.0, y);                    

        
        if(cdata->initial_equi==true)
        {
//          printf("myid:%d, RUNNNING INITIAL EQUI\n",myid);
          flag = CVode(cvode_mem, colrad_tequi, y, &t, CV_NORMAL);    

          int i_global,j_global,k_global;

#ifndef LOADBALANCE
          i_global = ((i - 1) + my_coord.x * (local_fd_dim.x - 2));
          j_global = ((j - 1) + my_coord.y * (local_fd_dim.y - 2));
          k_global =  ((k-1) + my_coord.z*(local_fd_dim.z-2));          
#else
          i_global = ((i - 1) + myid * (local_fd_dim.x - 2));
          j_global=k_global=1;
#endif
long int nje;
long int nfe;
long int nsetups;
long int nni;
long int nst;
long int ncfn;
long int netf;
CVodeGetNumJacEvals(cvode_mem, &nje);
CVodeGetNumRhsEvals(cvode_mem, &nfe);
CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
CVodeGetNumSteps(cvode_mem, &nst);
CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
CVodeGetNumErrTestFails(cvode_mem, &netf);

printf("myid:%d, COLRAD Cell %d was equilibrated, nfe:%ld, nje:%ld, nsetups:%ld, nni:%ld,nst:%ld,ncfn:%ld,netf:%ld\n",
      myid,i_global, nfe,nje,nsetups,nni,nst,ncfn,netf);


        } //initial equi
        else //NORMAL
        {
          cdata->dens=node.dens;
          cdata->Tinit=Te0; 

// printf("myid:%d, COLRAD Cell %d begin step, Te0:%f, ne0:%.2e\n", myid,(i - 1) + myid * (local_fd_dim.x - 2),Te0, ne0);          
          cdata->Heatrate=(Te0-Told)/dt; //TEST

          flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
          colrad_ptotal+=(fd_vol)*1e-30*cdata->P_TOTAL; // d.h. ptotal ist Gesamt-Leisung

          int i_global;
#ifndef LOADBALANCE          
          i_global = ((i - 1) + my_coord.x * (local_fd_dim.x - 2));
#else
          i_global = ((i - 1) + myid * (local_fd_dim.x - 2));          
#endif          

long int nje;
long int nfe;
long int nsetups;
long int nni;
long int nst;
long int ncfn;
long int netf;
CVodeGetNumJacEvals(cvode_mem, &nje);
CVodeGetNumRhsEvals(cvode_mem, &nfe);
CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
CVodeGetNumSteps(cvode_mem, &nst);
CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
CVodeGetNumErrTestFails(cvode_mem, &netf);
if(myid==1 && i_global==18) //DEBUG
printf("myid:%d, COLRAD Cell %d step done,HR:%.2e, Told:%f, T0:%f, T1:%f, ne0:%.2e, ne1:%.2e, nfe:%ld, nje:%ld, nsetups:%ld, nni:%ld,nst:%ld,ncfn:%ld,netf:%ld\n",
      myid,i_global, cdata->Heatrate, Told, node.temp*11604.5, Ith(y,0), node.ne, Ith(y,2),
      nfe,nje,nsetups,nni,nst,ncfn,netf);

          //ni0=cdata->ni;
        }
 

        //REASSIGN NEW TE AND NE
        node.temp=Ith(y,0)/11604.5;
        node.ne=Ith(y,2);
        node.Z=node.ne/ni0;
        node.P_EE=cdata->P_EE;
        node.P_EI=cdata->P_EI;
        node.P_MPI2=cdata->P_MPI2;
        node.P_MPI3=cdata->P_MPI3;
        node.P_RR=cdata->P_RR;

        node2.temp=node.temp; //auch in l2 speichern! wichtig!
        node2.ne=node.ne;     //sonst geht alles im diffusion-step vor 
        node2.Z=node.Z;       //ttm-writeout verloren!
        node2.P_EE=node.P_EE;
        node2.P_EI=node.P_EI;
        node2.P_MPI2=node.P_MPI2;
        node2.P_MPI3=node.P_MPI3;
        node2.P_RR=node.P_RR;


        if(node.temp <0 || isnan(node.temp) !=0 )        
        {
          char errstr[255];
          sprintf(errstr,"ERROR in COLRAD: Te became Nan or <0\n");
          error(errstr);
        }
        if(node.ne <0 || isnan(node.ne) !=0 )        
        {
          char errstr[255];
          sprintf(errstr,"ERROR in COLRAD: ne became Nan or <0\n");
          error(errstr);
        }
#ifndef LOADBALANCE
      } // for k
    }  //for j
#endif    
  }
 if(cdata->initial_equi==true)
 {
  colrad_write(0);
  cdata->initial_equi=false;
  if(myid==0)
    printf("Initial equi done\n");
 }
 else if(steps % ttm_int ==0)
 {
   // colrad_write(steps);
 } 

 #ifdef TIMING
 MPI_Barrier(cpugrid);
 // if(myid==0)
 // {

    gettimeofday(&end, NULL);
    realtype delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
             end.tv_usec - start.tv_usec) / 1.e6;

    // if(myid==0)
    printf("myid:%d, telaps:%f, COLRAD STEP DONE, kee:%.4e,keeb:%.4e, kei:%.4e, k3b:%.4e\n",
      myid,delta, k_EE_MAX, k_EE_REV_MAX, k_EI_MAX, k_EI_REV_MAX);
 // }
 #endif
  
}

// *********************************************************
//                      INIT FUNC
// *********************************************************
void colrad_init(void)
{
  //Init gsl-integration stuff
  winteg_inner= gsl_integration_workspace_alloc (integ_meshdim); //Integration workspace allocation  
  winteg_outer= gsl_integration_workspace_alloc (integ_meshdim); //Integration workspace allocation  
  winteg_fermi= gsl_integration_workspace_alloc (integ_meshdim); //Integration workspace allocation  
  winteg_exc=   gsl_integration_workspace_alloc (integ_meshdim); 

  // winteg_rb_inner=gsl_integration_romberg_alloc(20);
  // winteg_rb_outer=gsl_integration_romberg_alloc(20);

	HBAR=planck/2.0/pi;
	LASERFREQ=LIGHTSPEED/lambda;
	SUNMatrix A;
	N_Vector abstol;	
	int i,j,k;
	


	if(myid==0)
	{
    printf("*****************************************\n");
    printf("*      COLLISIONAL RADIATIVE MODEL      *\n");
    printf("*****************************************\n");
    printf(" READING ENERGY LEVELS \n");	
    printf(" reltol:%.4e\n",colrad_reltol);
    printf(" abstol:%.4e\n",colrad_abstol);
	}
	colrad_read_states();
  if(myid==0)
  {
    printf(" Nr. of Equations:%d                  *\n",total_species+3);    
  }

MPI_Barrier(cpugrid);
#ifdef OMP
    #pragma omp parallel
    {
      #pragma omp single
      {
        num_threads=omp_get_num_threads();
        printf("myid:%d, omp threads:%d\n",myid,num_threads);
      }
    }
#endif 
if(myid==0)
{  
  printf("*****************************************\n");
}

  cdata=(colrad_UserData) malloc(sizeof *cdata); 
  cdata->initial_equi=true;
  cdata->P_TOTAL=0.0;

	//total_species=z0_len+z1_len; //wird bereits in read-states reduced
	neq=total_species+3;

	for(i=1;i<local_fd_dim.x-1;i++)
	{

#ifndef LOADBALANCE    
		for(j=1;j<local_fd_dim.y-1;j++)
		{
			for(k=1;k<local_fd_dim.z-1;k++)
			{
#endif        
				node.y=N_VNew_Serial(neq);
        node.P_EE=0.0;
        node.P_EI=0.0;
        node.P_MPI2=0.0;
        node.P_MPI3=0.0;
        node.P_RR=0.0;

				//l2... <--brauche nur 1 mal speicher alloc'en
				//aber in DIFF-LOOP darauf achten, dass beim swappen
				//l2.y immer auf l1.y zeigt und nicht ins Leere
#ifndef LOADBALANCE        
			}
		}
#endif

	}
	abstol = N_VNew_Serial(neq);
	N_Vector Vdummy=N_VNew_Serial(neq); //wird bei re-init ersetzt

	// *************************  //
	// *   TOLERANCES          *
	// *************************     
     	for(i=0;i<neq;i++)
     	{
       	  Ith(abstol,i) = colrad_abstol;//fmax(Ith(colrad_y,i)*1e-6,10.0);

     	}
      Ith(abstol,0)=5.0; //Temp

	// *************************  //
	// *   SOLVER & CVODE INIT *
	// *************************
    	cvode_mem=NULL;
    	cvode_mem = CVodeCreate(CV_BDF);
    	CVodeInit(cvode_mem, colrad_ydot, 0, Vdummy);

    	CVodeSetUserData(cvode_mem, cdata);
    	CVodeSVtolerances(cvode_mem, colrad_reltol, abstol);

   	SUNLinearSolver LS;
    	A=SUNDenseMatrix(neq, neq);
	#ifdef LAPACK
	LS = SUNLinSol_LapackDense(Vdummy, A);
	#else
	LS=SUNLinSol_Dense(Vdummy, A);
	#endif
	CVodeSetLinearSolver(cvode_mem, LS, A); 


	//////////////////////////////////////////////////////////////////////
	//            OPTIONAL FUNCTIONS: SIEHE S.42, TABLE 4.2
	//////////////////////////////////////////////////////////////////////
	//CVodeSetIterType(cvode_mem,CV_NEWTON); //ist ja eh schon gesetzt
	//Keine gute idee mindt zu setzen
	//CVodeSetMaxOrd(cvode_mem,12);  //Bei adams default=12	
	//CVodeSetMaxErrTestFails(cvode_mem,7); //Default=7
	//CVodeSetMaxNonlinIters(cvode_mem,3); //Default=3
	//CVodeSetMaxConvFails(cvode_mem,10); //Default=10
	//
	//CVodeSetMaxStep(cvode_mem,1e-15);
	//CVodeSetNonlinConvCoef(cvode_mem,0.01); //Default=0.1
	//Max Steps muss sehr hoch sein bei großen steps,
	//sonst beschwert sich newton
	CVodeSetMaxNumSteps(cvode_mem,2500); //Default=500
	// CVodeSetInitStep(cvode_mem,1e-17);  //ACHTUNG:Wenn zu klein --> BS
	//CVodeSetEpsLin(cvode_mem,0.01); //Default=0.05;

  CVodeSetMaxStepsBetweenJac(cvode_mem,50); //default 50 --> guter performance boost
  // CVodeSetMaxStepsBetweenJac(cvode_mem,50); //default 50 --> guter performance boost
	
	//N_VDestroy(dummy);	
	
/*
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  CVodeFree(&cvode_mem);	
*/
  
  //colrad_read(0);
}


void colrad_Saha_init(int i,int j,int k)
{  
  // int i,j,k;
  // for(i=1;i<local_fd_dim.x-1;i++)
  // {
  //   for(j=1;i<local_fd_dim.y-1;j++)
  //   {
  //     for(k=1;i<local_fd_dim.z-1;k++)
  //     {
        realtype Te0,Ti0,ne0,ni0,rho0;
        N_Vector y;
        y=node.y;

        Te0=node.temp*11604.5;
        Ti0=node.md_temp*11604.5;        
        rho0=node.dens; ///1e10;
        ni0=rho0/AMU/26.9185; //1e28; //1e26/m^3 entspricht etwa 1e-4/Angtrom^3        
        ne0=node.ne;

        Ith(y,0)=Te0;
        Ith(y,1)=Ti0;
        Ith(y,2)=ne0;

        do_Saha(Te0,ni0,ne0,y);
  //     }
  //   } 
  // }
}

void colrad_read_states(void)
{
  // **************************************************************************
  // *            READ STATES FILES
  // **********************************************************************
  int lcnt = 1;
  int i, j;
  realtype *buf; //buffer 1d array for communication
  FILE* fin = NULL;
  char line[255];

  if (myid == 0)
  {
    fin = fopen("Al0_states.txt", "r");

    if (fin == NULL)
    {
      char errstr[255];
      sprintf(errstr, "ERROR in colrad_read_states: File %s not found\n", "Al0_states.txt");
      error(errstr);

    }
    while (1) {
      if (fgets(line, MAXLINE, fin) == NULL) break;
      lcnt++;
    }
    alloc2darr(realtype, STATES_z0, lcnt, 6);

    lcnt = 0;
    rewind(fin);
    while (1)
    {
      if (fgets(line, MAXLINE, fin) == NULL) break;
      sscanf(line, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
             &STATES_z0[lcnt][0], &STATES_z0[lcnt][1], &STATES_z0[lcnt][2],
             &STATES_z0[lcnt][3], &STATES_z0[lcnt][4], &STATES_z0[lcnt][5]);
      lcnt++;
    }
    z0_len = lcnt;
    fclose(fin);
    total_species += z0_len;
  }


  //NOW COMMUNICATE
  MPI_Bcast(&z0_len, 1, MPI_INT, 0, cpugrid);
  alloc1darr(realtype, buf, z0_len * 6);

  MPI_Barrier(cpugrid);
  if (myid == 0)
  {
    for (i = 0; i < z0_len * 6; i += 6) //fill 1d buff-array
    {
      buf[i] =   (realtype) STATES_z0[i / 6][0];
      buf[i + 1] = (realtype) STATES_z0[i / 6][1];
      buf[i + 2] = (realtype) STATES_z0[i / 6][2];
      buf[i + 3] = (realtype) STATES_z0[i / 6][3];
      buf[i + 4] = (realtype) STATES_z0[i / 6][4];
      buf[i + 5] = (realtype) STATES_z0[i / 6][5];

      //printf("E:%.4e,i/6:%d\n",STATES_z0[i/6][2], i/6);
    }
  }

  MPI_Barrier(cpugrid);
  MPI_Bcast(buf, z0_len * 6, REALTYPE, 0, cpugrid);


  //NOW RECONSTRUCT on other procs
  if (myid > 0)
  {
    alloc2darr(realtype, STATES_z0, z0_len, 6);
    for (i = 0; i < z0_len * 6; i += 6)
    {
      STATES_z0[i / 6][0] = buf[i];
      STATES_z0[i / 6][1] = buf[i + 1];
      STATES_z0[i / 6][2] = buf[i + 2];
      STATES_z0[i / 6][3] = buf[i + 3];
      STATES_z0[i / 6][4] = buf[i + 4];
      STATES_z0[i / 6][5] = buf[i + 5];


    }
  }
  free(buf);

  // ***********************************************
  //Read Al, Z=+1
  // **********************************************
  if (myid == 0)
  {
    lcnt = 1;
    fin = fopen("Al1_states.txt", "r");
    if (fin == NULL)
    {
      char errstr[255];
      sprintf(errstr, "ERROR in colrad_read_states: File %s not found\n", "Al1_states.txt");
      error(errstr);
    }
    while (1) {
      if (fgets(line, MAXLINE, fin) == NULL) break;
      lcnt++;
    }

    alloc2darr(realtype, STATES_z1, lcnt, 6);

    lcnt = 0;
    rewind(fin);
    while (1)
    {
      if (fgets(line, MAXLINE, fin) == NULL) break;
      sscanf(line, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
             &STATES_z1[lcnt][0], &STATES_z1[lcnt][1], &STATES_z1[lcnt][2],
             &STATES_z1[lcnt][3], &STATES_z1[lcnt][4], &STATES_z1[lcnt][5]);
      lcnt++;
    }
    z1_len = lcnt;
    fclose(fin);
    total_species += z1_len;

  }

  //NOW COMMUNICATE
  MPI_Bcast(&z1_len, 1, MPI_INT, 0, cpugrid);
  alloc1darr(realtype, buf, z1_len * 6);
  if (myid == 0)
  {
    for (i = 0; i < z1_len * 6; i += 6) //fill 1d buff-array
    {
      buf[i] =   (realtype) STATES_z1[i / 6][0];
      buf[i + 1] = (realtype) STATES_z1[i / 6][1];
      buf[i + 2] = (realtype) STATES_z1[i / 6][2];
      buf[i + 3] = (realtype) STATES_z1[i / 6][3];
      buf[i + 4] = (realtype) STATES_z1[i / 6][4];
      buf[i + 5] = (realtype) STATES_z1[i / 6][5];
    }
  }

  MPI_Bcast(buf, z1_len * 6, REALTYPE, 0, cpugrid);
  //NOW RECONSTRUCT on other procs
  if (myid > 0)
  {
    alloc2darr(realtype, STATES_z1, z1_len, 6);
    for (i = 0; i < z1_len * 6; i += 6)
    {
      STATES_z1[i / 6][0] = buf[i];
      STATES_z1[i / 6][1] = buf[i + 1];
      STATES_z1[i / 6][2] = buf[i + 2];
      STATES_z1[i / 6][3] = buf[i + 3];
      STATES_z1[i / 6][4] = buf[i + 4];
      STATES_z1[i / 6][5] = buf[i + 5];

    }
  }
  free(buf);


  // ***********************************************
  //Read Al, Z=+2
  // **********************************************
  #if MAXLEVEL > 1
  if (myid == 0)
  {
    lcnt = 1;
    fin = fopen("Al2_states.txt", "r");
    if (fin == NULL)
    {
      char errstr[255];
      sprintf(errstr, "ERROR in colrad_read_states: File %s not found\n", "Al2_states.txt");
      error(errstr);
    }
    while (1) {
      if (fgets(line, MAXLINE, fin) == NULL) break;
      lcnt++;
    }

    alloc2darr(realtype, STATES_z2, lcnt, 6);

    lcnt = 0;
    rewind(fin);
    while (1)
    {
      if (fgets(line, MAXLINE, fin) == NULL) break;
      sscanf(line, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
             &STATES_z2[lcnt][0], &STATES_z2[lcnt][1], &STATES_z2[lcnt][2],
             &STATES_z2[lcnt][3], &STATES_z2[lcnt][4], &STATES_z2[lcnt][5]);
      lcnt++;
    }
    z2_len = lcnt;
    fclose(fin);
    total_species += z2_len;
  }

  //NOW COMMUNICATE
  MPI_Bcast(&z2_len, 1, MPI_INT, 0, cpugrid);
  alloc1darr(realtype, buf, z2_len * 6);
  if (myid == 0)
  {
    for (i = 0; i < z2_len * 6; i += 6) //fill 1d buff-array
    {
      buf[i] =   (realtype) STATES_z2[i / 6][0];
      buf[i + 1] = (realtype) STATES_z2[i / 6][1];
      buf[i + 2] = (realtype) STATES_z2[i / 6][2];
      buf[i + 3] = (realtype) STATES_z2[i / 6][3];
      buf[i + 4] = (realtype) STATES_z2[i / 6][4];
      buf[i + 5] = (realtype) STATES_z2[i / 6][5];
    }
  }

  MPI_Bcast(buf, z2_len * 6, REALTYPE, 0, cpugrid);
  //NOW RECONSTRUCT on other procs
  if (myid > 0)
  {
    alloc2darr(realtype, STATES_z2, z2_len, 6);
    for (i = 0; i < z2_len * 6; i += 6)
    {
      STATES_z2[i / 6][0] = buf[i];
      STATES_z2[i / 6][1] = buf[i + 1];
      STATES_z2[i / 6][2] = buf[i + 2];
      STATES_z2[i / 6][3] = buf[i + 3];
      STATES_z2[i / 6][4] = buf[i + 4];
      STATES_z2[i / 6][5] = buf[i + 5];

    }
  }
  free(buf);
  #endif //MAXLEVEL > 1


  // ***********************************************
  //Read Al, Z=+3
  // **********************************************
  #if MAXLEVEL > 2
  if (myid == 0)
  {
    lcnt = 1;
    fin = fopen("Al3_states.txt", "r");
    if (fin == NULL)
    {
      char errstr[255];
      sprintf(errstr, "ERROR in colrad_read_states: File %s not found\n", "Al3_states.txt");
      error(errstr);
    }
    while (1) {
      if (fgets(line, MAXLINE, fin) == NULL) break;
      lcnt++;
    }

    alloc2darr(realtype, STATES_z3, lcnt, 6);

    lcnt = 0;
    rewind(fin);
    while (1)
    {
      if (fgets(line, MAXLINE, fin) == NULL) break;
      sscanf(line, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
             &STATES_z3[lcnt][0], &STATES_z3[lcnt][1], &STATES_z3[lcnt][2],
             &STATES_z3[lcnt][3], &STATES_z3[lcnt][4], &STATES_z3[lcnt][5]);



      lcnt++;

    }
    z3_len = lcnt;
    fclose(fin);
    total_species += z3_len;
  }

  //NOW COMMUNICATE
  MPI_Bcast(&z3_len, 1, MPI_INT, 0, cpugrid);
  alloc1darr(realtype, buf, z3_len * 6);
  if (myid == 0)
  {
    for (i = 0; i < z3_len * 6; i += 6) //fill 1d buff-array
    {
      buf[i] =   (realtype) STATES_z3[i / 6][0];
      buf[i + 1] = (realtype) STATES_z3[i / 6][1];
      buf[i + 2] = (realtype) STATES_z3[i / 6][2];
      buf[i + 3] = (realtype) STATES_z3[i / 6][3];
      buf[i + 4] = (realtype) STATES_z3[i / 6][4];
      buf[i + 5] = (realtype) STATES_z3[i / 6][5];
    }
  }

  MPI_Bcast(buf, z3_len * 6, REALTYPE, 0, cpugrid);
  //NOW RECONSTRUCT on other procs
  if (myid > 0)
  {
    alloc2darr(realtype, STATES_z3, z3_len, 6);
    for (i = 0; i < z3_len * 6; i += 6)
    {
      STATES_z3[i / 6][0] = buf[i];
      STATES_z3[i / 6][1] = buf[i + 1];
      STATES_z3[i / 6][2] = buf[i + 2];
      STATES_z3[i / 6][3] = buf[i + 3];
      STATES_z3[i / 6][4] = buf[i + 4];
      STATES_z3[i / 6][5] = buf[i + 5];

    }
  }
  free(buf);
  #endif //MAXLEVEL > 2

  // ***********************************************
  //Read Al, Z=+4
  // **********************************************
  #if MAXLEVEL > 3
  if (myid == 0)
  {
    lcnt = 1;
    fin = fopen("Al4_states.txt", "r");
    if (fin == NULL)
    {
      char errstr[255];
      sprintf(errstr, "ERROR in colrad_read_states: File %s not found\n", "Al4_states.txt");
      error(errstr);
    }
    while (1) {
      if (fgets(line, MAXLINE, fin) == NULL) break;
      lcnt++;
    }

    alloc2darr(realtype, STATES_z4, lcnt, 6);

    lcnt = 0;
    rewind(fin);
    while (1)
    {
      if (fgets(line, MAXLINE, fin) == NULL) break;
      sscanf(line, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
             &STATES_z4[lcnt][0], &STATES_z4[lcnt][1], &STATES_z4[lcnt][2],
             &STATES_z4[lcnt][3], &STATES_z4[lcnt][4], &STATES_z4[lcnt][5]);
      lcnt++;
    }
    z4_len = lcnt;
    fclose(fin);
    total_species += z4_len;
  }

  //NOW COMMUNICATE
  MPI_Bcast(&z4_len, 1, MPI_INT, 0, cpugrid);
  alloc1darr(realtype, buf, z4_len * 6);
  if (myid == 0)
  {
    for (i = 0; i < z4_len * 6; i += 6) //fill 1d buff-array
    {
      buf[i] =   (realtype) STATES_z4[i / 6][0];
      buf[i + 1] = (realtype) STATES_z4[i / 6][1];
      buf[i + 2] = (realtype) STATES_z4[i / 6][2];
      buf[i + 3] = (realtype) STATES_z4[i / 6][3];
      buf[i + 4] = (realtype) STATES_z4[i / 6][4];
      buf[i + 5] = (realtype) STATES_z4[i / 6][5];
    }
  }

  MPI_Bcast(buf, z4_len * 6, REALTYPE, 0, cpugrid);
  //NOW RECONSTRUCT on other procs
  if (myid > 0)
  {
    alloc2darr(realtype, STATES_z4, z4_len, 6);
    for (i = 0; i < z4_len * 6; i += 6)
    {
      STATES_z4[i / 6][0] = buf[i];
      STATES_z4[i / 6][1] = buf[i + 1];
      STATES_z4[i / 6][2] = buf[i + 2];
      STATES_z4[i / 6][3] = buf[i + 3];
      STATES_z4[i / 6][4] = buf[i + 4];
      STATES_z4[i / 6][5] = buf[i + 5];

    }
  }
  free(buf);
  #endif //MAXLEVEL > 3

  // **********************************
  // * ALLOC ARRAYS
  // **********************************
  alloc2darr(realtype, k_EE_z0_z0, z0_len, z0_len);
  alloc2darr(realtype, k_EE_z0_z0_b, z0_len, z0_len);

  alloc2darr(realtype, k_EE_z1_z1, z1_len, z1_len);
  alloc2darr(realtype, k_EE_z1_z1_b, z1_len, z1_len);

  #if MAXLEVEL > 1
  alloc2darr(realtype, k_EE_z2_z2, z2_len, z2_len);
  alloc2darr(realtype, k_EE_z2_z2_b, z2_len, z2_len);
  #endif

  #if MAXLEVEL > 2
  alloc2darr(realtype, k_EE_z3_z3, z3_len, z3_len);
  alloc2darr(realtype, k_EE_z3_z3_b, z3_len, z3_len);
  #endif

  #if MAXLEVEL > 3
  alloc2darr(realtype, k_EE_z4_z4, z4_len, z4_len);
  alloc2darr(realtype, k_EE_z4_z4_b, z4_len, z4_len);
  #endif
  // **********************************************
  // Now Thermal Ionization and Recomb. rate coeff. arrays

  // z0->z1
  alloc2darr(realtype, k_EI_z0_z1, z0_len, z1_len);
  alloc2darr(realtype, k_EI_z1_z0, z0_len, z1_len);

  #if MAXLEVEL > 1
  //z1->z2
  alloc2darr(realtype, k_EI_z1_z2, z1_len, z2_len);
  alloc2darr(realtype, k_EI_z2_z1, z1_len, z2_len);
  #endif

  #if MAXLEVEL > 2
  //z2->z3
  alloc2darr(realtype, k_EI_z2_z3, z2_len, z3_len);
  alloc2darr(realtype, k_EI_z3_z2, z2_len, z3_len);
  #endif

  #if MAXLEVEL > 3
  //z3->z4
  alloc2darr(realtype, k_EI_z3_z4, z3_len, z4_len);
  alloc2darr(realtype, k_EI_z4_z3, z3_len, z4_len);
  #endif
  // k_MPI arrays
  //z0->z1
  alloc3darr(realtype, k_MPI_z0_z1, z0_len, z1_len, 2);
  alloc3darr(realtype, k_MPI_z1_z0, z0_len, z1_len, 2);

  //z1->z2
  #if MAXLEVEL > 1
  alloc3darr(realtype, k_MPI_z1_z2, z1_len, z2_len, 2);
  alloc3darr(realtype, k_MPI_z2_z1, z1_len, z2_len, 2);
  #endif

  //z2->z3
  #if MAXLEVEL > 2
  alloc3darr(realtype, k_MPI_z2_z3, z2_len, z3_len, 2);
  alloc3darr(realtype, k_MPI_z3_z2, z2_len, z3_len, 2);
  #endif


  //z3->z4
  #if MAXLEVEL > 3
  alloc3darr(realtype, k_MPI_z3_z4, z3_len, z4_len, 2);
  alloc3darr(realtype, k_MPI_z4_z3, z3_len, z4_len, 2);
  #endif

  MPI_Bcast(&total_species, 1, MPI_INT, 0, cpugrid);
}





void do_Saha(realtype Te,realtype totalc,realtype ne,N_Vector y) //Bei init
{
  double tmp;
  double Zav=ne/totalc;
  // double Ti=Te;

  double Q_z0,Q_z1,Q_z2,Q_z3,Q_z4,Q_z5; //partition functions
  // double IPD_SP=0.0; //Stewart and Pyatt =(POWR(1.0+a/Debye,2.0/3.0)-1.0)/2.0/(Zav+1.0);
  // double IPD_IS=0.0; //Ion-sphere model (high dens,low T) =3.0*1.0*ECHARGE*ECHARGE/2.0/a/BOLTZMAN/Te;
  // double IPD_DH=0.0; //Debye-Hückel (high T,low dens) = 1.0*ECHARGE*ECHARGE/Debye/BOLTZMAN/Te;
  double r10,r21,r32,r43,r54; //ratios of ion. conc.
  double DeltaE;
  double IPD=0.0;
  double n0,n1,n2,n3,n4,n5;
  double p; //probability


  int i;


  r10=r21=r32=r43=r54=0;


  Zav=ne/totalc;
  double EF=fermi_E(ne);
  double mu=chempot(ne,Te);
  //Debye=SQRTR(BOLTZMAN*Te/4.0/pi/ECHARGE/ECHARGE/ne/(Zav+1));
  tmp=POWR(2.0*pi*EMASS*BOLTZMAN*Te,1.5)/planck/planck/planck; //ACHTUNG: Das ist wieder thermal-de brogilie Lambda
                                                              //Nutze richtiges chempot!!!!
  #ifdef DOIPD
  double IPD0,IPD1,IPD2,IPD3,IPD4;
  double z; //DOI after ionization (e.g.=1 for Al0)
  double r0; //Ion sphere radius
  r0=POWR(3.0/4.0/pi/totalc,1.0/3.0);
  double debye=SQRTR(BOLTZMAN*Te/4.0/pi/POWR(totalc+ne,2.0));
  // Atoms, solids, and plasmas in super-intense laser fields S.220
  IPD0=1.0*3.0/2.0/r0*ECHARGE*ECHARGE*(POWR(1.0+POWR(debye/r0,3.0),2.0/3.0)-POWR(debye/r0,2.0))/4.0/pi/ECONST;
  IPD1=2.0*3.0/2.0/r0*ECHARGE*ECHARGE*(POWR(1.0+POWR(debye/r0,3.0),2.0/3.0)-POWR(debye/r0,2.0))/4.0/pi/ECONST;
  IPD2=3.0*3.0/2.0/r0*ECHARGE*ECHARGE*(POWR(1.0+POWR(debye/r0,3.0),2.0/3.0)-POWR(debye/r0,2.0))/4.0/pi/ECONST;
  IPD3=4.0*3.0/2.0/r0*ECHARGE*ECHARGE*(POWR(1.0+POWR(debye/r0,3.0),2.0/3.0)-POWR(debye/r0,2.0))/4.0/pi/ECONST;  
  IPD4=5.0*3.0/2.0/r0*ECHARGE*ECHARGE*(POWR(1.0+POWR(debye/r0,3.0),2.0/3.0)-POWR(debye/r0,2.0))/4.0/pi/ECONST;  
  #endif  

  //compute partition functions
  Q_z0=0.0;
  Q_z1=0.0;
  Q_z2=0.0;
  Q_z3=0.0;
  Q_z4=0.0;
  Q_z5=0.0;
  for(i=0;i<z0_len;++i)
  {
    #ifdef DOIPD
    IPD=IPD0;
    #endif
    double Ei=STATES_z0[i][2]*eV2J-IPD+mu;
    double qi=STATES_z0[i][3]*exp(-(STATES_z0[i][2]-0.0)*eV2J/BOLTZMAN/Te); 
    if(Ei<0) //depressed state- > truncate partiation function 
      qi=0;

    Q_z0+=qi; //Mit Energien, relativ zum Grundzustand (level=0) dieses Ions

// if(myid==0)
// printf("Z0, qi:%.4e, exp:%.4e, gi:%.4e,Q:%.4e\n", qi, 
//         exp(-(STATES_z0[i][2]-0.0)*eV2J/BOLTZMAN/Te),
//         STATES_z0[i][3], Q_z0);

  }

  for(i=0;i<z1_len;++i)
  { 
    #ifdef DOIPD
    IPD=IPD0;
    #endif
    double Ei=STATES_z1[i][2]*eV2J-IPD+mu;
    double qi=STATES_z1[i][3]*exp(-(STATES_z1[i][2]-STATES_z1[0][2])*eV2J/BOLTZMAN/Te);     
    if(Ei<0) qi=0;

    Q_z1+=qi;

// if(myid==0)
// printf("Z1, qi:%.4e, exp:%.4e, gi:%.4e,Q:%.4e\n", qi, 
//         exp(-(STATES_z1[i][2]-STATES_z1[0][2])*eV2J/BOLTZMAN/Te),
//         STATES_z1[i][3],Q_z1);

  }

  #if MAXLEVEL > 1
  for(i=0;i<z2_len;++i)
  {
    #ifdef DOIPD
    IPD=IPD1;
    #endif
    double Ei=STATES_z2[i][2]*eV2J-IPD+mu;
    double qi=STATES_z2[i][3]*exp(-(STATES_z2[i][2]-STATES_z2[0][2])*eV2J/BOLTZMAN/Te);     
    if(Ei<0) qi=0;

    Q_z2+=qi;

// if(myid==0)
// printf("Z2, qi:%.4e, exp:%.4e, gi:%.4e,Q:%.4e\n", qi, 
//         exp(-(STATES_z2[i][2]-STATES_z2[0][2])*eV2J/BOLTZMAN/Te),
//         STATES_z2[i][3], Q_z2);    
  }
  #endif

  #if MAXLEVEL > 2
  for(i=0;i<z3_len;++i)
  {
    #ifdef DOIPD
    IPD=IPD2;
    #endif
    double Ei=STATES_z3[i][2]*eV2J-IPD+mu;
    double qi=STATES_z3[i][3]*exp(-(STATES_z3[i][2]-STATES_z3[0][2])*eV2J/BOLTZMAN/Te);     
    if(Ei<0) qi=0;

    Q_z3+=qi;

// if(myid==0)
// printf("Z3, qi:%.4e, exp:%.4e, gi:%.4e,Q:%.4e\n", qi, 
//         exp(-(STATES_z3[i][2]-STATES_z3[0][2])*eV2J/BOLTZMAN/Te),
//         STATES_z3[i][3],Q_z3);        
  }
  #endif

  #if MAXLEVEL > 3
   for(i=0;i<z4_len;++i)
   {
    #ifdef DOIPD
    IPD=IPD3;
    #endif
    double Ei=STATES_z4[i][2]*eV2J-IPD+mu;
    double qi=STATES_z4[i][3]*exp(-(STATES_z4[i][2]-STATES_z4[0][2])*eV2J/BOLTZMAN/Te);     
    if(Ei<0) qi=0;
 
     Q_z4+=qi;

// if(myid==0)
// printf("Z4, qi:%.4e, exp:%.4e, gi:%.4e,Q:%.4e\n", qi, 
//         exp(-(STATES_z4[i][2]-STATES_z4[0][2])*eV2J/BOLTZMAN/Te),
//         STATES_z4[i][3],Q_z4 );

   }
  #endif

  /*
  for(int i=0;i<z5_len;++i)
    Q_z5+=STATES_z5[i][3]*exp(-(STATES_z5[i][2]-STATES_z5[0][2])*eV2J/BOLTZMAN/Te);
  */
  // double tmp2,tmp3;
  ////////////
  // Z=0->1 //
  ////////////

  #ifdef DOIPD
  // z=1.0;
  // IPD=3.0*z*ECHARGE*ECHARGE/2.0/r0/4.0/pi/ECONST;  //Ion-Sphere model
  IPD=IPD0;
  #endif

  DeltaE=(STATES_z1[0][2]-0.0)*eV2J-IPD+mu;
  // DeltaE=fmax(0.0,DeltaE);
  p=exp(-DeltaE/BOLTZMAN/Te);
  // r10=2.0/ne*tmp*Q_z1/Q_z0*p; //r10= ratio of ion-concentrations n(Z=1)/n(Z=0); g_e=2 (statistical weight of electron)
  r10=Q_z1/Q_z0*p;
  // double r01=Q_z0/Q_z1*p;
  //printf("IPD0:%f\n",IPD*J2eV);
  // if(isnan(r10)!=0 || isinf(r10)!=0)
  // {
  //   char errstr[255];
  //   sprintf(errstr,"ERROR in Saha-Init: r10 is inf or nan:  p0:%.4e,Q_z1:%.4e,Q_z0:%.4e\n",p,Q_z1,Q_z0);
  //   error(errstr);
  // }
  // if(myid==0)
  //   printf("ERROR in Saha-Init: r10 is inf or nan:  p0:%.4e,Q_z1:%.4e,Q_z0:%.4e\n",p,Q_z1,Q_z0);


  ////////////
  // Z=1->2 //
  ////////////

  #if MAXLEVEL > 1

  #ifdef DOIPD
    // z=2.0;
    // IPD=3.0*z*ECHARGE*ECHARGE/2.0/r0/4.0/pi/ECONST;  //Ion-Sphere model
    IPD=IPD1;
  #endif  
    DeltaE=(STATES_z2[0][2]-STATES_z1[0][2])*eV2J-IPD+mu;
    // DeltaE=fmax(0.0,DeltaE);
    p=exp(-DeltaE/BOLTZMAN/Te);
    // p=exp(DeltaE/BOLTZMAN/Te);
    // r21=2.0/ne*tmp*Q_z2/Q_z1*p;
    r21=Q_z2/Q_z1*p;
    // double r12=Q_z1/Q_z2*p;
  // if(isnan(r21)!=0 || isinf(r21)!=0)
  // {
  //   char errstr[255];
  //   sprintf(errstr,"ERROR in Saha-Init: r21 is inf or nan:  p2:%.4e,Q_z2:%.4e,Q_z1:%.4e\n",p,Q_z2,Q_z1);
  //   error(errstr);
  // }
#endif
  // ////////////
  // // Z=2->3 //
  // ////////////
  #if MAXLEVEL > 2  
  #ifdef DOIPD  
    // z=3.0;
    // IPD=3.0*z*ECHARGE*ECHARGE/2.0/r0/4.0/pi/ECONST;  //Ion-Sphere model
    IPD=IPD2;
  #endif  
    DeltaE=(STATES_z3[0][2]-STATES_z2[0][2])*eV2J-IPD+mu;
    // DeltaE=fmax(0.0,DeltaE);
    p=exp(-DeltaE/BOLTZMAN/Te);
    // p=exp(DeltaE/BOLTZMAN/Te);
    // r32=2.0/ne*tmp*Q_z3/Q_z2*p;
    r32=Q_z3/Q_z2*p;
    // double r23=Q_z2/Q_z3*p;

  if(isnan(r32)!=0 || isinf(r32)!=0)
  {
    char errstr[255];
    sprintf(errstr,"ERROR in Saha-Init: r32 is inf or nan:  p3:%.4e,Q_z3:%.4e,Q_z2:%.4e\n",p,Q_z3,Q_z2);
    error(errstr);
  }
  #endif  
  // ////////////
  // // Z=3->4 //
  // ////////////
  #if MAXLEVEL > 3
  #ifdef DOIPD  
  // z=4.0;
  // IPD=3.0*z*ECHARGE*ECHARGE/2.0/r0/4.0/pi/ECONST;  //Ion-Sphere model
  IPD=IPD3;
  #endif  
  DeltaE=(STATES_z4[0][2]-STATES_z3[0][2])*eV2J-IPD+mu;
  // DeltaE=fmax(0.0,DeltaE);
  p=exp(-DeltaE/BOLTZMAN/Te);
  // p=exp(DeltaE/BOLTZMAN/Te);
  // r43=2.0/ne*tmp*Q_z4/Q_z3*p;
  r43=Q_z4/Q_z3*p;
  // double r34=Q_z3/Q_z4*p;
// printf("r34:%.4e,p:%.4e,DeltaE/kT:%.4e\n",r34,p,DeltaE/BOLTZMAN/Te);
  // if(isnan(r43)!=0 || isinf(r43)!=0)
  // {
  //   char errstr[255];
  //   sprintf(errstr,"ERROR in Saha-Init: r43 is inf or nan:  p4:%.4e,Q_z4:%.4e,Q_z3:%.4e\n",p,Q_z4,Q_z3);
  //   error(errstr);
  // }
#endif
  //concentrations from ratios and totalc

  n0=totalc/(r43*r32*r21*r10+r32*r21*r10+r21*r10+r10+1.0); // ?
  //n0=totalc*Zav/(4*r43*r32*r21*r10 + 3*r32*r21*r10+ 2*r21*r10+1.0);
  n1=r10*n0;
  
  //D.h. Neutrals komplett druck-ionisiert
  if(Q_z0==0)
  {
    n0=0.0;
    // n1=totalc*Zav/(4*r43*r32*r21 + 3*r32*r21+ 2*r21+1.0);
    n1=totalc/(r43*r32*r21 + r32*r21+ r21+1.0);
  }

  n2=r21*n1;
  n3=r32*n2;
  n4=r43*n3;
  n5=r54*n4;
  // Zav=(1*n1+4*n2+9*n3+16*n4+25*n5)/(1*n1+2*n2+3*n3+4*n4+5*n5);

  // n4=Zav*totalc/(r12*r23*r34+2*r23*r34+3*r34+4);
  // n3=r34*n4;
  // n2=r23*n3;
  // n1=r12*n2;
  // n0=r01*n1;

// if(myid==0)
// {
  // printf("n0:%.4e,n1:%.4e,n2:%.4e,n3:%.4e,n4:%.4e\n",n0,n1,n2,n3,n4);
  // printf("r10:%.4e,r21:%.4e,r32:%.4e,r43:%.4e\n",r10,r21,r32,r43);
  // printf("Qz4:%.4e,Qz3:%.4e,Qz2:%.4e,Qz1:%.4e,Qz0:%.4e\n",Q_z4, Q_z3,Q_z2,Q_z1,Q_z0);
// }


  /*
  printf("********************************************************************************************************************** \n");
  printf(" Initial distribution of concentration according to generalized Saha-Equation\n");
  printf(" Ti=Te=%f\n",Te);
  printf(" totalc=%.4e\n",totalc);
  printf(" Zmean=%.4e,ne:%.6e\n",ne/totalc,ne);
  printf(" [Al0]:%.2e,[Al1]:%.2e,[Al2]:%.2e,[Al3]:%.2e,[Al4]:%.2e,[Al5]:%.2e\n",n0,n1,n2,n3,n4,n5);
  printf(" I0=%.4e\n",I0);
  printf(" t_FWHM=%.4e\n",tFWHM);
  printf(" Lambda=%.4e\n",lambda);
  printf(" NEQ:%d\n",neq);
  printf("********************************************************************************************************************** \n");
  */


  // ************************
  // *  FILL Z0 STATES
  // ************************
  int ishift=3;
  //Now fill states
  for(i=0;i<z0_len;++i)
  {  
    DeltaE=(STATES_z0[i][2]-0)*eV2J; // bzgl. ground state
    double Ei=0.0;
    double prob=exp(-DeltaE/BOLTZMAN/Te);
    #ifdef DOIPD
    IPD=IPD0;
    Ei=STATES_z0[i][2]*eV2J-IPD+mu;    
    if(Ei<0) prob=0.0;    
    #endif
    if(Q_z0>0)
    {
      Ith(y,i+ishift)=n0/Q_z0*STATES_z0[i][3]*prob;
      if(Ith(y,i+ishift) < MINCONC)
        Ith(y,i+ishift)=0.0;

      // printf("ITH(%d):%.4e,prob:%.4e,Q_z0:%.4e,dE:%.4e\n",
      //   i+ishift,Ith(y,i+ishift), prob, Q_z0, DeltaE);
    }

    if(isnan(Ith(y,i+ishift))!=0 || isinf(Ith(y,i+ishift))!=0)
    {
      char errstr[255];
      sprintf(errstr,"ERROR in do_Saha z0: y is inf or nan! prob:%.4e, Ei:%.4e,DeltaE:%.4e\n", 
        prob,Ei,DeltaE);
      error(errstr);
    }
  }
    // ************************
    // *  FILL Z1 STATES
    // ************************
  for(i=0;i<z1_len;++i)
  {
    // DeltaE=(STATES_z1[i][2]-STATES_z1[0][2])*eV2J-IPD0+EF;
    // DeltaE=(STATES_z1[i][2]-STATES_z1[0][2])*eV2J-IPD1+EF;
    DeltaE=(STATES_z1[i][2]-STATES_z1[0][2])*eV2J;
    double Ei=0.0;
    double prob=exp(-DeltaE/BOLTZMAN/Te);
    #ifdef DOIPD
    IPD=IPD1;
    // Ei=(STATES_z1[i][2]-STATES_z1[0][2])*eV2J-IPD0+EF;    
    // Ei=(STATES_z1[i][2]-STATES_z1[0][2])*eV2J-IPD1+EF;   
    Ei=(STATES_z1[i][2])*eV2J-IPD0+mu;   
    if(Ei<0) prob=0.0;    
    #endif    
    if(Q_z1>0)
      Ith(y,i+ishift+z0_len)=n1/Q_z1*STATES_z1[i][3]*prob;

    if(Ith(y,i+ishift+z0_len) < MINCONC)
      Ith(y,i+ishift+z0_len)=0.0;

    if(isnan(Ith(y,i+ishift+z0_len))!=0 || isinf(Ith(y,i+ishift+z0_len))!=0)
    {
      char errstr[255];
      sprintf(errstr,"ERROR in do_Saha z1: y is inf or nan! prob:%.4e, Ei:%.4e,DeltaE:%.4e\n", 
        prob,Ei,DeltaE);
      error(errstr);
    }    
  }
  // ************************
  // *  FILL Z2 STATES
  // ************************
  #if MAXLEVEL > 1
  for(i=0;i<z2_len;++i)
  {
    DeltaE=(STATES_z2[i][2]-STATES_z2[0][2])*eV2J;//-IPD2+EF;
    double Ei=0.0;    
    double prob=exp(-DeltaE/BOLTZMAN/Te);
    #ifdef DOIPD
    IPD=IPD2;
    Ei=STATES_z2[i][2]*eV2J-IPD1+mu;    
    if(Ei<0) prob=0.0;    
    #endif 
    if(Q_z2>0)       
      Ith(y,i+ishift+z0_len+z1_len)=n2/Q_z2*STATES_z2[i][3]*prob;

    if(Ith(y,i+ishift+z0_len+z1_len) < MINCONC)
      Ith(y,i+ishift+z0_len+z1_len)=0.0;

    if(isnan(Ith(y,i+ishift+z0_len+z1_len))!=0 || isinf(Ith(y,i+ishift+z0_len+z1_len))!=0)
    {
      char errstr[255];
      sprintf(errstr,"ERROR in do_Saha z2: y is inf or nan! prob:%.4e, Ei:%.4e,DeltaE:%.4e\n", 
        prob,Ei,DeltaE);
      error(errstr);
    }

  }
  #endif
  // ************************
  // *  FILL Z3 STATES
  // ************************  
  #if MAXLEVEL > 2
  for(i=0;i<z3_len;++i)
  {
    DeltaE=(STATES_z3[i][2]-STATES_z3[0][2])*eV2J;//-IPD3+EF;
    double Ei=0.0;    
    double prob=exp(-DeltaE/BOLTZMAN/Te);
    #ifdef DOIPD
    IPD=IPD3;
    Ei=STATES_z3[i][2]*eV2J-IPD2+mu;    
    if(Ei<0) prob=0.0;    
    #endif            
    if(Q_z3>0)
      Ith(y,i+ishift+z0_len+z1_len+z2_len)=n3/Q_z3*STATES_z3[i][3]*prob;

    if(Ith(y,i+ishift+z0_len+z1_len+z2_len) < MINCONC)
      Ith(y,i+ishift+z0_len+z1_len+z2_len)=0.0;

    if(isnan(Ith(y,i+ishift+z0_len+z1_len+z2_len))!=0 || isinf(Ith(y,i+ishift+z0_len+z1_len+z2_len))!=0)
    {
      char errstr[255];
      sprintf(errstr,"ERROR in do_Saha z3: y is inf or nan! prob:%.4e, Ei:%.4e,DeltaE:%4.e\n", 
        prob,Ei,DeltaE);
      error(errstr);
    }    
  }
  #endif
  // ************************
  // *  FILL Z4 STATES
  // ************************
  #if MAXLEVEL > 3
  for(i=0;i<z4_len;++i)
  {
    DeltaE=(STATES_z4[i][2]-STATES_z4[0][2])*eV2J;//-IPD4+EF;
    double prob=exp(-DeltaE/BOLTZMAN/Te);
    double Ei=0.0;    
    #ifdef DOIPD
    IPD=IPD4;
    Ei=STATES_z4[i][2]*eV2J-IPD3+mu;    
    if(Ei<0) prob=0.0;    
    #endif                
    if(Q_z4>0)
      Ith(y,i+ishift+z0_len+z1_len+z2_len+z3_len)=n4/Q_z4*STATES_z4[i][3]*prob;

    if(Ith(y,i+ishift+z0_len+z1_len+z2_len+z3_len) < MINCONC)
      Ith(y,i+ishift+z0_len+z1_len+z2_len+z3_len)=0.0;

    if(isnan(Ith(y,i+ishift+z0_len+z1_len+z2_len+z3_len))!=0 || isinf(Ith(y,i+ishift+z0_len+z1_len+z2_len+z3_len))!=0)
    {
      char errstr[255];
      sprintf(errstr,"ERROR in do_Saha z4: y is inf or nan! prob:%.4e, Ei:%.4e,DeltaE:%4.e\n", 
        prob,Ei,DeltaE);
      error(errstr);
    }

  }
  #endif  

  double totalc_check=0.0;
  double n0_check=0.0;
  double n1_check=0.0;
  double n2_check=0.0;
  double n3_check=0.0;
  double n4_check=0.0;

  for(i=3;i<neq;i++)
  {
    totalc_check+=Ith(y,i);
    if(i-3 < z0_len)
      n0_check+=Ith(y,i);
    
    if(i-3 >= z0_len && i-3-z0_len < z1_len)
      n1_check+=Ith(y,i);

    if(i-3 >= z0_len+z1_len && i-3-z0_len-z1_len < z2_len)
      n2_check+=Ith(y,i);

    if(i-3 >= z0_len+z1_len+z2_len && i-3-z0_len-z1_len-z2_len < z3_len)
      n3_check+=Ith(y,i);    

    if(i-3 >= z0_len+z1_len+z2_len + z3_len && i-3-z0_len-z1_len-z2_len - z3_len < z4_len)
      n4_check+=Ith(y,i);      

    // if(myid==0) printf("y[%d]:%.4e\n",i,Ith(y,i));
  } 

/*
  if(ABS(totalc-totalc_check) > 0.01*totalc || ABS(n0_check-n0) > 0.01* n0 || 
     ABS(n1_check-n1) > 0.01* n1 || ABS(n2_check-n2) > 0.01* n2 || ABS(n3_check-n3) > 0.01* n3 ||
     ABS(n4_check-n4) > 0.01* n4 )
  if(myid==0)
  {
    char errstr[400];
    sprintf(errstr,"Inconsistency in do_Saha: totalc_check= %.4e, and totalc=%.4e,sum_n:%.4e,"
                    "n0:%.4e, n1:%.4e,n2:%.4e,n3:%.4e, n4:%.4e\n"
                    "n0_check:%.4e, n1_check:%.4e,n2_check:%.4e,n3_check:%.4e, n4_check:%.4e\n",
      totalc_check,totalc,n1+n2+n3+n4, n0, n1, n2, n3, n4,n0_check,n1_check,n2_check,n3_check,n4_check);
    error(errstr);
  }
*/

}





// ************************************************************************************************************
//                                      ACTION
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////
int colrad_ydot(realtype t, N_Vector y, N_Vector colrad_ydot, void *user_data)
{

/*
  double t0=1e-12;
  double I0=1e17;
  double tFWHM=100e-15;
  double sigmat=tFWHM/2.3548;
  double sigmatsq=sigmat*sigmat;
  */

  //It=I0*exp(-POWR((t-t0),2.0)/2.0/sigmatsq);

  colrad_UserData data;
  data = (colrad_UserData) user_data;

  realtype It=data->It;
  It=0; //It ist LOKALE Intesität!

  bool initial_equi=data->initial_equi; //falls ja, wird temperatur nicht variiert.

  realtype Eexc;

  P_E_EE=0.0;
  P_E_EI=0.0;
  P_E_MPI2=0.0;
  P_E_MPI3=0.0;
  P_E_RAD_RECOMB=0.0;

  realtype DeltaE;
  realtype kfwd,krev;
  realtype kfwd2; // für 3-Photon absorption


  int i,j;
  int ishift=3;  // UND NICHT =4 (colrad_ydot fängt bei 0 das zählen an!)
  int shift,shift1,shift2;
  realtype ne=Ith(y,2);
  realtype Te=Ith(y,0);
  realtype Ti=Ith(y,1);


  //FOR REABSORPTION
  realtype Ajk;
  realtype sigma_PI,abs_coeff,thickness,tau,wstark_freq,lambda0,wstark_len,escape_factor,groundstate_ioniz;
  thickness=1e-3;
  realtype sigma_tmp=64.0*POWR(pi,4.0)*POWR(ECHARGE,10.0)*EMASS/3.0/SQRTR(3.0)/POWR(4.0*pi*ECONST,5.0)/POWR(planck,6.0)/LIGHTSPEED/POWR(LASERFREQ,3.0)/POWR(13.6*eV2J,2.0);


  //Pre-Zero <--- MUI IMPORTANTE!!!!
  //and count totalc for ipd and stark effect    
  totalc=0.0;
  // #ifdef OMP
  // #pragma omp parallel for reduction(+: totalc)
  // #endif
  for(i=0;i<neq;i++)
  {
    Ith(colrad_ydot,i)=0.0;
    if(i >= 3) totalc+=Ith(y,i);
    if(isnan(Ith(y,i))!=0)
    {
      printf("myid:%d, WARNING Ith(y,%d) became NaN! z0len:%d,z1len:%d,z2len:%d\n",myid,i,z0_len,z1_len,z2_len);
      return 1;
    }
  }
  data->ni=totalc;  
  if(totalc<0)
    return 1; // <-- RHS fail

  // IPD KRAM
  realtype IPD0,IPD1,IPD2,IPD3,IPD4;
  IPD0=IPD1=IPD2=IPD3=0.0;
  realtype EF=fermi_E(ne);
  data->EF=EF;
#ifdef DOIPD
  realtype r0,debye;
  r0=POWR(3.0/4.0/pi/totalc,1.0/3.0);
  debye=SQRTR(BOLTZMAN*Te/4.0/pi/POWR(totalc+ne,2.0));
  // Atoms, solids, and plasmas in super-intense laser fields S.220
  IPD0=1.0*3.0/2.0/r0*ECHARGE*ECHARGE*(POWR(1.0+POWR(debye/r0,3.0),2.0/3.0)-POWR(debye/r0,2.0))/4.0/pi/ECONST;
  IPD1=2.0*3.0/2.0/r0*ECHARGE*ECHARGE*(POWR(1.0+POWR(debye/r0,3.0),2.0/3.0)-POWR(debye/r0,2.0))/4.0/pi/ECONST;
  IPD2=3.0*3.0/2.0/r0*ECHARGE*ECHARGE*(POWR(1.0+POWR(debye/r0,3.0),2.0/3.0)-POWR(debye/r0,2.0))/4.0/pi/ECONST;
  IPD3=4.0*3.0/2.0/r0*ECHARGE*ECHARGE*(POWR(1.0+POWR(debye/r0,3.0),2.0/3.0)-POWR(debye/r0,2.0))/4.0/pi/ECONST;
  IPD4=5.0*3.0/2.0/r0*ECHARGE*ECHARGE*(POWR(1.0+POWR(debye/r0,3.0),2.0/3.0)-POWR(debye/r0,2.0))/4.0/pi/ECONST;

  data->IPD0=IPD0;
  data->IPD1=IPD1;
  data->IPD2=IPD2;
  data->IPD3=IPD3;
  data->IPD4=IPD4;

#endif

  int retval=colrad_GetCoeffs(y,It,data);
  if(retval !=0 )
    return 1; //d.h. failure of RHS



//printf("IPD0:%.4e,IPD1:%.4e,IPD2:%.4e\n", IPD0*J2eV,IPD1*J2eV,IPD2*J2eV);
  //**********************************************
  //Z=0, Exec/De-Exec  + SPONTANEOUS EMISSION
  //**********************************************
#ifdef OMP  //OMP FUNZT AUF DIESE WEISE NICHT !!!
  //#pragma omp parallel for schedule(static) collapse(2) private(DeltaE,kfwd,krev,lambda0,Ajk,escape_factor,Eexc) num_threads(num_threads)
  #pragma omp for nowait schedule(static) private(j,DeltaE,kfwd,krev,lambda0,Ajk,escape_factor,Eexc)
#endif
  for(i=0;i<z0_len;++i)
  {
    for(j=0;j<z0_len;++j)
    {
      if(j<=i) continue; // MUI IMPORTANTE
      realtype engi=STATES_z0[i][2]*eV2J-IPD0+EF;
      if(engi < 0) continue; //depressed state is continuum

      DeltaE=(STATES_z0[j][2]-STATES_z0[i][2])*eV2J;
      kfwd=k_EE_z0_z0[i][j]*Ith(y,i+ishift)*ne;
      krev=k_EE_z0_z0_b[i][j]*Ith(y,j+ishift)*ne;

      //exec. reduces conc. of i state and increases conc. of j state                     
      Ith(colrad_ydot,i+ishift) -=kfwd;
      Ith(colrad_ydot,j+ishift) +=kfwd;

      //de-excec. increases conc. of i state & decr. conc. of j state
      Ith(colrad_ydot,i+ishift)+=krev;
      Ith(colrad_ydot,j+ishift)-=krev;


      Eexc= (-kfwd+krev)*DeltaE;
      P_E_EE+=Eexc;

#ifdef SPONT
      // ********SPONT EMISSION ********* //
      //Spont.emiss: delta-n !=>0 und delta-l=+-1 (nur opt.allowed trans.)
      if((STATES_z0[j][5]-STATES_z0[i][5])>0 && STATES_z0[j][4]-STATES_z0[i][4]==1)
      {
        Ajk=0.0;
        escape_factor=1.0;
        lambda0=planck*LIGHTSPEED/DeltaE;
        Ajk=EinsteinCoeff(STATES_z0[i][5],STATES_z0[j][5],STATES_z0[j][3],DeltaE);

        krev=Ith(y,ishift+j)*Ajk; //neues krev
        Ith(colrad_ydot,ishift+j) -= krev;
        Ith(colrad_ydot,ishift+i) += krev;
        //P_A_SE+=krev*DeltaE; ??               
      }
#endif

    }
  }

    // *************************************
  //Z=1, Exec/De-Exec  & SPONTANEOUS EMISSION
  // ***************************************
  shift2=z0_len;
#ifdef OMP
  //#pragma omp parallel for schedule(dynamic,1) collapse(2) private(DeltaE,kfwd,krev,lambda0,Ajk,escape_factor,Eexc) num_threads(num_threads)
  #pragma omp for nowait schedule(static) private(j,DeltaE,kfwd,krev,lambda0,Ajk,escape_factor,Eexc)
#endif
  for(i=0;i<z1_len;++i)
  {
    for(j=0;j<z1_len;++j)
    {
      if(i<=j) continue; // MUI IMPORTANTE
      realtype engi=STATES_z1[i][2]*eV2J-IPD0+EF;            
      if(engi < 0) continue; //depressed state is continuum

      DeltaE=(STATES_z1[j][2]-STATES_z1[i][2])*eV2J;
      kfwd=k_EE_z1_z1[i][j]*Ith(y,i+ishift+shift2)*ne;
      krev=k_EE_z1_z1_b[i][j]*Ith(y,j+ishift+shift2)*ne;

      //exec. reduces conc. of i state and increases conc. of j state
      Ith(colrad_ydot,i+ishift+shift2) -=kfwd;
      Ith(colrad_ydot,j+ishift+shift2) +=kfwd;

      //de-excec. increases conc. of i state & decr. conc. of j state
      Ith(colrad_ydot,i+ishift+shift2)+=krev;
      Ith(colrad_ydot,j+ishift+shift2)-=krev;

// if(myid==1)
//   printf("Z=1:kfwd:%.4e,krev:%.4e, yi:%.4e,yj:%.4e\n", kfwd,krev, Ith(y,i+ishift+shift2), Ith(y,j+ishift+shift2));

      Eexc= (-kfwd+krev)*DeltaE;
      P_E_EE+=Eexc;

#ifdef SPONT
      // ********SPONT EMISSION ********* //
      //Spont.emiss: delta-n !=>0 und delta-l=+-1 (nur opt.allowed trans.)
      if((STATES_z1[j][5]-STATES_z1[i][5])>0 && STATES_z1[j][4]-STATES_z1[i][4]==1)
      {
        escape_factor=1.0;
        lambda0=planck*LIGHTSPEED/DeltaE;
        Ajk=EinsteinCoeff(STATES_z1[i][5],STATES_z1[j][5],STATES_z1[j][3],DeltaE);

        krev=Ith(y,ishift+shift2+j)*Ajk;
        Ith(colrad_ydot,ishift+shift2+j) -= krev;
        Ith(colrad_ydot,ishift+shift2+i) += krev;
        //P_A_SE+=krev*DeltaE; ??
      }
#endif
    }
  }

//   // *************************************
//   //Z=2, Exec/De-Exec  & SPONTANEOUS EMISSION
//   // ***************************************
#if MAXLEVEL > 1
  shift2=z0_len+z1_len;
#ifdef OMP
  //#pragma omp parallel for schedule(dynamic,1) collapse(2) private(DeltaE,kfwd,krev,lambda0,Ajk,escape_factor,Eexc) num_threads(num_threads)
  #pragma omp for nowait schedule(static) private(j,DeltaE,kfwd,krev,lambda0,Ajk,escape_factor,Eexc)
#endif
  for(i=0;i<z2_len;++i)
  {
    for(j=0;j<z2_len;++j)
    {
      if(j<=i) continue; // MUI IMPORTANTE
      realtype engi=STATES_z2[i][2]*eV2J-IPD1+EF;
      if(engi < 0) continue; //depressed state is continuum

      DeltaE=(STATES_z2[j][2]-STATES_z2[i][2])*eV2J;
      kfwd=k_EE_z2_z2[i][j]*Ith(y,i+ishift+shift2)*ne;
      krev=k_EE_z2_z2_b[i][j]*Ith(y,j+ishift+shift2)*ne;

      //exec. reduces conc. of i state and increases conc. of j state
      Ith(colrad_ydot,i+ishift+shift2) -=kfwd;
      Ith(colrad_ydot,j+ishift+shift2) +=kfwd;

      //de-excec. increases conc. of i state & decr. conc. of j state
      Ith(colrad_ydot,i+ishift+shift2)+=krev;
      Ith(colrad_ydot,j+ishift+shift2)-=krev;

// if(myid==1)
//  printf("Z=2:kfwd:%.4e,krev:%.4e, yi:%.4e,yj:%.4e\n", kfwd,krev, Ith(y,i+ishift+shift2), Ith(y,j+ishift+shift2));

// printf("kfwd:%.5e, krev:%.5e,i:%d,j:%d,keefwd:%.4e,keeb:%.4e,T:%.4e\n",kfwd,krev,i,j,k_EE_z2_z2[i][j], k_EE_z2_z2_b[i][j],Te);

      Eexc= (-kfwd+krev)*DeltaE;
      P_E_EE+=Eexc;

#ifdef SPONT
      // ********SPONT EMISSION ********* //
      //Spont.emiss: delta-n !=>0 und delta-l=+-1 (nur opt.allowed trans.)
      if((STATES_z2[j][5]-STATES_z2[i][5])>0 && STATES_z2[j][4]-STATES_z2[i][4]==1)
      {
        escape_factor=1.0;
        lambda0=planck*LIGHTSPEED/DeltaE;
        Ajk=EinsteinCoeff(STATES_z2[i][5],STATES_z2[j][5],STATES_z2[j][3],DeltaE);

        krev=Ith(y,ishift+shift2+j)*Ajk;
        Ith(colrad_ydot,ishift+shift2+j) -= krev;
        Ith(colrad_ydot,ishift+shift2+i) += krev;
        //P_A_SE+=krev*DeltaE; ??
      }
#endif
    }
  }
#endif

  // *************************************
  //Z=3, Exec/De-Exec  & SPONTANEOUS EMISSION
  // ***************************************
#if MAXLEVEL > 2

  shift2=z0_len+z1_len+z2_len;
#ifdef OMP
//  #pragma omp parallel for schedule(dynamic,1) collapse(2) private(DeltaE,kfwd,krev,lambda0,Ajk,escape_factor,Eexc) num_threads(num_threads)
  #pragma omp for nowait schedule(static) private(j,DeltaE,kfwd,krev,lambda0,Ajk,escape_factor,Eexc)
#endif
  for(i=0;i<z3_len;++i)
  {
    for(j=0;j<z3_len;++j)
    {
      if(j<=i) continue; // MUI IMPORTANTE
      realtype engi=STATES_z3[i][2]*eV2J-IPD2+EF;
      if(engi < 0) continue; //depressed state is continuum

      DeltaE=(STATES_z3[j][2]-STATES_z3[i][2])*eV2J;
      kfwd=k_EE_z3_z3[i][j]*Ith(y,i+ishift+shift2)*ne;
      krev=k_EE_z3_z3_b[i][j]*Ith(y,j+ishift+shift2)*ne;

 // if(myid==1)
 //  printf("Z=3 kfwd:%.4e,krev:%.4e, yi:%.4e,yj:%.4e\n", kfwd,krev, Ith(y,i+ishift+shift2), Ith(y,j+ishift+shift2));

      //exec. reduces conc. of i state and increases conc. of j state
      Ith(colrad_ydot,i+ishift+shift2) -=kfwd;
      Ith(colrad_ydot,j+ishift+shift2) +=kfwd;

      //de-excec. increases conc. of i state & decr. conc. of j state
      Ith(colrad_ydot,i+ishift+shift2)+=krev;
      Ith(colrad_ydot,j+ishift+shift2)-=krev;


      Eexc= (-kfwd+krev)*DeltaE;
      P_E_EE+=Eexc;

// if(Ith(y,i+ishift+shift2) >0 || Ith(y,j+ishift+shift2)>0)
//   printf("myid:%d, i:%d,j:%d,kEE3:%.4e, kEE3b:%.4e,kfwd:%.4e,krev:%.4e,ni:%.4e,nj:%.4e,ne:%.4e\n",myid,i,j,
//   k_EE_z3_z3[i][j],k_EE_z3_z3_b[i][j],kfwd,krev,
//   Ith(y,i+ishift+shift2),Ith(y,j+ishift+shift2),ne );

#ifdef SPONT
      // ********SPONT EMISSION ********* //
      //Spont.emiss: delta-n !=>0 und delta-l=+-1 (nur opt.allowed trans.)
      if((STATES_z3[j][5]-STATES_z3[i][5])>0 && STATES_z3[j][4]-STATES_z3[i][4]==1)
      {
        escape_factor=1.0;
        lambda0=planck*LIGHTSPEED/DeltaE;
        Ajk=EinsteinCoeff(STATES_z3[i][5],STATES_z3[j][5],STATES_z3[j][3],DeltaE);

        krev=Ith(y,ishift+shift2+j)*Ajk;
        Ith(colrad_ydot,ishift+shift2+j) -= krev;
        Ith(colrad_ydot,ishift+shift2+i) += krev;
        //P_A_SE+=krev*DeltaE; ??
      }
#endif
    }
  }
#endif //MAXLEVEL > 2


  // *************************************
  //Z=4, Exec/De-Exec  & SPONTANEOUS EMISSION
  // ***************************************
#if MAXLEVEL > 3

  shift2=z0_len+z1_len+z2_len+z3_len;
#ifdef OMP
//  #pragma omp parallel for schedule(dynamic,1) collapse(2) private(DeltaE,kfwd,krev,lambda0,Ajk,escape_factor,Eexc) num_threads(num_threads)
  #pragma omp for nowait schedule(static) private(j,DeltaE,kfwd,krev,lambda0,Ajk,escape_factor,Eexc)
#endif
  for(i=0;i<z4_len;++i)
  {
    for(j=0;j<z4_len;++j)
    {
      if(j<=i) continue; // MUI IMPORTANTE      

      DeltaE=(STATES_z4[j][2]-STATES_z4[i][2])*eV2J-IPD3+EF;
      kfwd=k_EE_z4_z4[i][j]*Ith(y,i+ishift+shift2)*ne;
      krev=k_EE_z4_z4_b[i][j]*Ith(y,j+ishift+shift2)*ne;

      //exec. reduces conc. of i state and increases conc. of j state
      Ith(colrad_ydot,i+ishift+shift2) -=kfwd;
      Ith(colrad_ydot,j+ishift+shift2) +=kfwd;

      //de-excec. increases conc. of i state & decr. conc. of j state
      Ith(colrad_ydot,i+ishift+shift2)+=krev;
      Ith(colrad_ydot,j+ishift+shift2)-=krev;

// if(myid==1)
//  printf("Z=4:kfwd:%.4e,krev:%.4e, yi:%.4e,yj:%.4e\n", kfwd,krev, Ith(y,i+ishift+shift2), Ith(y,j+ishift+shift2));

      Eexc= (-kfwd+krev)*DeltaE;
      P_E_EE+=Eexc;

#ifdef SPONT
      // ********SPONT EMISSION ********* //
      //Spont.emiss: delta-n !=>0 und delta-l=+-1 (nur opt.allowed trans.)
      if((STATES_z4[j][5]-STATES_z4[i][5])>0 && STATES_z4[j][4]-STATES_z4[i][4]==1)
      {
        escape_factor=1.0;
        lambda0=planck*LIGHTSPEED/DeltaE;
        Ajk=EinsteinCoeff(STATES_z3[i][5],STATES_z3[j][5],STATES_z3[j][3],DeltaE);

        krev=Ith(y,ishift+shift2+j)*Ajk;
        Ith(colrad_ydot,ishift+shift2+j) -= krev;
        Ith(colrad_ydot,ishift+shift2+i) += krev;
        //P_A_SE+=krev*DeltaE; ??
      }
#endif
    }
  }
#endif // MAXLEVEL > 3  


  // ********************************************************
  //                    NOW IONIZATION RATES
  // ********************************************************
  // ***************************
  //Z=0->Z=1, Ioniz./Recomb.
  // ***************************
  shift1=ishift;
  shift2=ishift+z0_len;
#ifdef OMP
 // #pragma omp parallel for schedule(dynamic,1) collapse(2) private(DeltaE,kfwd,kfwd2,krev,\
                  escape_factor,groundstate_ioniz,sigma_PI,abs_coeff,tau,wstark_freq,lambda0,Eexc) num_threads(num_threads)
#pragma omp for nowait schedule(static) private(j,DeltaE,kfwd,kfwd2,krev,\
                  escape_factor,groundstate_ioniz,sigma_PI,abs_coeff,tau,wstark_freq,lambda0,Eexc)
#endif
  for(i=0;i<z0_len;++i)
  {
    for(j=0;j<z1_len;++j)
    {
      DeltaE=(STATES_z1[j][2]-STATES_z0[i][2])*eV2J-IPD0+EF;

#ifdef DOIPD
      if(DeltaE<0) continue;
      DeltaE=MAX(0.0,DeltaE);
#endif
      //COLL IONIZ
      Ith(colrad_ydot,i+shift1)          -= k_EI_z0_z1[i][j]*Ith(y,i+shift1)*ne;
      Ith(colrad_ydot,2)                 += k_EI_z0_z1[i][j]*Ith(y,i+shift1)*ne;   //Ne inc.
      Ith(colrad_ydot,j+shift2)          += k_EI_z0_z1[i][j]*Ith(y,i+shift1)*ne;

      //3-body recomb
      Ith(colrad_ydot,j+shift2)         -= k_EI_z1_z0[i][j]*Ith(y,j+shift2)*ne*ne;
      Ith(colrad_ydot,i+shift1)         += k_EI_z1_z0[i][j]*Ith(y,j+shift2)*ne*ne;
      Ith(colrad_ydot,2)                -= k_EI_z1_z0[i][j]*Ith(y,j+shift2)*ne*ne;

      Eexc  = -k_EI_z0_z1[i][j]*Ith(y,i+shift1)*ne*DeltaE;
      Eexc +=  k_EI_z1_z0[i][j]*Ith(y,j+shift2)*ne*ne*DeltaE;
      P_E_EI += Eexc;

      //////////
      // MPI
      //////////
      escape_factor=1.0;
#ifdef MULTIPHOTON
#ifdef STARK
      if(k_MPI_z1_z0[i][j][0]>0 && STATES_z1[j][5]>STATES_z0[i][5])
      {
        //außerdem geht die quasi-static. approx. nur mit nu-nl>0
        groundstate_ioniz=(STATES_z1[0][2]-STATES_z0[0][2])*eV2J;
        sigma_PI=sigma_tmp*POWR(DeltaE,2.5)/SQRTR(groundstate_ioniz);
        abs_coeff=sigma_PI*ne; //cross-section mit Pauli-blocking-faktor korrigieren!
        tau=abs_coeff*thickness;
        wstark_freq=StarkWidth(STATES_z1[j][5],STATES_z0[i][5],Te/11605,Ti/11605,Ti/11605,ne,totalc);
        lambda0=planck*LIGHTSPEED/DeltaE;
        wstark_len=wstark_freq*lambda0*lambda0/LIGHTSPEED;
        escape_factor=EscapeFactor(wstark_len*1e10,tau);
      }
#endif

//ACHTUNG: Diskussion ob MPI sich überhaupt für Potential-Lowering interessiert...?
if(DeltaE >0 ) //DeltaE ist bereits abzgl. IPD
{
      kfwd= k_MPI_z0_z1[i][j][0]*Ith(y,i+shift1); // Einheit: k_MPI = 1/s
      kfwd2=k_MPI_z0_z1[i][j][1]*Ith(y,i+shift1); // *wup;   // 2 photon- and 3-photon ionization!

      krev=k_MPI_z1_z0[i][j][0]*Ith(y,j+shift2)*ne;//*escape_factor; // *wlo;   // einheit: k_RADRECOMB= m^3/s

      //kfwd2=0.0; // 3-photon-ioniz. off
      Ith(colrad_ydot,shift1+i)        -=(kfwd +kfwd2);
      Ith(colrad_ydot,2)               +=(kfwd +kfwd2);   //Ne inc.
      Ith(colrad_ydot,shift2+j)        +=(kfwd +kfwd2);

      //krev=fmin(krev,kfwd);
      Ith(colrad_ydot,shift2+j)        -=krev;
      Ith(colrad_ydot,shift1+i)        +=krev;
      Ith(colrad_ydot,2)-=krev;

      // für 2 photonen und rad.recomb
      P_E_MPI2  += kfwd* (2.0*planck*LASERFREQ-(DeltaE)); //kein heating durch rad.recomb->Photon verschwindet (bisher ohne self-absorption)
      //jetzt für 3 photonen
      P_E_MPI3  += kfwd2*(3.0*planck*LASERFREQ-(DeltaE));
      //jetzt rad recomb
      P_E_RAD_RECOMB  -= krev*(DeltaE)*escape_factor;      //Radiative cooling, c.f. http://www.astronomy.ohio-state.edu/~dhw/A825/notes8.pdf S.2
}
#endif //MPI

    }
  }

  // ***************************
  //  Z=1->Z=2, Ioniz./Recomb.
  // ***************************
#if MAXLEVEL > 1

  shift1=ishift+z0_len;
  shift2=ishift+z0_len+z1_len;
#ifdef OMP
  //#pragma omp parallel for schedule(dynamic,1) collapse(2) private(DeltaE,kfwd,kfwd2,krev,\
                  escape_factor,groundstate_ioniz,sigma_PI,abs_coeff,tau,wstark_freq,lambda0,Eexc) num_threads(num_threads)
#pragma omp for nowait schedule(static) private(j,DeltaE,kfwd,kfwd2,krev,\
                  escape_factor,groundstate_ioniz,sigma_PI,abs_coeff,tau,wstark_freq,lambda0,Eexc)                  
#endif
  for(i=0;i<z1_len;++i)
  {
    for(j=0;j<z2_len;++j)
    {
      DeltaE=(STATES_z2[j][2]-STATES_z1[i][2])*eV2J-IPD1+EF;

#ifdef DOIPD
      if(DeltaE<0) continue;
      DeltaE=MAX(0.0,DeltaE);
#endif
      //COLL IONIZ
      Ith(colrad_ydot,i+shift1)          -= k_EI_z1_z2[i][j]*Ith(y,i+shift1)*ne;
      Ith(colrad_ydot,2)                 += k_EI_z1_z2[i][j]*Ith(y,i+shift1)*ne;   //Ne inc.
      Ith(colrad_ydot,j+shift2)          += k_EI_z1_z2[i][j]*Ith(y,i+shift1)*ne;

      //3-body recomb
      Ith(colrad_ydot,shift2+j)          -= k_EI_z2_z1[i][j]*Ith(y,j+shift2)*ne*ne;
      Ith(colrad_ydot,shift1+i)          += k_EI_z2_z1[i][j]*Ith(y,j+shift2)*ne*ne;
      Ith(colrad_ydot,2)                 -= k_EI_z2_z1[i][j]*Ith(y,j+shift2)*ne*ne;

      Eexc  = -k_EI_z1_z2[i][j]*Ith(y,i+shift1)*ne*DeltaE;
      Eexc +=  k_EI_z2_z1[i][j]*Ith(y,j+shift2)*ne*ne*DeltaE;
      P_E_EI += Eexc;
      //////////
      // MPI
      //////////
      escape_factor=1.0;
#ifdef MULTIPHOTON
#ifdef STARK
      if(k_MPI_z2_z1[i][j][0]>0 && STATES_z2[j][5]>STATES_z1[i][5])
      {
        //außerdem geht die quasi-static. approx. nur mit nu-nl>0
        groundstate_ioniz=(STATES_z2[0][2]-STATES_z1[0][2])*eV2J;
        sigma_PI=sigma_tmp*POWR(DeltaE,2.5)/SQRTR(groundstate_ioniz);
        abs_coeff=sigma_PI*ne; //cross-section mit Pauli-blocking-faktor korrigieren!
        tau=abs_coeff*thickness;
        wstark_freq=StarkWidth(STATES_z2[j][5],STATES_z1[i][5],Te/11605,Ti/11605,Ti/11605,ne,totalc);
        lambda0=planck*LIGHTSPEED/DeltaE;
        wstark_len=wstark_freq*lambda0*lambda0/LIGHTSPEED;
        escape_factor=EscapeFactor(wstark_len*1e10,tau);
      }
#endif

if(DeltaE >0)
{
      kfwd= k_MPI_z1_z2[i][j][0]*Ith(y,i+shift1); // *wup;
      kfwd2=k_MPI_z1_z2[i][j][1]*Ith(y,i+shift1); // *wup;   // 2 photon- and 3-photon ionization!

      krev=k_MPI_z2_z1[i][j][0]*Ith(y,shift2+j)*ne; // *wlo;   //eigentlich Rad-recomb.!=k_MPI_rev, einheit von k=m^3/s

      //kfwd2=0.0; // 3-photon-ioniz. off
      Ith(colrad_ydot,i+shift1)                -=(kfwd +kfwd2);
      Ith(colrad_ydot,2)                       +=(kfwd +kfwd2);   //Ne inc.
      Ith(colrad_ydot,shift2+j)                +=(kfwd +kfwd2);

      //krev=fmin(krev,kfwd);
      Ith(colrad_ydot,shift2+j)                -=krev;
      Ith(colrad_ydot,shift1+i)                +=krev;
      Ith(colrad_ydot,2)-=krev;

      // für 2 photonen und rad.recomb
      P_E_MPI2  += kfwd* (2.0*planck*LASERFREQ-(DeltaE)); //kein heating durch rad.recomb->Photon verschwindert (bisher ohne self-absorption)
      //jetzt für 3 photonen
      P_E_MPI3  += kfwd2*(3.0*planck*LASERFREQ-(DeltaE));
      //jetzt rad recomb
      P_E_RAD_RECOMB  -= krev*(DeltaE)*escape_factor;
}
#endif //MPI

    }
  }
#endif // MAXLEVEL > 1

  // ***************************
  //  Z=2->Z=3, Ioniz./Recomb.
  // ***************************
#if MAXLEVEL > 2

  shift1=ishift+z0_len+z1_len;
  shift2=ishift+z0_len+z1_len+z2_len;
#ifdef OMP
  //#pragma omp parallel for schedule(dynamic,1) collapse(2) private(DeltaE,kfwd,kfwd2,krev,\
                  escape_factor,groundstate_ioniz,sigma_PI,abs_coeff,tau,wstark_freq,lambda0,Eexc) num_threads(num_threads)
#pragma omp for nowait schedule(static) private(j,DeltaE,kfwd,kfwd2,krev,\
                  escape_factor,groundstate_ioniz,sigma_PI,abs_coeff,tau,wstark_freq,lambda0,Eexc)                  
#endif
  for(i=0;i<z2_len;++i)
  {
    for(j=0;j<z3_len;++j)
    {
      DeltaE=(STATES_z3[j][2]-STATES_z2[i][2])*eV2J-IPD2+EF;

#ifdef DOIPD
      if(DeltaE<0) continue;
      DeltaE=MAX(0.0,DeltaE);
#endif
      //COLL IONIZ
      Ith(colrad_ydot,i+shift1)          -= k_EI_z2_z3[i][j]*Ith(y,i+shift1)*ne;
      Ith(colrad_ydot,2)                 += k_EI_z2_z3[i][j]*Ith(y,i+shift1)*ne;   //Ne inc.
      Ith(colrad_ydot,j+shift2)          += k_EI_z2_z3[i][j]*Ith(y,i+shift1)*ne;

      //3-body recomb
      Ith(colrad_ydot,shift2+j)          -= k_EI_z3_z2[i][j]*Ith(y,j+shift2)*ne*ne;
      Ith(colrad_ydot,shift1+i)          += k_EI_z3_z2[i][j]*Ith(y,j+shift2)*ne*ne;
      Ith(colrad_ydot,2)                 -= k_EI_z3_z2[i][j]*Ith(y,j+shift2)*ne*ne;

      Eexc  = -k_EI_z2_z3[i][j]*Ith(y,i+shift1)*ne*DeltaE;
      Eexc +=  k_EI_z3_z2[i][j]*Ith(y,j+shift2)*ne*ne*DeltaE;
      P_E_EI += Eexc;
      //////////
      // MPI
      //////////
      escape_factor=1.0;
#ifdef MULTIPHOTON
#ifdef STARK
      if(k_MPI_z3_z2[i][j][0]>0 && STATES_z3[j][5]>STATES_z2[i][5])
      {
        //außerdem geht die quasi-static. approx. nur mit nu-nl>0
        groundstate_ioniz=(STATES_z3[0][2]-STATES_z2[0][2])*eV2J;
        sigma_PI=sigma_tmp*POWR(DeltaE,2.5)/SQRTR(groundstate_ioniz);
        abs_coeff=sigma_PI*ne; //cross-section mit Pauli-blocking-faktor korrigieren!
        tau=abs_coeff*thickness;
        wstark_freq=StarkWidth(STATES_z3[j][5],STATES_z2[i][5],Te/11605,Ti/11605,Ti/11605,ne,totalc);
        lambda0=planck*LIGHTSPEED/DeltaE;
        wstark_len=wstark_freq*lambda0*lambda0/LIGHTSPEED;
        escape_factor=EscapeFactor(wstark_len*1e10,tau);
      }
#endif

if(DeltaE >0)
{
      kfwd= k_MPI_z2_z3[i][j][0]*Ith(y,i+shift1); // *wup;
      kfwd2=k_MPI_z3_z2[i][j][1]*Ith(y,i+shift1); // *wup;   // 2 photon- and 3-photon ionization!

      krev=k_MPI_z3_z2[i][j][0]*Ith(y,shift2+j)*ne; // *wlo;   //eigentlich Rad-recomb.!=k_MPI_rev, einheit von k=m^3/s

      //kfwd2=0.0; // 3-photon-ioniz. off
      Ith(colrad_ydot,i+shift1)                -=(kfwd +kfwd2);
      Ith(colrad_ydot,2)                       +=(kfwd +kfwd2);   //Ne inc.
      Ith(colrad_ydot,shift2+j)                +=(kfwd +kfwd2);

      //krev=fmin(krev,kfwd);
      Ith(colrad_ydot,shift2+j)                -=krev;
      Ith(colrad_ydot,shift1+i)                +=krev;
      Ith(colrad_ydot,2)-=krev;

      // für 2 photonen und rad.recomb
      P_E_MPI2  += kfwd* (2.0*planck*LASERFREQ-(DeltaE)); //kein heating durch rad.recomb->Photon verschwindert (bisher ohne self-absorption)
      //jetzt für 3 photonen
      P_E_MPI3  += kfwd2*(3.0*planck*LASERFREQ-(DeltaE));
      //jetzt rad recomb
      P_E_RAD_RECOMB  -= krev*(DeltaE)*escape_factor;
}
#endif //MPI

    }
  }
#endif // MAXLEVEL > 2

  // ***************************
  //  Z=3->Z=4, Ioniz./Recomb.
  // ***************************
#if MAXLEVEL > 3

  shift1=ishift+z0_len+z1_len+z2_len;
  shift2=ishift+z0_len+z1_len+z2_len+z3_len;
#ifdef OMP
 //#pragma omp parallel for schedule(dynamic,1) collapse(2) private(DeltaE,kfwd,kfwd2,krev,\
                   escape_factor,groundstate_ioniz,sigma_PI,abs_coeff,tau,wstark_freq,lambda0,Eexc) num_threads(num_threads)
#pragma omp for nowait schedule(static) private(j,DeltaE,kfwd,kfwd2,krev,\
                  escape_factor,groundstate_ioniz,sigma_PI,abs_coeff,tau,wstark_freq,lambda0,Eexc)                   
#endif
  for(i=0;i<z3_len;++i)
  {
    for(j=0;j<z4_len;++j)
    {
      DeltaE=(STATES_z4[j][2]-STATES_z3[i][2])*eV2J-IPD3+EF;
#ifdef DOIPD
      if(DeltaE<0) continue;
      DeltaE=MAX(0.0,DeltaE);
#endif
      //COLL IONIZ
      Ith(colrad_ydot,i+shift1)          -= k_EI_z3_z4[i][j]*Ith(y,i+shift1)*ne;
      Ith(colrad_ydot,2)                 += k_EI_z3_z4[i][j]*Ith(y,i+shift1)*ne;   //Ne inc.
      Ith(colrad_ydot,j+shift2)          += k_EI_z3_z4[i][j]*Ith(y,i+shift1)*ne;

      //3-body recomb
      Ith(colrad_ydot,shift2+j)          -= k_EI_z4_z3[i][j]*Ith(y,j+shift2)*ne*ne;
      Ith(colrad_ydot,shift1+i)          += k_EI_z4_z3[i][j]*Ith(y,j+shift2)*ne*ne;
      Ith(colrad_ydot,2)                 -= k_EI_z4_z3[i][j]*Ith(y,j+shift2)*ne*ne;

      Eexc  = -k_EI_z3_z4[i][j]*Ith(y,i+shift1)*ne*DeltaE;
      Eexc +=  k_EI_z4_z3[i][j]*Ith(y,j+shift2)*ne*ne*DeltaE;
      P_E_EI += Eexc;

      //////////
      // MPI
      //////////
      escape_factor=1.0;
#ifdef MULTIPHOTON
#ifdef STARK
      if(k_MPI_z4_z3[i][j][0]>0 && STATES_z4[j][5]>STATES_z3[i][5])
      {
        //außerdem geht die quasi-static. approx. nur mit nu-nl>0
        groundstate_ioniz=(STATES_z4[0][2]-STATES_z3[0][2])*eV2J;
        sigma_PI=sigma_tmp*POWR(DeltaE,2.5)/SQRTR(groundstate_ioniz);
        abs_coeff=sigma_PI*ne; //cross-section mit Pauli-blocking-faktor korrigieren!
        tau=abs_coeff*thickness;
        wstark_freq=StarkWidth(STATES_z4[j][5],STATES_z3[i][5],Te/11605,Ti/11605,Ti/11605,ne,totalc);
        lambda0=planck*LIGHTSPEED/DeltaE;
        wstark_len=wstark_freq*lambda0*lambda0/LIGHTSPEED;
        escape_factor=EscapeFactor(wstark_len*1e10,tau);
      }
#endif

if(DeltaE >0)
{
      kfwd= k_MPI_z3_z4[i][j][0]*Ith(y,i+shift1); // *wup;
      kfwd2=k_MPI_z4_z3[i][j][1]*Ith(y,i+shift1); // *wup;   // 2 photon- and 3-photon ionization!

      krev=k_MPI_z4_z3[i][j][0]*Ith(y,shift2+j)*ne; // *wlo;   //eigentlich Rad-recomb.!=k_MPI_rev, einheit von k=m^3/s

      //kfwd2=0.0; // 3-photon-ioniz. off
      Ith(colrad_ydot,i+shift1)                -=(kfwd +kfwd2);
      Ith(colrad_ydot,2)                       +=(kfwd +kfwd2);   //Ne inc.
      Ith(colrad_ydot,shift2+j)                +=(kfwd +kfwd2);

      //krev=fmin(krev,kfwd);
      Ith(colrad_ydot,shift2+j)                -=krev;
      Ith(colrad_ydot,shift1+i)                +=krev;
      Ith(colrad_ydot,2)-=krev;

      // für 2 photonen und rad.recomb
      P_E_MPI2  += kfwd* (2.0*planck*LASERFREQ-(DeltaE)); //kein heating durch rad.recomb->Photon verschwindert (bisher ohne self-absorption)
      //jetzt für 3 photonen
      P_E_MPI3  += kfwd2*(3.0*planck*LASERFREQ-(DeltaE));
      //jetzt rad recomb
      P_E_RAD_RECOMB  -= krev*(DeltaE)*escape_factor;
}
#endif //MPI

    }
  }
#endif //MAXLEVEL > 3  

  // ********************** THERMO ******************************************

  // double cvinv=1.0/(1.5*BOLTZMAN*ne);

  realtype P_E_TOTAL=P_E_EI+P_E_EE+P_E_MPI2+P_E_MPI3+P_E_RAD_RECOMB;
  data->P_TOTAL=P_E_TOTAL;
  data->P_EE=P_E_EE;
  data->P_EI=P_E_EI;
  data->P_MPI2=P_E_MPI2;
  data->P_MPI3=P_E_MPI3;
  data->P_RR=P_E_RAD_RECOMB;

  // printf("myid:%d, PEI:%.4e, PEE:%.4e,MPI2:%.4e, MPI3:%.4e, RADREC:%.4e\n",
  //         myid,
  //         P_E_EI, 
  //         P_E_EE, 
  //         P_E_MPI2, 
  //         P_E_MPI3, 
  //         P_E_RAD_RECOMB);

  //BEI PRE-EQUILIBRIERUNG T=CONST !
  if(initial_equi==false)
  {
    realtype cv=EOS_cve_from_r_te(data->dens, Te);  

    cv *= 1e30/11604.5*eV2J;   // von eV/(eV*A^3) to  J/(K*m^3) 
    realtype cvinv=1.0/cv;    

    // double cvinv=  1.0/Cv(Te/11604.5, ne);
    //double cvinv=1.0/(1.5*BOLTZMAN*Te);
    Ith(colrad_ydot,0) =  cvinv*P_E_TOTAL;
    Ith(colrad_ydot,0) += data->Heatrate;

// if(myid==1)
//   for(i=0;i<neq;i++)
//   {
//     printf("theta:%.4e\n", cvinv*P_E_TOTAL); 
//   }

  }
  else
  {
    Ith(colrad_ydot,0)=0.0;
  }




  return 0; // 0 heisst alles ok
}


// ********************************************************************+
// *                COMPUTE RATE COEFFS
// ********************************************************************+
int colrad_GetCoeffs(N_Vector y,realtype It,void *user_data)
{
  int i,j,k;
  realtype kronecker;
  realtype a;
  realtype DeltaE;
  realtype Te,ne;
  realtype expint;
  realtype G2;
  realtype I_1,I_2;

  Te=Ith(y,0);
  if(Te <0 || isnan(Te)!=0) return -1;
  ne=Ith(y,2);
  if(ne <0 || isnan(ne)!=0) return -1;

  // double v_e=SQRTR(8.0*BOLTZMAN*Te/pi/EMASS);
  realtype E_ion_H=13.6*eV2J;
  // double alpha_i=0.05;
  // double alpha_e=0.05;
  // double beta_i=4.0;
  // double four_pi_a0_sq=4.0*pi*POWR(bohr_radius,2.0);
  

  // double E_ion_H_div_kTe_sq=POWR((E_ion_H/BOLTZMAN/Te),2.0);
  // double two_pi_me_kT_hsq=2.0*pi*EMASS*BOLTZMAN*Te/POWR(planck,2.0);
  // double log54_beta_i=log(5.0/4.0*beta_i);
  realtype kbTe=(BOLTZMAN*Te);
  kbTe=Te/11604.5;

  realtype tmp0,tmp1,tmp2;

  colrad_UserData data;
  data = (colrad_UserData) user_data;
  realtype IPD0=data->IPD0*J2eV;
  realtype IPD1=data->IPD1*J2eV;
  realtype IPD2=data->IPD2*J2eV;
  realtype IPD3=data->IPD3*J2eV;
  realtype IPD4=data->IPD4*J2eV;
  realtype EF=data->EF*J2eV;

  realtype Tinit=data->Tinit;
  bool  initial_equi = data->initial_equi;

  if(initial_equi)
  {
    if(ABS(Te-Tinit) > Tinit*0.03)  // cvode war zu eifrig und probiert extreme temp. aus
      return -1;
  }
  else
  {
    // if(ABS(Te-Tinit) > Tinit*0.5)  // cvode war zu eifrig und probiert extreme temp. aus
    //   return -1;    
    if(Te>1e4)  //ACHTUNG: Keine gute lösung...eher hotfix
      return -1;
  }

  //MPI
  // realtype k_RR_fact1=32*pi*POWR(bohr_radius,2.0)/3.0/175700.00067;
  // realtype k_RR_fact2=POWR((E_ion_H/BOLTZMAN/Te),2.0);
  // realtype sigma_MPI_2;//=sigma1/LASERFREQ/POWR(planck*LASERFREQ,2.0); //MPI-cross-sect. (2-photon)
  // realtype sigma_MPI_3;
  // realtype sigma1;

  //const. zur berechnung von sigma1
  // realtype sigma_tmp=64.0*POWR(pi,4.0)*POWR(ECHARGE,10.0)*EMASS/3.0/SQRTR(3.0)/POWR(4.0*pi*ECONST,5.0)/POWR(planck,6.0)/LIGHTSPEED/POWR(LASERFREQ,3.0)/POWR(13.6*eV2J,2.0);
  // realtype I_sq=It*It; //for 2-photon-ioniz.
  // realtype I_cu=I_sq*It; // for 3-photon-ioniz.
  int fail=0;

  // realtype pow_two_pi_me_kT_hsq_tmp1= POWR(two_pi_me_kT_hsq,1.5); //ACHTUNG: Das ist thermal De-Broglie Lambda
                                                               //Ich muss das für die rückwärts-raten irgendie 
                                                               //durch chempot ersetzen sonst pass das net
                                                               //Ebenso in Saha
#ifdef MULTIPHOTON
  realtype twophoton_energy=2.0*planck*LASERFREQ*J2eV;
  realtype threephoton_energy=3.0*planck*LASERFREQ*J2eV;
  realtype nu_div_hnu_sq=LASERFREQ/POWR(planck*LASERFREQ,2.0);
  realtype nu_div_nu_div_hnu_cub=LASERFREQ/LASERFREQ/POWR(planck*LASERFREQ,3.0);
#endif


  realtype mu=chempot(ne,Te);
  realtype mu_eV=mu*J2eV;
  //realtype fermi_factor=eval_fermi_integrand(ne,Te,mu);
  realtype fermi_factor=1.0;
  if(fermi_factor==-1)
    return -1;

  // realtype kmax_estim=1e4; 
  realtype nesq=gsl_pow_2(ne);


  //PREZERO RATE-COEFFS CODEBLOCK 
  {
    //pre-zero k_EE's
    for(i=0;i<z0_len;i++)
    {
            for(j=0;j<z0_len;j++)
            {
                    k_EE_z0_z0[i][j]=0.0;
                    k_EE_z0_z0_b[i][j]=0.0;
            }
    }

    for(i=0;i<z1_len;i++)
    {
            for(j=0;j<z1_len;j++)
            {
                    k_EE_z1_z1[i][j]=0.0;
                    k_EE_z1_z1_b[i][j]=0.0;
            }
    }

    #if MAXLEVEL > 1
    for(i=0;i<z2_len;i++)
    {
            for(j=0;j<z2_len;j++)
            {
                    k_EE_z2_z2[i][j]=0.0;
                    k_EE_z2_z2_b[i][j]=0.0;
            }
    }
    #endif

    #if MAXLEVEL > 2 
    for(i=0;i<z3_len;i++)
    {
            for(j=0;j<z3_len;j++)
            {
                    k_EE_z3_z3[i][j]=0.0;
                    k_EE_z3_z3_b[i][j]=0.0;
            }
    }
    #endif

    #if MAXLEVEL > 3
    for(i=0;i<z4_len;i++)
    {
            for(j=0;j<z4_len;j++)
            {
                    k_EE_z4_z4[i][j]=0.0;
                    k_EE_z4_z4_b[i][j]=0.0;
            }
    }
    #endif


    //pre-zero k_MPI's
    for(i=0;i<z0_len;i++)
    {
            for(j=0;j<z1_len;j++)
            {
                    k_EI_z0_z1[i][j]=0.0;
                    k_EE_z1_z1[i][j]=0.0;
                    for(k=0;k<2;k++)
                    {
                      k_MPI_z0_z1[i][j][k]=0.0;
                      k_MPI_z1_z0[i][j][k]=0.0;
                    }
            }
    }

    #if MAXLEVEL > 1  
    for(i=0;i<z1_len;i++)
    {
            for(j=0;j<z2_len;j++)
            {
                    k_EI_z1_z2[i][j]=0.0;
                    k_EI_z2_z1[i][j]=0.0;
                    for(k=0;k<2;k++)
                    {
                      k_MPI_z1_z2[i][j][k]=0.0;
                      k_MPI_z2_z1[i][j][k]=0.0;
                    }
            }
    }
    #endif

    #if MAXLEVEL > 2 
     for(i=0;i<z2_len;i++)
    {
            for(j=0;j<z3_len;j++)
            {
                    k_EI_z2_z3[i][j]=0.0;
                    k_EI_z3_z2[i][j]=0.0;
                    for(k=0;k<2;k++)
                    {
                      k_MPI_z2_z3[i][j][k]=0.0;
                      k_MPI_z3_z2[i][j][k]=0.0;
                    }
            }
    }
    #endif  

    #if MAXLEVEL > 3  
    for(i=0;i<z3_len;i++)
    {
            for(j=0;j<z4_len;j++)
            {
                    k_EI_z3_z4[i][j]=0.0;
                    k_EI_z4_z3[i][j]=0.0;
                    for(k=0;k<2;k++)
                    {
                      k_MPI_z3_z4[i][j][k]=0.0;
                      k_MPI_z4_z3[i][j][k]=0.0;
                    }
            }
    }
    #endif
  }

//ACHTUNG: Nur um initial equi zu überspringen!
 // if(initial_equi==true)
 // return 0; 
  ///////////////////////////////
  // Elec. excitation for Z=0
  ///////////////////////////////
  fail=0;
#ifdef OMP
//#pragma omp parallel for schedule(dynamic,1) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2) num_threads(num_threads)
  #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,fparams_exc, winteg_exc)
  // #pragma omp for simd schedule(static) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2)
#endif
  
  for(i=0;i<z0_len;++i)
  {
    for(j=0;j<z0_len;++j)
    {
        if(j<=i) continue;
        kronecker=0.0; //optically allowed transition
        if(STATES_z0[i][4]==STATES_z0[j][4])  // l_j==l_i ?
          kronecker=1.0;//optically forbidden transition

          DeltaE=(STATES_z0[j][2]-STATES_z0[i][2]);

#ifdef DOIPD
          realtype Ei=(STATES_z0[i][2])-IPD0+mu_eV;//+EF;
          if(Ei<0) continue;
#endif          
        
          k_EE_z0_z0[i][j]=eval_excitation_integral(ne,Te,mu,DeltaE*eV2J,kronecker); // in m^3/s    
          k_EE_z0_z0_b[i][j]=eval_dexcitation_integral(ne,Te,mu,DeltaE*eV2J,kronecker)*STATES_z0[i][3]/STATES_z0[j][3];

          k_EE_MAX=MAX(k_EE_z0_z0[i][j],k_EE_MAX);
          k_EE_REV_MAX=MAX(k_EE_z0_z0_b[i][j],k_EE_REV_MAX);          

      }
    }

if(fail==1)
  return -1;

 
  ////////////////////////
  // Elec. exc. for Z=1
  ///////////////////////
  fail=0;
#ifdef OMP
  #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,fparams_exc,winteg_exc)
#endif
  for(i=0;i<z1_len;++i)
  {
    for(j=0;j<z1_len;++j)
    {
        if(j<=i) continue;
        kronecker=0.0;
        
        if(STATES_z1[i][4]==STATES_z1[j][4])
          kronecker=1.0;
#ifdef DOIPD
        //realtype Ei=(STATES_z1[i][2]-STATES_z1[0][2])-IPD1+EF;
        realtype Ei=STATES_z1[i][2]-IPD0+mu_eV;
        if(Ei<0) continue;
#endif

        DeltaE=(STATES_z1[j][2]-STATES_z1[i][2]);
        k_EE_z1_z1[i][j]=eval_excitation_integral(ne,Te,mu,DeltaE*eV2J,kronecker); // in m^3/s
        k_EE_z1_z1_b[i][j]=eval_dexcitation_integral(ne,Te,mu,DeltaE*eV2J,kronecker)*STATES_z1[i][3]/STATES_z1[j][3];

        k_EE_MAX=MAX(k_EE_z1_z1[i][j],k_EE_MAX);
        k_EE_REV_MAX=MAX(k_EE_z1_z1_b[i][j],k_EE_REV_MAX);        
      }
    }
  if(fail==1) return -1;


  ////////////////////////
  // Elec. exc. for Z=2
  ///////////////////////
#if MAXLEVEL > 1

  fail=0;
#ifdef OMP
  #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,fparams_exc,winteg_exc)
#endif
  for(i=0;i<z2_len;++i)
  {
    for(j=0;j<z2_len;++j)
    {
        if(j<=i) continue;
        kronecker=0.0;
        
        if(STATES_z2[i][4]==STATES_z2[j][4])
          kronecker=1.0;

#ifdef DOIPD
        // realtype Ei=(STATES_z2[i][2]-STATES_z2[0][2])-IPD2+EF;
        realtype Ei=STATES_z2[i][2]-IPD1+mu_eV;
        if(Ei<0) continue;
#endif

        DeltaE=(STATES_z2[j][2]-STATES_z2[i][2]);
        k_EE_z2_z2[i][j]=eval_excitation_integral(ne,Te,mu,DeltaE*eV2J,kronecker);         
        k_EE_z2_z2_b[i][j]=eval_dexcitation_integral(ne,Te,mu,DeltaE*eV2J,kronecker)*STATES_z2[i][3]/STATES_z2[j][3];

        k_EE_MAX=MAX(k_EE_z2_z2[i][j],k_EE_MAX);
        k_EE_REV_MAX=MAX(k_EE_z2_z2_b[i][j],k_EE_REV_MAX);        
        
      }
    }
  if(fail==1) return -1;

#endif //MAXLEVEL > 1

  ////////////////////////
  // Elec. exc. for Z=3
  ///////////////////////
#if MAXLEVEL > 2

  fail=0;
  #ifdef OMP
  #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,fparams_exc,winteg_exc)
  #endif
  for(i=0;i<z3_len;++i)
  {
    for(j=0;j<z3_len;++j)
    {
        if(j<=i) continue;
        kronecker=0.0;
        
        if(STATES_z3[i][4]==STATES_z3[j][4])
          kronecker=1.0;

#ifdef DOIPD
        // realtype Ei=(STATES_z3[i][2]-STATES_z3[0][2])-IPD3+EF;
        realtype Ei=STATES_z3[i][2]-IPD2+mu_eV;
        if(Ei<0) continue;
#endif
        DeltaE=(STATES_z3[j][2]-STATES_z3[i][2]);
        k_EE_z3_z3[i][j]=eval_excitation_integral(ne,Te,mu,DeltaE*eV2J,kronecker); 
        k_EE_z3_z3_b[i][j]=eval_dexcitation_integral(ne,Te,mu,DeltaE*eV2J,kronecker)*STATES_z3[i][3]/STATES_z3[j][3];  

        k_EE_MAX=MAX(k_EE_z3_z3[i][j],k_EE_MAX);
        k_EE_REV_MAX=MAX(k_EE_z3_z3_b[i][j],k_EE_REV_MAX);

// printf("myid:%d,kfw:%.4e,krev:%.4e,i:%d,j:%d\n",myid,k_EE_z3_z3[i][j],k_EE_z3_z3_b[i][j],i,j);  
      }
    }
  if(fail==1)
          return -1;  
#endif // MAXLEVEL > 2


  ////////////////////////
  // Elec. exc. for Z=4
  ///////////////////////
#if MAXLEVEL > 3

  fail=0;
#ifdef OMP
  #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,fparams_exc,winteg_exc)
#endif
  for(i=0;i<z4_len;++i)
  {
    for(j=0;j<z4_len;++j)
    {
        if(j<=i) continue;
        kronecker=0.0;
        
        if(STATES_z4[i][4]==STATES_z4[j][4])
          kronecker=1.0;
#ifdef DOIPD
        // realtype Ei=(STATES_z4[i][2]-STATES_z4[0][2])-IPD4+EF;
        realtype Ei=STATES_z4[i][2]-IPD3+mu_eV;
        if(Ei<0) continue;
#endif          
        DeltaE=(STATES_z4[j][2]-STATES_z4[i][2]);
        k_EE_z4_z4[i][j]=eval_excitation_integral(ne,Te,mu,DeltaE*eV2J,kronecker); 
        k_EE_z4_z4_b[i][j]=eval_dexcitation_integral(ne,Te,mu,DeltaE*eV2J,kronecker);        

        k_EE_MAX=MAX(k_EE_z4_z4[i][j],k_EE_MAX);
        k_EE_REV_MAX=MAX(k_EE_z4_z4_b[i][j],k_EE_REV_MAX);        
      }
    }
  if(fail==1)
          return -1;
#endif // MAXLEVEL > 3


  // *************************************************
  // *           NOW IONIZATION COEFFS
  // *************************************************
  /////////////////
  // Ioniz. 0->1
  /////////////////
        
  fail=0;
#ifdef OMP
   // #pragma omp parallel for schedule(static) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2,sigma1,sigma_MPI_2,sigma_MPI_3,tmp0,tmp1,tmp2) num_threads(num_threads)
   // #pragma omp for nowait schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2,sigma1,sigma_MPI_2,sigma_MPI_3,tmp0,tmp1,tmp2)
  // #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,fparams_inner, fparams_outer, winteg_inner, winteg_outer,winteg_rb_inner,winteg_rb_outer)
#endif
  for(i=0;i<z0_len;++i)
  {
    for(j=0;j<z1_len;++j)
    {

      DeltaE=(STATES_z1[j][2]-STATES_z0[i][2])-IPD0;//+EF;
if(DeltaE <0 )
  continue;


if(RATEIONIZMAX*ne*fermi_factor*Ith(y,i+3)> MINRATE)
      k_EI_z0_z1[i][j]=MAX(0.0,double_integral_ionization(ne,Te, mu, DeltaE*eV2J));

if(RATERECOMBMAX*nesq*fermi_factor*Ith(y,j+z0_len+3) > MINRATE)
      k_EI_z1_z0[i][j]=STATES_z0[i][3]/STATES_z1[j][3]*double_integral_recombination(ne,Te, mu, DeltaE*eV2J);      

k_EI_MAX=MAX(k_EI_MAX,k_EI_z0_z1[i][j]);
k_EI_REV_MAX=MAX(k_EI_REV_MAX,k_EI_z1_z0[i][j]);
  
      k_EI_z0_z1[i][j]*=fermi_factor;
      k_EI_z1_z0[i][j]*=fermi_factor;

#ifdef MULTIPHOTON
      // *******************
      // * MPI 2 PHOTONS   *
      // *******************
      realtype dE_SI=DeltaE*eV2J; 
      if(twophoton_energy >= DeltaE-IPD0 && DeltaE-IPD0 > 0.0 )
      {
        sigma1=sigma_tmp*POWR(dE_SI,2.5)/SQRTR(dE_SI);
        sigma_MPI_2=sigma1*sigma1/nu_div_hnu_sq;
        k_MPI_z0_z1[i][j][0]=sigma_MPI_2*I_sq;
      }
      // *******************
      // * MPI 3 PHOTONS   *
      // *******************   
      if(threephoton_energy >= DeltaE-IPD0 && DeltaE-IPD0 >0.0 ) 
      {
        sigma1=sigma_tmp*POWR(dE_SI,2.5)/SQRTR(dE_SI);
        sigma_MPI_3=sigma1*sigma1*sigma1/2.0/nu_div_nu_div_hnu_cub;
        k_MPI_z0_z1[i][j][1]=sigma_MPI_3*I_cu;//*prob*beta_pi(Te,mu,DeltaE);
      }
#endif      
      // ***************
      // * RAD RECOMB  *
      // ***************
      // if(DeltaE>0)
      // {
      //   if(expint > 0 )
      //   {
      //     //k_MPI_z1_z0[i][j][0]=v_e*k_RR_fact1*1.0*k_RR_fact2*POWR((DeltaE)*J2eV/STATES_z1[j][2],1.5)*expint*exp(a);
      //     k_MPI_z1_z0[i][j][0]=v_e*k_RR_fact1*1.0*k_RR_fact2*POWR((DeltaE)/STATES_z1[j][2],1.5)*expint*exp(a);
      //   }
      //   else 
      //   {
      //     k_MPI_z1_z0[i][j][0]=0.0; //exp(a) kann +Inf werden--> unsinnige rate coeff.. wenn einer der faktoren=0 --> rest egal
      //   }
      // }

    }
  }
  if(fail==1) return -1;


   /////////////////
  // Ioniz. 1->2
  ///////////////// 
#if MAXLEVEL > 1
  fail=0;
#ifdef OMP
 // #pragma omp parallel for schedule(static) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2,sigma1,sigma_MPI_2,sigma_MPI_3,tmp0,tmp1,tmp2) num_threads(num_threads)
  // #pragma omp for nowait schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2,sigma1,sigma_MPI_2,sigma_MPI_3,tmp0,tmp1,tmp2)
  // #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,fparams_inner, fparams_outer, winteg_inner, winteg_outer,winteg_rb_inner,winteg_rb_outer)
#endif
  for(i=0;i<z1_len;++i)
  {
    for(j=0;j<z2_len;++j)
    {

      DeltaE=(STATES_z2[j][2]-STATES_z1[i][2])-IPD1; //+EF;

if(DeltaE <0 )
  continue;

if(RATEIONIZMAX*ne*fermi_factor*Ith(y,i+z0_len+3)> MINRATE)
      k_EI_z1_z2[i][j]=MAX(0.0,double_integral_ionization(ne,Te, mu, DeltaE*eV2J));

// if(kmax_estim*ne*ne*Ith(y,j+z0_len+z1_len+3)>MINRATE)
if(RATERECOMBMAX*nesq*fermi_factor*Ith(y,j+z0_len+z1_len+3) > MINRATE)      
      k_EI_z2_z1[i][j]=STATES_z1[i][3]/STATES_z2[j][3]*double_integral_recombination(ne,Te, mu, DeltaE*eV2J);

k_EI_MAX=MAX(k_EI_MAX,k_EI_z1_z2[i][j]);
k_EI_REV_MAX=MAX(k_EI_REV_MAX,k_EI_z2_z1[i][j]);

      k_EI_z1_z2[i][j]*=fermi_factor;
      k_EI_z2_z1[i][j]*=fermi_factor;



#ifdef MULTIPHOTON
      realtype dE_SI=DeltaE*eV2J; 
      // *******************
      // * MPI 2 PHOTONS   *
      // *******************
      if(twophoton_energy >= DeltaE - IPD1)
      {
        sigma1=sigma_tmp*POWR(dE_SI,2.5)/SQRTR(dE_SI);
        sigma_MPI_2=sigma1*sigma1/nu_div_hnu_sq;
        k_MPI_z1_z2[i][j][0]=sigma_MPI_2*I_sq;
      }
      // *******************
      // * MPI 3 PHOTONS   *
      // *******************  
      if(threephoton_energy > DeltaE - IPD1)
      {
        sigma1=sigma_tmp*POWR(dE_SI,2.5)/SQRTR(dE_SI);
        sigma_MPI_3=sigma1*sigma1*sigma1/2.0/nu_div_nu_div_hnu_cub;
        k_MPI_z1_z2[i][j][1]=sigma_MPI_3*I_cu;//*prob*beta_pi(Te,mu,DeltaE);
      }
#endif      
      // **************
      // * RAD RECOMB *
      // **************
      // if(DeltaE>0)
      // {

      //   if(expint > 0 )
      //   {
      //     k_MPI_z2_z1[i][j][0]=v_e*k_RR_fact1*4.0*k_RR_fact2*POWR((DeltaE)/STATES_z2[j][2],1.5)*expint*exp(a);          
      //   }
      //   else 
      //   {
      //     k_MPI_z2_z1[i][j][0]=0.0; //exp(a) kann +Inf werden--> unsinnige rate coeff.. wenn einer der faktoren=0 --> rest egal
      //   }      
      // }      

    }
  }

  if(fail==1) return -1;

#endif //MAXLEVEL > 1

  /////////////////
  // Ioniz. 2->3
  ///////////////// 
#if MAXLEVEL > 2

  fail=0;
#ifdef OMP
  // #pragma omp parallel for   schedule(static) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2,sigma1,sigma_MPI_2,sigma_MPI_3,tmp0,tmp1,tmp2) num_threads(num_threads)
  // #pragma omp for nowait schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2,sigma1,sigma_MPI_2,sigma_MPI_3,tmp0,tmp1,tmp2)
  // #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,fparams_inner, fparams_outer, winteg_inner, winteg_outer,winteg_rb_inner,winteg_rb_outer)
#endif
  for(i=0;i<z2_len;++i)
  {
    for(j=0;j<z3_len;++j)
    {

      DeltaE=(STATES_z3[j][2]-STATES_z2[i][2])-IPD2;//+EF;

if(DeltaE <0 )
  continue;

if(RATEIONIZMAX*ne*fermi_factor*Ith(y,i+z0_len+z1_len+3)> MINRATE)
      k_EI_z2_z3[i][j]=MAX(0.0,double_integral_ionization(ne,Te, mu, DeltaE*eV2J));

// if(kmax_estim*ne*ne*Ith(y,j+z0_len+z1_len+z2_len+3)>MINRATE)    
if(RATERECOMBMAX*nesq*fermi_factor*Ith(y,j+z0_len+z1_len+z2_len+3) > MINRATE)
      k_EI_z3_z2[i][j]=STATES_z2[i][3]/STATES_z3[j][3]*double_integral_recombination(ne,Te, mu, DeltaE*eV2J);

k_EI_MAX=MAX(k_EI_MAX,k_EI_z1_z2[i][j]);
k_EI_REV_MAX=MAX(k_EI_REV_MAX,k_EI_z2_z1[i][j]);

      k_EI_z2_z3[i][j]*=fermi_factor;
      k_EI_z3_z2[i][j]*=fermi_factor;

#ifdef MULTIPHOTON
      realtype dE_SI=DeltaE*eV2J;
      // *******************
      // * MPI 2 PHOTONS   *
      // *******************
      if(twophoton_energy > DeltaE - IPD2)
      {
        sigma1=sigma_tmp*POWR(dE_SI,2.5)/SQRTR(dE_SI);
        sigma_MPI_2=sigma1*sigma1/nu_div_hnu_sq;
        k_MPI_z2_z3[i][j][0]=sigma_MPI_2*I_sq;
      }
      // *******************
      // * MPI 3 PHOTONS   *
      // *******************  
      if(threephoton_energy > DeltaE -IPD2)
      {
        sigma1=sigma_tmp*POWR(dE_SI,2.5)/SQRTR(dE_SI);
        sigma_MPI_3=sigma1*sigma1*sigma1/2.0/nu_div_nu_div_hnu_cub;
        k_MPI_z2_z3[i][j][1]=sigma_MPI_3*I_cu;//*prob*beta_pi(Te,mu,DeltaE);
      }
#endif            
      
    }
  }
  if(fail==1) return -1;

  #endif // MAXLEVEL > 2

   /////////////////
  // Ioniz. 3->4
  ///////////////// 
#if MAXLEVEL > 3  

  fail=0;
#ifdef OMP
 // #pragma omp parallel for schedule(static) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2,sigma1,sigma_MPI_2,sigma_MPI_3,tmp0,tmp1,tmp2) num_threads(num_threads)
  // #pragma omp for nowait schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2,sigma1,sigma_MPI_2,sigma_MPI_3,tmp0,tmp1,tmp2)
  // #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,fparams_inner, fparams_outer, winteg_inner, winteg_outer,winteg_rb_inner,winteg_rb_outer)
#endif
  for(i=0;i<z3_len;++i)
  {
    for(j=0;j<z4_len;++j)
    {

      DeltaE=(STATES_z4[j][2]-STATES_z3[i][2])-IPD3;//+EF;

if(DeltaE <0 )
  continue;      


if(RATEIONIZMAX*ne*fermi_factor*Ith(y,i+z0_len+z1_len+z2_len+3)> MINRATE)
      k_EI_z3_z4[i][j]=MAX(0.0,double_integral_ionization(ne,Te, mu, DeltaE*eV2J));

// if(kmax_estim*ne*ne*Ith(y,j+z0_len+z1_len+z2_len+z3_len+3)>MINRATE)    
if(RATERECOMBMAX*nesq*fermi_factor*Ith(y,j+z0_len+z1_len+z2_len+z3_len+3) > MINRATE)      
      k_EI_z4_z3[i][j]=STATES_z3[i][3]/STATES_z4[j][3]*double_integral_recombination(ne,Te, mu, DeltaE*eV2J);

k_EI_MAX=MAX(k_EI_MAX,k_EI_z3_z4[i][j]);
k_EI_REV_MAX=MAX(k_EI_REV_MAX,k_EI_z4_z3[i][j]);      

      k_EI_z3_z4[i][j]*=fermi_factor;
      k_EI_z4_z3[i][j]*=fermi_factor;

#ifdef MULTIPHOTON
      realtype dE_SI=DeltaE*eV2J;
      // *******************
      // * MPI 2 PHOTONS   *
      // *******************
      if(twophoton_energy> DeltaE - IPD3)
      {
        sigma1=sigma_tmp*POWR(dE_SI,2.5)/SQRTR(dE_SI);
        sigma_MPI_2=sigma1*sigma1/nu_div_hnu_sq;
        k_MPI_z3_z4[i][j][0]=sigma_MPI_2*I_sq;
      }
      // *******************
      // * MPI 3 PHOTONS   *
      // *******************  
      if(threephoton_energy > DeltaE - IPD3)
      {        
        sigma1=sigma_tmp*POWR(dE_SI,2.5)/SQRTR(dE_SI);
        sigma_MPI_3=sigma1*sigma1*sigma1/2.0/nu_div_nu_div_hnu_cub;
        k_MPI_z3_z4[i][j][1]=sigma_MPI_3*I_cu;//*prob*beta_pi(Te,mu,DeltaE);
      }
#endif
      // **************
      // * RAD RECOMB *
      // **************
      // if(DeltaE>0)
      // {
      //   if(expint>0)
      //   {
      //     k_MPI_z4_z3[i][j][0]=v_e*k_RR_fact1*4.0*k_RR_fact2*POWR((DeltaE)/STATES_z4[j][2],1.5)*expint*exp(a);  
      //   }
      //   else
      //   {
      //     k_MPI_z4_z3[i][j][0]=0.0;
      //   }       
        
      // }      
      
    }
  }

  if(fail==1) return -1; 
#endif // MAXLEVEL > 3
// if(k_EI_MAX >0 || k_EI_REV_MAX > 0)
  // printf("myid:%d, k_EI_MAX:%.4e, k_EI_REV_MAX:%.4e,fac:%.4e\n",myid,k_EI_MAX,k_EI_REV_MAX,fermi_factor);
    return 0;

}

// ***************************************************************************************



// ********************************************************************
// * WRITE COLRAD CONCENTRATIONS TO FILE FOR RESTART
// ********************************************************************
int colrad_write(int number)
{
    FILE *outfile;
    char fname[255];
    sprintf(fname, "%s.%d.%d.colrad", outfilename, myid,number);
    outfile = fopen(fname, "w");
    if (NULL == outfile)
    {
      char errstr[255];
      sprintf(errstr,"ERROR: cannot open colrad outfile %s\n",fname);
      error(errstr);
    }

    int i,j,k,l;    
    for(i=1; i < local_fd_dim.x-1; i++)
    {
      for(j=1; j < local_fd_dim.y-1; j++)
      {
        for(k=1; k < local_fd_dim.z-1; k++) 
        {
          fprintf(outfile,"%d %d %d",i,j,k);
          for(l=0;l<neq;l++)
          {
            fprintf(outfile, " %.4e ", Ith(node.y,l));
          }
          fprintf(outfile,"\n");
        }        
      }

    }
    fclose(outfile);
  return 0; //alles ok
}

// ********************************************************************
// * READ COLRAD CONCENTRATIONS FOR RESTART
// ********************************************************************
int colrad_read(int number)
{
  FILE *infile;
  char fname[255];
  int i,j,k,l;

  sprintf(fname, "%s.%d.%d.colrad", outfilename, myid,number);
  infile=fopen(fname,"r");
  if(infile==NULL)
  {
    char errstr[255];
    sprintf(errstr,"ERROR: Colrad infile %s not found\n", fname);
    error(errstr);

  }
  realtype tmp;


    char **tokens;
    size_t numtokens;

    char line[MAX_LINE_LENGTH];
    int linenr=1;
    for(i=1; i < local_fd_dim.x-1; i++)
    {
      for(j=1; j < local_fd_dim.y-1; j++)
      {
        for(k=1; k < local_fd_dim.z-1; k++) 
        {

          //read data
          if (fgets (line, MAX_LINE_LENGTH, infile) == NULL) {
            char errstr[255];
            sprintf(errstr,"Error Reading colrad-input-file: %s in line %d.\n", fname,linenr);
            error(errstr);
          }
          tokens = strsplit(line, ", \t\n", &numtokens); //strsplit in imd_ttm.char
          for(l=0;l<neq;l++)
          {

            sscanf(tokens[l+3], "%lf",  &tmp);                              
            Ith(node.y,l)=tmp;
            
          }
          linenr++;

        for (l = 0; l < numtokens; l++) {
          free(tokens[l]);
        }
        if (tokens != NULL)
          free(tokens);          

        }        
      }

    }  
  fclose(infile);
  return 0;
}



// ***************************************
// *    INTEGRATION STUFF
// *************************************

realtype inner_integrand_ionization(realtype x, void *p) // x=E_strich
{
  struct my_f_params * params = (struct my_f_params *)p;
  realtype E_prime=x;
  realtype ne=params->ne;
  realtype T=params->T;
  realtype mu=params->mu;  
  realtype DeltaE = params->DeltaE;  
  realtype E=params->E; //brauche ich nachher für E''=E-E'-DELTAE
  realtype E_prime_prime=E-E_prime-DeltaE;

  realtype Pauli_E_prime=1.0-1.0/(1.0+EXPR((E_prime-mu)/BOLTZMAN/T));
  realtype Pauli_E_prime_prime=1.0-1.0/(1.0+EXPR((E_prime_prime-mu)/BOLTZMAN/T));

  realtype f=Pauli_E_prime*Pauli_E_prime_prime; //F(E), sigmaderiv und SQRTR(2*eng1/emass) im outer integrand
  return f;
}


realtype inner_integrand_recombination(realtype x, void *p) // x=incoming elec. energy : ACHTUNG: Für x=0 --> NaN
{  
  struct my_f_params * params = (struct my_f_params *)p;
  realtype E_prime=x;
  realtype ne=params->ne;
  realtype T=params->T;
  realtype mu=params->mu;  
  realtype DeltaE = params->DeltaE;    
  realtype E=params->E; 
  
  // Anmerkung: E'' und E' sind die enregien der einfallenden elektronen
  //            während E=DeltaE+E'+E'' energie des sekunddären Elektrons!
  // realtype E=DeltaE+E_prime+E_prime_prime; 
  realtype E_prime_prime=E-DeltaE-E_prime;
  realtype fermi_fun_E_prime=      1.0/(1.0+EXPR((E_prime-mu)/BOLTZMAN/T));
  realtype fermi_fun_E_prime_prime=1.0/(1.0+EXPR((E_prime_prime-mu)/BOLTZMAN/T));
  realtype f= fermi_fun_E_prime * fermi_fun_E_prime_prime;
  return f;
}


realtype outer_integrand_recombination(realtype x,void *p)
{
  struct my_f_params * params = (struct my_f_params *)p;
  realtype E=x; //other incoming electron (inneres Integral behandelt E_prime)  
                          //Das sekundäre Elektron hat Energie E,
  realtype DeltaE = params->DeltaE;  
  realtype ne=params->ne;
  realtype T=params->T;
  realtype mu=params->mu;  

    
  if(E <= DeltaE)  
    return 0.0; //Cross section wird =0

  realtype fermi_fun=1.0/(1.0+EXPR((E-mu)/BOLTZMAN/T));
  realtype Pauli_E = 1.0-fermi_fun;

  fparams_inner.T=T;
  fparams_inner.ne=ne;
  fparams_inner.mu=mu;
  fparams_inner.DeltaE=DeltaE;
  fparams_inner.E=E;

  gsl_function gslfun_inner;
  gslfun_inner.function=&inner_integrand_recombination;
  gslfun_inner.params=&fparams_inner;  

  realtype integ_inner;
  realtype integ_err;

  realtype y=E/DeltaE;
  // realtype sigma_deriv = 4.0*M_PI*bohr_radius_sq*E_ion_H_sq_J / gsl_pow_2(DeltaE) * alpha_i*(y-1.0)/gsl_pow_2(y)*LOGR(5*beta_i*y/4) 
                       // /2.0/(E-DeltaE);   
                   

  //konstanten nach double_integral_recombination verfrachtet                       
  realtype sigma_deriv = (y-1.0)/POWR(y,2.0)*LOGR(beta_i*1.25*y)/(E-DeltaE);                       


  gsl_integration_qag(&gslfun_inner, 0, E-DeltaE, integ_abstol, integ_reltol, integ_meshdim,1,
                      winteg_inner, &integ_inner, &integ_err);  

  // size_t evals;
  // gsl_integration_romberg(&gslfun_inner, 0, E-DeltaE, 1e-4, integ_reltol, &integ_inner,
                         // &evals, winteg_rb_inner);

  return sigma_deriv*Pauli_E*E*integ_inner;

}

realtype double_integral_recombination(realtype ne,realtype T, realtype mu, realtype DeltaE)
{

// return 0;
  gsl_function gslfun_outer;
  gslfun_outer.function = &outer_integrand_recombination;

  fparams_outer.T=T;
  fparams_outer.ne=ne;
  fparams_outer.mu=mu;
  fparams_outer.DeltaE=DeltaE;

  gslfun_outer.params = &fparams_outer;

  realtype integ_outer=0;
  realtype integ_err=0;

   //angepasste obere integr. grenze
   realtype eupper=0.0;
   if(mu> 0) eupper=POWR(3*T,0.33) * eV2J + mu +DeltaE;
   else      eupper=10.0*T/11604* eV2J+ DeltaE;

 // integrand_inner=@(eng1,eng2) sigma_deriv(eng1).*SQRTR(2*eng1/emass) ...
                  // .*F(eng1).*Pauli(eng2).*Pauli(eng1-DeltaE-eng2);


  // gsl_integration_qagiu(&gslfun_outer, 0.0, integ_abstol, integ_reltol, integ_meshdim,
                        // winteg_outer, &integ_outer, &integ_err); 

  gsl_integration_qag(&gslfun_outer, DeltaE, eupper, integ_abstol_recomb, integ_reltol, integ_meshdim,1,
                       winteg_outer, &integ_outer, &integ_err);  


  // integ_outer = integral_simpson(&outer_integrand_recombination, DeltaE, eupper, 1000, &fparams_inner);
  integ_outer *= 2.0*M_PI*bohr_radius_sq*E_ion_H_sq_J * gsl_pow_2(1.0/DeltaE)*alpha_i; //konstanten aus sigma_deriv herausgezogen

  // size_t evals;
  // gsl_integration_romberg(&gslfun_outer, DeltaE, muINF, 1e-6, integ_reltol, &integ_outer,
                         // &evals, winteg_rb_outer);


  //NICHT VERGESSEN: RATIO DER STATISTICAL WEIGHTS
  integ_outer *= recomb_const/ne/ne; //Später ne^2 entfernen und in ydot entfernen!
  // if(integ_outer < 1e-100) integ_outer=0.0;
  return MAX(integ_outer,0.0);

}



realtype outer_integrand_ionization(realtype x,void *p)
{
  struct my_f_params * params = (struct my_f_params *)p;
  realtype eng=x;
  realtype DeltaE = params->DeltaE;  
  realtype ne=params->ne;
  realtype T=params->T;
  realtype mu=params->mu;  


  realtype fermi_fun=1.0/(1.0+EXPR((eng-mu)/BOLTZMAN/T));

  // if(fermi_fun < 1e-100) return 0;

  if(eng <= DeltaE)
    return 0.0; //Cross section wird zu 0

  realtype y=eng/DeltaE;

  // realtype sigma_deriv = 4.0*M_PI*bohr_radius_sq*E_ion_H_sq_J * gsl_pow_2(1.0/DeltaE)*alpha_i*(y-1.0)/gsl_pow_2(y)*LOGR(5*beta_i*y/4) 
                       // /2.0/(eng-DeltaE); 

  //sigma_deriv: konstanten rausgezogen und nach double_integral_ionization verfrachtet
  realtype sigma_deriv = (y-1.0)/POWR(y,2.0)*LOGR(beta_i*1.25*y)/(eng-DeltaE);


  fparams_inner.T=T;
  fparams_inner.ne=ne;
  fparams_inner.mu=mu;
  fparams_inner.DeltaE=DeltaE;
  fparams_inner.E=eng;

  gsl_function gslfun_inner;
  gslfun_inner.function=&inner_integrand_ionization;
  gslfun_inner.params=&fparams_inner;  

  realtype integ_inner;
  realtype integ_err;
 

  // gsl_integration_qag(&gslfun_inner, 0.0, eng-DeltaE, integ_abstol, integ_reltol, integ_meshdim,1,
                      // winteg_inner, &integ_inner, &integ_err);  

  // size_t evals;
  // gsl_integration_romberg(&gslfun_inner, 0.0, eng-DeltaE, 1e-4, integ_reltol, &integ_inner,
                         // &evals, winteg_rb_inner);
  
  // integ_inner = integral_simpson(&inner_integrand_ionization, 1e-40, eng-DeltaE, 5000, &fparams_inner);

  realtype fa=inner_integrand_ionization(1e-20,p);
  realtype fb=inner_integrand_ionization(eng-DeltaE,p);
  realtype fm=inner_integrand_ionization((eng-DeltaE)/2,p);
  realtype h=eng-DeltaE-0.0;
  realtype Sinit=h/6*(fa+4*fm+fb);
  #pragma omp parallel
  {
    integ_inner = integral_simpson_par(inner_integrand_ionization, 1e-20, eng-DeltaE, 1e-4,Sinit,fa,fb,fm,10,p);
  } 

  return eng*fermi_fun*sigma_deriv*integ_inner;

}

realtype double_integral_ionization(realtype ne,realtype T, realtype mu, realtype DeltaE)
{
  // return 0;

  gsl_function gslfun_outer;
  gslfun_outer.function = &outer_integrand_ionization;

  fparams_outer.T=T;
  fparams_outer.ne=ne;
  fparams_outer.mu=mu;
  fparams_outer.DeltaE=DeltaE;

  gslfun_outer.params = &fparams_outer;

  realtype integ_outer=0;
  realtype integ_err=0;


  realtype eupper=0.0;
  if(mu> 0) eupper=POWR(3*T,0.33) * eV2J + mu +DeltaE; //entartet
  else      eupper=10.0*T/11604* eV2J+ DeltaE;         //nicht-entartet


  gsl_integration_qag(&gslfun_outer, DeltaE*1.001, eupper, integ_abstol_ioniz, integ_reltol, integ_meshdim,1,
                      winteg_outer, &integ_outer, &integ_err);  

  // integ_outer = integral_simpson(&outer_integrand_ionization, DeltaE*1.001, eupper, 1000, &fparams_outer);

  integ_outer *= 2.0*M_PI*bohr_radius_sq*E_ion_H_sq_J * gsl_pow_2(1.0/DeltaE)*alpha_i; //konstanten aus sigma_deriv herausgezogen

  //Romberg hängt sich hier auf
  // size_t evals;
  // gsl_integration_romberg(&gslfun_outer, DeltaE*1.001, muINF,0.0, integ_reltol, &integ_outer,
                          // &evals, winteg_rb_outer);

  // if(integ_outer<MINRATE) integ_outer=0.0;
  integ_outer *= ioniz_const / ne; //ACHTUNG: Später ne entfernenu und in ydot nicht mehr multiplizieren!
  return MAX(integ_outer,0.0);

}

realtype fermi_integrand(realtype x, void *p)
{
  struct my_f_params * params = (struct my_f_params *)p;
  realtype eng=x;
  realtype ne=params->ne;
  realtype T=params->T;
  realtype mu=params->mu;  

  realtype vel=SQRTR(2.0*eng/EMASS);
  realtype fermi_fun=1.0/(1.0+EXPR((eng-mu)/BOLTZMAN/T));
  realtype F=double_emass_pow_3_2/2.0/ne/hbar_cub/pi_sq*SQRTR(eng)*fermi_fun;  // DOS * f_fermi

  return F;
}
realtype eval_fermi_integrand(realtype ne,realtype T, realtype mu)
{
  fparams_fermi.T=T;
  fparams_fermi.ne=ne;
  fparams_fermi.mu=mu; 

  gsl_function fun;
  fun.function=&fermi_integrand;
  fun.params=&fparams_fermi;

  realtype integ_err=0;
  realtype integ_result=0;


  gsl_error_handler_t *old_error_handler=gsl_set_error_handler_off ();


  // int code= gsl_integration_qags (&fun, mu, muINF, integ_abstol, integ_reltol, integ_meshdim,
  //                       winteg_fermi, &integ_result, &integ_err);   

  // int code=gsl_integration_qagiu(&fun, mu, integ_abstol, integ_reltol, integ_meshdim,
                        // winteg_fermi, &integ_result, &integ_err); 

  int code= gsl_integration_qag(&fun, 0, 3.0*muINF, integ_abstol, integ_reltol, integ_meshdim,1,
                                winteg_fermi, &integ_result, &integ_err);    


  gsl_set_error_handler(old_error_handler); //reset the error handler 

  return (code==GSL_SUCCESS ? integ_result : -1);  // RHS abbrechen -> neuer versuch

  //return integ_result;
}




realtype eval_excitation_integral(realtype ne,realtype T,realtype mu, realtype DeltaE, int allowed)
{
  gsl_function fun;
  fun.function = &integrand_excitation;

  fparams_exc.T=T;
  fparams_exc.ne=ne;
  fparams_exc.mu=mu;
  fparams_exc.DeltaE=DeltaE;
  fparams_exc.allowed=allowed;

  fun.params = &fparams_exc;

  realtype integ_result=0;
  realtype integ_err=0;





// integ_result = integral_simpson(&integrand_excitation, DeltaE, 20*DeltaE, 5000, &fparams_exc);


  gsl_error_handler_t *old_error_handler=gsl_set_error_handler_off ();

  // int code= gsl_integration_qags (&fun, DeltaE, muINF, integ_abstol, integ_reltol, integ_meshdim,
  //                         winteg_exc, &integ_result, &integ_err);    

  int code= gsl_integration_qag(&fun, DeltaE, muINF, 1e-12, integ_reltol, integ_meshdim,1,
                          winteg_exc, &integ_result, &integ_err);    


 gsl_set_error_handler(old_error_handler); //reset the error handler 

  if (code != GSL_SUCCESS)
  {
    //print integrand
    int i=0;
    realtype dx=(muINF-DeltaE)/250;    
    for(i=0;i<250;i++)
    {
      printf("x:%.4e,T:%f, ne:%.2e, dE:%.2e, integ:%.4e\n",
              dx*i,fparams_exc.T, fparams_exc.ne, fparams_exc.DeltaE, integrand_excitation(dx*i,&fparams_exc));
    }

    error("ERROR in eval_excitation_integral\n");
  }

  //return integ_result;
  // if(myid==1) printf("integ:%.4e\n",integ_result);
  // if(integ_result < MINRATE) integ_result=0.0;

  
  return MAX(integ_result,0.0);

}


realtype eval_dexcitation_integral(realtype ne,realtype T,realtype mu, realtype DeltaE, int allowed)
{  
  gsl_function fun;
  fun.function = &integrand_deexcitation;

  fparams_exc.T=T;
  fparams_exc.ne=ne;
  fparams_exc.mu=mu;
  fparams_exc.DeltaE=DeltaE;
  fparams_exc.allowed=allowed;

  fun.params = &fparams_exc;

  realtype integ_result=0;
  realtype integ_err=0;


  // gsl_integration_qagiu(&fun, DeltaE, integ_abstol, integ_reltol, integ_meshdim,
                       // winteg_exc, &integ_result, &integ_err);    

  //Wird sonst immer zu null
  // gsl_integration_qags (&fun, DeltaE, muINF, integ_abstol, integ_reltol, integ_meshdim,
  //                         winteg_exc, &integ_result, &integ_err);    

  //Variante Aslan
  fparams_exc.mu=mu+DeltaE;
  fun.function = &integrand_excitation;


   // gsl_integration_qags (&fun, DeltaE, muINF, integ_abstol, integ_reltol, integ_meshdim,
                         // winteg_exc, &integ_result, &integ_err);    

   gsl_integration_qag(&fun, DeltaE, muINF, integ_abstol, integ_reltol, integ_meshdim,1,
                          winteg_exc, &integ_result, &integ_err);    


    // size_t neval;
    // gsl_integration_qng(&fun, DeltaE, muINF, integ_abstol, integ_reltol,  //FAST?
    //                             &integ_result, &integ_err, &neval);

  // if(integ_result < MINRATE) integ_result=0.0;
  return MAX(integ_result,0.0);

}

realtype integrand_deexcitation(realtype x,void *p)
{
  struct my_f_params * params = (struct my_f_params *)p;
  realtype eng=x;
  realtype DeltaE = params->DeltaE;  
  realtype ne=params->ne;
  realtype T=params->T;
  realtype mu=params->mu;  
  int allowed=params->allowed;

  if(eng <= DeltaE)
    return 0.0;

  realtype vel=SQRTR(2.0*eng/EMASS);
  // if(vel< 1e-100) return 0;

  realtype fermi_fun=1.0/(1.0+EXPR((eng-DeltaE-mu)/BOLTZMAN/T)); //mod
  // if(fermi_fun < 1e-100) return 0;

  realtype sigma=0.0;
  realtype y=eng/DeltaE;

  realtype Pauli=1.0-1.0/(1.0+EXPR((eng+mu)/BOLTZMAN/T)); //mod
  // if(Pauli<1e-100) return 0;

  if(allowed==1)
    sigma=4.0*M_PI*bohr_radius_sq*E_ion_H_sq_J* gsl_pow_2(1.0/DeltaE)*alpha_i*(y-1.0)/gsl_pow_2(y)*LOGR(5*beta_i*y/4);
  else
    sigma=4.0*M_PI*bohr_radius_sq*alpha_i*(y-1.0)/gsl_pow_2(y);

  realtype F=double_emass_pow_3_2/2.0/ne/hbar_cub/pi_sq*SQRTR(eng)*fermi_fun*SQRTR(eng/(eng-DeltaE));  //ACHTUNG: Letzter Term --> divergent
  return vel*sigma*F*Pauli;
}



realtype  integral_simpson(realtype (*f)(realtype, void*), realtype a, realtype b, int n, void* p)
{
  // struct my_f_params * params = (struct my_f_params *)p;
  realtype nd=(realtype) n;
  realtype h=(b-a)/nd;
  realtype h_half=h*0.5;

  realtype f0=f(a,p);
  realtype fn=f(b,p);
  
  int i;
  //summation 1
  realtype sum1=0.0;
  realtype sum2=0.0;

  sum2+=f((2*b-h)*0.5,p); // weil die loop nur von 1 bis N-1 geht, aber summation 2 muss von 1 bis N gehen.
                                // deswegen wird das letzte Elemn. bereits hier addiert

  #ifdef OMP
   #pragma omp parallel for reduction(+: sum1, sum2)
  #endif
  for(i=1;i<n-1;i++)
  {
    realtype xk=a+i*h;
    realtype f_xk=f(xk,p);
    sum1+=f_xk;
    //realtype xhalf=(x_k-1 + x_k)/2 = (x_k-h + x_k)/2
    //realtype xhalf=(2*xk-h)*0.5;
    realtype xhalf=0.5*xk-h_half;
    realtype f_xhalf=f(xhalf,p);
    sum2+=f_xhalf;

    // printf("a:%.4e,b:%.4e, xk:%.4e, f_xk:%.5e, fhalf:%.4e\n", a,b, xk, f_xk, f_xhalf);
  }
  realtype result= h/6.0*(f0+fn+2.0*sum1+4.0*sum2); 
  return MAX(0.0,result);
  
}

realtype integrand_excitation(realtype x,void *p)
{
  struct my_f_params * params = (struct my_f_params *)p;
  realtype eng=x;
  realtype DeltaE = params->DeltaE;  
  realtype ne=params->ne;
  realtype T=params->T;
  realtype mu=params->mu;  
  int allowed=params->allowed;
  
  if(eng <= DeltaE)
    return 0.0; //siehe flychk manual

  realtype vel=SQRTR(2.0*eng/EMASS);
  realtype fermi_fun=1.0/(1.0+EXPR((eng-mu)/BOLTZMAN/T));

  // if(fermi_fun < 1e-100) return 0;

  realtype sigma=0.0;
  realtype y=eng/DeltaE;
  realtype Pauli=1.0-1.0/(1.0+EXPR((eng-DeltaE+mu)/BOLTZMAN/T));
  // if(Pauli < 1e-100) return 0;
  
  if(allowed==1)
    sigma=4.0*M_PI*bohr_radius_sq*E_ion_H_sq_J* gsl_pow_2(1.0/DeltaE)*alpha_i*(y-1.0)/gsl_pow_2(y)*logf(5*beta_i*y/4);  
  else
    sigma=4.0*M_PI*bohr_radius_sq*alpha_i*(y-1.0)/gsl_pow_2(y);

#ifdef USEFLOAT
  float F=double_emass_pow_3_2/2.0/ne/hbar_cub/pi_sq*sqrtf(eng)*fermi_fun;  
#else
  realtype F=double_emass_pow_3_2/2.0/ne/hbar_cub/pi_sq*SQRTR(eng)*fermi_fun;  
#endif  

  // printf("F:%.4e, Pauli:%.4e, eng:%.4e, sigma:%.4e, y:%.4e, de:%.4e,ne:%.4e,T:%.4e\n", 
         // F, Pauli, eng, sigma, y, DeltaE, ne, T);
  return vel*sigma*F*Pauli;
}



// ****************************************************************************
// *  PARALLEL SIMPSON INTEGERAL STUFF
// ****************************************************************************


realtype  integral_simpson_par(
     realtype (*f)(realtype,void*), /* function to integrate */
     realtype ah,                   /* left interval boundary  */
     realtype bh,                   /* right interval boundary */
     realtype tolh,
     realtype S,
     realtype fah,
     realtype fbh,
     realtype fmid,
     int recmax, 
     void* p)   //param. struct

{
  simpson_error=0;

  realtype a, b, tolerance;       /* vars poped from stack */
  realtype integral_result;       /* value to return */
  realtype h;                     /* interval size */
  realtype mid;                   /* center of interval */
  // realtype one_trapezoid_area;    /* area of one trapezoid */
  // realtype two_trapezoid_area;    /* area of two trapezoids */
  realtype left_area;             /* integral of left half interval */
  realtype right_area;            /*     "       right    "         */

  stack_t stack;
  work_t work;

  int ready, idle, busy;

  /* prepare stack */
  work.a = ah;
  work.b = bh;
  work.tol = tolh;
  work.S=S;
  work.fa=fah;
  work.fb=fbh;
  work.rec=recmax;
  work.fm=fmid;

  create_stack(&stack, sizeof(work_t));
  push_stack(stack, &work);

  integral_result = 0.0;

  busy = 0;

#pragma omp parallel default(none) \
    shared(stack, integral_result,f,busy,simpson_error,p) \
    private(a,b,tolerance, work, h, mid, \
            idle, ready)
  {

    ready = 0;
    idle = 1;

    while(!ready && !simpson_error)
    {

#pragma omp critical (stack)
      {
        if (!empty_stack(stack)){
          /* we have new work */ 
          pop_stack(stack, &work);

          if (idle){
            /* say others i'm busy */
            busy += 1;
            idle = 0;
          }

        }else{
          /* no new work on stack */
          if (!idle){
            busy -= 1;
            idle = 1;
          }

          /* nobody has anything to do; let us leave the loop */
          if (busy == 0)
            ready = 1;

        }
      } /* end critical(stack) */

      if (idle)
        continue;

      b = work.b;
      a = work.a;
      tolerance = work.tol;
      realtype S=work.S; // vorheriges TOTAL integral
      realtype fa=work.fa;
      realtype fb=work.fb;
      realtype fm=work.fm;
      int rec=work.rec;

      h = (b - a)/2;
      mid = (a+b)/2;
      realtype lm=(a+mid)/2;
      realtype rm=(mid+b)/2;    

      realtype flm=f(lm,p);
      realtype frm=f(rm,p);
      
      realtype Sl=h/6*(fa+4*flm+fm);
      realtype Sr=h/6*(fm+4*frm+fb);

      realtype delta=Sl+Sr-S;

      if(rec <= 0 || fabs(delta) <= 15*tolerance)
      {
#pragma omp critical (integral_result) //hier stand davor nur result (war ein fehler)
        integral_result += Sl+Sr+delta/15;
      }
      else // error not acceptable
      {
        // serious numerical trouble: it won't converge
        if ((tolerance/2 == tolerance) || fabs(a-lm) <= tolerance || tolerance < tolmax) 
        { 
            simpson_error = 1; 
#pragma omp critical (integral_result)            
            integral_result = S;      
        }
        else
        {
          //push new subintervals to stack
          work.a = a;
          work.b = mid;
          work.tol = tolerance/2;
          work.S = Sl;
          work.fa=fa;
          work.fb=fm;
          work.fm=flm;
          work.rec=rec-1;


    #pragma omp critical (stack)
          {
            //Linke Seite
            push_stack(stack, &work);      
            //Rechte seite
            work.a = mid;
            work.b = b;
            work.tol = tolerance/2;
            work.S=Sr;
            work.fa=fm;
            work.fb=fb;
            work.fm=frm;
            work.rec=rec-1;

            push_stack(stack, &work);

          }
        }
      }
    } /* while */
  } /* end omp parallel */

  return integral_result;
}

/******************************************
 * create new stack
 ******************************************/
void 
create_stack(
             stack_t* stack,     /* stack to create */
             int element_size)   /* size of a stack element */
{
  int initial_size = INITIAL_STACK_SIZE;

  /* allocate memory for new stack struct */
  (*stack) = (stack_t) malloc(sizeof(struct stack_s));
  if (!(*stack)){
    char errstr[255];
    sprintf(errstr, "error: could not allocate memory for stack.. Abort.\n"); 
    error(errstr);
    // exit(1);
  }    

  /* allocate memory for stack elements */
  (*stack)->elements = (void*) malloc(element_size * initial_size);
  (*stack)->mem_reserve = initial_size; 
  if (!(*stack)->elements){
    char errstr[255];
    sprintf(errstr, "error: could not allocate memory for stack.. Abort.\n");
    error(errstr);
  }

  (*stack)->el_size = element_size;
  (*stack)->el_count = 0;

}

/*****************************************
 * check if the stack is empty 
 *****************************************/
int 
empty_stack
(stack_t stack)
{
  return stack->el_count <= 0;
}


/*****************************************
 * push a element on stack
 *****************************************/
void 
push_stack(
           stack_t stack,    /* target stack */
           void* element)    /* element to push */
{
  int i, new_reserve;
  int log2_count;

  /* check if we need more memory for stack */    
  if (stack->el_count >= stack->mem_reserve){

      /* calculate new size for the stack
         it should be a power of two */
      for (i = stack->el_count, log2_count = 0; 
           i > 0; 
           i>>1, log2_count++);
      new_reserve = 1 << log2_count;
 
      /* reallocate memory for phase thread tables 
         and nullify new values */
      stack->elements = (void *) realloc(stack->elements, 
                      stack->el_size * new_reserve);
      if (!stack->elements){
        char errstr [255];
        sprintf(errstr, "error: can't reallocate stack.. Aborting\n");
        error(errstr);
        // exit(1);
      }

      stack->mem_reserve = new_reserve;
  }
  
  /* now push the element on top of the stack */
  memcpy((char*)stack->elements + stack->el_count*stack->el_size, 
            element, stack->el_size);
  stack->el_count++;

}


/*****************************************
 * pop an element from stack
 *****************************************/
void 
pop_stack(
          stack_t stack,    /* target stack */
          void* element)    /* where poped el. should be stored */
{
  if (stack->el_count <= 0){
    char errstr[255];
    sprintf(errstr, "error: trying to pop from empty stack.\n");
    error(errstr);
    // exit(2);
  }

  stack->el_count--;
  memcpy(element, 
         (char*)stack->elements + stack->el_count*stack->el_size, 
         stack->el_size);
}



