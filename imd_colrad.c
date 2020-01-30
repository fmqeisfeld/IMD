#include "imd.h"


// ***************************************************
// * TODO:
// * mpi2 und mpi3 anpassen weil energies delE und kbTe in eV
// * SPEEDUP pragma omp   auch in ydot
// * floats wo möglich
// ****************************************************


#define USEFLOAT  // hauptsächlich in der funktion genexptint. Profiling zeigte, dass
                  // hier die meiste zeit verbraucht wird -> float verdoppelt performance
#define OMP
#define LAPACK
//#define MULTIPHOTON
// #define SPONT  //<-- spontante emission, Kaum effekt 
//#define STARK  //<-- reabsorption via stark effect
#define DOIPD    //

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
#define EXPR(x) exp(x)      //<--NIEMALS expf draus machen!

#ifdef USEFLOAT
float fak(float t, float x, float j,float s); //aux. function for genexpint
float genexpint(float x,float ss,float j);
#else
double fak(double t, double x, double j,double s); //aux. function for genexpint
double genexpint(double x,double ss,double j);
#endif
double ExpInt(double x);        //gsl hat keine float-variante für expint


const int MAX_LINE_LENGTH=3000; //wird von colrad_read benuthzt
                                //Hängt davon ab wieviele states betrachtet werden
                                //sollte eigentlich besser dynamisch berechnet werden



// *********************************************************
// PHYSICAL CONSTANTS
// *********************************************************
// const double eV2J=1.6021766E-19;
const double eV2H=0.03674932; //eV to Hartree
const double colrad_reltol=1e-6;
const double colrad_abstol=10;

// const double J2eV=6.2415091E18;
const double planck=6.62607004E-34; // J/s
const double bohr_radius=0.52917721067E-10; // m
const int MAXLINE = 255;
const double  pi=3.141592653589793;

const double colrad_tequi=1e-3;//TEST// 1e-12; //bei initial equi ohne Temperatur-variation erst einmal 
                                //die Saha-besetzungsdichten equilibrieren

//const double  LIGHTSPEED=2.997925458e8; // m/s
double  HBAR;
double  LASERFREQ;


int num_threads;
//const double  EMASS=9.10938356e-31;     // kg
//const double  ECONST=8.854187817e-12;   // As/Vm
//const double  BOLTZMAN=1.38064852e-23;  // J/K
//const double  ECHARGE=1.60217662e-19;  // C
//const double  AMU=1.66053904020e-27;   // atomic mass unit




typedef struct {
  realtype It; //Intesity
  realtype IPD0,IPD1,IPD2,IPD3,IPD4;
  double EF;
  realtype P_EE,P_EI,P_MPI2,P_MPI3,P_RR;
  double P_TOTAL; //komplette colrad-Leistungsdichte, für eng-filge
  double dens; //weil cv F(dens,temp) 
  double ni; //für Z=ne/ni
  bool initial_equi;
} *colrad_UserData;
colrad_UserData  cdata;

// *********************************************************
//                    MAIN
// *********************************************************
void do_colrad(double dt)
{
  int flag;
  double t;
  double tout=dt;
  int i,j,k;
  N_Vector y;
  double Te0,Ti0,rho0,ni0,ne0;
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
    for(j=1;j<local_fd_dim.y-1;j++)
    {
      for(k=1;k<local_fd_dim.z-1;k++)
      {
        if(l1[i][j][k].natoms < fd_min_atoms) continue;
        y=l1[i][j][k].y;

        Te0=l1[i][j][k].temp*11604.5;
        Ti0=l1[i][j][k].md_temp*11604.5;        
        rho0=l1[i][j][k].dens;
        ni0=rho0/AMU/26.9185; //1e28; //1e26/m^3 entspricht etwa 1e-4/Angtrom^3   

        if(cdata->initial_equi==true)     
        {
          double Zmean=MeanCharge(Te0, rho0, atomic_charge, atomic_weight,i,j,k);          
          ne0= Zmean* l1[i][j][k].dens / (atomic_weight * AMU);          
          l1[i][j][k].ne=ne0; //Saha init greift darauf zu
          colrad_Saha_init(i, j, k);          
        }
        else //NORMAL
        {
          ne0=l1[i][j][k].ne;
        }

        
        Ith(y,0)=Te0;
        Ith(y,1)=Ti0;
        Ith(y,2)=ne0;

        flag = CVodeReInit(cvode_mem, 0.0, y);                    

        
        if(cdata->initial_equi==true)
        {
//          printf("myid:%d, RUNNNING INITIAL EQUI\n",myid);
          flag = CVode(cvode_mem, colrad_tequi, y, &t, CV_NORMAL);          
          int i_global,j_global,k_global;

          i_global = ((i - 1) + my_coord.x * (local_fd_dim.x - 2));
          j_global = ((j - 1) + my_coord.y * (local_fd_dim.y - 2));
          k_global =  ((k-1) + my_coord.z*(local_fd_dim.z-2));          
  // printf("COLRAD Cell %d,%d,%d was equilibrated\n",i_global,j_global,k_global);
        }
        else //NORMAL
        {
          cdata->dens=l1[i][j][k].dens;
          flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
          colrad_ptotal+=(fd_vol)*1e-30*cdata->P_TOTAL; // d.h. ptotal ist Gesamt-Leisung
          //ni0=cdata->ni;
        }

// printf("myid:%d, after i:%d,j:%d,k:%d,Tefin:%.4e,Zfin:%.4e,Zvgl:%.4e,pcell:%.4e,ptot:%.4e\n",
//     myid,i,j,k, Ith(y,0), Ith(y,2)/ni0,
//     MeanCharge(Ith(y,0), rho0, atomic_charge, atomic_weight,i,j,k),
//     cdata->P_TOTAL, colrad_ptotal);
//printf("myid:%d, i:%d, Delta(Te):%.4e,Zmean:%.4e\n",myid,i, Ith(y,0)-Te0,Ith(y,2)/ni0);      

        //REASSIGN NEW TE AND NE
        l1[i][j][k].temp=Ith(y,0)/11604.5;
        l1[i][j][k].ne=Ith(y,2);
        l1[i][j][k].Z=l1[i][j][k].ne/ni0;
        l1[i][j][k].P_EE=cdata->P_EE;
        l1[i][j][k].P_EI=cdata->P_EI;
        l1[i][j][k].P_MPI2=cdata->P_MPI2;
        l1[i][j][k].P_MPI3=cdata->P_MPI3;
        l1[i][j][k].P_RR=cdata->P_RR;

        l2[i][j][k].temp=l1[i][j][k].temp; //auch in l2 speichern! wichtig!
        l2[i][j][k].ne=l1[i][j][k].ne;     //sonst geht alles im diffusion-step vor 
        l2[i][j][k].Z=l1[i][j][k].Z;       //ttm-writeout verloren!
        l2[i][j][k].P_EE=l1[i][j][k].P_EE;
        l2[i][j][k].P_EI=l1[i][j][k].P_EI;
        l2[i][j][k].P_MPI2=l1[i][j][k].P_MPI2;
        l2[i][j][k].P_MPI3=l1[i][j][k].P_MPI3;
        l2[i][j][k].P_RR=l1[i][j][k].P_RR;

if(myid==0)
{
  int i_global = ((i - 1) + my_coord.x * (local_fd_dim.x - 2));
  if (i_global==16 || i_global==17)
  {
    printf("ig:%d, ne:%.4e, Z:%.4e,natoms:%d, ni0:%.4e,ni:%.4e, rho:%.4e,mpi2:%.4e,mpi3:%.4e\n", 
      i_global, l1[i][j][k].ne, l1[i][j][k].Z, l1[i][j][k].natoms, ni0, cdata->ni, l1[i][j][k].dens,
      l1[i][j][k].P_MPI2,l1[i][j][k].P_MPI3);
  }
}
// printf("myid:%d,COLRAD,i:%d,ne:%.4e,Z:%.4e,temp:%.4e\n", myid,i,l1[i][j][k].ne, l1[i][j][k].Z,l1[i][j][k].temp);

        if(l1[i][j][k].temp <0 || isnan(l1[i][j][k].temp) !=0 )        
        {
          char errstr[255];
          sprintf(errstr,"ERROR in COLRAD: Te became Nan or <0\n");
          error(errstr);
        }
        if(l1[i][j][k].ne <0 || isnan(l1[i][j][k].ne) !=0 )        
        {
          char errstr[255];
          sprintf(errstr,"ERROR in COLRAD: ne became Nan or <0\n");
          error(errstr);
        }

      }
    } 
  }
 if(cdata->initial_equi==true)
 {
  colrad_write(0);
  cdata->initial_equi=false;
 }
 else if(steps % ttm_int ==0)
 {
  colrad_write(steps);
 } 

 #ifdef TIMING
 MPI_Barrier(cpugrid);
 // if(myid==0)
 // {

    gettimeofday(&end, NULL);
    double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
             end.tv_usec - start.tv_usec) / 1.e6;

    if(myid==0)
      printf("t-elapsed:%f\n",delta);
 // }
 #endif
  
}

// *********************************************************
//                      INIT FUNC
// *********************************************************
void colrad_init(void)
{
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
		for(j=1;j<local_fd_dim.y-1;j++)
		{
			for(k=1;k<local_fd_dim.z-1;k++)
			{
				l1[i][j][k].y=N_VNew_Serial(neq);
        l1[i][j][k].P_EE=0.0;
        l1[i][j][k].P_EI=0.0;
        l1[i][j][k].P_MPI2=0.0;
        l1[i][j][k].P_MPI3=0.0;
        l1[i][j][k].P_RR=0.0;

				//l2... <--brauche nur 1 mal speicher alloc'en
				//aber in DIFF-LOOP darauf achten, dass beim swappen
				//l2.y immer auf l1.y zeigt und nicht ins Leere
			}

		}
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
      Ith(abstol,0)=2.0; //Temp

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
	CVodeSetMaxNumSteps(cvode_mem,1500); //Default=500
	//CVodeSetInitStep(cvode_mem,1e-20);  //ACHTUNG:Wenn zu klein --> BS
	//CVodeSetEpsLin(cvode_mem,0.01); //Default=0.05;

  CVodeSetMaxStepsBetweenJac(cvode_mem,200); //default 50 --> guter performance boost
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
        double Te0,Ti0,ne0,ni0,rho0;
        N_Vector y;
        y=l1[i][j][k].y;

        Te0=l1[i][j][k].temp*11604.5;
        Ti0=l1[i][j][k].md_temp*11604.5;        
        rho0=l1[i][j][k].dens; ///1e10;
        ni0=rho0/AMU/26.9185; //1e28; //1e26/m^3 entspricht etwa 1e-4/Angtrom^3        
        ne0=l1[i][j][k].ne;

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
  double *buf; //buffer 1d array for communication
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
    alloc2darr(double, STATES_z0, lcnt, 6);

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
  alloc1darr(double, buf, z0_len * 6);

  MPI_Barrier(cpugrid);
  if (myid == 0)
  {
    for (i = 0; i < z0_len * 6; i += 6) //fill 1d buff-array
    {
      buf[i] =   (double) STATES_z0[i / 6][0];
      buf[i + 1] = (double) STATES_z0[i / 6][1];
      buf[i + 2] = (double) STATES_z0[i / 6][2];
      buf[i + 3] = (double) STATES_z0[i / 6][3];
      buf[i + 4] = (double) STATES_z0[i / 6][4];
      buf[i + 5] = (double) STATES_z0[i / 6][5];

      //printf("E:%.4e,i/6:%d\n",STATES_z0[i/6][2], i/6);
    }
  }

  MPI_Barrier(cpugrid);
  MPI_Bcast(buf, z0_len * 6, MPI_DOUBLE, 0, cpugrid);


  //NOW RECONSTRUCT on other procs
  if (myid > 0)
  {
    alloc2darr(double, STATES_z0, z0_len, 6);
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

    alloc2darr(double, STATES_z1, lcnt, 6);

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
  alloc1darr(double, buf, z1_len * 6);
  if (myid == 0)
  {
    for (i = 0; i < z1_len * 6; i += 6) //fill 1d buff-array
    {
      buf[i] =   (double) STATES_z1[i / 6][0];
      buf[i + 1] = (double) STATES_z1[i / 6][1];
      buf[i + 2] = (double) STATES_z1[i / 6][2];
      buf[i + 3] = (double) STATES_z1[i / 6][3];
      buf[i + 4] = (double) STATES_z1[i / 6][4];
      buf[i + 5] = (double) STATES_z1[i / 6][5];
    }
  }

  MPI_Bcast(buf, z1_len * 6, MPI_DOUBLE, 0, cpugrid);
  //NOW RECONSTRUCT on other procs
  if (myid > 0)
  {
    alloc2darr(double, STATES_z1, z1_len, 6);
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

    alloc2darr(double, STATES_z2, lcnt, 6);

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
  alloc1darr(double, buf, z2_len * 6);
  if (myid == 0)
  {
    for (i = 0; i < z2_len * 6; i += 6) //fill 1d buff-array
    {
      buf[i] =   (double) STATES_z2[i / 6][0];
      buf[i + 1] = (double) STATES_z2[i / 6][1];
      buf[i + 2] = (double) STATES_z2[i / 6][2];
      buf[i + 3] = (double) STATES_z2[i / 6][3];
      buf[i + 4] = (double) STATES_z2[i / 6][4];
      buf[i + 5] = (double) STATES_z2[i / 6][5];
    }
  }

  MPI_Bcast(buf, z2_len * 6, MPI_DOUBLE, 0, cpugrid);
  //NOW RECONSTRUCT on other procs
  if (myid > 0)
  {
    alloc2darr(double, STATES_z2, z2_len, 6);
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

    alloc2darr(double, STATES_z3, lcnt, 6);

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
  alloc1darr(double, buf, z3_len * 6);
  if (myid == 0)
  {
    for (i = 0; i < z3_len * 6; i += 6) //fill 1d buff-array
    {
      buf[i] =   (double) STATES_z3[i / 6][0];
      buf[i + 1] = (double) STATES_z3[i / 6][1];
      buf[i + 2] = (double) STATES_z3[i / 6][2];
      buf[i + 3] = (double) STATES_z3[i / 6][3];
      buf[i + 4] = (double) STATES_z3[i / 6][4];
      buf[i + 5] = (double) STATES_z3[i / 6][5];
    }
  }

  MPI_Bcast(buf, z3_len * 6, MPI_DOUBLE, 0, cpugrid);
  //NOW RECONSTRUCT on other procs
  if (myid > 0)
  {
    alloc2darr(double, STATES_z3, z3_len, 6);
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

    alloc2darr(double, STATES_z4, lcnt, 6);

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
  alloc1darr(double, buf, z4_len * 6);
  if (myid == 0)
  {
    for (i = 0; i < z4_len * 6; i += 6) //fill 1d buff-array
    {
      buf[i] =   (double) STATES_z4[i / 6][0];
      buf[i + 1] = (double) STATES_z4[i / 6][1];
      buf[i + 2] = (double) STATES_z4[i / 6][2];
      buf[i + 3] = (double) STATES_z4[i / 6][3];
      buf[i + 4] = (double) STATES_z4[i / 6][4];
      buf[i + 5] = (double) STATES_z4[i / 6][5];
    }
  }

  MPI_Bcast(buf, z4_len * 6, MPI_DOUBLE, 0, cpugrid);
  //NOW RECONSTRUCT on other procs
  if (myid > 0)
  {
    alloc2darr(double, STATES_z4, z4_len, 6);
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
  alloc2darr(double, k_EE_z0_z0, z0_len, z0_len);
  alloc2darr(double, k_EE_z0_z0_b, z0_len, z0_len);

  alloc2darr(double, k_EE_z1_z1, z1_len, z1_len);
  alloc2darr(double, k_EE_z1_z1_b, z1_len, z1_len);

  #if MAXLEVEL > 1
  alloc2darr(double, k_EE_z2_z2, z2_len, z2_len);
  alloc2darr(double, k_EE_z2_z2_b, z2_len, z2_len);
  #endif

  #if MAXLEVEL > 2
  alloc2darr(double, k_EE_z3_z3, z3_len, z3_len);
  alloc2darr(double, k_EE_z3_z3_b, z3_len, z3_len);
  #endif

  #if MAXLEVEL > 3
  alloc2darr(double, k_EE_z4_z4, z4_len, z4_len);
  alloc2darr(double, k_EE_z4_z4_b, z4_len, z4_len);
  #endif
  // **********************************************
  // Now Thermal Ionization and Recomb. rate coeff. arrays

  // z0->z1
  alloc2darr(double, k_EI_z0_z1, z0_len, z1_len);
  alloc2darr(double, k_EI_z1_z0, z0_len, z1_len);

  #if MAXLEVEL > 1
  //z1->z2
  alloc2darr(double, k_EI_z1_z2, z1_len, z2_len);
  alloc2darr(double, k_EI_z2_z1, z1_len, z2_len);
  #endif

  #if MAXLEVEL > 2
  //z2->z3
  alloc2darr(double, k_EI_z2_z3, z2_len, z3_len);
  alloc2darr(double, k_EI_z3_z2, z2_len, z3_len);
  #endif

  #if MAXLEVEL > 3
  //z3->z4
  alloc2darr(double, k_EI_z3_z4, z3_len, z4_len);
  alloc2darr(double, k_EI_z4_z3, z3_len, z4_len);
  #endif
  // k_MPI arrays
  //z0->z1
  alloc3darr(double, k_MPI_z0_z1, z0_len, z1_len, 2);
  alloc3darr(double, k_MPI_z1_z0, z0_len, z1_len, 2);

  //z1->z2
  #if MAXLEVEL > 1
  alloc3darr(double, k_MPI_z1_z2, z1_len, z2_len, 2);
  alloc3darr(double, k_MPI_z2_z1, z1_len, z2_len, 2);
  #endif

  //z2->z3
  #if MAXLEVEL > 2
  alloc3darr(double, k_MPI_z2_z3, z2_len, z3_len, 2);
  alloc3darr(double, k_MPI_z3_z2, z2_len, z3_len, 2);
  #endif


  //z3->z4
  #if MAXLEVEL > 3
  alloc3darr(double, k_MPI_z3_z4, z3_len, z4_len, 2);
  alloc3darr(double, k_MPI_z4_z3, z3_len, z4_len, 2);
  #endif

  MPI_Bcast(&total_species, 1, MPI_INT, 0, cpugrid);
}





void do_Saha(double Te,double totalc,double ne,N_Vector y) //Bei init
{
  double tmp;
  double Zav=ne/totalc;
  // double Ti=Te;

  double Q_z0,Q_z1,Q_z2,Q_z3,Q_z4,Q_z5; //partition functions
  // double IPD_SP=0.0; //Stewart and Pyatt =(pow(1.0+a/Debye,2.0/3.0)-1.0)/2.0/(Zav+1.0);
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
  //Debye=sqrt(BOLTZMAN*Te/4.0/pi/ECHARGE/ECHARGE/ne/(Zav+1));
  tmp=pow(2.0*pi*EMASS*BOLTZMAN*Te,1.5)/planck/planck/planck; //ACHTUNG: Das ist wieder thermal-de brogilie Lambda
                                                              //Nutze richtiges chempot!!!!
  #ifdef DOIPD
  double IPD0,IPD1,IPD2,IPD3,IPD4;
  double z; //DOI after ionization (e.g.=1 for Al0)
  double r0; //Ion sphere radius
  r0=pow(3.0/4.0/pi/totalc,1.0/3.0);
  double debye=sqrt(BOLTZMAN*Te/4.0/pi/pow(totalc+ne,2.0));
  // Atoms, solids, and plasmas in super-intense laser fields S.220
  IPD0=1.0*3.0/2.0/r0*ECHARGE*ECHARGE*(pow(1.0+pow(debye/r0,3.0),2.0/3.0)-pow(debye/r0,2.0))/4.0/pi/ECONST;
  IPD1=2.0*3.0/2.0/r0*ECHARGE*ECHARGE*(pow(1.0+pow(debye/r0,3.0),2.0/3.0)-pow(debye/r0,2.0))/4.0/pi/ECONST;
  IPD2=3.0*3.0/2.0/r0*ECHARGE*ECHARGE*(pow(1.0+pow(debye/r0,3.0),2.0/3.0)-pow(debye/r0,2.0))/4.0/pi/ECONST;
  IPD3=4.0*3.0/2.0/r0*ECHARGE*ECHARGE*(pow(1.0+pow(debye/r0,3.0),2.0/3.0)-pow(debye/r0,2.0))/4.0/pi/ECONST;  
  IPD4=5.0*3.0/2.0/r0*ECHARGE*ECHARGE*(pow(1.0+pow(debye/r0,3.0),2.0/3.0)-pow(debye/r0,2.0))/4.0/pi/ECONST;  
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
    double Ei=STATES_z0[i][2]*eV2J-IPD+EF;
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
    double Ei=STATES_z1[i][2]*eV2J-IPD+EF;
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
    double Ei=STATES_z2[i][2]*eV2J-IPD+EF;
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
    double Ei=STATES_z3[i][2]*eV2J-IPD+EF;
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
    double Ei=STATES_z4[i][2]*eV2J-IPD+EF;
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

  DeltaE=(STATES_z1[0][2]-0.0)*eV2J-IPD+EF;
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
    DeltaE=(STATES_z2[0][2]-STATES_z1[0][2])*eV2J-IPD+EF;
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
    DeltaE=(STATES_z3[0][2]-STATES_z2[0][2])*eV2J-IPD+EF;
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
  DeltaE=(STATES_z4[0][2]-STATES_z3[0][2])*eV2J-IPD+EF;
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
    Ei=STATES_z0[i][2]*eV2J-IPD+EF;    
    if(Ei<0) prob=0.0;    
    #endif
    if(Q_z0>0)
    {
      Ith(y,i+ishift)=n0/Q_z0*STATES_z0[i][3]*prob;
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
    Ei=(STATES_z1[i][2])*eV2J-IPD0+EF;   
    if(Ei<0) prob=0.0;    
    #endif    
    if(Q_z1>0)
      Ith(y,i+ishift+z0_len)=n1/Q_z1*STATES_z1[i][3]*prob;

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
    Ei=STATES_z2[i][2]*eV2J-IPD1+EF;    
    if(Ei<0) prob=0.0;    
    #endif 
    if(Q_z2>0)       
      Ith(y,i+ishift+z0_len+z1_len)=n2/Q_z2*STATES_z2[i][3]*prob;

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
    Ei=STATES_z3[i][2]*eV2J-IPD2+EF;    
    if(Ei<0) prob=0.0;    
    #endif            
    if(Q_z3>0)
      Ith(y,i+ishift+z0_len+z1_len+z2_len)=n3/Q_z3*STATES_z3[i][3]*prob;

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
    Ei=STATES_z4[i][2]*eV2J-IPD3+EF;    
    if(Ei<0) prob=0.0;    
    #endif                
    if(Q_z4>0)
      Ith(y,i+ishift+z0_len+z1_len+z2_len+z3_len)=n4/Q_z4*STATES_z4[i][3]*prob;

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

}





// ************************************************************************************************************
//                                      ACTION
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////
int colrad_ydot(double t, N_Vector y, N_Vector colrad_ydot, void *user_data)
{

/*
  double t0=1e-12;
  double I0=1e17;
  double tFWHM=100e-15;
  double sigmat=tFWHM/2.3548;
  double sigmatsq=sigmat*sigmat;
  */

  //It=I0*exp(-pow((t-t0),2.0)/2.0/sigmatsq);

  colrad_UserData data;
  data = (colrad_UserData) user_data;

  double It=data->It;
  It=0; //It ist LOKALE Intesität!

  bool initial_equi=data->initial_equi; //falls ja, wird temperatur nicht variiert.

  double Eexc;

  P_E_EE=0;
  P_E_EI=0;
  P_E_MPI2=0.0;
  P_E_MPI3=0.0;
  P_E_RAD_RECOMB=.0;

  double DeltaE;
  double kfwd,krev;
  double kfwd2; // für 3-Photon absorption


  int i,j;
  int ishift=3;  // UND NICHT =4 (colrad_ydot fängt bei 0 das zählen an!)
  int shift,shift1,shift2;
  double ne=Ith(y,2);
  double Te=Ith(y,0);
  double Ti=Ith(y,1);


  //FOR REABSORPTION
  double Ajk;
  double sigma_PI,abs_coeff,thickness,tau,wstark_freq,lambda0,wstark_len,escape_factor,groundstate_ioniz;
  thickness=1e-3;
  double sigma_tmp=64.0*pow(pi,4.0)*pow(ECHARGE,10.0)*EMASS/3.0/sqrt(3.0)/pow(4.0*pi*ECONST,5.0)/pow(planck,6.0)/LIGHTSPEED/pow(LASERFREQ,3.0)/pow(13.6*eV2J,2.0);


  //Pre-Zero <--- MUI IMPORTANTE!!!!
  //and count totalc for ipd and stark effect    
  totalc=0.0;
  // #ifdef OMP
  // #pragma omp parallel for reduction(+: totalc)
  // #endif
  for(i=0;i<neq;i++)
  {
    Ith(colrad_ydot,i)=0.0;
    if(i>=3) totalc+=Ith(y,i);
  }
  data->ni=totalc;  
  if(totalc<0)
    return 1; // <-- RHS fail

  // IPD KRAM
  double IPD0,IPD1,IPD2,IPD3,IPD4;
  IPD0=IPD1=IPD2=IPD3=0.0;
  double EF=fermi_E(ne);
  data->EF=EF;
#ifdef DOIPD
  double r0,debye;
  r0=pow(3.0/4.0/pi/totalc,1.0/3.0);
  debye=sqrt(BOLTZMAN*Te/4.0/pi/pow(totalc+ne,2.0));
  // Atoms, solids, and plasmas in super-intense laser fields S.220
  IPD0=1.0*3.0/2.0/r0*ECHARGE*ECHARGE*(pow(1.0+pow(debye/r0,3.0),2.0/3.0)-pow(debye/r0,2.0))/4.0/pi/ECONST;
  IPD1=2.0*3.0/2.0/r0*ECHARGE*ECHARGE*(pow(1.0+pow(debye/r0,3.0),2.0/3.0)-pow(debye/r0,2.0))/4.0/pi/ECONST;
  IPD2=3.0*3.0/2.0/r0*ECHARGE*ECHARGE*(pow(1.0+pow(debye/r0,3.0),2.0/3.0)-pow(debye/r0,2.0))/4.0/pi/ECONST;
  IPD3=4.0*3.0/2.0/r0*ECHARGE*ECHARGE*(pow(1.0+pow(debye/r0,3.0),2.0/3.0)-pow(debye/r0,2.0))/4.0/pi/ECONST;
  IPD4=5.0*3.0/2.0/r0*ECHARGE*ECHARGE*(pow(1.0+pow(debye/r0,3.0),2.0/3.0)-pow(debye/r0,2.0))/4.0/pi/ECONST;

  data->IPD0=IPD0;
  data->IPD1=IPD1;
  data->IPD2=IPD2;
  data->IPD3=IPD3;
  data->IPD4=IPD4;
if(isnan(IPD0)!=0)
printf("IPD0 is NaN! r0:%.4e, totalc:%.4e, ne:%.4e\n",r0,totalc,ne);

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
      double engi=STATES_z0[i][2]*eV2J-IPD0+EF;
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
      double engi=STATES_z1[i][2]*eV2J-IPD0+EF;            
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
      double engi=STATES_z2[i][2]*eV2J-IPD1+EF;
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
  for(i=0;i<z2_len;++i)
  {
    for(j=0;j<z2_len;++j)
    {
      if(j<=i) continue; // MUI IMPORTANTE
      double engi=STATES_z3[i][2]*eV2J-IPD2+EF;
      if(engi < 0) continue; //depressed state is continuum

      DeltaE=(STATES_z3[j][2]-STATES_z3[i][2])*eV2J;
      kfwd=k_EE_z3_z3[i][j]*Ith(y,i+ishift+shift2)*ne;
      krev=k_EE_z3_z3_b[i][j]*Ith(y,j+ishift+shift2)*ne;

      //exec. reduces conc. of i state and increases conc. of j state
      Ith(colrad_ydot,i+ishift+shift2) -=kfwd;
      Ith(colrad_ydot,j+ishift+shift2) +=kfwd;

      //de-excec. increases conc. of i state & decr. conc. of j state
      Ith(colrad_ydot,i+ishift+shift2)+=krev;
      Ith(colrad_ydot,j+ishift+shift2)-=krev;


      Eexc= (-kfwd+krev)*DeltaE;
      P_E_EE+=Eexc;

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
        sigma_PI=sigma_tmp*pow(DeltaE,2.5)/sqrt(groundstate_ioniz);
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
        sigma_PI=sigma_tmp*pow(DeltaE,2.5)/sqrt(groundstate_ioniz);
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
        sigma_PI=sigma_tmp*pow(DeltaE,2.5)/sqrt(groundstate_ioniz);
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
        sigma_PI=sigma_tmp*pow(DeltaE,2.5)/sqrt(groundstate_ioniz);
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

  double P_E_TOTAL=P_E_EI+P_E_EE+P_E_MPI2+P_E_MPI3+P_E_RAD_RECOMB;
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
    //double cvinv=  1.0/EOS_cve_from_r_te(data->dens, Te);
    double cvinv=  1.0/Cv(Te/11604.5, ne);
    Ith(colrad_ydot,0) =  cvinv*P_E_TOTAL;
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
int colrad_GetCoeffs(N_Vector y,double It,void *user_data)
{
  int i,j,k;
  double kronecker;
  double a;
  double DeltaE;
  double Te,ne;
  double expint;
  double G2;
  double I_1,I_2;

  Te=Ith(y,0);
  if(Te <0 || isnan(Te)!=0) return -1;
  ne=Ith(y,2);
  if(ne <0 || isnan(ne)!=0) return -1;

  double v_e=sqrt(8.0*BOLTZMAN*Te/pi/EMASS);
  double E_ion_H=13.6*eV2J;
  double alpha_i=0.05;
  double alpha_e=0.05;
  double beta_i=4.0;
  double four_pi_a0_sq=4.0*pi*pow(bohr_radius,2.0);
  

  double E_ion_H_div_kTe_sq=pow((E_ion_H/BOLTZMAN/Te),2.0);
  double two_pi_me_kT_hsq=2.0*pi*EMASS*BOLTZMAN*Te/pow(planck,2.0);
  double log54_beta_i=log(5.0/4.0*beta_i);
  double kbTe=(BOLTZMAN*Te);
  kbTe=Te/11604.5;

  double tmp0,tmp1,tmp2;

  colrad_UserData data;
  data = (colrad_UserData) user_data;
  double IPD0=data->IPD0*J2eV;
  double IPD1=data->IPD1*J2eV;
  double IPD2=data->IPD2*J2eV;
  double IPD3=data->IPD3*J2eV;
  double IPD4=data->IPD4*J2eV;
  double EF=data->EF*J2eV;

  //MPI
  double k_RR_fact1=32*pi*pow(bohr_radius,2.0)/3.0/175700.00067;
  double k_RR_fact2=pow((E_ion_H/BOLTZMAN/Te),2.0);
  double sigma_MPI_2;//=sigma1/LASERFREQ/pow(planck*LASERFREQ,2.0); //MPI-cross-sect. (2-photon)
  double sigma_MPI_3;
  double sigma1;

  //const. zur berechnung von sigma1
  double sigma_tmp=64.0*pow(pi,4.0)*pow(ECHARGE,10.0)*EMASS/3.0/sqrt(3.0)/pow(4.0*pi*ECONST,5.0)/pow(planck,6.0)/LIGHTSPEED/pow(LASERFREQ,3.0)/pow(13.6*eV2J,2.0);
  double I_sq=It*It; //for 2-photon-ioniz.
  double I_cu=I_sq*It; // for 3-photon-ioniz.
  int fail=0;

  double pow_two_pi_me_kT_hsq_tmp1= pow(two_pi_me_kT_hsq,1.5); //ACHTUNG: Das ist thermal De-Broglie Lambda
                                                               //Ich muss das für die rückwärts-raten irgendie 
                                                               //durch chempot ersetzen sonst pass das net
                                                               //Ebenso in Saha
#ifdef MULTIPHOTON
  double twophoton_energy=2.0*planck*LASERFREQ*J2eV;
  double threephoton_energy=3.0*planck*LASERFREQ*J2eV;
  double nu_div_hnu_sq=LASERFREQ/pow(planck*LASERFREQ,2.0);
  double nu_div_nu_div_hnu_cub=LASERFREQ/LASERFREQ/pow(planck*LASERFREQ,3.0);
#endif
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
  ///////////////////////////////
  // Elec. excitation for Z=0
  ///////////////////////////////
  fail=0;
#ifdef OMP
//#pragma omp parallel for schedule(dynamic,1) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2) num_threads(num_threads)
  #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2)
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
          double Ei=(STATES_z0[i][2])-IPD0+EF;
          if(Ei<0) continue;
#endif          

          a=DeltaE/kbTe;
          expint=ExpInt(a);

if(expint==-1)
#ifndef OMP
        return -1;
#else
        fail=1;
#endif

          G2=genexpint(a,1.0,1.0);
          I_1=exp(-a)/a - expint;
          I_2=I_1*log54_beta_i+expint/a-G2;
          k_EE_z0_z0[i][j]=v_e*four_pi_a0_sq*(kronecker*alpha_e*a*a*I_1 +  (1.0-kronecker)*alpha_i*E_ion_H_div_kTe_sq*I_2);

      }
    }

if(fail==1)
  return -1;

  // *********************     
  // * NOW REVERSE RATE  *
  // *********************
#ifdef OMP
//#pragma omp parallel for   schedule(dynamic,1) collapse(2) private(kronecker,tmp0,tmp1,tmp2) num_threads(num_threads)
#pragma omp for schedule(static) private(j,kronecker,tmp0,tmp2)
#endif
  for(i=0;i<z0_len;++i)
  {
    for(j=0;j<z0_len;++j)
    {
        if(j<=i ) continue;

        kronecker=1.0; //Achtung:Im Reverse process ist dies ein anderes kronecker delta
        //tmp0=pow(2,(1-kronecker))*STATES_z0[j][3]/STATES_z0[i][3];
        tmp0=STATES_z0[j][3]/STATES_z0[i][3];

        //tmp1=3.0*(1.0-kronecker)/2.0;
        tmp1=0.0; //3.0*(1.0-kronecker)/2.0;
        tmp2=EXPR(-(STATES_z0[j][2]-STATES_z0[i][2])/kbTe);

        //k_EE_z0_z0_b[i][j]=k_EE_z0_z0[i][j]/tmp0/(pow(two_pi_me_kT_hsq,tmp1))/tmp2;
        if(k_EE_z0_z0[i][j]> 0)
        {
          k_EE_z0_z0_b[i][j]=k_EE_z0_z0[i][j]/tmp0/tmp2;
        }
        else
        {
          k_EE_z0_z0_b[i][j]=0.0;
        }
    }
  }
  ////////////////////////
  // Elec. exc. for Z=1
  ///////////////////////
  fail=0;
#ifdef OMP
// #pragma omp parallel for   schedule(dynamic,1) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2) num_threads(num_threads)
  #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2)
  // #pragma omp for simd schedule(static) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2)
#endif
  for(i=0;i<z1_len;++i)
  {
    for(j=0;j<z1_len;++j)
    {
        if(j<=i) continue;
        kronecker=0.0;
        {
          if(STATES_z1[i][4]==STATES_z1[j][4])
            kronecker=1.0;
#ifdef DOIPD
          //double Ei=(STATES_z1[i][2]-STATES_z1[0][2])-IPD1+EF;
          double Ei=STATES_z1[i][2]-IPD0+EF;
          if(Ei<0) continue;
#endif


          DeltaE=(STATES_z1[j][2]-STATES_z1[i][2]);
          a=DeltaE/kbTe;
          expint=ExpInt(a);

if(expint==-1)
#ifndef OMP
        return -1;
#else
        fail=1;
#endif
          G2=genexpint(a,1.0,1.0);
          I_1=EXPR(-a)/a - expint;
          I_2=I_1*log54_beta_i+expint/a-G2;
          k_EE_z1_z1[i][j]=v_e*four_pi_a0_sq*(kronecker*alpha_e*a*a*I_1 +(1-kronecker)*alpha_i*E_ion_H_div_kTe_sq*I_2);
        }
      }
    }
  if(fail==1) return -1;

  // *********************     
  // * NOW REVERSE RATE  *
  // *********************
#ifdef OMP
//#pragma omp parallel for   schedule(dynamic,1) collapse(2) private(kronecker,tmp0,tmp1,tmp2) num_threads(num_threads)
#pragma omp for schedule(static) private(j,kronecker,tmp0,tmp2)
#endif
  for(i=0;i<z1_len;++i)
  {
    for(j=0;j<z1_len;++j)
    {
        if(i==j || j<i ) continue;
        tmp0=STATES_z1[j][3]/STATES_z1[i][3];

        //tmp1=3.0*(1.0-kronecker)/2.0;
        tmp1=0.0; //3.0*(1.0-kronecker)/2.0;

        tmp2=EXPR(-(STATES_z1[j][2]-STATES_z1[i][2])/kbTe);

        //k_EE_z0_z0_b[i][j]=k_EE_z0_z0[i][j]/tmp0/(pow(two_pi_me_kT_hsq,tmp1))/tmp2;
        if(k_EE_z1_z1[i][j]>0)
        {
          k_EE_z1_z1_b[i][j]=k_EE_z1_z1[i][j]/tmp0/tmp2;
        }
        else
        {
          k_EE_z1_z1_b[i][j]=0.0;
        }

    }
  }


  ////////////////////////
  // Elec. exc. for Z=2
  ///////////////////////
#if MAXLEVEL > 1

  fail=0;
#ifdef OMP
// #pragma omp parallel for   schedule(dynamic,1) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2) num_threads(num_threads)
  #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2)
  // #pragma omp for simd schedule(static) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2)
#endif
  for(i=0;i<z2_len;++i)
  {
    for(j=0;j<z2_len;++j)
    {
        if(j<=i) continue;
        kronecker=0.0;
        {
          if(STATES_z2[i][4]==STATES_z2[j][4])
            kronecker=1.0;

#ifdef DOIPD
          // double Ei=(STATES_z2[i][2]-STATES_z2[0][2])-IPD2+EF;
          double Ei=STATES_z2[i][2]-IPD1+EF;
          if(Ei<0) continue;
#endif

          DeltaE=(STATES_z2[j][2]-STATES_z2[i][2]);
          a=DeltaE/kbTe;
          expint=ExpInt(a);

if(expint==-1)
#ifndef OMP
        return -1;
#else
        fail=1;
#endif
          G2=genexpint(a,1.0,1.0);
          I_1=EXPR(-a)/a - expint;
          I_2=I_1*log54_beta_i+expint/a-G2;
          k_EE_z2_z2[i][j]=v_e*four_pi_a0_sq*(kronecker*alpha_e*a*a*I_1 +(1-kronecker)*alpha_i*E_ion_H_div_kTe_sq*I_2);
        }
      }
    }
  if(fail==1) return -1;

  // *********************     
  // * NOW REVERSE RATE  *
  // *********************
#ifdef OMP
//#pragma omp parallel for   schedule(dynamic,1) collapse(2) private(kronecker,tmp0,tmp1,tmp2) num_threads(num_threads)
#pragma omp for schedule(static) private(j,kronecker,tmp0,tmp2)  
#endif
  for(i=0;i<z2_len;++i)
  {
    for(j=0;j<z2_len;++j)
    {
        if(j<=i ) continue;
        tmp0=STATES_z2[j][3]/STATES_z2[i][3];

        //tmp1=3.0*(1.0-kronecker)/2.0;
        tmp1=0.0; //3.0*(1.0-kronecker)/2.0;

        tmp2=EXPR(-(STATES_z2[j][2]-STATES_z2[i][2])/kbTe);
        //ACHTUNG: tmp2 kann -->0 gehen
        // --> krev = NaN
        // --> prüfe ob fwd-rate nicht sowieso = 0 ist

        //k_EE_z0_z0_b[i][j]=k_EE_z0_z0[i][j]/tmp0/(pow(two_pi_me_kT_hsq,tmp1))/tmp2;
        if(k_EE_z2_z2[i][j] > 0)
        {
          k_EE_z2_z2_b[i][j]=k_EE_z2_z2[i][j]/tmp0/tmp2;
        }
        else
        {
          k_EE_z2_z2_b[i][j]=0.0;
        }
    }
  }
#endif //MAXLEVEL > 1

  ////////////////////////
  // Elec. exc. for Z=3
  ///////////////////////
#if MAXLEVEL > 2

  fail=0;
  #ifdef OMP
  // #pragma omp parallel for   schedule(dynamic,1) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2) num_threads(num_threads)
  #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2)
  // #pragma omp for simd schedule(static) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2)
  #endif
  for(i=0;i<z3_len;++i)
  {
    for(j=0;j<z3_len;++j)
    {
        if(j<=i) continue;
        kronecker=0.0;
        {
          if(STATES_z3[i][4]==STATES_z3[j][4])
            kronecker=1.0;

#ifdef DOIPD
          // double Ei=(STATES_z3[i][2]-STATES_z3[0][2])-IPD3+EF;
          double Ei=STATES_z3[i][2]-IPD2+EF;
          if(Ei<0) continue;
#endif

          DeltaE=(STATES_z3[j][2]-STATES_z3[i][2]);
          a=DeltaE/kbTe;
          expint=ExpInt(a);

if(expint==-1)
#ifndef OMP
        return -1;
#else
        fail=1;
#endif

          G2=genexpint(a,1.0,1.0);
          I_1=EXPR(-a)/a - expint;
          I_2=I_1*log54_beta_i+expint/a-G2;
          k_EE_z3_z3[i][j]=v_e*four_pi_a0_sq*(kronecker*alpha_e*a*a*I_1 +(1-kronecker)*alpha_i*E_ion_H_div_kTe_sq*I_2);
        }
      }
    }
  if(fail==1)
          return -1;
  // *********************     
  // * NOW REVERSE RATE  *
  // *********************
  #ifdef OMP
  //#pragma omp parallel for   schedule(dynamic,1) collapse(2) private(kronecker,tmp0,tmp1,tmp2) num_threads(num_threads)
  #pragma omp for schedule(static) private(j,kronecker,tmp0,tmp2)        
  #endif
  for(i=0;i<z3_len;++i)
  {
    for(j=0;j<z3_len;++j)
    {
        if(j<=i ) continue;
        tmp0=STATES_z3[j][3]/STATES_z3[i][3];

        //tmp1=3.0*(1.0-kronecker)/2.0;
        tmp1=0.0; //3.0*(1.0-kronecker)/2.0;

        tmp2=EXPR(-(STATES_z3[j][2]-STATES_z3[i][2])/kbTe);
        //k_EE_z0_z0_b[i][j]=k_EE_z0_z0[i][j]/tmp0/(pow(two_pi_me_kT_hsq,tmp1))/tmp2;
        if(k_EE_z3_z3[i][j] > 0)
        {
          k_EE_z3_z3_b[i][j]=k_EE_z3_z3[i][j]/tmp0/tmp2;
        }
        else
        {
          k_EE_z3_z3_b[i][j]=0.0;
        }
        

    }
  }
#endif // MAXLEVEL > 2
  ////////////////////////
  // Elec. exc. for Z=4
  ///////////////////////
#if MAXLEVEL > 3

  fail=0;
#ifdef OMP
// #pragma omp parallel for   schedule(dynamic,1) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2) num_threads(num_threads)
  #pragma omp for simd schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2)
  // #pragma omp for simd schedule(static) collapse(2) private(kronecker,DeltaE,a,expint,G2,I_1,I_2)
#endif
  for(i=0;i<z4_len;++i)
  {
    for(j=0;j<z4_len;++j)
    {
        if(j<=i) continue;
        kronecker=0.0;
        {
          if(STATES_z4[i][4]==STATES_z4[j][4])
            kronecker=1.0;
#ifdef DOIPD
          // double Ei=(STATES_z4[i][2]-STATES_z4[0][2])-IPD4+EF;
          double Ei=STATES_z4[i][2]-IPD3+EF;
          if(Ei<0) continue;
#endif          
          DeltaE=(STATES_z4[j][2]-STATES_z4[i][2]);
          a=DeltaE/kbTe;
          expint=ExpInt(a);

if(expint==-1)
#ifndef OMP
        return -1;
#else
        fail=1;
#endif

          G2=genexpint(a,1.0,1.0);
          I_1=EXPR(-a)/a - expint;
          I_2=I_1*log54_beta_i+expint/a-G2;
          k_EE_z4_z4[i][j]=v_e*four_pi_a0_sq*(kronecker*alpha_e*a*a*I_1 +(1-kronecker)*alpha_i*E_ion_H_div_kTe_sq*I_2);
        }
      }
    }
  if(fail==1)
          return -1;

  // *********************     
  // * NOW REVERSE RATE  *
  // *********************
#ifdef OMP
//#pragma omp parallel for   schedule(dynamic,1) collapse(2) private(kronecker,tmp0,tmp1,tmp2) num_threads(num_threads)
#pragma omp for schedule(static) private(j,kronecker,tmp0,tmp2)        
#endif        
  for(i=0;i<z4_len;++i)
  {
    for(j=0;j<z4_len;++j)
    {
        if(j<=i ) continue;
        tmp0=STATES_z4[j][3]/STATES_z4[i][3];

        //tmp1=3.0*(1.0-kronecker)/2.0;
        tmp1=0.0; //3.0*(1.0-kronecker)/2.0;

        tmp2=EXPR(-(STATES_z4[j][2]-STATES_z4[i][2])/kbTe);
        //k_EE_z0_z0_b[i][j]=k_EE_z0_z0[i][j]/tmp0/(pow(two_pi_me_kT_hsq,tmp1))/tmp2;
        if(k_EE_z4_z4[i][j] > 0)
        {
          k_EE_z4_z4_b[i][j]=k_EE_z4_z4[i][j]/tmp0/tmp2;
        }
        else
        {
          k_EE_z4_z4_b[i][j]=0.0;
        }        
    }
  }
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
   #pragma omp for nowait schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2,sigma1,sigma_MPI_2,sigma_MPI_3,tmp0,tmp1,tmp2)
#endif
  for(i=0;i<z0_len;++i)
  {
    for(j=0;j<z1_len;++j)
    {
      kronecker=0.0;
      if(STATES_z0[i][4]==STATES_z1[j][4])
        kronecker=1.0;

      DeltaE=(STATES_z1[j][2]-STATES_z0[i][2])-IPD0+EF;
if(DeltaE <0 )
  continue;

#ifdef DOIPD
      DeltaE=MAX(0.0,DeltaE);
#endif      
      a=DeltaE/kbTe;

//printf("myid:%d, dE:%.4e, a:%.4e\n",myid,DeltaE,a);

      expint=ExpInt(a);


if(expint==-1)
#ifndef OMP
        return -1;
#else
        fail=1;
#endif
      G2=genexpint(a,1.0,1.0);
      I_1=EXPR(-a)/a - expint;
      I_2=I_1*log54_beta_i+expint/a-G2;
      k_EI_z0_z1[i][j]=v_e*four_pi_a0_sq*alpha_i*E_ion_H_div_kTe_sq*I_2;

#ifdef MULTIPHOTON
      // *******************
      // * MPI 2 PHOTONS   *
      // *******************
      double dE_SI=DeltaE*eV2J; 
      if(twophoton_energy >= DeltaE-IPD0 && DeltaE-IPD0 > 0.0 )
      {
        sigma1=sigma_tmp*pow(dE_SI,2.5)/sqrt(dE_SI);
        sigma_MPI_2=sigma1*sigma1/nu_div_hnu_sq;
        k_MPI_z0_z1[i][j][0]=sigma_MPI_2*I_sq;
      }
      // *******************
      // * MPI 3 PHOTONS   *
      // *******************   
      if(threephoton_energy >= DeltaE-IPD0 && DeltaE-IPD0 >0.0 ) 
      {
        sigma1=sigma_tmp*pow(dE_SI,2.5)/sqrt(dE_SI);
        sigma_MPI_3=sigma1*sigma1*sigma1/2.0/nu_div_nu_div_hnu_cub;
        k_MPI_z0_z1[i][j][1]=sigma_MPI_3*I_cu;//*prob*beta_pi(Te,mu,DeltaE);
      }
#endif      
      // ***************
      // * RAD RECOMB  *
      // ***************
      if(DeltaE>0)
      {
        if(expint > 0 )
        {
          //k_MPI_z1_z0[i][j][0]=v_e*k_RR_fact1*1.0*k_RR_fact2*pow((DeltaE)*J2eV/STATES_z1[j][2],1.5)*expint*EXPR(a);
          k_MPI_z1_z0[i][j][0]=v_e*k_RR_fact1*1.0*k_RR_fact2*pow((DeltaE)/STATES_z1[j][2],1.5)*expint*EXPR(a);
        }
        else 
        {
          k_MPI_z1_z0[i][j][0]=0.0; //EXPR(a) kann +Inf werden--> unsinnige rate coeff.. wenn einer der faktoren=0 --> rest egal
        }
      }
      // ******************
      // * 3-body recomb  *
      // ******************
      kronecker=0;
      tmp0=2*STATES_z1[j][3]/STATES_z0[i][3];
      // tmp1=3.0*(1.0-kronecker)/2.0;
      tmp2=EXPR(-(DeltaE-IPD0)/kbTe);
      if(k_EI_z0_z1[i][j]==0) //wohl wichtig
        k_EI_z1_z0[i][j]=0; 
      else       
        k_EI_z1_z0[i][j]=k_EI_z0_z1[i][j]/tmp0/pow_two_pi_me_kT_hsq_tmp1/tmp2;

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
  #pragma omp for nowait schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2,sigma1,sigma_MPI_2,sigma_MPI_3,tmp0,tmp1,tmp2)
#endif
  for(i=0;i<z1_len;++i)
  {
    for(j=0;j<z2_len;++j)
    {
      kronecker=0.0;
      if(STATES_z1[i][4]==STATES_z2[j][4])
        kronecker=1.0;

      DeltaE=(STATES_z2[j][2]-STATES_z1[i][2])-IPD1+EF;

if(DeltaE <0 )
  continue;

#ifdef DOIPD
      DeltaE=MAX(0.0,DeltaE);
#endif            
      a=DeltaE/kbTe;
      expint=ExpInt(a);

if(expint==-1)
#ifndef OMP
        return -1;
#else
        fail=1;
#endif

      G2=genexpint(a,1.0,1.0);
      I_1=EXPR(-a)/a - expint;
      I_2=I_1*log54_beta_i+expint/a-G2;
      k_EI_z1_z2[i][j]=v_e*four_pi_a0_sq*alpha_i*E_ion_H_div_kTe_sq*I_2;

#ifdef MULTIPHOTON
      double dE_SI=DeltaE*eV2J; 
      // *******************
      // * MPI 2 PHOTONS   *
      // *******************
      if(twophoton_energy >= DeltaE - IPD1)
      {
        sigma1=sigma_tmp*pow(dE_SI,2.5)/sqrt(dE_SI);
        sigma_MPI_2=sigma1*sigma1/nu_div_hnu_sq;
        k_MPI_z1_z2[i][j][0]=sigma_MPI_2*I_sq;
      }
      // *******************
      // * MPI 3 PHOTONS   *
      // *******************  
      if(threephoton_energy > DeltaE - IPD1)
      {
        sigma1=sigma_tmp*pow(dE_SI,2.5)/sqrt(dE_SI);
        sigma_MPI_3=sigma1*sigma1*sigma1/2.0/nu_div_nu_div_hnu_cub;
        k_MPI_z1_z2[i][j][1]=sigma_MPI_3*I_cu;//*prob*beta_pi(Te,mu,DeltaE);
      }
#endif      
      // **************
      // * RAD RECOMB *
      // **************
      if(DeltaE>0)
      {
        if(expint > 0 )
        {
          k_MPI_z2_z1[i][j][0]=v_e*k_RR_fact1*4.0*k_RR_fact2*pow((DeltaE)/STATES_z2[j][2],1.5)*expint*EXPR(a);          
        }
        else 
        {
          k_MPI_z2_z1[i][j][0]=0.0; //EXPR(a) kann +Inf werden--> unsinnige rate coeff.. wenn einer der faktoren=0 --> rest egal
        }      
      }
      // ******************
      // * 3-body recomb  *
      // ******************        
      kronecker=0;
      tmp0=2*STATES_z2[j][3]/STATES_z1[i][3];
      // tmp1=3.0*(1.0-kronecker)/2.0;
      tmp2=EXPR(-(DeltaE-IPD1)/kbTe);
      if(k_EI_z1_z2[i][j]==0) //wohl wichtig
        k_EI_z2_z1[i][j]=0; 
      else             
        k_EI_z2_z1[i][j]=k_EI_z1_z2[i][j]/tmp0/pow_two_pi_me_kT_hsq_tmp1/tmp2;      

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
  #pragma omp for nowait schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2,sigma1,sigma_MPI_2,sigma_MPI_3,tmp0,tmp1,tmp2)
#endif
  for(i=0;i<z2_len;++i)
  {
    for(j=0;j<z3_len;++j)
    {
      kronecker=0.0;
      if(STATES_z2[i][4]==STATES_z3[j][4])
      kronecker=1.0;

      DeltaE=(STATES_z3[j][2]-STATES_z2[i][2])-IPD2+EF;

if(DeltaE <0 )
  continue;

#ifdef DOIPD
      DeltaE=MAX(0.0,DeltaE);
#endif            
      a=DeltaE/kbTe;
         expint=ExpInt(a);

if(expint==-1)
#ifndef OMP
        return -1;
#else
        fail=1;
#endif

      G2=genexpint(a,1.0,1.0);
      I_1=EXPR(-a)/a - expint;
      I_2=I_1*log54_beta_i+expint/a-G2;
      k_EI_z2_z3[i][j]=v_e*four_pi_a0_sq*alpha_i*E_ion_H_div_kTe_sq*I_2;

#ifdef MULTIPHOTON
      double dE_SI=DeltaE*eV2J;
      // *******************
      // * MPI 2 PHOTONS   *
      // *******************
      if(twophoton_energy > DeltaE - IPD2)
      {
        sigma1=sigma_tmp*pow(dE_SI,2.5)/sqrt(dE_SI);
        sigma_MPI_2=sigma1*sigma1/nu_div_hnu_sq;
        k_MPI_z2_z3[i][j][0]=sigma_MPI_2*I_sq;
      }
      // *******************
      // * MPI 3 PHOTONS   *
      // *******************  
      if(threephoton_energy > DeltaE -IPD2)
      {
        sigma1=sigma_tmp*pow(dE_SI,2.5)/sqrt(dE_SI);
        sigma_MPI_3=sigma1*sigma1*sigma1/2.0/nu_div_nu_div_hnu_cub;
        k_MPI_z2_z3[i][j][1]=sigma_MPI_3*I_cu;//*prob*beta_pi(Te,mu,DeltaE);
      }
#endif      
      // **************
      // * RAD RECOMB *
      // **************
      if(DeltaE>0)
      {
        if(expint > 0)
        {
          k_MPI_z3_z2[i][j][0]=v_e*k_RR_fact1*4.0*k_RR_fact2*pow((DeltaE)/STATES_z3[j][2],1.5)*expint*EXPR(a);  
        }
        else
        {
          k_MPI_z3_z2[i][j][0]=0.0;
        }
        
      }
      
      // ******************
      // * 3-body recomb  *
      // ******************         
      kronecker=0;
      tmp0=2*STATES_z3[j][3]/STATES_z2[i][3];
      // tmp1=3.0*(1.0-kronecker)/2.0;
      tmp2=EXPR(-(DeltaE-IPD2)/kbTe);

      if(k_EI_z2_z3[i][j] > 0)
      {
        k_EI_z3_z2[i][j]=k_EI_z2_z3[i][j]/tmp0/pow_two_pi_me_kT_hsq_tmp1/tmp2;
      }
      else
      {
        k_EI_z3_z2[i][j]=0.0;
      }
      
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
  #pragma omp for nowait schedule(static) private(j,kronecker,DeltaE,a,expint,G2,I_1,I_2,sigma1,sigma_MPI_2,sigma_MPI_3,tmp0,tmp1,tmp2)
#endif
  for(i=0;i<z3_len;++i)
  {
    for(j=0;j<z4_len;++j)
    {
      kronecker=0.0;
      if(STATES_z3[i][4]==STATES_z4[j][4])
      kronecker=1.0;

      DeltaE=(STATES_z4[j][2]-STATES_z3[i][2])-IPD3+EF;

if(DeltaE <0 )
  continue;      

#ifdef DOIPD
      DeltaE=MAX(0.0,DeltaE);
#endif            
      a=DeltaE/kbTe;
      expint=ExpInt(a);

if(expint==-1)
#ifndef OMP
        return -1;
#else
        fail=1;
#endif

      G2=genexpint(a,1.0,1.0);
      I_1=EXPR(-a)/a - expint;
      I_2=I_1*log54_beta_i+expint/a-G2;
      k_EI_z3_z4[i][j]=v_e*four_pi_a0_sq*alpha_i*E_ion_H_div_kTe_sq*I_2;

#ifdef MULTIPHOTON
      double dE_SI=DeltaE*eV2J;
      // *******************
      // * MPI 2 PHOTONS   *
      // *******************
      if(twophoton_energy> DeltaE - IPD3)
      {
        sigma1=sigma_tmp*pow(dE_SI,2.5)/sqrt(dE_SI);
        sigma_MPI_2=sigma1*sigma1/nu_div_hnu_sq;
        k_MPI_z3_z4[i][j][0]=sigma_MPI_2*I_sq;
      }
      // *******************
      // * MPI 3 PHOTONS   *
      // *******************  
      if(threephoton_energy > DeltaE - IPD3)
      {        
        sigma1=sigma_tmp*pow(dE_SI,2.5)/sqrt(dE_SI);
        sigma_MPI_3=sigma1*sigma1*sigma1/2.0/nu_div_nu_div_hnu_cub;
        k_MPI_z3_z4[i][j][1]=sigma_MPI_3*I_cu;//*prob*beta_pi(Te,mu,DeltaE);
      }
#endif
      // **************
      // * RAD RECOMB *
      // **************
      if(DeltaE>0)
      {
        if(expint>0)
        {
          k_MPI_z4_z3[i][j][0]=v_e*k_RR_fact1*4.0*k_RR_fact2*pow((DeltaE)/STATES_z4[j][2],1.5)*expint*EXPR(a);  
        }
        else
        {
          k_MPI_z4_z3[i][j][0]=0.0;
        }       
        
      }
      // ******************
      // * 3-body recomb  *
      // ******************
      //if(DeltaE>0.0)
      kronecker=0;
      tmp0=2*STATES_z4[j][3]/STATES_z3[i][3];
      // tmp1=3.0*(1.0-kronecker)/2.0;
      tmp2=EXPR(-(DeltaE-IPD3)/kbTe);
      if(k_EI_z3_z4[i][j] >9 )
      {
        k_EI_z4_z3[i][j]=k_EI_z3_z4[i][j]/tmp0/pow_two_pi_me_kT_hsq_tmp1/tmp2;
      }
      else
      {
        k_EI_z4_z3[i][j]= 0.0;
      }
      
      
    }
  }
  if(fail==1) return -1; 
#endif // MAXLEVEL > 3

  // printf("MYID:%d,geschaftt!\n",myid);
    return 0;

}

// ***************************************************************************************
double ExpInt(double x)
{
 if(x>350) 
 {
      return 0; //wird extrem klein -> inderflow
 }  
  
   
   return gsl_sf_expint_E1(x);

   //Approx by Swamee and Ohija
   double A=log((0.56146/x+0.65)*(1+x));
   //Numeric limit
   double B;
   B=x*x*x*x*exp(7.7*x)*pow(2.0+x,3.7);
   return pow(pow(A,-7.7)+B,-0.13);
}


#ifdef USEFLOAT
float fak(float t, float x, float j,float s) //aux. function for genexpint
{
  return powf(t,x-1)*powf( powf(logf(-logf(t)) , j) / -logf(t),s);
}
#else
double fak(double t, double x, double j,double s) //aux. function for genexpint
{
  return pow(t,x-1)*pow( pow(log(-log(t)) , j) / -log(t),s);
}
#endif



#ifdef USEFLOAT
float genexpint(float x,float ss,float j) 
{
  int maks=5;
  float eps=1E-6;

  float b=0.36787944117f; //exp(-1)
  float s_old=0;

  float i=0.0;
  float s=0.0;
  float t=0.0;
  float d=0.0;
  float dd=0.0;
  float sum=0.0;
  float m=0.0;
  int n,k;

  for(n=1;n<maks+1;++n)
  {
    if(n==1)
    {
        i=1.0;
        float logt=logf(0.5*b);
        float simplefak=-powf(0.5*b,x-1.0)*logf(-logt)/logt;
        s=b*simplefak;

    }
    else
    {
      m=i;
      d=b/(3*m);
      dd=2*d;
      t=0.5*d;
      sum=0.0;

      for(k=1;k<m+1;++k)
      {
        float logt=logf(t);
        float simplefak=-powf(t,x-1.0)*logf(-logt)/logt;        
        sum+=simplefak;
        t=t+dd;
        logt=logf(t);
        simplefak=-powf(t,x-1.0)*logf(-logt)/logt;
        sum+=simplefak;
        t=t+d;
      }
      i=i*3.0;
      s=(s+b*sum/m)/3.0;
    }
    s=s/tgammaf(j+1.0);
    if(fabs(s_old-s)<=eps)
      break;
    else
       s_old=s;
  }
  return s; ///tgamma(j+1.0);
  //return s/tgamma1pm1(j+1.0); //,policy<digits10<3> >());
}

#else
double genexpint(double x,double ss,double j) 
// G_k(x)= E_1^(k-1) (x) wobei  
// E_ss^j= 1/(Gamma(j+1)) * int_1^inf (ln(t))^j * t^(-ss) e^(-z*t) dt
{
  //Nuri Ozlap, Elgiz Bairamov, "uniform convergence and computation of the generalized exponential integrals"
  //J.Math.Chem. 49:520-530 (2011)
  // Siehe auch. The effcient computation of some generalisedexponential integrals
  //

  //  double j=1.0; //j=k-1 (k=2)
  //  double ss=1.0;

  // für x=60 ---> 4.7243e-45

  // if(x>=60)
  //  return 5e-45;

  int maks=5;
  double eps=1E-12;

  eps=1E-6;  

  double b=exp(-1.);
  double s_old=0;

  double i=0.0;
  double s=0.0;
  double t=0.0;
  double d=0.0;
  double dd=0.0;
  double sum=0.0;
  double m=0.0;
  int n,k;
  for(n=1;n<maks+1;++n)
  {
    if(n==1)
    {
        i=1.0;
        //s=b*fak(0.5*b,x,j,ss);
        double logt=log(0.5*b);
        double simplefak=-pow(0.5*b,x-1.0)*log(-logt)/logt;
        s=b*simplefak;

    }
    else
    {
      m=i;
      d=b/(3*m);
      dd=2*d;
      t=0.5*d;
      sum=0.0;

      for(k=1;k<m+1;++k)
      {
        double logt=log(t);
        double simplefak=-pow(t,x-1.0)*log(-logt)/logt;        

        //sum=sum+fak(t,x,j,ss);
        sum+=simplefak;

        t=t+dd;

        logt=log(t);
        simplefak=-pow(t,x-1.0)*log(-logt)/logt;

        //sum=sum+fak(t,x,j,ss);

        sum+=simplefak;
        t=t+d;
      }
      i=i*3.0;
      s=(s+b*sum/m)/3.0;
    }
    s=s/tgamma(j+1.0);

    if(fabs(s_old-s)<=eps)
      break;
    else
       s_old=s;
  }
  return s; ///tgamma(j+1.0);
  //return s/tgamma1pm1(j+1.0); //,policy<digits10<3> >());
}
#endif

// ************************************************
// *            RE-ABSORPTION KRAM
// ************************************************
double StarkWidth(double nu,double nl,double Te,double Ti,double Tr,double Ne,double Ni)
{
  //ACHTUNG (nu-nl) muss größer 0 sein !!! sonst w_qs=0 und => escape_fac=0

  //Te-input, Ti-input, Tr-input in eV !
  //Zi=degree of ionization

  double Zcore=13.0;
  double Zr=Zcore-1.0;

  Ni=Ni/(1/pow(bohr_radius,3.0)); // m^-3 to bohr_radius^-3
  Ne=Ne/(1/pow(bohr_radius,3.0));
  double Zi=Ne/Ni;

  Te=Te*eV2H;
  Ti=Ti*eV2H;
  Tr=Tr*eV2H; //radiators temp.
  double mi=26.981*AMU/EMASS; //mass in units of emass
  double me=1.0;
  double mr=1.0*AMU/EMASS;
  double ve=sqrt(Te/me);
  double vi=sqrt(Ti/mi);
  double vr=sqrt(Tr/mi);

  double re=pow(3.0/4.0/pi/Ne,1.0/3.0); //mean distance
  double ri=pow(3.0/4.0/pi/Ni,1.0/3.0);
  double qi=Zi;
  double qe=1.0;
  double qr=-Zr; //mui importante
  double F=2.0*pi*pow(4.0/15.0,2.0/3.0)*qe*pow(Ne,2.0/3.0); //Holstmarkian field
  double debye=sqrt(Te/4.0/pi/Ne/qe/qe);
  double eta_pert=(1.0+re/debye)*exp(-re/debye);
  double eta_rad=exp(-qr*qe/re/Te);
  double omega_mean=1.4385;
  double kronecker=0;
  if(nu-nl==1.0)
        kronecker=0.5;
  else
     kronecker=1.0;
  double w_qs=omega_mean*3.0*(nu*nu-nl*nl)/Zcore*F*eta_pert*eta_rad*kronecker;
  double f_F=sqrt(ve*ve+vr*vr)/re; //microfield freq.
  double R_s=w_qs/f_F/eta_rad; // static/dynamic stark ratio
  double f_s=R_s/(R_s+0.5); //quasi-staticity factor
  double w_ss=f_s*w_qs;  //Full stark width in hartree
  double w_J=w_ss*4.359745e-18; //hartree to Joule

  if(w_ss==0) //escape factor will be NaN
  {
        printf("WARNING:w_ss =0\n");
        printf("R_s:%f,w_qs:%f,kron:%f,f_F:%f\n",R_s,w_qs,kronecker,f_F);
  }
  //bandwidth conversion to meters
  //double delta_lambda=w_J*lambda0*lambda0/planck/LIGHTSPEED; //in m
  return w_J/planck; /// in Hz
  //return delta_lambda;
}

//eigentl. tansmission factor (c.f. 
double EscapeFactor(double w,double tau0)
{
  //w in Angstrom
  /*
  double a=5.441*exp(-18.55*tau0)+0.668*exp(-0.0515*tau0);
  double b=1.259e-5*pow(tau0,4.0)-2.694e-4*pow(tau0,3.0)+1.877e-3*pow(tau0,2.0)-1.625e-3*tau0+6.887e-3;
  double c=5.773*exp(-18.63*tau0)+0.7008*exp(-0.05994*tau0);
  double d=(2.245*tau0*tau0-0.08527*tau0+0.1311)/(tau0*tau0+8.936*tau0+2.649);
  
  double T=(a*w*w+b*w)/(c*w*w+d);
  */
  double lambda_lo=-10.0*w;
  double lambda_hi=+10.0*w;
  double T=trapz(lambda_lo,lambda_hi,100,w,tau0);
  return T;
}

double Lorentzian(double w,double lambda) //both params in angs!
{
  //Area-normalized
  //https://magicplot.com/wiki/fit_equations
  //return 1.0/pi*0.5*w/(pow(lambda,2.0)+pow(0.5*w,2.0));
  return 1.0/(w/2.0)/pi*1.0/(1.0+pow(lambda/(w/2.0),2.0));
}
double integrand(double w,double lambda,double tau,double phi0)
{
  double exponent=tau*phi0*Lorentzian(w,lambda);
  //printf("tau:%.4e,phi0:%.4e,Lor:%.4e,e:%.4e\n",
  //      tau,phi0,Lorentzian(w,lambda),exponent);
  return exp(-exponent)*Lorentzian(w,lambda);
}
double trapz(double a,double b,int n,double w,double tau)
{
  double phi0=Lorentzian(w,0.0);
  double area=0.0;
  double sum=0.0;
  double h=(b-a)/n;
  int x;
  for(x=1;x<n;++x)
  {
    sum+=integrand(w,((double) x)*h+a,tau,phi0);
  }
  area=h/2.0*(integrand(w,a,tau,phi0)+integrand(w,b,tau,phi0)+2*sum);
  return area;
}



// Springer Handbook of atomic,molecular and optical physics  S.838 (Part D)
// c.f. https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas-notation-data-and-formulas/atomic-spectroscopy
// Achtung: Geht nur wenn DeltaN > 0, also übergänge in selbe schalen nicht
// möglich. Grund: Näherung des Gaunt-Faktors
double EinsteinCoeff(double n1,double n2,double g2,double DeltaE)
{
  double Delta_Lambda=planck*LIGHTSPEED/DeltaE;
  double DeltaN=n2-n1;
  double eps1=1.0/n1/n1;
  double eps2=1.0/n2/n2;
  double Gaunt=1.0-0.25/fabs(DeltaN);
  double S; //line strength
  double A21;
  double z=13.0;
  double nominator=32.0/pi/sqrt(3)*pow(ECHARGE*bohr_radius/z,2.0)*pow(eps1*eps2,1.5)*Gaunt;
  double denom=pow(eps1-eps2,4.0);
  //S=32.0/pi/sqrt(3)*pow(ECHARGE*bohr_radius/z,2.0)*pow(eps1*eps2,1.5)/pow(eps1-eps2,4.0)*Gaunt;
  S=nominator/denom;
  A21=16.0*pow(pi,3.0)/3.0/planck/ECONST/pow(Delta_Lambda,3.0)/g2*S;

  /*
  {
        //printf("Error:A21 is NaN\n");
        printf("A21:%.4e,Gaunt:%f,deltaN:%f,S:%.4e,eps1:%f,eps2:%f,Nominator:%.4e,Denom:%.4e,g2:%f,Delta_Lambda:%.4e\n",
                A21,Gaunt,DeltaN,S,eps1,eps2,nominator,denom,g2,Delta_Lambda);
    
  }
  */
  return A21;
}

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
            fprintf(outfile, " %.4e ", Ith(l1[i][j][k].y,l));
          }
          fprintf(outfile,"\n");
        }        
      }

    }
    fclose(outfile);
  return 0; //alles ok
}

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
  double tmp;


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
            Ith(l1[i][j][k].y,l)=tmp;
            
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


