#include "imd.h"
/******************************************************************
*
*  ACHTUNG: IMD nutzt links-händiges koordinatensystem!
*  	    -->Die Gleichungen aus den Lehrbüchern dürfen nicht einfach so übernommen werden, da
*	       Standartmäßig ein rechtshändiges Koordsys. benutzt wird
*	       Im Detail: y=-x* , x=y*, z=z*
*	      ( F_x )     (  F_y  ) 
*	      ( F_y )  =  ( -F_x  )
*	      ( F_z )     (  F_z  )
*		       
* 
*  2D-und 1D- Maxwell Solver 
*  Bedingung: global_fd_dim.z = 1
*	      bzw. selbes in y-Richtung für 1D
*  Bisher nur TMZ-MODE für 2D und Z-Polarisation für 1D-Mode
*
*  Fuer kurze Tests bietet sich matlab-script an: fdtd_2d_suite.m
*
*  Dispersion wird gehandhabt mittels
*  Auxiliaray differential equation method (ADE)
*  Bisher nur für Drude-Medium
* 
*  Schema für Field-update equations aus: 
*  Taflove, "Advances in FDTD Computational Electrodynamics: Photonics and Nanotechnology" S.170, S.271, S.506
*  Genaueres zur ADE-Methode in  Taflove, "Computational Electrodynamics" Kapitel 9.4, S.213 oder S. 398
*  Genaueres zur Berenger-Split-Field Methode für PMLs ebenfalls im Taflove
*  In Taflove S.180 gibt es Drude-Lorenz Parameter füer Gold,Silber,Alu,Chrom
* 
*  Die update-equations sind auch zu finden in 
*  "Electromagnetics of Time-Varying Complex Media", S. 302
*
*  ACHTUNG: Man sollte darauf achten, dass an den Rändern, dort wo die EM-Wellen absorbiert werden vakuum ist!
*  ACHTUNG: Bei NanS mal die Optimierungs-Flags weglassen...Weiss noch nicht wo genau BS gemacht wird...
*  
*
* Der Waist-Radius w_0 ist der Radius,bei dem die optische Intensität auf 1/e^2 abfällt (bzw. E auf 1/e*E0)
* dabei geht die Größe in folgender Form in die Abhängigkeit der
* Intensität: I(x,y) =I_0 exp(-(x^2+y^2)/w_0^2)       (d.h. Zentrum x0=0,y0=0)
* bzw. E(x,y) = E_0*exp(-(x^2+y^2)/w_0^2/2)
* ***********************************************************************************************
* CODE LOG:
*************************************************************************************************
* 25.12.18: in fdtd_comm_ghost_cells die Felder an den nicht-pbc randzellen Nullen!
* 26.12.18: 1D- Maxwellsolver hinzugefügtmit ABC (2nd order) oder PML  + Total-Field/Scattered-Field
*	    BUGFIX: Bei Softsource muss Courant-Zahl Sc mit einfließen: E+=Sc*E0 
*
* 27.12.18: Es stehen nun 5 Methoden für Laser-Energie-Absorption zur Verfügung (siehe FDTD1D-do_fdtd routine)
*	    Fazit: Je genauer die Näherung, umso geringer E_abs (nur wenig), ausser für slowly-varying envelope mit 
*	           monochromatischem Laser (keine Dispersion)==> Hier ist E_abs ca. faktor 2 Kleiner als
*		   bei allen anderen
*
***********************************************************/
/*********************
*   AFTER TTM-INIT   *
**********************/
#define LORENTZ
//#define FDTD2D
#define FDTD1D
// *******************
//SOFT SOURCE
int fdtd_softsource(double t) //returns -1 if Intensity is too low --> no need for do_fdtd
{
  
  // E(x0)+= Hinc*dt/eps0/dx
  // Hinc=srcfun(n+1)/Imp0

  double nd,sine,spatial;

  //double t0,dur;
  int ilocal;
  double disty; //dist from srcy_int
//  double srcw;  // 1/e -width
  int j,jglobal;
 
  real Hinc,Einc; 
  real ezx2=(dt_FDTD/fd_dy)/(eps0);
  real ezy2=(dt_FDTD/fd_dx)/(eps0);

  sine=sin(2.0*M_PI*(c0/(fd_dx*Nlambdax))*(t+dt_FDTD));
  if(containssrc==1) 
  {
    ilocal=srcx_int+1-my_coord.x*(local_fd_dim.x-2); 

    // ************************
    // *    2D SOFTSOURCE
    // *************************    
#ifdef FDTD2D    
    for(j=1;j<local_fd_dim.y-1;++j)
    {
    	jglobal=(j-1) + my_coord.y*(local_fd_dim.y-2);
    	disty=(double) ABS(jglobal-srcy_int);
    	spatial=exp(-0.5* ((double) disty*disty/srcw_int/srcw_int));

    	if(jglobal>bw && jglobal<global_fd_dim.y-1-bw)
    	{
    	  Hinc=laser_timefun*sine*spatial/Imp0; //mal 0.5 weil energie sich auf TMZ und TEZ aufspaltet
    	  Einc=Hinc*Imp0;
        //TMZ
     	  l1[ilocal][j][1].Ezx+=dt_FDTD/eps0/fd_dx*Hinc;
    	  l1[ilocal][j][1].Ezy+=dt_FDTD/eps0/fd_dy*Hinc;

        //TEZ
    	  l1[ilocal][j][1].Hzx-=dt_FDTD/mu0/fd_dx*Einc; 
    	  l1[ilocal][j][1].Hzy-=dt_FDTD/mu0/fd_dy*Einc;
    	}
    }    
#endif

    // ************************
    // *    1D SOFTSOURCE
    // ************************* 
#ifdef FDTD1D
    //l1[ilocal][1][1].Ezx+=laser_timefun*sine*Sc;
    Hinc=laser_timefun*sine/Imp0; //*Sc
    l1[ilocal][1][1].Ezx+=ezy2*Hinc; //ezy2 ist korrekt
#endif
  }
//  iglobal = (i-1) + my_coord.x*(local_fd_dim.x-2);
  return 1; //alles ok
} //softsource
// *******************************************************************************************************

// *****************************************************
// *          1D VERSION OF UPDATE EQUATIONS           *
// *****************************************************
#ifdef FDTD1D
void do_fdtd(double t) //do a single fdtd-step
{ 
  int i;
  int iglobal;
  double ezj;
  double hy1,hy2,hx1,hx2;
  double ezx1,ezx2,ezy1,ezy2;
  double mu,eps,sigmax,sigmay,sigmastarx,sigmastary;
  double beta_p,k_p,gamma_p;
  double Ez_present,Jz_present;
  double qe;    //instantaneous electrical dissipation power density

  //double  atmp,btmp,integ_const; //variables for integration of qe
  //double   qtn1,qtn2; //qe(t(n+1)) und qe(t(n))

  double Ez_old,Ez_old_old,Ez_deriv;

///////////////// vars for harmonic approx. of joule heating (dielectric hysteresis) ///////////////////////
  double gamma_plasma;
  double omega_plasma;
  double sigma_omega;

  //Koeffs for drude-lorentz dispersion
  double OmegaL,ZetaL,TauL,RhoL,ezp;
  //Lorentz - Params
  double deltaeps_L=1.0;
  double Gamma_L=9.84485242e14;
  double Omega0_L=2.47640578e15; //interband freq  
  double omegapl_L=9.6929256835e15; //entspricht 800 nm
  double epsinf=2.73;
  double dtsq=dt_FDTD*dt_FDTD;
  double wlsq=Omega0_L*Omega0_L;

  gamma_p=1.7927354e14;
  omega_plasma=2.0312604e16;

  deltaeps_L=omegapl_L*omegapl_L/Omega0_L/Omega0_L;

  int skip_fdtd=0;
  skip_fdtd=fdtd_softsource(t);

  if(skip_fdtd==-1)  // Achtung:Resonanz-effekte können EM-feld auch lange nach
    return;	     // laser-strahl aufrecht erhalten
		     // dies gilt es zu pruefen
 // ***************************************************
  // *  Ezx  & Ezy update: von n --> n+1 (gegenwart)
  // **************************************************
  fdtd_comm_ghost_cells();
  int j=1;
  for(i=1;i<local_fd_dim.x-1;++i)
  {
  	iglobal =  ((i-1) + my_coord.x*(local_fd_dim.x-2));

  	if(l1[i][j][1].natoms>=fd_min_atoms) //mit dispersion
  	{

  	  // ezx1=(2.0*eps0-dt_FDTD*beta_p)/(2.0*eps0+dt_FDTD*beta_p);
  	  // ezx2=2.0*dt_FDTD/(2.0*eps0+dt_FDTD*beta_p)/fd_dx;
    	// zj=ezx2*fd_dx;
      // l1[i][j][1].Ezx=ezx1*l1[i][j][1].Ezx+ezx2*(l1[i][j][1].Hy-l1[i-1][j][1].Hy);

     	  //save "old" Ezx/Ezy for update of Jzx/Jzy befor Ezx/Ezy update to n+1
   	  l2[i][j][1].Ezx=l1[i][j][1].Ezx;     //@n 
      l2[i][j][1].Jzx=l1[i][j][1].Jzx;
      l2[i][j][1].Jzy=l1[i][j][1].Jzy;
      l2[i][j][1].Jlzx=l1[i][j][1].Jlzx;
      l2[i][j][1].Jlzy=l1[i][j][1].Jlzy;

      epsinf=l1[i][j][1].DL[0];
      omega_plasma=l1[i][j][1].DL[5];
      gamma_p=l1[i][j][1].DL[1];
      omegapl_L=l1[i][j][1].DL[2];
      Omega0_L=l1[i][j][1].DL[3];
      Gamma_L=l1[i][j][1].DL[4];

      wlsq=Omega0_L*Omega0_L;
      double wplsq=omegapl_L*omegapl_L;
      deltaeps_L=omegapl_L*omegapl_L/wlsq;
      double Cb=dt_FDTD/eps0/epsinf;
      //Update 
      l1[i][j][1].Ezx=l1[i][j][1].Ezx
                      +Cb*(l1[i][j][1].Hy-l1[i-1][j][1].Hy)/fd_dx
                      -Cb*l1[i][j][1].Jzx
                      -Cb*l1[i][j][1].Jlzx;

      l1[i][j][1].Pzx=l1[i][j][1].Pzx+dt_FDTD*l1[i][j][1].Jlzx;


  }
  else //standard für vakuum & pml
  {
      eps=eps0;
      sigmay=l1[i][j][1].sigmay;
      ezx1=(eps-0.5*dt_FDTD*sigmay)/(eps+0.5*dt_FDTD*sigmay);
      ezx2=(dt_FDTD/fd_dx)/(eps+0.5*dt_FDTD*sigmay);

      l1[i][j][1].Ezx=ezx1*l1[i][j][1].Ezx
                     +ezx2*(l1[i][j][1].Hy-l1[i-1][j][1].Hy);	  
  }	
}
// ********************************************
// *  Hy  update (zukunfts Hy @ n+3/2  *
// ********************************************
fdtd_comm_ghost_cells();
j=1;
for(i=1;i<local_fd_dim.x-1;++i)
{
  mu=mu0;//l1[i][j][1].mu;
  eps=l1[i][j][1].eps;

  sigmay=l1[i][j][1].sigmay;
  sigmastary=sigmay*mu/eps;

  
  hy1=(mu-0.5*dt_FDTD*sigmastary)/(mu+0.5*dt_FDTD*sigmastary);
  hy2=(dt_FDTD/fd_dx)/(mu+0.5*dt_FDTD*sigmastary);
  //Hy-Updt.
  l1[i][j][1].Hy=hy1*l1[i][j][1].Hy+hy2*(l1[i+1][j][1].Ezx-l1[i][j][1].Ezx);

  if(l1[i][j][1].natoms>=fd_min_atoms) //Für vakuum ist der rest nicht noetig
  {

    epsinf=l1[i][j][1].DL[0];
    omega_plasma=l1[i][j][1].DL[5];
    gamma_p=l1[i][j][1].DL[1];
    omegapl_L=l1[i][j][1].DL[2];
    Omega0_L=l1[i][j][1].DL[3];
    Gamma_L=l1[i][j][1].DL[4];

    double CD1=(2-dt_FDTD*gamma_p)/(2+dt_FDTD*gamma_p);
    double CL1=(2-dt_FDTD*Gamma_L)/(2+dt_FDTD*Gamma_L);
    double CD2=2*dt_FDTD/(2+dt_FDTD*gamma_p);
    double CL2=2*dt_FDTD/(2+dt_FDTD*Gamma_L);

    double omegapl_D_sq=omega_plasma*omega_plasma;
    double omegapl_L_sq=omegapl_L*omegapl_L;
    double omega0_L_sq=Omega0_L*Omega0_L;

    //J-Drude Updt.
    l1[i][j][1].Jzx=CD1*l1[i][j][1].Jzx+CD2*(omegapl_D_sq*eps0*l1[i][j][1].Ezx);
    //J-Lorentz Updt.
    l1[i][j][1].Jlzx=CL1*l1[i][j][1].Jlzx+CL2*(omegapl_L_sq*eps0*l1[i][j][1].Ezx-omega0_L_sq*l1[i][j][1].Pzx);

    //Leistungsdichte nur für aktive Zellen berechnen
    Ez_present=l1[i][j][1].Ezx; //@n+1
    Jz_present=(l1[i][j][1].Jzx+l2[i][j][1].Jzx)*0.5; //@n+1
                                                                                                                              
    double Jlz_present=0.5*(l1[i][j][1].Jlzx+l2[i][j][1].Jlzx);
    
    //Comp. power density
    double qe_Drude= gamma_p*(Jz_present*Jz_present)/eps0/omegapl_D_sq;
    //double qe_Lorentz=Gamma_L/deltaeps_L/eps0/wlsq*(Jlz_present*Jlz_present+Jlx_present*Jlx_present+Jly_present*Jly_present);
    double qe_Lorentz=Gamma_L/eps0/omegapl_L_sq*(Jlz_present*Jlz_present);
    qe=qe_Drude+qe_Lorentz;
    l1[i][j][1].source=6.3538562638e-26*qe;//convert to imd-units       

  }
}
  // ****************************************************************
  //  Hz source for TF/SF: Um Aufspaltung der Welle in beide Richtungen zu unterbinden! Siehe Matlab-skript fdtd.m
  //  Leider nicht 100% vermeidbar. Siehe fdtd2_4_6_18_serial.m --> etwa 1% von E0 wandert auch nach -x --> pml kümmert sich darum  
  //  Außerdem gibt es ja auch noch die reflektierte Welle. 
  // **************************************************************
  
if(containssrc)
{
  int ilocal=srcx_int+1-my_coord.x*(local_fd_dim.x-2);
  real sine=sin(2.0*M_PI*(c0/(fd_dx*Nlambdax))*t);
  real temporal=exp(-0.5*pow(t-laser_t_0,2.0)/laser_sigma_t_squared);
  real Einc=sine*temporal*E0;//*Sc;
  //H kompensieren, damit keine Welle Richtung -x ensteht!
  hy2=(dt_FDTD/fd_dx)/(mu0);
  l1[ilocal-1][1][1].Hy-=hy2*Einc;

  //real tf2=2.0*E0*exp(-0.5*pow(t-dt_FDTD-laser_t_0,2.0)/laser_sigma_t_squared)/Imp0;
  //l1[ilocal-1][1][1].Hy-=tf2*sine*Sc;
}

}
#endif //FDTD1D

// **********************************************
// *  2D VERSION OF MAXWELL SOLVER (TMZ-MODE)   *
// **********************************************
#ifdef FDTD2D
void do_fdtd(double t) //do a single fdtd-step
{ 
  int i,j; 
  int istart,iend;
  int jstart, jend; //for Ez-update. Shift by +1 if first global i/j
  int iglobal,jglobal;

  //medium props
  double mu,eps,sigmax,sigmay,sigmastarx,sigmastary;
  double beta_p,k_p;
 
  double omega_plasma=2.2955e+16; 
  double gamma_p=1.1174e+15;
  double Gamma_L=9.84485242e14;
  double Omega0_L=2.47640578e15; //interband freq  
  double omegapl_L=9.6929256835e15; //entspricht 800 nm
  double epsinf=2.73;


//Aus imd_tmm.c
/*
  gamma_p=1.7927354e14; //Drude-Lorentz
  omega_plasma=2.0312604e16;
  omegapl_L=9.6929256835e15;
  Gamma_L=9.84485242e14;
  Omega0_L=2.47640578e15;
  epsinf=2.73;
*/

  double dtsq=dt_FDTD*dt_FDTD;
  double wlsq=Omega0_L*Omega0_L;
  double deltaeps_L=omegapl_L*omegapl_L/wlsq;


  double qe;    //instantaneous electrical dissipation power density
  double Ez_present,Jz_present;
  
  //coeffs for update eqn.s
  double e1,e2,ej;
  double ezx1,ezx2, ezy1,ezy2,ex1,ex2,ey1,ey2; //wegen sigmax,sigmay etc.
  double hy1,hy2,hx1,hx2; 
  double hzx1,hzx2,hzy1,hzy2;

  //DRUDE-LORENTZ coeffs
  double b1,b2,b3;
  double c1,c2,c3,c4,c5; 
  double denom,denom2;
  

  int skip_fdtd=0;
  skip_fdtd=fdtd_softsource(t);

  if(skip_fdtd==-1)  // Achtung:Resonanz-effekte können EM-feld auch lange nach
    return;	     // laser-strahl aufrecht erhalten
		     // dies gilt es zu pruefen
 // ***************************************************
  // *  Ezx  & Ezy update: von n --> n+1 (gegenwart)
  // **************************************************
  fdtd_comm_ghost_cells();

  //Save old field values (incl. buffer cells)  
  //before update equations
  for(i=0;i<local_fd_dim.x;++i)    
  {
    for(j=0;j<local_fd_dim.y;++j)  
    {
          l2[i][j][1].Ezx=l1[i][j][1].Ezx; //@n
          l2[i][j][1].Ezy=l1[i][j][1].Ezy; //@n

          l2[i][j][1].Ex=l1[i][j][1].Ex; //@n
          l2[i][j][1].Ey=l1[i][j][1].Ey; //@n

	  l2[i][j][1].Jzx=l1[i][j][1].Jzx;
	  l2[i][j][1].Jzy=l1[i][j][1].Jzy;

	  l2[i][j][1].Jx=l1[i][j][1].Jx;
	  l2[i][j][1].Jy=l1[i][j][1].Jy;

	  
	  l2[i][j][1].Jlzx=l1[i][j][1].Jlzx;
	  l2[i][j][1].Jlzy=l1[i][j][1].Jlzy;

	  l2[i][j][1].Jlx=l1[i][j][1].Jlx;
	  l2[i][j][1].Jly=l1[i][j][1].Jly;

    }
  }


  istart=1;
  iend=local_fd_dim.x-1;
  jstart=1;
  jend=local_fd_dim.y-1;

  if(my_coord.y==0) jstart=2;
  if(my_coord.y==cpu_dim.y-1) jend=local_fd_dim.y-2;
  if(my_coord.x==cpu_dim.x-1) iend=local_fd_dim.x-2;

//TESTCASE
int jmin=10;
int jmax=global_fd_dim.y-10;
 
 
  for(i=istart;i<local_fd_dim.x-1;++i)    //in dieser loop gibts kein l[i+1], dafuer aber l[i-1] 
  {
    iglobal = (i-1) + my_coord.x*(local_fd_dim.x-2);

//TESTCASE
//if(iglobal<srcx_int) continue;

    for(j=jstart;j<local_fd_dim.y-1;++j)  //in dieser loop gibts kein l[j+1], dafuer aber l[j-1] 
    {
      jglobal=(j-1) + my_coord.y*(local_fd_dim.y-2);

	if(l1[i][j][1].natoms>=fd_min_atoms) //mit Jp
	{
	  //if(steps>1)
	  {
	  epsinf=l1[i][j][1].DL[0];
	  omega_plasma=l1[i][j][1].DL[5];
  	  gamma_p=l1[i][j][1].DL[1];
	  omegapl_L=l1[i][j][1].DL[2];
	  Omega0_L=l1[i][j][1].DL[3];
	  Gamma_L=l1[i][j][1].DL[4];
	  }

	  wlsq=Omega0_L*Omega0_L;
	  double wplsq=omegapl_L*omegapl_L;
	  deltaeps_L=omegapl_L*omegapl_L/wlsq;

	  double Cb=dt_FDTD/eps0/epsinf;

	  l1[i][j][1].Ezx=l1[i][j][1].Ezx+Cb*(l1[i][j][1].Hy-l1[i-1][j][1].Hy)/fd_dx-Cb*l1[i][j][1].Jzx-Cb*l1[i][j][1].Jlzx; 
	  l1[i][j][1].Ezx=l1[i][j][1].Ezx-Cb*(l1[i][j][1].Hx-l1[i][j-1][1].Hx)/fd_dy-Cb*l1[i][j][1].Jzy-Cb*l1[i][j][1].Jlzy; 

          l1[i][j][1].Ex=l1[i][j][1].Ex + Cb*(l1[i][j][1].Hzx+l1[i][j][1].Hzy - l1[i][j-1][1].Hzx-l1[i][j-1][1].Hzy)/fd_dy - Cb*l1[i][j][1].Jx-Cb*l1[i][j][1].Jlx;
          l1[i][j][1].Ey=l1[i][j][1].Ey - Cb*(l1[i][j][1].Hzx+l1[i][j][1].Hzy - l1[i-1][j][1].Hzx-l1[i-1][j][1].Hzy)/fd_dx - Cb*l1[i][j][1].Jy-Cb*l1[i][j][1].Jly;;

	 
          l1[i][j][1].Pzx=l1[i][j][1].Pzx+dt_FDTD*l1[i][j][1].Jlzx;
          l1[i][j][1].Pzy=l1[i][j][1].Pzy+dt_FDTD*l1[i][j][1].Jlzy;
          l1[i][j][1].Px=l1[i][j][1].Px+dt_FDTD*l1[i][j][1].Jlx;
          l1[i][j][1].Py=l1[i][j][1].Py+dt_FDTD*l1[i][j][1].Jly;
	}
	else //vakuum & pml
	{
          eps=eps0;
          sigmay=l1[i][j][1].sigmax;
          sigmax=l1[i][j][1].sigmay;

          ezx1=(eps-0.5*dt_FDTD*sigmax)/(eps+0.5*dt_FDTD*sigmax);
          ezx2=(dt_FDTD)/(eps+0.5*dt_FDTD*sigmax);

          ezy1=(eps-0.5*dt_FDTD*sigmay)/(eps+0.5*dt_FDTD*sigmay);
          ezy2=(dt_FDTD)/(eps+0.5*dt_FDTD*sigmay);

	  ex1=ezy1; //(eps-0.5*dt_FDTD*sigmay)/(eps+0.5*dt_FDTD*sigmay);
	  ex2=ezy2; //(dt_FDTD)/(eps+0.5*dt_FDTD*sigmay);
		
          ey1=ezx1; //(eps-0.5*dt_FDTD*sigmax)/(eps+0.5*dt_FDTD*sigmax);
          ey2=ezx2; //(dt_FDTD)/(eps+0.5*dt_FDTD*sigmax);

          l1[i][j][1].Ezx=ezx1*l1[i][j][1].Ezx+ezx2*(l1[i][j][1].Hy-l1[i-1][j][1].Hy)/fd_dx;
          l1[i][j][1].Ezy=ezy1*l1[i][j][1].Ezy-ezy2*(l1[i][j][1].Hx-l1[i][j-1][1].Hx)/fd_dy;

          l1[i][j][1].Ex=ex1*l1[i][j][1].Ex + ex2*(l1[i][j][1].Hzx+l1[i][j][1].Hzy - l1[i][j-1][1].Hzx-l1[i][j-1][1].Hzy)/fd_dy;
          l1[i][j][1].Ey=ey1*l1[i][j][1].Ey - ey2*(l1[i][j][1].Hzx+l1[i][j][1].Hzy - l1[i-1][j][1].Hzx-l1[i-1][j][1].Hzy)/fd_dx;

	}	
    }
  }
  // ********************************************
  // *  Hx  & Hy update (zukunfts Hx,y @ n+1  *
  // ********************************************
  fdtd_comm_ghost_cells();
  for(i=1;i<local_fd_dim.x-1;++i) // iend)
  {
    iglobal = (i-1) + my_coord.x*(local_fd_dim.x-2);

//TESTCASE
//if(iglobal<srcx_int) continue;

    for(j=1;j<local_fd_dim.y-1;++j)
    {
	jglobal=(j-1) + my_coord.y*(local_fd_dim.y-2);

        mu=l1[i][j][1].mu;
        eps=l1[i][j][1].eps;
        sigmay=l1[i][j][1].sigmax;
        sigmax=l1[i][j][1].sigmay;

        sigmastarx=sigmax*mu/eps;
        sigmastary=sigmay*mu/eps;

        hy1=(mu-0.5*dt_FDTD*sigmastarx)/(mu+0.5*dt_FDTD*sigmastarx);
        hy2=(dt_FDTD)/(mu+0.5*dt_FDTD*sigmastarx);
        hx1=(mu-0.5*dt_FDTD*sigmastary)/(mu+0.5*dt_FDTD*sigmastary);
        hx2=(dt_FDTD)/(mu+0.5*dt_FDTD*sigmastary);
	
	hzx1=hy1;
	hzx2=hy2;
	hzy1=hx1;
	hzy2=hx2;

        l2[i][j][1].Hy=l1[i][j][1].Hy;
        l2[i][j][1].Hx=l1[i][j][1].Hx;

        l2[i][j][1].Hzx=l1[i][j][1].Hzx;
        l2[i][j][1].Hzy=l1[i][j][1].Hzy;

        l1[i][j][1].Hy=hy1*l1[i][j][1].Hy+hy2*(l1[i+1][j][1].Ezx-l1[i][j][1].Ezx+l1[i+1][j][1].Ezy-l1[i][j][1].Ezy)/fd_dx;
        l1[i][j][1].Hx=hx1*l1[i][j][1].Hx-hx2*(l1[i][j+1][1].Ezx-l1[i][j][1].Ezx+l1[i][j+1][1].Ezy-l1[i][j][1].Ezy)/fd_dy;

	l1[i][j][1].Hzx=hzx1*l1[i][j][1].Hzx - hzx2*(l1[i+1][j][1].Ey-l1[i][j][1].Ey)/fd_dx;
	l1[i][j][1].Hzy=hzy1*l1[i][j][1].Hzy + hzy2*(l1[i][j+1][1].Ex-l1[i][j][1].Ex)/fd_dy;


        if(l1[i][j][1].natoms<fd_min_atoms) continue; //Für vakuum ist der rest nicht noetig

        epsinf=l1[i][j][1].DL[0];
        omega_plasma=l1[i][j][1].DL[5];
        gamma_p=l1[i][j][1].DL[1];
        omegapl_L=l1[i][j][1].DL[2];
        Omega0_L=l1[i][j][1].DL[3];
        Gamma_L=l1[i][j][1].DL[4];

	double CD1=(2-dt_FDTD*gamma_p)/(2+dt_FDTD*gamma_p);
	double CL1=(2-dt_FDTD*Gamma_L)/(2+dt_FDTD*Gamma_L);
	double CD2=2*dt_FDTD/(2+dt_FDTD*gamma_p);
	double CL2=2*dt_FDTD/(2+dt_FDTD*Gamma_L);

	double omegapl_D_sq=omega_plasma*omega_plasma;
	double omegapl_L_sq=omegapl_L*omegapl_L;
	double omega0_L_sq=Omega0_L*Omega0_L;

	//J-Drude
	l1[i][j][1].Jzx=CD1*l1[i][j][1].Jzx+CD2*(omegapl_D_sq*eps0*l1[i][j][1].Ezx);
	l1[i][j][1].Jzy=CD1*l1[i][j][1].Jzy+CD2*(omegapl_D_sq*eps0*l1[i][j][1].Ezy);
	l1[i][j][1].Jx= CD1*l1[i][j][1].Jx+ CD2*(omegapl_D_sq*eps0*l1[i][j][1].Ex);
	l1[i][j][1].Jy= CD1*l1[i][j][1].Jy+ CD2*(omegapl_D_sq*eps0*l1[i][j][1].Ey);

	//J-Lorentz
	l1[i][j][1].Jlzx=CL1*l1[i][j][1].Jlzx+CL2*(omegapl_L_sq*eps0*l1[i][j][1].Ezx-omega0_L_sq*l1[i][j][1].Pzx);
	l1[i][j][1].Jlzy=CL1*l1[i][j][1].Jlzy+CL2*(omegapl_L_sq*eps0*l1[i][j][1].Ezy-omega0_L_sq*l1[i][j][1].Pzy);	
	l1[i][j][1].Jlx=CL1*l1[i][j][1].Jlx  +CL2*(omegapl_L_sq*eps0*l1[i][j][1].Ex- omega0_L_sq*l1[i][j][1].Px);
	l1[i][j][1].Jly=CL1*l1[i][j][1].Jly  +CL2*(omegapl_L_sq*eps0*l1[i][j][1].Ey- omega0_L_sq*l1[i][j][1].Py);

	//Leistungsdichte nur für aktive Zellen
	//räuml. und zeitl. Interpolieren, da sich Jx und Jy 
        //versetzt von Jz befinden und außerdem an n+3/2 sowie n+1/2 existieren
        //Da sich die Nachbarkomponente an der Stelle i-1, bzw. j-1
        //befindet, ist auch sicher, dass diese auch bereits berechnet 
        //wurde und im Speicher liegt
        //Aufgrund dieser Interpolation müssen jetzt aber auch
        //Jx und Jy sowie kommuniziert werden!


        Ez_present=l1[i][j][1].Ezx+l1[i][j][1].Ezy; //@n+1
        Jz_present=(l1[i][j][1].Jzx+l1[i][j][1].Jzy+l2[i][j][1].Jzx+l2[i][j][1].Jzy)*0.5; //@n+1
	
        double Jx_present=0.5*((l1[i][j][1].Jx+l1[i][j-1][1].Jx)*0.5+(l2[i][j][1].Jx+l2[i][j-1][1].Jx)*0.5); 
        double Jy_present=0.5*((l1[i][j][1].Jy+l1[i-1][j][1].Jy)*0.5+(l1[i][j][1].Jy+l1[i-1][j][1].Jy)*0.5);; 
                                                                  
                                                                  
        double Jlz_present=0.5*(l1[i][j][1].Jlzx+l1[i][j][1].Jlzy+l2[i][j][1].Jlzx+l2[i][j][1].Jlzy);
        double Jlx_present=0.5*((l1[i][j][1].Jlx+l1[i][j-1][1].Jlx)*0.5+(l2[i][j][1].Jlx+l2[i][j-1][1].Jlx)*0.5);
        double Jly_present=0.5*((l1[i][j][1].Jly+l1[i-1][j][1].Jly)*0.5+(l2[i][j][1].Jly+l2[i-1][j][1].Jly)*0.5);


        //Comp. power density
        double qe_Drude= gamma_p*(Jz_present*Jz_present+Jx_present*Jx_present+Jy_present*Jy_present)/eps0/omegapl_D_sq;
        //double qe_Lorentz=Gamma_L/deltaeps_L/eps0/wlsq*(Jlz_present*Jlz_present+Jlx_present*Jlx_present+Jly_present*Jly_present);
        double qe_Lorentz=Gamma_L/eps0/omegapl_L_sq*(Jlz_present*Jlz_present+Jlx_present*Jlx_present+Jly_present*Jly_present);
        qe=qe_Drude+qe_Lorentz;
        l1[i][j][1].source=6.3538562638e-26*qe;//convert to imd-units 
    }
  }
}
#endif //FDTD2D

/******************************************************************************************************************************************************************************/
int init_fdtd(void)
{
	double freq;
        int i,j;
        //double fd_dz; //<-global
	double Sc_tmp;
	Imp0=sqrt(mu0/eps0);

if(global_fd_dim.z>1)
  error("global_fd_dim.z>1 not supported by FDTD yet\n");
/*
	if(myid==0 && global_fd_dim.z>1){
	  printf("ERROR: global_fd_dim.z>1 not supported by FDTD yet.\n");
 	  error("global_fd_dim.z>1 not supported by FDTD yet\n");
	  return -1;
	}
*/
#ifdef FDTD1D
	if(myid==0 && global_fd_dim.y>1){
	  printf("ERROR: global_fd_dim.y>1 not supported by FDTD1D.\n");
	  return -1;
	}
#endif
	
	fd_dx=fd_h.x/1e10;
	fd_dy=fd_h.y/1e10;
	fd_dz=fd_h.z/1e10;
	
#ifdef FDTD2D
	Sc_tmp=1.0/sqrt(2); //CFL number for Basic 2d-case (preliminary)	
	Sc_FDTD=MIN(Sc_tmp,Sc); //user kann kleineres Sc wählen
	dt_FDTD=Sc_FDTD*MIN(fd_dx,fd_dy)/c0;
#endif
#ifdef FDTD1D
	Sc_tmp=1.0;
	Sc_FDTD=MIN(Sc_tmp,Sc); //user kann kleineres Sc wählen
	dt_FDTD=Sc_FDTD*fd_dx/c0;
#endif

	//ACHTUNG: die dispersiven gleichungen haben selbes stabilitätskriterium wie standard fdtd. Kein adaptive timestep nötig
 	freq=c0/lambda;
        omega_laser=2.0*M_PI*freq;
	Nlambdax=c0/(freq*fd_dx);
	Nlambday=c0/(freq*fd_dy);

	
 	E0=sqrt(2.0*I0*Imp0); // V/m, E-feldstärke aus max intensity via vacuum impedanz
        //weil I_0=0.5*eps0*c0*E_0^2 (Intensität ist ein Zeitmittel ! )
 	E0_init=E0; //E0_init wird benutzt für laser_active ja/nein
	
	laser_sigma_t_squared=laser_sigma_t*laser_sigma_t;
	laser_sigma_t1_squared=laser_sigma_t1*laser_sigma_t1;

        //check if proc. contains source at all (ydim is theoretically inf)
	srcx_int=round(srcx/fd_dx);
	srcw_int=round(srcw/fd_dy);

	srcx_int=MAX(srcx_int,bw+10);
        if(srcx_int<=local_fd_dim.x-1-1+my_coord.x*(local_fd_dim.x-2) &&    
           srcx_int>=(my_coord.x*(local_fd_dim.x-2)))
	{
	  containssrc=1;
//	printf("myid:%d,my_coord.x:%d,local.dim.x:%d,srcx_int:%d\n",
//		myid,my_coord.x,local_fd_dim.x,srcx_int);
	}

#ifdef FDTD2D
	srcy_int=(int)(((double)global_fd_dim.y+0.5)/2.0);
	if(myid==0){
	  printf("*****************************************************\n");
  	  printf("*         FINITE-DIFFERENCE-TIME-DOMAIN 2D          *\n");
	  printf("*****************************************************\n");
	  printf("Courant Number:%.2f\n",Sc_FDTD);
  	  printf("dt_FDTD:%.4e s,fd_dx:%.4e m,fd_dy:%.4e m\n",dt_FDTD,fd_dx,fd_dy);
	  printf("lambda:%.4e\n",lambda);
	  printf("Nlambdax:%.4e.Nlambday:%.4e\n",Nlambdax,Nlambday);
	  printf("srcx:%.4e Angs,srcy:%.4e Angs\n",srcx*1e10,srcy_int*fd_dy*1e10);
	  printf("srcx_int:%d,srcy_int:%d (in lattice spacings)\n",srcx_int,srcy_int);
          printf("srcw_int:%d (lattice spacings),srcw (Angs):%.4e\n",(int) srcw_int,srcw*1e10);
	  printf("I0:%.4e W/m^2,E0:%.4e V/m\n",I0,E0);
	  printf("t0:%.4e,sigma_t:%.4e\n",laser_t_0,laser_sigma_t);
          printf("t1:%.4e,sigma_t1:%.4e\n",laser_t_1,laser_sigma_t1);
	  printf("pml-thickness:%d,bwx:%.4e,bwy:%.4e\n",bw,bw*fd_dx,bw*fd_dy);
  	  printf("***************************\n");
	}
	//laser_spot_area=(double) 2.0*M_PI*pow(srcw_int*fd_dy,2.0);
	
	//E0 muss im 2D-Fall angepasst werden, da 2 Modi existieren mit je 2 Aufspaltungen.
	//Außerdem wird im 2D-Fall NICHT die Total Field/Scattered Field Methode genutzt, sodass E0 verdoppelt werden muss,
	//da die Welle sich aufspaltet und nicht nur in 1 Richtung fortbewegt

	//Es muss gelten sqrt( (E_zx,0 +E_zy,0)^2 + Ex,0^2+Ey_0^2 ) = E_0
	//Für den Fall dass nur TEZ exisitert folgt:      E_0'=E_zx,0=E_zy,0 = E_0/2
        //Für den Fall dass TEZ und TMZ exisiteren folgt: E_0'=E_zy,0=E_zx,0=E_zy,0=E_x,0=E_y,0= E_0/sqrt(6)
	//Vorausgesetzt jeder Modus und jede Komponente soll gleich stark angeregt werden
	E0=E0*2.0;	 //Da sich Welle in +x/-y Richtung aufspaltet
 	//E0=E0*0.5;     //Korrektur für TEZ only
 	E0=E0/sqrt(6); //Korrektur für TMZ + TEZ
	
	laser_spot_area=(double) 2.0*srcw_int*fd_dy*fd_dz; //Die Bezugsfläche zur Abschätzung der Fluenz
	if(MIN(Nlambdax,Nlambday)<20)
	{
	  if(myid==0)
	    printf("WARNING: Too few points-per-lambda in x or y-direction!Nlambdax:%f,Nlambday:%f\n",Nlambdax,Nlambday);
	}
#endif
#ifdef FDTD1D
        if(myid==0){
          printf("***************************************************\n");
          printf("*         FINITE-DIFFERENCE-TIME-DOMAIN 1D        *\n");
          printf("***************************************************\n");
          printf("Z-Polarized\n");
	  printf("Courant Number:%.2f\n",Sc_FDTD);
          printf("dt_FDTD:%.4e s,fd_dx:%.4e m,fd_dy:%.4e m\n",dt_FDTD,fd_dx,fd_dy);
          printf("lambda:%.4e\n",lambda);
          printf("Nlambdax:%.4e\n",Nlambdax);
          printf("srcx:%.4e Angs\n",srcx*1e10);
          printf("srcx_int:%d (in lattice spacings)\n",srcx_int);
          printf("I0:%.4e W/m^2,E0:%.4e V/m\n",I0,E0);
          printf("t0:%.4e,sigma_t:%.4e\n",laser_t_0,laser_sigma_t);
          printf("t1:%.4e,sigma_t1:%.4e\n",laser_t_1,laser_sigma_t1);
	  printf("pml-thickness:%d,bwx:%.4e\n",bw,bw*fd_dx);
          printf("***************************\n");
        }
        //laser_spot_area=(double) 2.0*M_PI*pow(srcw_int*fd_dy,2.0);
        laser_spot_area=(double) fd_dy*fd_dz; //Die Bezugsfläche zur Abschätzung der Fluenz
        if(Nlambdax<20)
        {
          if(myid==0)
            printf("WARNING: Too few points-per-lambda in x-direction!! Nlambdax:%f\n",Nlambdax);
        }

#endif


	for(i=0;i<local_fd_dim.x;++i)
	  for(j=0;j<local_fd_dim.y;++j)
	  {
	    l1[i][j][1].eps=eps0;	//vacuum
	    l1[i][j][1].mu=mu0;
	    l1[i][j][1].sigmax=0;
	    l1[i][j][1].sigmay=0;
  	    l1[i][j][1].Ezx=l1[i][j][1].Ezy=l1[i][j][1].Ex=l1[i][j][1].Ey=0.0;
	    l1[i][j][1].Hx=l1[i][j][1].Hy=l1[i][j][1].Hzx=l1[i][j][1].Hzy=0.0;
	    l1[i][j][1].Jzx=l1[i][j][1].Jzy=l1[i][j][1].Jx=l1[i][j][1].Jy=0.0;
	    l1[i][j][1].Pzx=l1[i][j][1].Pzy=l1[i][j][1].Px=l1[i][j][1].Py=0.0;
            //l1[i][j][1].Pzxold=l1[i][j][1].Pzyold=l1[i][j][1].Pxold=l1[i][j][1].Pyold=0.0;  
	    l1[i][j][1].Jlzx=l1[i][j][1].Jlzy=l1[i][j][1].Jlx=l1[i][j][1].Jly=0.0;
	    //l1[i][j][1].Ezold=0.0;

            l2[i][j][1].eps=eps0;       //vacuum
            l2[i][j][1].mu=mu0;
            l2[i][j][1].sigmax=0;
            l2[i][j][1].sigmay=0;
            l2[i][j][1].Ezx=l2[i][j][1].Ezy=l2[i][j][1].Ex=l2[i][j][1].Ey=0.0;
            l2[i][j][1].Hx=l2[i][j][1].Hy=l2[i][j][1].Hzx=l2[i][j][1].Hzy=0.0;
            l2[i][j][1].Jzx=l2[i][j][1].Jzy=l2[i][j][1].Jx=l2[i][j][1].Jy=0.0;
            l2[i][j][1].Pzx=l2[i][j][1].Pzy=l2[i][j][1].Px=l2[i][j][1].Py=0.0;
            //l2[i][j][1].Pzxold=l2[i][j][1].Pzyold=l2[i][j][1].Pxold=l2[i][j][1].Pyold=0.0;
	    l2[i][j][1].Jlzx=l2[i][j][1].Jlzy=l2[i][j][1].Jlx=l2[i][j][1].Jly=0.0;

//            l2[i][j][1].Ezold=0.0;
	  }

	init_pml(); 
	fdtd_create_mpi_datatypes();
	return 0;
}


/**************************************  //perfectly matched layer (split-field berenger) for absorbtion of E/H-Fields at boundary ***********************************************/
void init_pml(void)	//perfectly matched layer (split-field berenger)
{
	int i,j,k;
	int iglobal,jglobal;
//TESTCASE
	bw=MAX(bw,5); //mindestens 15 layer dick sollte schon sein

   	real refl_coeff=1.0e-12; //1e-6
	real gradingorder=8.0; //6.0;
	real bf1,bf2,bf3,bf4; //boundary factos
	real jdist,idist;      //distance from boundary in units of lattice-spacings
	
	real sigma_max_x, sigma_max_y;
	sigma_max_x= (-log10(refl_coeff)*(gradingorder+1.0)*eps0*c0)/(2.0*(real) bw*fd_dx);
	sigma_max_y= (-log10(refl_coeff)*(gradingorder+1.0)*eps0*c0)/(2.0*(real) bw*fd_dy);
	
	//Assumed vacuum @pml-layers
	bf1=(1.0*sigma_max_x)/(pow(bw,(real)gradingorder)*((real)gradingorder+1.0));
   	bf2=bf1;

#ifdef FDTD1D
  sigma_max_y=sigma_max_x;
#endif
	bf3=(1.0*sigma_max_y)/(pow(bw,(real)gradingorder)*((real)gradingorder+1.0));
	bf4=bf3;

#ifdef FDTD2D
        //sigma_x matrix	
	for(j=1;j<local_fd_dim.y-1;++j)
	{
	  jglobal =  ((j-1) + my_coord.y*(local_fd_dim.y-2));	 
	  for(i=1;i<local_fd_dim.x-1;++i)
	  {		
	    iglobal =  ((i-1) + my_coord.x*(local_fd_dim.x-2));
	    if(jglobal>=global_fd_dim.y-bw)
	    {	     
	      jdist=((real) bw-((global_fd_dim.y-1)-jglobal));
	      l1[i][j][1].sigmax=bf1*(pow((jdist+0.5),(real)(gradingorder+1.0)) -
				      pow((jdist-0.5*(jglobal>global_fd_dim.y-1-bw)),(real)(gradingorder+1.0))); 
	    } 
	    else if(jglobal<=bw)
	    {
	      jdist=(real)(bw-jglobal);
              l1[i][j][1].sigmax=bf2*(pow((jdist+0.5),(real)(gradingorder+1.0)) -
                                      pow((jdist-0.5*(jglobal<bw)),(real)(gradingorder+1.0)));
	    }
	  }
	}
	//sigma_y matrix
        for(i=1;i<local_fd_dim.x-1;++i)
        {
          iglobal =  ((i-1) + my_coord.x*(local_fd_dim.x-2));
          for(j=1;j<local_fd_dim.y-1;++j)
          {
            iglobal =  ((i-1) + my_coord.x*(local_fd_dim.x-2));
            if(iglobal>=global_fd_dim.x-bw)
            {
              idist=(real)(bw- ((global_fd_dim.x-1)-iglobal));
              l1[i][j][1].sigmay=bf3*(pow((idist+0.5),(real)(gradingorder+1.0)) -
                                      pow((idist-0.5*(iglobal>global_fd_dim.x-1-bw)),(real)(gradingorder+1.0)));	      
            }
            else if(iglobal<=bw)
            {
              idist=(real) (bw-iglobal);
              l1[i][j][1].sigmay=bf4*(pow((idist+0.5),(real)(gradingorder+1.0)) -
                                      pow((idist-0.5*(iglobal<bw)),(real)(gradingorder+1.0)));
            }
          } 
        }
#endif

#ifdef FDTD1D
	//sigma_y matrix
 	j=1;
        for(i=1;i<local_fd_dim.x-1;++i)
        {
            iglobal =  ((i-1) + my_coord.x*(local_fd_dim.x-2));
            if(iglobal>=global_fd_dim.x-bw)
            {
              idist=(real)(bw- ((global_fd_dim.x-1)-iglobal));
              l1[i][j][1].sigmay=bf3*(pow((idist+0.5),(real)(gradingorder+1.0)) -
                                      pow((idist-0.5*(iglobal>global_fd_dim.x-1-bw)),(real)(gradingorder+1.0)));	      
            }
            else if(iglobal<=bw)
            {
              idist=(real) (bw-iglobal);
              l1[i][j][1].sigmay=bf4*(pow((idist+0.5),(real)(gradingorder+1.0)) -
                                      pow((idist-0.5*(iglobal<bw)),(real)(gradingorder+1.0)));
            }
        }
#endif
}


/**************************************** MPI-DERIVED DATATYPE IN ANALOGY to imd_ttm.c ********************************************************/
int fdtd_create_mpi_datatypes(void)
{
    ttm_Element tmpelement;
    ttm_Element * tmpelement_pointer;
    tmpelement_pointer=&tmpelement;
    MPI_Aint tmpaddr;

    /* elements to be sent:        Ezx,Ezy,Hx,Hy,
				   Hzx,Hzy,Ex,Ey
     *                             (MPI_UB to set upper bound and skip the rest) */


    int blockcounts[13]={1,1,1,1,1,1,1,1,1,1,1,1,1};
    MPI_Datatype types[13]={MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
			   MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
			   MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
			   MPI_UB};
    MPI_Aint displs[13];

    MPI_Address(tmpelement_pointer, &tmpaddr);
    MPI_Address(&tmpelement.Ezx, &displs[0]);
    MPI_Address(&tmpelement.Ezy, &displs[1]);
    MPI_Address(&tmpelement.Hx, &displs[2]);
    MPI_Address(&tmpelement.Hy, &displs[3]);

    MPI_Address(&tmpelement.Hzx, &displs[4]);
    MPI_Address(&tmpelement.Hzy, &displs[5]);
    MPI_Address(&tmpelement.Ex, &displs[6]);
    MPI_Address(&tmpelement.Ey, &displs[7]);

    MPI_Address(&tmpelement.Hzx, &displs[4]);
    MPI_Address(&tmpelement.Hzy, &displs[5]);
    MPI_Address(&tmpelement.Ex, &displs[6]);
    MPI_Address(&tmpelement.Ey, &displs[7]);

    MPI_Address(&tmpelement.Jx, &displs[8]);
    MPI_Address(&tmpelement.Jy, &displs[9]);
    MPI_Address(&tmpelement.Jlx, &displs[10]);
    MPI_Address(&tmpelement.Jly, &displs[11]);

    tmpelement_pointer++;
    MPI_Address(tmpelement_pointer, &displs[12]);

    displs[12]-=tmpaddr;
    displs[11]-=tmpaddr;
    displs[10]-=tmpaddr;
    displs[9]-=tmpaddr;

    displs[8]-=tmpaddr;
    displs[7]-=tmpaddr;
    displs[6]-=tmpaddr;
    displs[5]-=tmpaddr;

    displs[4]-=tmpaddr;
    displs[3]-=tmpaddr;
    displs[2]-=tmpaddr;
    displs[1]-=tmpaddr;
    displs[0]-=tmpaddr;

    MPI_Type_struct(13,blockcounts,displs,types,&mpi_element3);
    MPI_Type_commit(&mpi_element3);

 /* datatype for one string of elements along z (short of 2 lattice points) */
  MPI_Type_contiguous(local_fd_dim.z-2, mpi_element3, &mpi_zrow3);
  MPI_Type_commit(&mpi_zrow3);

  { /* add displacements to skip ghost layers (mpi_zrow_block) */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore=&tmpel;

    int blockcounts[3]={1,1,1};
    MPI_Datatype types[3]={MPI_LB, mpi_zrow3, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore,&displs[0]);
    tmpstore += 1;
    MPI_Address(tmpstore,&displs[1]);
    tmpstore += 1 + (local_fd_dim.z-2);
    MPI_Address(tmpstore,&displs[2]);

    displs[2]-=displs[0];
    displs[1]-=displs[0];
    displs[0]-=displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_zrow_block3);
    MPI_Type_commit(&mpi_zrow_block3);
  }


  /* datatype for one layer of elements perpendicular to x
   * (short of 2 strings along z-axis)                    */
  MPI_Type_contiguous(local_fd_dim.y-2, mpi_zrow_block3, &mpi_xplane3);
  MPI_Type_commit(&mpi_xplane3);
  { /* add displacements to skip ghost layers (mpi_xplane_block) */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore=&tmpel;

    int blockcounts[3]={1,1,1};
    MPI_Datatype types[3]={MPI_LB, mpi_xplane3, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore,&displs[0]);
    tmpstore += 2 + (local_fd_dim.z-2);
    MPI_Address(tmpstore,&displs[1]);
    tmpstore += (2 + (local_fd_dim.z-2))*((local_fd_dim.y-2) + 1);
    MPI_Address(tmpstore,&displs[2]);

    displs[2]-=displs[0];
    displs[1]-=displs[0];
    displs[0]-=displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_xplane_block3);
    MPI_Type_commit(&mpi_xplane_block3);
  }

  /* datatype for one layer of elements perpendicular to y (short of 2) */
  MPI_Type_vector((local_fd_dim.x-2), 1, 2+(local_fd_dim.y-2),
                  mpi_zrow_block3, &mpi_yplane3);
  MPI_Type_commit(&mpi_yplane3);

  { /* add displacements to skip ghost layers (mpi_yplane_block) */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore=&tmpel;

    int blockcounts[3]={1,1,1};
    MPI_Datatype types[3]={MPI_LB, mpi_yplane3, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore,&displs[0]);
    tmpstore += (2 + (local_fd_dim.z-2))*
                (2 + (local_fd_dim.y-2));
    MPI_Address(tmpstore,&displs[1]);
    tmpstore += (1 + (local_fd_dim.x-2))*
                (2 + (local_fd_dim.z-2))*
                (2 + (local_fd_dim.y-2));
    MPI_Address(tmpstore,&displs[2]);

    displs[2]-= displs[0];
    displs[1]-= displs[0];
    displs[0]-= displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_yplane_block3);
    MPI_Type_commit(&mpi_yplane_block3);
  }


  /* datatype for one string of elements along y (short of 2) */
  MPI_Type_vector((local_fd_dim.y-2), 1, 2+(local_fd_dim.z-2),
                  mpi_element3, &mpi_yrow3);
  MPI_Type_commit(&mpi_yrow3);
  /* add displacements to create datatype from which mpi_zplane will be built
   * (again, skipping ghost layers) */
  { /* mpi_yrow_block */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore=&tmpel;

    int blockcounts[3] = {1,1,1};
    MPI_Datatype types[3]={MPI_LB, mpi_yrow3, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore,&displs[0]);
    tmpstore += 2 + (local_fd_dim.z-2);
    MPI_Address(tmpstore,&displs[1]);
    tmpstore += (1 + (local_fd_dim.y-2)) *
                (2 + (local_fd_dim.z-2));
    MPI_Address(tmpstore,&displs[2]);

    displs[2]-=displs[0];
    displs[1]-=displs[0];
    displs[0]-=displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_yrow_block3);
    MPI_Type_commit(&mpi_yrow_block3);
  }
  /* datatype for one layer of elements perpendicular to z
   * (short of 2 strings) */
  MPI_Type_contiguous((local_fd_dim.x-2), mpi_yrow_block3, &mpi_zplane3);
  MPI_Type_commit(&mpi_zplane3);

  { /* add displacements to skip ghost layers (mpi_zplane_block) */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore=&tmpel;

    int blockcounts[3]={1,1,1};
    MPI_Datatype types[3]={MPI_LB, mpi_zplane3, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore,&displs[0]);
    tmpstore += (2 + (local_fd_dim.z-2)) *
                (2 + (local_fd_dim.y-2));
    MPI_Address(tmpstore,&displs[1]);
    tmpstore += (1 + (local_fd_dim.x-2)) *
                (2 + (local_fd_dim.z-2)) *
                (2 + (local_fd_dim.y-2));
    MPI_Address(tmpstore,&displs[2]);

    displs[2]-=displs[0];
    displs[1]-=displs[0];
    displs[0]-=displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_zplane_block3);
    MPI_Type_commit(&mpi_zplane_block3);
  }


   return 0;
}

/************************************************  MPI COMMUNICATION IN ANALOGY TO imd_ttm.c ********************************************************/
void fdtd_comm_ghost_cells(void)
{
  /** MPI communication (can occur before and/or during MD calculations?) */
  /* Remember:
   * east -> -x
   * west -> +x
   * north-> -y
   * south-> +y
   * up   -> -z
   * down -> +z
   * *************/
  int i,j,k;

  /* x direction */
  if(pbc_dirs.x==1 || (my_coord.x != 0 && my_coord.x != cpu_dim.x-1) )
  {
    /* send left slice to left neighbor. */
    /* Simultaneously receive slice from right neighbor */
    MPI_Sendrecv(&l1[1][0][0],1,mpi_xplane_block3,nbeast,9100,
        &l1[(local_fd_dim.x-2)+1][0][0],1,mpi_xplane_block3,nbwest,9100,
        cpugrid,&stati2[0]);


    /* send right slice to right neighbor. */
    /* Simultaneously receive slice from left neighbor */
    MPI_Sendrecv(&l1[(local_fd_dim.x-2)][0][0],1,mpi_xplane_block3,nbwest,9200,
        &l1[0][0][0],1,mpi_xplane_block3,nbeast,9200,
        cpugrid,&stati2[1]);
  }
  else /* no pbc and we are at the surface */
  {
    if (my_coord.x==0 && my_coord.x!=cpu_dim.x-1) /* left surface */
    {
      /* only receive from right */
      MPI_Recv(&l1[(local_fd_dim.x-2)+1][0][0],1,mpi_xplane_block3,nbwest,9100,
          cpugrid,&stati2[0]);

      /* only send to right */
      MPI_Send(&l1[(local_fd_dim.x-2)][0][0],1,mpi_xplane_block3,nbwest,9200,
          cpugrid);

      for (j=1;j<=(local_fd_dim.y-2);j++)
      {
          l1[0][j][1].Hx=0;
          l1[0][j][1].Hy=0;
          l1[0][j][1].Ezx=0;
          l1[0][j][1].Ezy=0;
	  l1[0][j][1].Ex=l1[0][j][1].Ey=0;
	  l1[0][j][1].Hzx=l1[0][j][1].Hzy=0;
	 
//	  l1[0][j][1].Ezold=0.0;
      }

    }
    else if (my_coord.x==cpu_dim.x-1 && my_coord.x!=0) /* right surface */
    {
      /* only send to left */
      MPI_Send(&l1[1][0][0],1,mpi_xplane_block3,nbeast,9100,
          cpugrid);

      /* only receive from left */
      MPI_Recv(&l1[0][0][0],1,mpi_xplane_block3,nbeast,9200,
          cpugrid,&stati2[1]);

      for (j=1;j<=(local_fd_dim.y-2);j++)
      {
          l1[local_fd_dim.x-1][j][1].Hx=0;
          l1[local_fd_dim.x-1][j][1].Hy=0;
          l1[local_fd_dim.x-1][j][1].Ezx=0;
          l1[local_fd_dim.x-1][j][1].Ezy=0;
      }

    } 
  }

  /* y direction */
if(global_fd_dim.y>1)
{
  if(pbc_dirs.y==1 || (my_coord.y != 0 && my_coord.y != cpu_dim.y-1) )
  {
    /* send left slice to left neighbor. */
    /* Simultaneously receive slice from right neighbor */
    MPI_Sendrecv(&l1[0][1][0],1,mpi_yplane_block3,nbnorth,910,
        &l1[0][(local_fd_dim.y-2)+1][0],1,mpi_yplane_block3,nbsouth,910,
        cpugrid,&stati[2]);
      /* send right slice to right neighbor. */
      /* Simultaneously receive slice from left neighbor */
    MPI_Sendrecv(&l1[0][(local_fd_dim.y-2)][0],1,mpi_yplane_block3,nbsouth,920,
        &l1[0][0][0],1,mpi_yplane_block3,nbnorth,920,
        cpugrid,&stati2[3]);
  }
  else /* no pbc and we are at the surface */
  {
    if (my_coord.y==0 && my_coord.y!=cpu_dim.y-1) /* left surface */
    {
        MPI_Recv(&l1[0][(local_fd_dim.y-2)+1][0],1,mpi_yplane_block3,nbsouth,910,
            cpugrid,&stati2[2]);
        MPI_Send(&l1[0][(local_fd_dim.y-2)][0],1,mpi_yplane_block3,nbsouth,920,
            cpugrid);

        for (i=1;i<=(local_fd_dim.x-2);i++)
        {
          l1[i][0][1].Hx=0;
          l1[i][0][1].Hy=0;
          l1[i][0][1].Ezx=0;
          l1[i][0][1].Ezy=0;
        }
    }
    else if (my_coord.y==cpu_dim.y-1 && my_coord.y!=0) /* right surface */
    {
        MPI_Send(&l1[0][1][0],1,mpi_yplane_block3,nbnorth,910,
            cpugrid);
        MPI_Recv(&l1[0][0][0],1,mpi_yplane_block3,nbnorth,920,
            cpugrid,&stati2[3]);

        for (i=1;i<=(local_fd_dim.x-2);i++)
        {
          l1[i][local_fd_dim.y-1][1].Hx=0;
          l1[i][local_fd_dim.y-1][1].Hy=0;
          l1[i][local_fd_dim.y-1][1].Ezx=0;
          l1[i][local_fd_dim.y-1][1].Ezy=0;
        }
    } 
  }
}
  //No communication in z-dir
}




/************************************ MEDIUM PROPERTIES *************************************************************/
double getgammap(int i,int j)
{
  //return 5e13;

  return 1e14; //TESTCASE
  double nueff=NueffInterpol(l1[i][j][1].dens,l1[i][j][1].temp,l1[i][j][1].md_temp);
  if(nueff==-1)
  {
      printf("steps:%d,proc:%d,i:%d,j:%d,k:1, ERROR during KappaInterpol in FILLMESH, Te:%f (eV), Ti:%f (eV) ,dens:%f (kg/m^3),atoms:%d\n",
	      steps,myid,i,j,l1[i][j][1].temp,l1[i][j][1].md_temp,l1[i][j][1].dens,l1[i][j][1].natoms);
      error("ERROR during NueffInterpol in getgammap (imd_fdtd.c)");
 
  }
  return nueff;
  // ~ 5e13
}
double getomegap(int i,int j) //plasma-freq
{
  //ACHTUNG: omega_p = cont. nur fuer testfall!! TESTCASE

  //return sqrt(l1[i][j][1].ne*ECHARGE*ECHARGE/EMASS/eps0);
  //return sqrt(l1[i][j][1].ne*3.187015417814758e+03);
  return 2.160294047429742e+16;
  //return 2.5243e+16;
}


