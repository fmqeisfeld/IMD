#define RHOMIN 51


#define ADVMODE 2  // 0=NO ADVECTION, 1=GODUNOV SOLVER VIA VCOM, 2=DISCRETE FLUX SOLVER (PREDICT ATOMIC FLUXES)
//#define ADV2D

#define DEBUG_LEVEL 1

#include <sys/time.h>
#include "imd.h"
#include <complex.h>
//#include "tricub_coeffmat.h"  //A[64][64] matrix fuer tricubic interpolation

#ifdef DEBUG
#include <assert.h>
#endif

#ifdef BUFCELLS
#define NBUFFC 2
#else
#define NBUFFC 0
#endif /*BUFCELLS*/

#define TTMOUTBUFLEN 400
// ******************************************************************************/
/* Macro to allow electr. heat capacity to vary
 * with electronic temperature (C_e=gamma*T_e) */
// EIN RELIKT AUS DER VERGANGENHEIT
// gut fuer den festkörper aber unbrauchbar fuer schmelze & gas
// Hier ist sogar die Freie-Elektronen-Näherung besser (Fermi-integral berechnen)
// Hierfür gibt es weiter unten im Code die entsprechende wide-range routine
#define FD_C ((fd_c==0)?(fd_gamma*l1[i][j][k].temp):(fd_c))

#define node l1[i][j][k]  //bequemer, muss oft geschrieben werden
// ****************************************************
// *         MAIN FUNC
// *****************************************************
void calc_ttm()
{
  int i, j, k;
  double fd_subtime;
  double tau_DIFF; //integr. timestep-size in IMD-time units
  double tau_FDTD; //integr. timestep-size in SI-units

  int fdtd_substeps;
  //int diff_substeps;
  int fdtd_div_diff; // fdtd_substeps/diff_substeps
  // ACHTUNG: integr.zeitschritt fuer Diffusion ist maximal so groß
  // wie integr.zeitschritt fuer maxwell-solver
  // Außerdem kann der benutzer via fd_n_timesteps eine untere grenze bestimmen

  int telaps_md;     //fuer timing
  int telaps_ttm;

  fd_n_timesteps=100;
  diff_substeps=fd_n_timesteps;

  tau_DIFF=(timestep)/((double) fd_n_timesteps);

  update_fd();
  //do_ADV(1.0);
  do_cell_activation();
  do_FILLMESH();
  ttm_fill_ghost_layers();  

    for (i = 0; i < diff_substeps; i++)
    {
      //do_FILLMESH();
      do_tmm(tau_DIFF);
      tmm_time += tau_DIFF * 10.18 * 1.0e-15; //in sek
      do_DIFF(tau_DIFF);
      do_FILLMESH();
      ttm_fill_ghost_layers();            
    }

//Calc internal eng and updt new U after diff
tot_elec_energy_local =0;
for (i=1; i<local_fd_dim.x-1; ++i)
{
for (j=1; j<local_fd_dim.y-1; ++j)
{
for (k=1; k<local_fd_dim.z-1; ++k)
{
      if(l1[i][j][k].natoms>=1)
      {
        //l1[i][j][k].U =EOS_ee_from_r_te(l1[i][j][k].dens, l1[i][j][k].temp * 11604.5) * 26.9815 * AMU * J2eV; //eV/Atom
        tot_elec_energy_local += l1[i][j][k].U*((double) l1[i][j][k].natoms);
      }      
      else
      {
        l1[i][j][k].U=0.0; //l2[i][j][k].U=0.0;
      }
}
}
}

MPI_Reduce(&Eabs_local, &Eabs_global, 1, MPI_DOUBLE, MPI_SUM, 0, cpugrid);

if(myid==0)
  printf("step:%d, Finc:%.4e, t-t0:%.4e, Refl:%.4e \n",steps,Eabs_global * eV2J / laser_spot_area,(tmm_time - laser_t_0) * 1e15,tmm_refl);



/*  
#if DEBUG_LEVEL>1
  printf("proc:%d,steps:%d,entered calc_ttm\n", myid);
#endif


  if (myid == 0) // In der Entwicklungsphase macht es Sinn IMMER die Performance im Auge zu behalten
  { // Hier wird er Timer für den MD-Teil ausgewertet
    gettimeofday(&eetime, NULL);
    telaps_md = ((eetime.tv_sec - sstime.tv_sec) * 1000000) + (eetime.tv_usec - sstime.tv_usec);
  }

#ifdef FDTD
  //dont waste comp.time on maxwell solver, if laser is inactive
  laser_timefun = E0 * exp(-0.5 * pow(t_SI + dt_FDTD - laser_t_0, 2.0) / laser_sigma_t_squared);

  if (laser_t_1 > 0.0)
    laser_timefun = laser_timefun + E0 * exp(-0.5 * pow(t_SI + dt_FDTD - laser_t_1, 2.0) / laser_sigma_t1_squared);

  if (laser_timefun < E0_init * 1e-4) // switch off source if intensity very low //TESTCASE
  {
    laser_active = false;
  }
  else laser_active = true;


  if (laser_active)
  {
    fdtd_substeps = (int) ((timestep * 10.18 * 1.0e-15) / (dt_FDTD) + 1.0); //dt_FDTD  ist in sek
    diff_substeps = MAX(fdtd_substeps, (int) (timestep / max_dt_ttm)); //max_dt_TTM ist in IMD-zeiteinheit (ACHTUNG:=inf bei mdstep=0)

    diff_substeps = MAX(diff_substeps, fd_n_timesteps);

    diff_substeps = (int) (((double)diff_substeps) / ((double) fdtd_substeps)); //weil es in diesem fall eine äußere fdtd-loop gibt
    diff_substeps = MAX(diff_substeps, 1); //damit wenigstens 1 einziger step gemacht wird

    tau_FDTD = timestep * 10.18 * 1.0e-15 / ((double) fdtd_substeps);
    tau_DIFF = tau_FDTD / ((double) diff_substeps) / (10.18 * 1.0e-15);
  }
  else //Kein laser aktiv --> keine äußere fdtd-loop noetig
  {
    fdtd_substeps = MAX(fd_n_timesteps, (int) (timestep / max_dt_ttm));
    diff_substeps = 1;
    tau_FDTD = timestep * 10.18 * 1.0e-15 / ((double) fdtd_substeps);
    tau_DIFF = tau_FDTD / (10.18 * 1.0e-15);
  }
#else // ifdef FDTD
  diff_substeps = MAX(fd_n_timesteps, (int) (timestep / max_dt_ttm));
  tau_DIFF = timestep / ((double) diff_substeps);
  //laser active? (ohne maxwell solver) <-- in rescale_ttm
#endif // ifdef FDTD

  // ************************************ UNABHÄNGIG VON LASER-ENERGIE ABSORPTIONS-METHODE ******************************************
  update_fd();   //calc md-temp,fluxes and clear old arrays

  if (steps > 0)
  {
    do_ADV(timestep);      //Advection-step fuer die elektronen-temperatur
    ttm_fill_ghost_layers(); //new Te?
    do_FILLMESH(); //calc electronic optical and transport properties
    ttm_fill_ghost_layers(); //new kappa?
  }
  CFL_maxdt(); //<--koennte sich geändert haben durch advection

  ///////////////////////////////////////////////////
  // FINITE-DIFFERENCE-subloop  mit MAXWELL-SOLVER //
  // In diesem Fall gibts noch eine äußere Loop    //
  // d.h. diff_substeps muss bezogen werden auf    //
  // maxwell_substeps            //
  ///////////////////////////////////////////////////
  if (myid == 0)
  {
    gettimeofday(&sstime, NULL); //start timer for FD-stuff
    //printf("fdtd_substeps:%d,lact:%d\n",fdtd_substeps,laser_active);
  }
  if (steps > 0)
  {
    if (laser_active)
      MPI_Reduce(&Eabs_local, &Eabs_global, 1, MPI_DOUBLE, MPI_SUM, 0, cpugrid);
#ifdef FDTD
    for (i = 0; i < fdtd_substeps; i++)
    {
      if (laser_active) do_fdtd(t_SI); //braucht die aktuelle sim-Zeit in sek.
//ttm_writeout(i); //DEBUG
//if(myid==0) printf("substeps:%d\n",i);
      for (k = 1; k <= diff_substeps; k++)
      {
        do_DIFF(tau_DIFF);
        ttm_fill_ghost_layers();
        do_FILLMESH();    //ACHTUNG:Eigentlich müsste nach FILLMESH
        //Wieder kommuniziert werden, da kappa geupdatet ist
        //d.h. in der naechsten iteration ist für 1 FD-step
        //kappa falsch. Nehme an,dass der Unterschied eher gering ist
      }
      t_SI += tau_FDTD; //in Si-units
    }
// *****************************************
// * PERFORMANCE MESSEN und SCREEN-OUTPUT
// ******************************************
    if (myid == 0)
    {
      gettimeofday(&eetime, NULL);
      telaps_ttm = ((eetime.tv_sec - sstime.tv_sec) * 1000000) + (eetime.tv_usec - sstime.tv_usec);

      printf("mdsteps:%d,laser_act:%d,FDTD-steps:%d, DIFF-steps:%d, t-t0:%.4e fs,E(t):%.4e V/m,t_SI:%.4e (s), elapsMD:%d,elapsTTM:%d,Eabs:%.4e (J/m^2)\n",
             steps, (int) laser_active, fdtd_substeps, diff_substeps, (t_SI - laser_t_0) * 1e15, laser_timefun, t_SI,
             telaps_md, telaps_ttm, Eabs_global * 1.6021766e-19 / laser_spot_area);

      gettimeofday(&sstime, NULL); //start timer for all non-ttm stuff
    }

#else //ifdef FDTD
    //////////////////////////////////////////////////
    //FINITE-DIFFERENCE-subloop  ohne MAXWELL-SOLVER//
    //////////////////////////////////////////////////
    for (i = 0; i < diff_substeps; i++)
    {
#ifdef LASER
      laser_rescale_ttm();
#endif
#ifdef TMM
      do_tmm(tau_DIFF);
      tmm_time += tau_DIFF * 10.18 * 1.0e-15; //in sek
#endif
      do_DIFF(tau_DIFF);
      ttm_fill_ghost_layers();
      do_FILLMESH();
    }
// ************************************************************
// PERFORMANCE MESSEN und SCREEN OUTPUT (OHNE MAXWELL SOLVER)
// ************************************************************
//MPI_Reduce(&tot_elec_energy_local,&tot_elec_energy_global,1,MPI_DOUBLE,MPI_SUM,0,cpugrid); //DEBUG PURPOSE

    if (myid == 0)
    {
      gettimeofday(&eetime, NULL);
      telaps_ttm = ((eetime.tv_sec - sstime.tv_sec) * 1000000) + (eetime.tv_usec - sstime.tv_usec);

#ifdef LASER
      printf("mdsteps:%d,laser:%d,DIFF-steps:%d, elapsMD:%d,elapsTTM:%d,Eabs:%.2e\n",
             steps, laser_active, diff_substeps, telaps_md, telaps_ttm, Eabs_global * 1.6021766e-19 / laser_spot_area);
#endif
#ifdef TMM
      printf("mdsteps:%d,laser:%d,DIFF-steps:%d, elapsMD:%d,elapsTTM:%d,t-t0:%.2f,refl:%.2e,trans:%.2e,abs:%.2e,Eabs:%.4e J/m^2\n",
             steps, laser_active, diff_substeps, telaps_md, telaps_ttm, (tmm_time - laser_t_0) * 1e15,
             tmm_refl, tmm_trans, tmm_abs, Eabs_global * 1.6021766e-19 / laser_spot_area);
#endif

      gettimeofday(&sstime, NULL); //start timer for all non-ttm stuff
    }
#endif  //#ifdef FDTD --> else
  } //if step>0

  */
}

// *********************************************************************************************************************************************
// * UPDATE_FD COMPUTES MD-TEMPERATURES AND ATOMIC FLUXES ACROSS FD-CELLS  AND CLEARS OLD ARRAYS
// *********************************************************************************************************************************************
void update_fd()
{
  natoms_local = 0;
  int i, j, k, i_global, j_global, k_global;
  //int fd_tags=(local_fd_dim.x-2)*(local_fd_dim.y-2)*(local_fd_dim.z-2)+1; //every MPI-message needs to be unique
  int loopvar, l;

  int tot_neighs; //for density calc.

  loopvar = 0;
  l = 0;
  cellptr p;
  real tot_mass = 0;
#if DEBUG_LEVEL>1
  printf("steps:%d,proc:%d,entered update_fd\n", steps, myid);
#endif

  //Clear dirichlet-maxx-array
#ifdef DIRICHLET
  for (j = 0; j < global_fd_dim.y; ++j)
    dirichlet_maxx_local[j] = -1;
#endif

#if ADVMODE==2
// *********************************************************
// Clear edges & corners of fluxes &temp (erstmal nur 2D)  *
// and precalc some params
// *********************************************************
  for (k = 0; k < 8; k++)
  {
    l1[0][0][1].flux[k] = l1[local_fd_dim.x - 1][local_fd_dim.y - 1][1].flux[k] = 0;
    l1[local_fd_dim.x - 1][0][1].flux[k] = l1[0][local_fd_dim.y - 1][1].flux[k] = 0;
  }
  l1[0][0][1].temp = l1[local_fd_dim.x - 1][local_fd_dim.y - 1][1].temp = 0;
  l1[local_fd_dim.x - 1][0][1].temp = l1[0][local_fd_dim.y - 1][1].temp = 0;

#endif

  // **********************************************************************
  // * FIRST LOOP: Clear arrays,get md-temps and detect
  // * activation/deactivation of FD-Cells
  // ************************************************************************
  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {
    i_global =  ((i - 1) + my_coord.x * (local_fd_dim.x - 2));

    //Clear dirichlet arrays
#ifdef DIRICHLET
    dirichlet_maxy_local[i_global] = -999999;
    dirichlet_miny_local[i_global] = 99999;
#endif

    for (j = 1; j < local_fd_dim.y - 1; ++j)
    {
      j_global = ((j - 1) + my_coord.y * (local_fd_dim.y - 2));
      for (k = 1; k < local_fd_dim.z - 1; ++k)
      {
        k_global =  ((k - 1) + my_coord.z * (local_fd_dim.z - 2));

        tot_neighs = 0;
        //CLEAR ARRAYS
        l1[i][j][k].vcomx = 0.0;
        l1[i][j][k].vcomy = 0.0;
        l1[i][j][k].vcomz = 0.0;

        l1[i][j][k].source = l2[i][j][k].source = 0.0;
#if ADVMODE==2
        l1[i][j][k].flux[0] = l1[i][j][k].flux[1] = l1[i][j][k].flux[2] = l1[i][j][k].flux[3] = 0;
        l1[i][j][k].flux[4] = l1[i][j][k].flux[5] = l1[i][j][k].flux[6] = l1[i][j][k].flux[7] = 0;
#endif

        l1[i][j][k].natoms_old = l2[i][j][k].natoms_old = l1[i][j][k].natoms;
l1[i][j][k].natoms = l2[i][j][k].natoms = 0; // <-- nach imd_forces_nbl.c verschoben
        tot_mass = 0.0;

        l1[i][j][k].xi = l2[i][j][k].xi = 0.0;
        l1[i][j][k].source = l2[i][j][k].source = 0.0;

        /* loop over encompassed MD cells */
        for (loopvar = 0; loopvar < n_md_cells; loopvar++)
        {
          p = l1[i][j][k].md_cellptrs[loopvar];

          /* add number of atoms of this MD cell*/
l1[i][j][k].natoms += p->n; // <-- nach imd_forces_nbl.c verschoben
          natoms_local += p->n;

          for (l = 0; l < p->n; l++) //loop over atoms
          {
            tot_neighs += NUMNEIGHS(p, l);


            //NUMNEIGHS(p, l) = 0; // nach jedem step clearen.wird in imd_forces_nbl.c wieder aufsummiert // <-- nach imd_forces_nbl.c verschoben
            l1[i][j][k].vcomx += IMPULS(p, l, X);
            l1[i][j][k].vcomy += IMPULS(p, l, Y);
            l1[i][j][k].vcomz += IMPULS(p, l, Z);
            tot_mass += MASSE(p, l);
// ***********************************************
// ATOMIC FLUXES FOR ADVECTION, erstmal nur 2D
// *************************************************
#if ADVMODE==2

//REMINDER
// flux[0] : teilchen erhalten von +x,y
// flux[1] : teilchen erhalten von -x,y
// flux[2] : teilchen erhalten von x,+y
// flux[3] : teilchen erhalten von +x,-y
// flux[4] : teilchen erhalten von +x,+y
// flux[5] : teilchen erhalten von x,-y
// flux[6] : teilchen erhalten von -x,+y
// flux[7] : teilchen erhalten von -x,-y
            if (steps > 0)
              if (p->fdi[l] != i_global || p->fdj[l] != j_global)
              {
                // +x
                if (p->fdi[l] > i_global && p->fdj[l] > j_global)
                  l1[i][j][k].flux[4]++; //aus +x,+y
                else if (p->fdi[l] > i_global && p->fdj[l] < j_global)
                  l1[i][j][k].flux[3]++;  // aus +x,-y
                else if (p->fdi[l] > i_global && p->fdj[l] == j_global)
                  l1[i][j][k].flux[0]++; // aus +x,y=y
                //-x
                else if (p->fdi[l] < i_global && p->fdj[l] < j_global) // aus -x,-y
                  l1[i][j][k].flux[7]++;
                else if (p->fdi[l] < i_global && p->fdj[l] > j_global)  // aus -x,+y
                  l1[i][j][k].flux[6]++;
                else if (p->fdi[l] < i_global && p->fdj[l] == j_global) // aus -x,y=y
                  l1[i][j][k].flux[1]++;
                //x=x
                else if (p->fdi[l] == i_global && p->fdj[l] > j_global) //aus x=x,+y
                  l1[i][j][k].flux[2]++;
                else if (p->fdi[l] == i_global && p->fdj[l] < j_global) //aus x=x,-y
                  l1[i][j][k].flux[5]++;
              }
            p->fdi[l] = i_global;
            p->fdj[l] = j_global;
            p->fdk[l] = k_global;
#endif
          } //loop through atoms
        } //loop through md-cells

        int foo;
        for (foo = 0; foo < 8; foo++)
          l2[i][j][k].flux[foo] = l1[i][j][k].flux[foo];

        l2[i][j][k].natoms = l1[i][j][k].natoms;
        //COMPUTE DENSITY OF TTM-CELLS FROM NEIGHBORS PER ATOM in kg/m^3
        if (l1[i][j][k].natoms > 0)
        {
          l1[i][j][k].dens = l2[i][j][k].dens = (double) tot_neighs / ((double)l1[i][j][k].natoms) * atomic_weight / neighvol * 1660.53907; //kg/m^3


if (l1[i][j][k].dens == 0) //in step 0...noch keine neigh list?
          {
            l1[i][j][k].dens = l2[i][j][k].dens = (double) l1[i][j][k].natoms * atomic_weight / fd_vol * 1660.53907;
          }



l1[i][j][k].dens= l2[i][j][k].dens = (double) l1[i][j][k].natoms * atomic_weight / fd_vol * 1660.53907;       


          l1[i][j][k].vcomx /= tot_mass;
          l1[i][j][k].vcomy /= tot_mass;
          l1[i][j][k].vcomz /= tot_mass;
          // ****************
          // *   DIRICHLET  *
          // ****************
#ifdef DIRICHLET
          //find outermost cells for each iglobal-row
          //es muss eine zelle sein, die auch in diff-loop als aktiv betrachtet wird!
          if (l1[i][j][k].natoms >= fd_min_atoms)
          {
            if (j_global > dirichlet_maxy_local[i_global]) dirichlet_maxy_local[i_global] = j_global;
            if (j_global < dirichlet_miny_local[i_global]) dirichlet_miny_local[i_global] = j_global;
            if (i_global > dirichlet_maxx_local[j_global]) dirichlet_maxx_local[j_global] = i_global;
          }
#endif
        } //end if natoms > 0
        else
        {
          l1[i][j][k].dens = l2[i][j][k].dens = 0.0;
          l1[i][j][k].vcomx = 0.0;
          l1[i][j][k].vcomy = 0.0;
          l1[i][j][k].vcomz = 0.0;
        }

        l2[i][j][k].vcomx = l1[i][j][k].vcomx;
        l2[i][j][k].vcomy = l1[i][j][k].vcomy;
        l2[i][j][k].vcomz = l1[i][j][k].vcomz;

        l1[i][j][k].md_temp = l2[i][j][k].md_temp = 0;
        // **********************************************************************************************************
        // * 2nd loop over MD cells and atoms For MD-temperature from kinetic energies of particles                 *        
        // **********************************************************************************************************
        for (loopvar = 0; loopvar < n_md_cells; loopvar++) //loop over md-cells
        {
          p = l1[i][j][k].md_cellptrs[loopvar];
          for (l = 0; l < p->n; l++) //loop over atoms
          {
            l1[i][j][k].md_temp += MASSE(p, l) * SQR(IMPULS(p, l, X) / MASSE(p, l)
                                   - l1[i][j][k].vcomx);
            l1[i][j][k].md_temp += MASSE(p, l) * SQR(IMPULS(p, l, Y) / MASSE(p, l)
                                   - l1[i][j][k].vcomy);
            l1[i][j][k].md_temp += MASSE(p, l) * SQR(IMPULS(p, l, Z) / MASSE(p, l)
                                   - l1[i][j][k].vcomz);
          }
        }


        if (l1[i][j][k].natoms >= fd_min_atoms)
        {


l1[i][j][k].md_temp /= 3.0 * l1[i][j][k].natoms;          
//l1[i][j][k].md_temp/=3 * l1[i][j][k].dens*fd_vol*1e-30/26.9815/AMU;          //SMOOTH TEMP?



        }
        else l1[i][j][k].md_temp = 0.0;
        l2[i][j][k].md_temp = l1[i][j][k].md_temp;

#if DEBUG_LEVEL>0
        if (l1[i][j][k].natoms >= fd_min_atoms)
        {
          if (isnan(l1[i][j][k].md_temp) != 0 || l1[i][j][k].md_temp <= 0.0) //warum steps>0? --> fuer "on the fly"-probe (imd_generate),
          { //hat do_maxwell noch keine temp zugewiesen!
            if (steps > 0)
            {
              printf("steps:%d,proc:%d,i:%d,j:%d,k:%d,mdtemp:%f,atoms:%d\n", steps, myid, i, j, k, l1[i][j][k].md_temp, l1[i][j][k].natoms);
              error("md_temp is NaN or <=0");
            }
          }
        }
#endif
        // ******************************************************************
        // *  INITIALIZE ELECTRON TEMPERATURE AND CHECK EOS PLAUSIBILITY
        // *****************************************************************
        if (steps < 1)
        {
          l1[i][j][k].natoms_old = l2[i][j][k].natoms_old = l1[i][j][k].natoms;
          if (l1[i][j][k].natoms >= fd_min_atoms)
          {
            l1[i][j][k].temp = l2[i][j][k].temp = l1[i][j][k].md_temp;

l1[i][j][k].U =EOS_ee_from_r_te(l1[i][j][k].dens, l1[i][j][k].temp * 11604.5) * 26.9815 * AMU * J2eV; // eV/Atom
l2[i][j][k].U=l1[i][j][k].U;

  //PLAUSIBILITY EOS CHECK:
  double echeck= EOS_ee_from_r_te(l1[i][j][k].dens, l1[i][j][k].temp * 11604.5)*26.9815 * AMU * J2eV;
  double tcheck = EOS_te_from_r_ee(l1[i][j][k].dens, echeck/26.9815/AMU * eV2J) / 11604.5;
  double tinit=l1[i][j][k].temp;
  if(ABS(tcheck -tinit) > tinit*0.01) // 1% unterschied
  {
    char errstr[255];

    sprintf(errstr,"ERROR: EOS Plausibility check failed, TfromU != Tinit. Tinit:%.4e, TfromU:%.4e\n"
                    "Maybe Interpolation table too sparse or increase tolerance",tinit,tcheck); 
    error(errstr);
  }
  

          }
        }
      } // for k
    } // for j
  } //for i
#if DEBUG_LEVEL>1
  printf("steps:%d,proc:%d,update_fd complete\n", steps, myid);
#endif
  
//NOW COMMUNICATE DIRICHLET ARRAYS
#ifdef DIRICHLET
  MPI_Allreduce(&dirichlet_maxy_local[0], &dirichlet_maxy_global[0], global_fd_dim.x, MPI_INT, MPI_MAX, cpugrid);
  MPI_Allreduce(&dirichlet_miny_local[0], &dirichlet_miny_global[0], global_fd_dim.x, MPI_INT, MPI_MIN, cpugrid);
  MPI_Allreduce(&dirichlet_maxx_local[0], &dirichlet_maxx_global[0], global_fd_dim.y, MPI_INT, MPI_MAX, cpugrid);
#endif

}
// ***************************************************************************************************************************
// * do_FILLMESH computes wide-range electronic transport and optical properties of fd-cells
// ****************************************************************************************************************************
void do_FILLMESH(void)
{
  int i, j, k, i_global, j_global, k_global;

#if DEBUG_LEVEL>1
  printf("steps:%d,proc:%d,entered FILLMESH\n", steps, myid);
#endif

  // ********************************************************************************************
  // * Dont skip this loop, even if no laser is active, since here fd_k,fd_g and Z is computed
  // *********************************************************************************************
#if DEBUG_LEVEL>1
  printf("proc:%d,steps:%d,entered 2nd loop\n", myid, steps);
#endif

  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {
    i_global =  ((i - 1) + my_coord.x * (local_fd_dim.x - 2));
    for (j = 1; j < local_fd_dim.y - 1; ++j)
    {
      j_global =  ((j - 1) + my_coord.y * (local_fd_dim.y - 2));
      for (k = 1; k < local_fd_dim.z - 1; ++k)
      {
        k_global =  ((k - 1) + my_coord.z * (local_fd_dim.z - 2));


        if (l1[i][j][k].natoms >= fd_min_atoms)
        {
          //////////////////////////////////////////////////////////////          
          //           IONISATIONSGRAD UND ELEK.DICHTE
          ////////////////////////////////////////////////////////////// 
          l1[i][j][k].Z=MeanCharge(l1[i][j][k].temp*11604.5, l1[i][j][k].dens, atomic_charge, atomic_weight,i,j,k);

          // l1[i][j][k].Z=QfromT(l1[i][j][k].temp,l1[i][j][k].dens);

          if (l1[i][j][k].Z == -1.0)
          {
            printf("steps:%d,proc:%d,i:%d,j:%d,k:%d, ERROR during QfromT in FILLMESH, Te:%f (K), dens:%f (kg/m^3),atoms:%d\n",
                   steps, myid, i, j, k, l1[i][j][k].temp * 11604.5, l1[i][j][k].dens, l1[i][j][k].natoms);
            error("ERROR during QfromT in FILLMESH");
          }

          l1[i][j][k].ne = l1[i][j][k].Z * l1[i][j][k].dens / (atomic_weight * AMU); //Assumption: Quasi-neutrality condition
          l2[i][j][k].ne = l1[i][j][k].ne;

          //////////////////////////////////////////////////////////////          
          //                      Wärmekapazität
          //////////////////////////////////////////////////////////////      
          l1[i][j][k].Ce = EOS_cve_from_r_te(l1[i][j][k].natoms, l1[i][j][k].temp * 11604.5); //Interpol.Tabelle          
          l2[i][j][k].Ce=l1[i][j][k].Ce;

#if DEBUG_LEVEL>0
          if (l1[i][j][k].Ce < 0)
          {
            printf("steps:%d,proc:%d,i:%d,j:%d,k:%d, ERROR during CfromT in FILLMESH, Ce:%.4e, Te:%f (K), dens:%f (kg/m^3)\n",
                   steps, myid, i, j, k, l1[i][j][k].Ce, l1[i][j][k].temp * 11604.5, l1[i][j][k].dens);
            error("ERROR during CfromT in FILLMESH");
          }
#endif
          //////////////////////////////////////////////////////////////          
          //                      KAPPA 
          //////////////////////////////////////////////////////////////          
          l1[i][j][k].fd_k = getKappa(l1[i][j][k].temp, l1[i][j][k].md_temp, l1[i][j][k].ne, l1[i][j][k].Z); //Hardcoding ist faster
// node.fd_k=fd_k;          
          l2[i][j][k].fd_k = l1[i][j][k].fd_k;

          /* //TRIKUBISCHE INTERPOLATION AUS TABELLE // trikub. interpol ist recht langsam 
                     l1[i][j][k].fd_k=KappaInterpol(l1[i][j][k].dens,l1[i][j][k].temp,l1[i][j][k].md_temp);
                     if(l1[i][j][k].fd_k==-1.0)
                     {
                        printf("steps:%d,proc:%d,i:%d,j:%d,k:%d, ERROR during KappaInterpol in FILLMESH, Te:%f (eV), Ti:%f (eV) ,dens:%f (kg/m^3),atoms:%d\n",
                                steps,myid,i,j,k,l1[i][j][k].temp,l1[i][j][k].md_temp,l1[i][j][k].dens,l1[i][j][k].natoms);
                        error("ERROR during KappaInterpol in FILLMESH");
                     }
          */

#if DEBUG_LEVEL>0  //ACHTUNG: Das ist nur bei Hardcoding noetig. Bei Interpolation, sollte die Interpol-Funktion error-check machen
          //NaN kappa
          if (isnan(l1[i][j][k].fd_k) != 0 || l1[i][j][k].fd_k < 0) //&& steps>0 || l1[i][j][k].fd_k<0 && steps>0)
          {
            printf("proc:%d,i:%d,j:%d,k:%d,steps:%d, atoms:%d,fd_k is NaN %f\n", myid, i, j, k, steps, l1[i][j][k].natoms, l1[i][j][k].fd_k);
            error("fd_k is NaN.");
          }
#endif
          //////////////////////////////////////////////////////////////   
          //                   GAMMA  (KOPPL.CONST)
          //////////////////////////////////////////////////////////////          
          
          l1[i][j][k].fd_g = getGamma(l1[i][j][k].temp, l1[i][j][k].md_temp, l1[i][j][k].ne, l1[i][j][k].Z);
// node.fd_g=fd_g;          
          l2[i][j][k].fd_g = l1[i][j][k].fd_g;
          //NaN gamma
#if DEBUG_LEVEL>0
          if (isnan(l1[i][j][k].fd_g) != 0 || l1[i][j][k].fd_g < 0)
          {
            printf("proc:%d,i:%d,j:%d,k:%d,steps:%d, fd_g is NaN %f\n", myid, i, j, k, steps, l1[i][j][k].fd_g);
            error("fd_g is NaN.");
          }
#endif
          //////////////////////////////////////////////////////////////          
          //        DRUDE - LORENTZ PARAMS FUER MAXWELL SOLVER
          //////////////////////////////////////////////////////////////   
#ifdef FDTD
          
          int fitresult = fitDL(i, j, k);
#if DEBUG_LEVEL>0
          if (fitresult == -1.0)
          {
            printf("steps:%d,proc:%d,i:%d,j:%d,k:%d, ERROR during fitDL in FILLMESH, Te:%f (K), dens:%f (kg/m^3)\n",
                   steps, myid, i, j, k, l1[i][j][k].temp * 11604.5, l1[i][j][k].dens);
            error("ERROR during fitDL in FILLMESH");
          }
#endif //DEBUG LEVEL
#endif  //FDTD     
//
        }// if >= min_atoms ....
        else
        {
          l1[i][j][k].fd_k = 0.0;
          l1[i][j][k].fd_g = 0.0;
          l1[i][j][k].Z = 0.0;
          l1[i][j][k].ne = 0.0;
          l1[i][j][k].Ce = 0.0;
          //l1[i][j][k].temp=0.0;  //nicht nullen, sonst funktioniert advection nicht
        }
      } // for k ....
    }   // for j
  }     // for i

#if DEBUG_LEVEL>1
  printf("steps:%d,proc:%d,FILLMESH complete\n", steps, myid);
#endif
}
// **********************************************************************************
// * ROUTINE TO COMMUNICATE ONCE EVERY MD-STEP
// * USING POINT-2-POINT MPI_Sendrecv(...) with local copies
// ***********************************************************************************
void do_COMMFLUX(void)
{
#if ADVMODE==0
  return;
#endif
  int i, j, k, i_global, j_global, k_global;
#if DEBUG_LEVEL>1
//  if(DEBUG_LEVEL>0)
  printf("steps:%d,proc:%d,entered do_COMMFLUX\n", steps, myid);
#endif

//   // //Ersst mal interne energie berechnen in jeder Zelle
//   for (i=1;i<local_fd_dim.x-1;++i)
//   {
//     for (j=1; j<local_fd_dim.y-1; ++j)
//     {
//       for (k=1; k<local_fd_dim.z-1; ++k)
//       {
//         //if(l1[i][j][k].natoms < fd_min_atoms) continue;
// /*        
//         if(l1[i][j][k].dens < RHOMIN)
//         {
//           l1[i][j][k].U=l2[i][j][k].U=0.0;
//           continue;
//         } 
// */        
//         //Brauche nicht in eV/Atom umrechnen...
//         //l1[i][j][k].U=EOS_ee_from_r_te(l1[i][j][k].dens, l1[i][j][k].temp*11604.5);//*26.9815*AMU*J2eV; //J/kg -> eV/Atom        
//         //l2[i][j][k].U=l1[i][j][k].U; //wegen diff-step
// if(l1[i][j][k].natoms >= 1)
//         //U : temp * eV / Vol / temp ----> eV/ mass / temp *(atoms*mass(atom)) ---> eV/Atom
//         l1[i][j][k].U=l2[i][j][k].U=l1[i][j][k].temp*0.00291669/(l1[i][j][k].dens*6.02214076e-4)*26.9815;
//         //l1[i][j][k].U=l2[i][j][k].U=l1[i][j][k].temp*0.00291669*fd_vol/node.natoms;
//       else
//         l1[i][j][k].U=l2[i][j][k].U=0.0;

//         // l1[i][j][k].U=EOS_ee_from_r_te(l1[i][j][k].dens, l1[i][j][k].temp*11604.5)*l1[i][j][k].dens*(fd_vol*1e-30)*J2eV;
//         // l2[i][j][k].U=l1[i][j][k].U; //wegen diff-step

//         //CHECK
//         //double tmp=EOS_te_from_r_ee(l1[i][j][k].dens, l1[i][j][k].U/(26.9815*AMU*J2eV))/11604.5;
//         //printf("r:%.4e,t:%.4e,ufromt:%.4e,tfromu:%.4e\n",
//         //        l1[i][j][k].dens, l1[i][j][k].temp, l1[i][j][k].U,tmp);

//       }
//     }
//   }

  ///////////////////////////////////////////////////////////////////////////////////////////
  // Kommunikation von diversen größen, die nur jeden MD-step gebraucht werden
  // statt jeden fd-step (dafür gibts ja fill_ghost_cells)
  ///////////////////////////////////////////////////////////////////////////////////////////
  /** MPI communication (can occur before and/or during MD calculations?) */
  /* Remember:
   * east -> -x
   * west -> +x
   * north-> -y
   * south-> +y
   * up   -> -z
   * down -> +z
   * *************/
#if DEBUG_LEVEL>1
  printf("proc:%d,steps:%d,entered comm loop\n", steps, myid);
#endif

  //first comm. yz-plane (+x/-x)
  for (j = 1; j < local_fd_dim.y - 1; ++j)
  {
    for (k = 1; k < local_fd_dim.z - 1; ++k)
    {
      if (my_coord.x != 0 && my_coord.x != cpu_dim.x - 1) //BULK
      {
        //flux
        MPI_Sendrecv(l1[local_fd_dim.x - 2][j][k].flux, 8, MPI_INT, nbwest, 7302, //+x
                     l1[0][j][k].flux, 8, MPI_INT, nbeast, 7302,                  //-x
                     cpugrid, &stati[0]);

        MPI_Sendrecv(l1[1][j][k].flux, 8, MPI_INT, nbeast, 7402,                  //-x
                     l1[local_fd_dim.x - 1][j][k].flux, 8, MPI_INT, nbwest, 7402, //+x
                     cpugrid, &stati[1]);
        //selbes fuer interne eng.
        MPI_Sendrecv(&l1[local_fd_dim.x - 2][j][k].U, 1, MPI_DOUBLE, nbwest, 5302, //+x
                     &l1[0][j][k].U, 1, MPI_DOUBLE, nbeast, 5302,                  //-x
                     cpugrid, &stati[0]);

        MPI_Sendrecv(&l1[1][j][k].U, 1, MPI_DOUBLE, nbeast, 5402,                  //-x
                     &l1[local_fd_dim.x - 1][j][k].U, 1, MPI_DOUBLE, nbwest, 5402, //+x
                     cpugrid, &stati[1]);
      }
      else if (my_coord.x == 0) //SURFACE:only comm. with nbwest
      {
        //flux
        MPI_Sendrecv(l1[local_fd_dim.x - 2][j][k].flux, 8, MPI_INT, nbwest, 7302, //+x
                     l1[local_fd_dim.x - 1][j][k].flux, 8, MPI_INT, nbwest, 7402, //+x
                     cpugrid, &stati[1]);
        //selbes fuer u
        MPI_Sendrecv(&l1[local_fd_dim.x - 2][j][k].U, 1, MPI_DOUBLE, nbwest, 5302, //+x
                     &l1[local_fd_dim.x - 1][j][k].U, 1, MPI_DOUBLE, nbwest, 5402, //+x
                     cpugrid, &stati[1]);
      }
      else if (my_coord.x == cpu_dim.x - 1) //SURFACE:only comm with east
      {
        //flux
        MPI_Sendrecv(l1[1][j][k].flux, 8, MPI_INT, nbeast, 7402,         //-x
                     l1[0][j][k].flux, 8, MPI_INT, nbeast, 7302,         //-x
                     cpugrid, &stati[0]);
        //selbes fuer u
        MPI_Sendrecv(&l1[1][j][k].U, 1, MPI_DOUBLE, nbeast, 5402,         //-x
                     &l1[0][j][k].U, 1, MPI_DOUBLE, nbeast, 5302,         //-x
                     cpugrid, &stati[0]);
      }
    }
  }

  // 2nd: comm. xz-plane (+y/-y)
  // north-> -y
  // south-> +y
#ifdef ADV2D
  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {
    for (k = 1; k < local_fd_dim.z - 1; ++k)
    {
      if (pbc_dirs.y == 1 || (my_coord.y != 0 && my_coord.y != cpu_dim.y - 1)) //BULK
      {
        //flux
        MPI_Sendrecv(l1[i][local_fd_dim.y - 2][k].flux, 8, MPI_INT, nbsouth, 732,  //+y
                     l1[i][0][k].flux, 8, MPI_INT, nbnorth, 732,                  //-y
                     cpugrid, &stati[2]);

        MPI_Sendrecv(l1[i][1][k].flux, 8, MPI_INT, nbnorth, 742,                  //-y
                     l1[i][local_fd_dim.y - 1][k].flux, 8, MPI_INT, nbsouth, 742, //+y
                     cpugrid, &stati[3]);
        //selbes fuer u
        MPI_Sendrecv(&l1[i][local_fd_dim.y - 2][k].U, 1, MPI_DOUBLE, nbsouth, 730,  //+y
                     &l1[i][0][k].U, 1, MPI_DOUBLE, nbnorth, 730,                  //-y
                     cpugrid, &stati[2]);

        MPI_Sendrecv(&l1[i][1][k].U, 1, MPI_DOUBLE, nbnorth, 740,                  //-y
                     &l1[i][local_fd_dim.y - 1][k].U, 1, MPI_DOUBLE, nbsouth, 740, //+y
                     cpugrid, &stati[3]);
      }
      else //no pbc and surface
      {
        if (my_coord.y == 0 && my_coord.y != cpu_dim.y - 1) //only communicate with +y
        {
          //flux
          MPI_Sendrecv(l1[i][local_fd_dim.y - 2][k].flux, 8, MPI_INT, nbsouth, 732,  //+y
                       l1[i][local_fd_dim.y - 1][k].flux, 8, MPI_INT, nbsouth, 742,  //+y
                       cpugrid, &stati[3]);
          MPI_Sendrecv(&l1[i][local_fd_dim.y - 2][k].U, 1, MPI_DOUBLE, nbsouth, 730,  //+y
                       &l1[i][local_fd_dim.y - 1][k].U, 1, MPI_DOUBLE, nbsouth, 740,  //+y
                       cpugrid, &stati[3]);
        }
        else if (my_coord.y == cpu_dim.y - 1 && my_coord.y != 0) // only communicate with -y
        {
          //flux
          MPI_Sendrecv(l1[i][1][k].flux, 8, MPI_INT, nbnorth, 742,      //-y
                       l1[i][0][k].flux, 8, MPI_INT, nbnorth, 732,       //-y
                       cpugrid, &stati[2]);
          MPI_Sendrecv(&l1[i][1][k].U, 1, MPI_DOUBLE, nbnorth, 740,      //-y
                       &l1[i][0][k].U, 1, MPI_DOUBLE, nbnorth, 730,       //-y
                       cpugrid, &stati[2]);
        }
      }
    }
  }
#endif //ADV2D

  // Zusätzliche Kommunikation für ADV_MODE==2
  // ********************************
  // EDGES COMM (erstmal nur 2D)
  //  3rd comm +x/+y <---> -x/-y
  //  +x/+y = nbws,  -x/-y = nben
  //
  // 4th comm +x/-y <---> -x/+y
  //  +x/-y = nbnw, <---> -x/+y=nbse
  //
  // *********************************
#ifdef ADV2D
  if ((my_coord.y != 0 && my_coord.y != cpu_dim.y - 1) && (my_coord.x != 0 && my_coord.x != cpu_dim.x - 1)) //BULK
  {
    MPI_Sendrecv(l1[local_fd_dim.x - 2][local_fd_dim.y - 2][1].flux, 8, MPI_INT, nbws, 632, //nach +x/+y
                 l1[0][0][1].flux, 8, MPI_INT, nben, 632,                       //von -x/-y
                 cpugrid, &stati[2]);

    MPI_Sendrecv(l1[1][1][1].flux, 8, MPI_INT, nben, 631,                       //nach -x/-y
                 l1[local_fd_dim.x - 1][local_fd_dim.y - 1][1].flux, 8, MPI_INT, nbws, 631, //von +x/+y
                 cpugrid, &stati[3]);

    /////////////////////

    MPI_Sendrecv(l1[local_fd_dim.x - 2][1][1].flux, 8, MPI_INT, nbnw, 622,    //nach +x/-y
                 l1[0][local_fd_dim.y - 1][1].flux, 8, MPI_INT, nbse, 622,                //von -x/+y
                 cpugrid, &stati[2]);

    MPI_Sendrecv(l1[1][local_fd_dim.y - 2][1].flux, 8, MPI_INT, nbse, 621,    //nach -x/+y
                 l1[local_fd_dim.x - 1][0][1].flux, 8, MPI_INT, nbnw, 621,   //von +x/-y
                 cpugrid, &stati[3]);

    ///////////////////////
    //Selbes für interne eng.
    ///////////////////////
    MPI_Sendrecv(&l1[local_fd_dim.x - 2][local_fd_dim.y - 2][1].U, 1, MPI_DOUBLE, nbws, 532, //nach +x/+y
                 &l1[0][0][1].U, 1, MPI_DOUBLE, nben, 532,                                   //von -x/-y
                 cpugrid, &stati[2]);

    MPI_Sendrecv(&l1[1][1][1].U, 1, MPI_DOUBLE, nben, 531,                                   //nach -x/-y
                 &l1[local_fd_dim.x - 1][local_fd_dim.y - 1][1].U, 1, MPI_DOUBLE, nbws, 531,  //von +x/+y
                 cpugrid, &stati[3]);
    /////////////////////
    MPI_Sendrecv(&l1[local_fd_dim.x - 2][1][1].U, 1, MPI_DOUBLE, nbnw, 522,             //nach +x/-y
                 &l1[0][local_fd_dim.y - 1][1].U, 1, MPI_DOUBLE, nbse, 522,                //von -x/+y
                 cpugrid, &stati[2]);

    MPI_Sendrecv(&l1[1][local_fd_dim.y - 2][1].U, 1, MPI_DOUBLE, nbse, 521,                //nach -x/+y
                 &l1[local_fd_dim.x - 1][0][1].U, 1, MPI_DOUBLE, nbnw, 521,                //von +x/-y
                 cpugrid, &stati[3]);
  }
  //Jetzt noch fuer Ränder
  else if (my_coord.x == 0)
  {
    if (my_coord.y == 0)
    {
      MPI_Sendrecv(l1[local_fd_dim.x - 2][local_fd_dim.y - 2][1].flux, 8, MPI_INT, nbws, 632, //nach +x/+y
                   l1[local_fd_dim.x - 1][local_fd_dim.y - 1][1].flux, 8, MPI_INT, nbws, 631, cpugrid, &stati[2]); //von +x/+y
      //selbes fuer u
      MPI_Sendrecv(&l1[local_fd_dim.x - 2][local_fd_dim.y - 2][1].U, 1, MPI_DOUBLE, nbws, 532,
                   &l1[local_fd_dim.x - 1][local_fd_dim.y - 1][1].U, 1, MPI_DOUBLE, nbws, 531, cpugrid, &stati[2]);

    }
    else if (my_coord.y == cpu_dim.y - 1)
    {
      MPI_Sendrecv(l1[local_fd_dim.x - 2][0][1].flux, 8, MPI_INT, nbnw, 622,
                   l1[local_fd_dim.x - 1][1][1].flux, 8, MPI_INT, nbnw, 621, cpugrid, &stati[3]);
      //selbes fuer u
      MPI_Sendrecv(&l1[local_fd_dim.x - 2][0][1].U, 1, MPI_DOUBLE, nbnw, 522,
                   &l1[local_fd_dim.x - 1][1][1].U, 1, MPI_DOUBLE, nbnw, 521, cpugrid, &stati[3]);
    }
    //else if(my_coord.y!=cpu_dim.y-1 && my_coord.y!=0) //In diesem Fall beides
    else
    {
      MPI_Sendrecv(l1[local_fd_dim.x - 2][local_fd_dim.y - 2][1].flux, 8, MPI_INT, nbws, 632,
                   l1[local_fd_dim.x - 1][local_fd_dim.y - 1][1].flux, 8, MPI_INT, nbws, 631, cpugrid, &stati[2]);

      MPI_Sendrecv(l1[local_fd_dim.x - 2][0][1].flux, 8, MPI_INT, nbnw, 622,
                   l1[local_fd_dim.x - 1][1][1].flux, 8, MPI_INT, nbnw, 621, cpugrid, &stati[3]);
      //selbes fuer u
      MPI_Sendrecv(&l1[local_fd_dim.x - 2][local_fd_dim.y - 2][1].U, 1, MPI_DOUBLE, nbws, 532,
                   &l1[local_fd_dim.x - 1][local_fd_dim.y - 1][1].U, 1, MPI_DOUBLE, nbws, 531, cpugrid, &stati[2]);

      MPI_Sendrecv(&l1[local_fd_dim.x - 2][0][1].U, 1, MPI_DOUBLE, nbnw, 522,
                   &l1[local_fd_dim.x - 1][1][1].U, 1, MPI_DOUBLE, nbnw, 521, cpugrid, &stati[3]);
    }
  }
  else if (my_coord.x == cpu_dim.x - 1)
  {
    if (my_coord.y == 0)
    {
      MPI_Sendrecv(l1[1][local_fd_dim.y - 2][1].flux, 8, MPI_INT, nbse, 621,
                   l1[0][local_fd_dim.y - 1][1].flux, 8, MPI_INT, nbse, 622, cpugrid, &stati[2]);

      //selbes fuer u
      MPI_Sendrecv(&l1[1][local_fd_dim.y - 2][1].U, 1, MPI_DOUBLE, nbse, 521,
                   &l1[0][local_fd_dim.y - 1][1].U, 1, MPI_DOUBLE, nbse, 522, cpugrid, &stati[2]);

    }
    else if (my_coord.y == cpu_dim.y - 1)
    {
      MPI_Sendrecv(l1[1][1][1].flux, 8, MPI_INT, nben, 631,
                   l1[0][0][1].flux, 8, MPI_INT, nben, 632, cpugrid, &stati[3]);
      //selbes fuer u
      MPI_Sendrecv(&l1[1][1][1].U, 1, MPI_DOUBLE, nben, 531,
                   &l1[0][0][1].U, 1, MPI_DOUBLE, nben, 532, cpugrid, &stati[3]);
    }
    else //beides
    {
      MPI_Sendrecv(l1[1][local_fd_dim.y - 2][1].flux, 8, MPI_INT, nbse, 621,
                   l1[0][local_fd_dim.y - 1][1].flux, 8, MPI_INT, nbse, 622, cpugrid, &stati[2]);

      MPI_Sendrecv(l1[1][1][1].flux, 8, MPI_INT, nben, 631,
                   l1[0][0][1].flux, 8, MPI_INT, nben, 632, cpugrid, &stati[3]);
      //selbes fuer u
      MPI_Sendrecv(&l1[1][local_fd_dim.y - 2][1].U, 1, MPI_DOUBLE, nbse, 521,
                   &l1[0][local_fd_dim.y - 1][1].U, 1, MPI_DOUBLE, nbse, 522, cpugrid, &stati[2]);

      MPI_Sendrecv(&l1[1][1][1].U, 1, MPI_DOUBLE, nben, 531,
                   &l1[0][0][1].U, 1, MPI_DOUBLE, nben, 532, cpugrid, &stati[3]);
    }
  }
  else  //d.h. my_coord.x!=0 und my_coord.x!=cpu_dim.x-1 und my_coord.y ist entweder 0 oder cpu_dim.y-1
  {
    if (my_coord.y == 0) //comm. mit +x,+y und -x,+y
    {
      MPI_Sendrecv(l1[local_fd_dim.x - 2][local_fd_dim.y - 2][1].flux, 8, MPI_INT, nbws, 632,
                   l1[local_fd_dim.x - 1][local_fd_dim.y - 1][1].flux, 8, MPI_INT, nbws, 631, cpugrid, &stati[2]);

      MPI_Sendrecv(l1[1][local_fd_dim.y - 2][1].flux, 8, MPI_INT, nbse, 621,
                   l1[0][local_fd_dim.y - 1][1].flux, 8, MPI_INT, nbse, 622, cpugrid, &stati[3]);

      //selbes fuer u
      MPI_Sendrecv(&l1[local_fd_dim.x - 2][local_fd_dim.y - 2][1].U, 1, MPI_DOUBLE, nbws, 532,
                   &l1[local_fd_dim.x - 1][local_fd_dim.y - 1][1].U, 1, MPI_DOUBLE, nbws, 531, cpugrid, &stati[2]);

      MPI_Sendrecv(&l1[1][local_fd_dim.y - 2][1].U, 1, MPI_DOUBLE, nbse, 521,
                   &l1[0][local_fd_dim.y - 1][1].U, 1, MPI_DOUBLE, nbse, 522, cpugrid, &stati[3]);
    }
    else if (my_coord.y == cpu_dim.y - 1) // comm. mit -x,-y und +x,-y
    {
      MPI_Sendrecv(l1[1][1][1].flux, 8, MPI_INT, nben, 631,
                   l1[0][0][1].flux, 8, MPI_INT, nben, 632, cpugrid, &stati[2]);

      MPI_Sendrecv(l1[local_fd_dim.x - 2][1][1].flux, 8, MPI_INT, nbnw, 622,
                   l1[local_fd_dim.x - 1][0][1].flux, 8, MPI_INT, nbnw, 621, cpugrid, &stati[3]);

      //selbes fuer u
      MPI_Sendrecv(&l1[1][1][1].U, 1, MPI_DOUBLE, nben, 531,
                   &l1[0][0][1].U, 1, MPI_DOUBLE, nben, 532, cpugrid, &stati[2]);

      MPI_Sendrecv(&l1[local_fd_dim.x - 2][1][1].U, 1, MPI_DOUBLE, nbnw, 522,
                   &l1[local_fd_dim.x - 1][0][1].U, 1, MPI_DOUBLE, nbnw, 521, cpugrid, &stati[3]);

    }
  }
#endif //ADVMODE2d
#if DEBUG_LEVEL>1
//  if(DEBUG_LEVEL>0)
  printf("proc:%d,steps:%d,do_COMMFLUX complete\n", myid, steps);
#endif

// // Wozu brauche ich E_new noch???  -> bald entfernen, erstmal lassen,wollen ja keine unnoetigen bugs
// #ifdef MPI
//   {
//     double E_new_reduced;
//     MPI_Reduce(&E_new_local, &E_new_reduced, 1, MPI_DOUBLE, MPI_SUM, 0, cpugrid);
//     if (myid == 0) E_new += E_new_reduced * fd_h.x * fd_h.y * fd_h.z / natoms;
//   }
// #else
//   E_new += E_new_local * fd_h.x * fd_h.y * fd_h.z / natoms;
// #endif // MPI

//   E_new_local = 0.0; // ebenfalls unnoetig
}

/****************************************************************************
 *        TTM INIT                  *
*****************************************************************************/
/* init_ttm(): initialize FD stuff, ttm_init */
void init_ttm()
{
  int i, j, k;

  if (myid == 0)
    printf("imdrestart=%d\n", imdrestart);
#if DEBUG_LEVEL>1
  printf("proc:%d,entered init_ttm()\n", myid);
#endif
  // *****************************************
  // * READ AND BCAST INTERPOLATION TABLES
  // ******************************************
  // read_bc_interp(&QfromT_interp,"EOS_QfromT.txt");
  // read_bc_interp(&CfromT_interp,"EOS_CfromT.txt"); //für CFL maxdt

nn_read_table(&intp_cve_from_r_te, "EOS_cve_from_r_te.txt");
nn_read_table(&intp_ee_from_r_tesqrt, "EOS_ee_from_r_tesqrt.txt");

  //read_tricub_interp(&kappa_interp,"kappa.txt"); //Hardcoding ist schneller
  //Lese Drude-Lorentz Interpolationstabellen


#ifdef FDTD
  read_tricub_interp(&Lop1i, "DL1.txt");
  read_tricub_interp(&Lop2i, "DL2.txt");
  read_tricub_interp(&Lop3i, "DL3.txt");
  read_tricub_interp(&Lop4i, "DL4.txt");
  read_tricub_interp(&Lop5i, "DL5.txt");
#endif


  fd_vol = fd_h.x * fd_h.y * fd_h.z;
  natoms_local = 0;
  n_md_cells = fd_ext.x * fd_ext.y * fd_ext.z;
  max_dt_ttm = timestep / ((double) fd_n_timesteps); //nur zu beginn...im weiteren verlauf adaptiv

  //neighvol needed for per-atom density calc. only for single-species simulations
// #ifdef NBL
//   neighvol = pow(cellsz, 1.5) * 4.0 / 3.0 * M_PI;
// #else
//   neighvol = pow(sqrt(pair_pot.end[0]), 3.0) * 4.0 / 3.0 * M_PI;
// #endif
neighvol = pow(sqrt(pair_pot.end[0]), 3.0) * 4.0 / 3.0 * M_PI;

  /* Check if cell_dim and fd_ext are commensurate */
  if ( fd_one_d == 1 || fd_one_d == 0 )
    if ( (cell_dim.x - NBUFFC) % fd_ext.x != 0 )
    {
      char buf[255];
      sprintf(buf,
              "cell_dim and fd_ext are not commensurate:\ncell_dim.x=%d; fd_ext.x=%d\n",
              cell_dim.x, fd_ext.x );
      error (buf);
    }
  if ( fd_one_d == 2 || fd_one_d == 0 )
    if ( (cell_dim.y - NBUFFC) % fd_ext.y != 0 )
    {
      char buf[255];
      sprintf(buf,
              "cell_dim and fd_ext are not commensurate:\ncell_dim.y=%d; fd_ext.y=%d\n",
              cell_dim.y, fd_ext.y );
      error (buf);
    }
  if ( fd_one_d == 3 || fd_one_d == 0 )
    if ( (cell_dim.z - NBUFFC) % fd_ext.z != 0 )
    {
      char buf[255];
      sprintf(buf,
              "cell_dim and fd_ext are not commensurate:\ncell_dim.z=%d; fd_ext.z=%d\n",
              cell_dim.z, fd_ext.z );
      error (buf);
    }

  /* local size of FD lattice, will add 2 layers later for boundaries */
  local_fd_dim.x = (cell_dim.x - NBUFFC) / fd_ext.x;
  local_fd_dim.y = (cell_dim.y - NBUFFC) / fd_ext.y;
  local_fd_dim.z = (cell_dim.z - NBUFFC) / fd_ext.z;

  /* get these to the right sizes */
  global_fd_dim.x = local_fd_dim.x * cpu_dim.x;
  global_fd_dim.y = local_fd_dim.y * cpu_dim.y;
  global_fd_dim.z = local_fd_dim.z * cpu_dim.z;
  local_fd_dim.x += 2; /* 2 for ghost layers */
  local_fd_dim.y += 2; /* 2 for ghost layers */
  local_fd_dim.z += 2; /* 2 for ghost layers */

  /* Time to initialize our FD lattice... */

  /* Allocate memory for two lattices */
  lattice1 = (ttm_Element*) malloc(
               (local_fd_dim.x) * (local_fd_dim.y) * (local_fd_dim.z) * sizeof(ttm_Element));
  lattice2 = (ttm_Element*) malloc(
               (local_fd_dim.x) * (local_fd_dim.y) * (local_fd_dim.z) * sizeof(ttm_Element));
  l1 = (ttm_Element***) malloc( local_fd_dim.x * sizeof(ttm_Element**) );
  l2 = (ttm_Element***) malloc( local_fd_dim.x * sizeof(ttm_Element**) );
  for (i = 0; i < local_fd_dim.x; i++)
  {
    l1[i] = (ttm_Element**) malloc( local_fd_dim.y * sizeof(ttm_Element*) );
    l2[i] = (ttm_Element**) malloc( local_fd_dim.y * sizeof(ttm_Element*) );

    for (j = 0; j < local_fd_dim.y; j++)
    {
      l1[i][j] = lattice1 + i * local_fd_dim.y * local_fd_dim.z + j * local_fd_dim.z;
      l2[i][j] = lattice2 + i * local_fd_dim.y * local_fd_dim.z + j * local_fd_dim.z;

//      if (i!=0 && j!=0 && i!=local_fd_dim.x-1 && j!=local_fd_dim.y-1)
      { /* we are not in a ghost layer*/
        for (k = 0; k < local_fd_dim.z; k++)
        {
          int tmpindex = 0;
          int xc, yc, zc;

          /* Initialize md_cellptrs, temp,
           * and source of this FE Cell.
           * Also set fd_cell_idx in encompassed MD cells.
           **********/

          /* allocate MD-cell pointer arrays */
          if (i != 0 && j != 0 && i != local_fd_dim.x - 1 && j != local_fd_dim.y - 1 && k != 0 && k != local_fd_dim.z - 1)
          {
            l1[i][j][k].md_cellptrs = (cellptr*)malloc(n_md_cells * sizeof(cellptr));
            l2[i][j][k].md_cellptrs = (cellptr*)malloc(n_md_cells * sizeof(cellptr));

            /* loop over encompassed MD cells */
            for (xc = (i - 1) * fd_ext.x + NBUFFC / 2; xc < i * fd_ext.x + NBUFFC / 2; xc++)
            {
              for (yc = (j - 1) * fd_ext.y + NBUFFC / 2; yc < j * fd_ext.y + NBUFFC / 2; yc++)
              {
                for (zc = (k - 1) * fd_ext.z + NBUFFC / 2; zc < k * fd_ext.z + NBUFFC / 2; zc++)
                {
                  cellptr p;

                  /* pointer to this MD cell */
                  p = l1[i][j][k].md_cellptrs[tmpindex] =
                        l2[i][j][k].md_cellptrs[tmpindex] =
                          PTR_3D_V(cell_array, xc, yc, zc, cell_dim);

                  /* write array indices of our FD cell to this MD cell*/
                  p->fd_cell_idx.x = i;
                  p->fd_cell_idx.y = j;
                  p->fd_cell_idx.z = k;

                  natoms_local += p->n;
                  tmpindex++;
                }
              }
            }
          }

#ifdef DEBUG
          assert(tmpindex == n_md_cells);
#endif

          /* no incoming thermal power per default */
          l2[i][j][k].source = l1[i][j][k].source = 0.0;
          l2[i][j][k].fd_k = l1[i][j][k].fd_k = 0;
          l1[i][j][k].fd_g = l1[i][j][k].fd_g = 0;
          l2[i][j][k].proc = l1[i][j][k].proc = myid;
          l2[i][j][k].xi = l1[i][j][k].xi = 0.0;
          l2[i][j][k].ne = l1[i][j][k].ne = 0.0;
          l2[i][j][k].temp = l1[i][j][k].temp = 0.0;
          l2[i][j][k].dens = l1[i][j][k].dens = 0.0;
          l2[i][j][k].natoms = l1[i][j][k].natoms = 0;
          l2[i][j][k].natoms_old = l1[i][j][k].natoms_old = 0;
          l2[i][j][k].U = l1[i][j][k].U = 0.0;
          int bar;
#ifdef FDTD
          for (bar = 0; bar < 6; bar++)
            l2[i][j][k].DL[bar] = l1[i][j][k].DL[bar] = 0.0;
#endif
#if ADVMODE==2
          int foo;
          for (foo = 0; foo < 8; foo++)
          {
            l1[i][j][k].flux[foo] = l2[i][j][k].flux[foo] = 0;
          }
#endif
        }
      }
    }
  }

  //NOW INITIALIZE DIRICHLET BOUNDARY ARRAYS
#ifdef DIRICHLET
  alloc1darr(int, dirichlet_maxy_local, global_fd_dim.x);
  alloc1darr(int, dirichlet_maxy_global, global_fd_dim.x);

  alloc1darr(int, dirichlet_miny_local, global_fd_dim.x);
  alloc1darr(int, dirichlet_miny_global, global_fd_dim.x);

  alloc1darr(int, dirichlet_maxx_local, global_fd_dim.y);
  alloc1darr(int, dirichlet_maxx_global, global_fd_dim.y);

  for (i = 0; i < global_fd_dim.x; ++i)
  {
    dirichlet_maxy_local[i] = dirichlet_maxy_global[i] = -99999;
    dirichlet_miny_local[i] = dirichlet_miny_global[i] = 999999;
  }
  for (i = 0; i < global_fd_dim.y; ++i)
    dirichlet_maxx_local[i] = dirichlet_maxx_global[i] = -99999;

  dirichlet_surfx_int = round(dirichlet_surfx / fd_h.x);
#endif

#ifdef LASERYZ
  laser_p_peak /= (global_fd_dim.y * global_fd_dim.z);
  laser_p_peak1 /= (global_fd_dim.y * global_fd_dim.z);
#endif

  //MPI_Barrier(cpugrid);
  if (myid == 0)
  {
    printf("***************************************************\n");
    printf("*            TWO-TEMPERATURE MODEL   	      *\n");
    printf("***************************************************\n");
    printf("Global FD cell array dimensions: %d x %d x %d\n",
           global_fd_dim.x, global_fd_dim.y, global_fd_dim.z);
    printf("Local FD cell array dimensions: %d x%d x %d\n",
           local_fd_dim.x, local_fd_dim.y, local_fd_dim.z);
    printf("fd_h.x:%f A, fd_h.y:%f A,fd_h.z:%f A\n", fd_h.x, fd_h.y, fd_h.z);
    printf("Volume of one FD cell: %e [cubic Angstroms]\n",
           fd_h.x * fd_h.y * fd_h.z );
    printf("Volume of whole sample: %e [cubic Angstroms]\n",
           fd_h.x * fd_h.y * fd_h.z *
           global_fd_dim.x * global_fd_dim.y * global_fd_dim.z );
#if DEBUG_LEVEL>0
    printf("DEBUG_LEVEL>0\n");
#endif
#ifdef DIRICHLET
    printf("dirichlet_surfx=%.2f, dirichlet_surfx_int=%d\n", dirichlet_surfx, dirichlet_surfx_int);
#endif
  }

#ifdef DEBUG
  printf(
    "Found %d atoms initializing FD lattice in process number %d.\n",
    natoms_local, myid );
#endif

#ifdef MPI
  // create MPI datatypes
  ttm_create_mpi_datatypes();
#endif

  MPI_Bcast(&imdrestart, 1, MPI_INT, 0, cpugrid);
  if (imdrestart == 0)
  {
    update_fd(); /* get md_temp and v_com etc. */
  }

//DEBUG
//ttm_writeout(9999);
  /***********************
  *  MY MOD: restart ttm *
  ************************/
  if (imdrestart > 0)
  {

    int readstep = imdrestart * checkpt_int;
    int readttm = readstep / ttm_int;
    ttm_read(readttm);
    /*
        if(myid==0)
      printf("ttm_read finished\n");
    */
#ifdef FDTD
    t_SI = (double) imdrestart * (double) checkpt_int * timestep * 10.18 / 1.0e15;
    if (myid == 0)
      printf("t_SI:%.4e s\n", t_SI);
#endif
    MPI_Barrier(cpugrid);
    ttm_writeout(100000);  //only for debug
#ifdef PDECAY
    if (myid == 0)
      printf("Using ramp_start:%f and ramp_end:%f\n", ramp_start, ramp_end);
#endif
    //write_config(100000,steps);
    update_fd();
    do_FILLMESH(); //calc electronic optical and transport properties
    ttm_fill_ghost_layers();
    CFL_maxdt(); //braucht kappa von nachbarzellen,deswegen vorher fill_ghost_layers
  }

  /* time to contact neighbors and fill ghost layers */
  ttm_fill_ghost_layers();
}



// ******************************************************
// *              ADVECTION STEP      *
// *              *
// *  ADVMODE=1: --> VIA VCOM USING GODUNOV SCHEME  *
// *  ADVMODE=2: --> VIA ATOMIX FLUXES      *
// ******************************************************
void do_ADV(double tau)
{
  if (steps < 1) return;
  int i, j, k;
  int i_global, j_global, k_global;

  k = 1; //erstmal nur 2D

  do_COMMFLUX();

  //NOW CALC advection
  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {
    for (j = 1; j < local_fd_dim.y - 1; ++j)
    {
#ifdef FDTD
      SwapTTM(i, j, k);
#endif

      i_global = ((i - 1) + my_coord.x * (local_fd_dim.x - 2));
      j_global = ((j - 1) + my_coord.y * (local_fd_dim.y - 2));
      //k_global =  ((k-1) + my_coord.z*(local_fd_dim.z-2));

//REMINDER
// flux[0] : teilchen erhalten von +x,y
// flux[1] : teilchen erhalten von -x,y
// flux[2] : teilchen erhalten von x,+y
// flux[3] : teilchen erhalten von +x,-y
// flux[4] : teilchen erhalten von +x,+y
// flux[5] : teilchen erhalten von x,-y
// flux[6] : teilchen erhalten von -x,+y
// flux[7] : teilchen erhalten von -x,-y

      double Nold = (double) l1[i][j][k].natoms_old;
      double Nnew = (double) l1[i][j][k].natoms; //tmp;


      if (Nnew > 0)
      {
l2[i][j][k].U= l1[i][j][k].U * Nold / Nnew
              + tau*(
                          // +x/-x
+(double) l1[i][j][k].flux[0] * l1[i + 1][j][k].U //erhalten von +x,y
- (double) l1[i + 1][j][k].flux[1] * l1[i][j][k].U //nach +x,y abgeflossen

+ (double) l1[i][j][k].flux[1] * l1[i - 1][j][k].U //erhalten von -x,y
- (double) l1[i - 1][j][k].flux[0] * l1[i][j][k].U //nach -x,y abgeflossen
                          // +y/-y
#ifdef ADV2D
                          + (double) l1[i][j][k].flux[2] * l1[i][j + 1][k].U //erhalten von x,+y
                          - (double) l1[i][j + 1][k].flux[5] * l1[i][j][k].U //nach x,+y abgeflossen

                          + (double) l1[i][j][k].flux[5] * l1[i][j - 1][k].U //erhalten von -y,x
                          - (double) l1[i][j - 1][k].flux[2] * l1[i][j][k].U //nach -y,x abgeflossen
                          //aus +x, +y  und nach +x,+y
                          + (double) l1[i][j][k].flux[4] * l1[i + 1][j + 1][k].U //erhalten von +x,+y
                          - (double) l1[i - 1][j - 1][k].flux[7] * l1[i][j][k].U //nach +x,-y abgeflossen
                          //aus +x, -y nach +x,-y
                          + (double) l1[i][j][k].flux[3] * l1[i + 1][j - 1][k].U //erhalten von +x,-y
                          - (double) l1[i - 1][j + 1][k].flux[6] * l1[i][j][k].U //nach +x,-y abgeflossen
                          //aus -x,+y und nach -x,+y
                          + (double) l1[i][j][k].flux[6] * l1[i - 1][j + 1][k].U //erhalten von -x,+y
                          - (double) l1[i - 1][j + 1][k].flux[3] * l1[i][j][k].U //nach -x,+y abgeflossen
                          //aus -x,-y nach -x,-y
                          + (double) l1[i][j][k].flux[7] * l1[i - 1][j - 1][k].U //erhalten von -x,-y
                          - (double) l1[i - 1][j - 1][k].flux[4] * l1[i][j][k].U //abgeflossen nach -x,-y
#endif
                        ) / Nnew;


        
        l2[i][j][k].temp = EOS_te_from_r_ee(l1[i][j][k].dens, l2[i][j][k].U / (26.9815 * AMU * J2eV)) / 11604.5;
        
        //IDEE: Evtl. statt te_from_re mittels Cv und DeltaU?

        if(l2[i][j][k].temp<=0.0 )        
        {
          printf("\nmyid:%d,steps:%d,i:%d,j:%d Temp is <Tmin:%.4e in do_ADV, nnew:%.f,nold:%f,n"
           "dens: %.6e , dens_old: %.6e natoms:%d,  atoms_old:%d\n"
           "from +x,y: %d , atold:%d , t: %.4e , to +x,y: %d\n"
           "from -x,y: %d , atold:%d , t: %.4e , to -x,y: %d\n"
           "from x,+y: %d , atold:%d , t: %.4e , to x,+y: %d\n"
           "from x,-y: %d , atold:%d , t: %.4e , to x,-y: %d\n"
           "from +x,+y: %d , atold:%d , t: %.4e , to +x,+y: %d\n"
           "from +x,-y: %d , atold:%d , t: %.4e , to +x,-y: %d\n"
           "from -x,+y: %d , atold:%d , t: %.4e , to -x,+y: %d\n"
           "from -x,-y: %d , atld:%d , t: %.4e , to -x,-y: %d\n\n",
            myid,steps,i,j, l2[i][j][k].temp,
            Nnew,Nold,
           l1[i][j][k].dens, l2[i][j][k].dens, l1[i][j][k].natoms, l2[i][j][k].natoms,
           l1[i][j][k].flux[0],l2[i+1][j][k].natoms,   l1[i+1][j][k].U,   l1[i+1][j][k].flux[1],
           l1[i][j][k].flux[1],l2[i-1][j][k].natoms,   l1[i-1][j][k].U,   l1[i-1][j][k].flux[0],
           l1[i][j][k].flux[2],l2[i][j+1][k].natoms,   l1[i][j+1][k].U,   l1[i][j+1][k].flux[5],
           l1[i][j][k].flux[5],l2[i][j-1][k].natoms,   l1[i][j-1][k].U,   l1[i][j-1][k].flux[2],
           l1[i][j][k].flux[4],l2[i+1][j+1][k].natoms, l1[i+1][j+1][k].U, l1[i-1][j-1][k].flux[7],
           l1[i][j][k].flux[3],l2[i+1][j-1][k].natoms, l1[i+1][j-1][k].U, l1[i-1][j+1][k].flux[6],
           l1[i][j][k].flux[6],l2[i-1][j+1][k].natoms, l1[i-1][j+1][k].U, l1[i-1][j+1][k].flux[3],
           l1[i][j][k].flux[7],l2[i-1][j-1][k].natoms, l1[i-1][j-1][k].U, l1[i-1][j-1][k].flux[4]);
        }
        
      }
      else if (Nnew < 1)
      {
        l2[i][j][k].U = 0.0;
        l2[i][j][k].temp = 0.0;
      }
    } //for j
  } //for i

  l3 = l1;
  l1 = l2;
  l2 = l3;

}

void do_cell_activation(void)
{
  int i,j,k;
  int i_global,j_global,k_global;
    //NOW CHECK IF NEW CELL ACTIVATED!!
  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {
    i_global = ((i - 1) + my_coord.x * (local_fd_dim.x - 2));
    for (j = 1; j < local_fd_dim.y - 1; ++j)
    {
      j_global = ((j - 1) + my_coord.y * (local_fd_dim.y - 2));
      for (k = 1; k < local_fd_dim.z - 1; ++k)
      {        
        k_global =  ((k-1) + my_coord.z*(local_fd_dim.z-2));

        if (l1[i][j][k].natoms_old >= fd_min_atoms && l1[i][j][k].natoms < fd_min_atoms)
        {
          // ZELLE DEAKTIVIERT
//#if DEBUG_LEVEL>1
          printf("Warning:FD cell deactivated on proc %d on step %d at i:%d,j:%d,k%d with %d atoms and temp:%.4e\n", myid, steps,
                 i_global, j_global, k_global, l1[i][j][k].natoms, l1[i][j][k].temp);
// #endif
          // Cell deactivated. Deduce its electronic energy from E_new_local
          l1[i][j][k].xi = 0.0;
        }
        // ZELLE AKTIVIERT
        else if (l1[i][j][k].natoms_old < fd_min_atoms && l1[i][j][k].natoms >= fd_min_atoms)
        {

// #if DEBUG_LEVEL>1
          // ZELLE AKTIVIERT
          printf("Warning:New FD cell activated on proc %d at ig:%d,jg:%d,kg:%d with %d atoms on step:%d and T=%.4e, atoms_old:%d\n",
                 myid, i, j, k, l1[i][j][k].natoms, steps, l1[i][j][k].temp, l1[i][j][k].natoms_old);
// #endif
          // *****************************************************
          // * NEU AKTIVIERTE ZELLE MIT UNSINNIGER TEMPERATUR,   *
          // * d.h. ADVECTION HAT NICHT FUNKTIONIERT             *
          // * ---> WENDE ALTES SCHEMA AN UND BERECHNE MITTEL    *
          // * AUS NACHBARZELLEN               *
          // *****************************************************
          if (isnan(l1[i][j][k].temp) != 0 || l1[i][j][k].temp <= 0.003) //Temp zu klein (etwa 35K) ->wide-range props werden bullshit --> Diffusion instabil
          {

// #if DEBUG_LEVEL>0
            printf("proc:%d,steps:%d,ig:%d,jg:%d,kg:%d WARNING: Freshly activated cell with Te is NaN or < Tmin:%.4e,atoms:%d, dens: %.6e ,"
                   "using neighbor cells or mdtemp\n",
                   myid, steps, i_global, j_global, 0, l1[i][j][k].temp, l1[i][j][k].natoms, l1[i][j][k].dens);
// #endif

            // Freshly activated cell. Gets avg. electron energy of active
            // neighbor cells, the created energy is added to E_new_local
            int n_neighbors = 0;
            double E_el_neighbors = 0.0;
            // 6 indices: -x,x,-y,y,-z,z //

            // Was folgt, ist das ursprüngliche schema, mit den neu aktivierten Zellen
            // eine elec-temp. aus dem mittel der nachbarzelen zugewiesen wird ---> verletzt Energieerhaltung
            if (l1[i + 1][j][k].natoms >= fd_min_atoms)
            {
              E_el_neighbors += SQR(l1[i + 1][j][k].temp);
              n_neighbors++;
            }
            if (l1[i - 1][j][k].natoms >= fd_min_atoms)
            {
              E_el_neighbors += SQR(l1[i - 1][j][k].temp);
              n_neighbors++;
            }
            if (l1[i][j + 1][k].natoms >= fd_min_atoms)
            {
              E_el_neighbors += SQR(l1[i][j + 1][k].temp);
              n_neighbors++;
            }
            if (l1[i][j - 1][k].natoms >= fd_min_atoms)
            {
              E_el_neighbors += SQR(l1[i][j - 1][k].temp);
              n_neighbors++;
            }
            if (l1[i][j][k + 1].natoms >= fd_min_atoms)
            {
              E_el_neighbors += SQR(l1[i][j][k + 1].temp);
              n_neighbors++;
            }
            if (l1[i][j][k - 1].natoms >= fd_min_atoms)
            {
              E_el_neighbors += SQR(l1[i][j][k - 1].temp);
              n_neighbors++;
            }
            {
              if (n_neighbors != 0)
              {
                l1[i][j][k].temp = sqrt(E_el_neighbors / ((double)n_neighbors));
                l2[i][j][k].temp = l1[i][j][k].temp;

//MIT ADV Variante    
//l1[i][j][k].U = l2[i][j][k].U= EOS_ee_from_r_te(l1[i][j][k].dens, l1[i][j][k].temp * 11604.5) * 26.9815 * AMU * J2eV; // eV/Atom                

//NO-ADV Variante                
l1[i][j][k].U = l2[i][j][k].U= l1[i+1][j][k].U;
if(isnan(node.U)!= 0)
{
  printf("ERROR: U is NaN, step:%d, u+1:%.4e,i:%d\n",steps, l1[i+1][j][k].U,i);
  error("U is Nan");
}


                {
// #if DEBUG_LEVEL>0
                  printf("proc:%d,steps:%d,i:%d,j:%d,k:%d, Te is NaN or <=0, using neighbor cells=>Te=%f\n",
                         myid, steps, i_global, j_global, 0, l1[i][j][k].temp);
// #endif
//HOTFIX: still < Tmin? --> use md-temp
                  if (l1[i][j][k].temp < 0.003)
                  {
// #if DEBUG_LEVEL>0
                    printf("proc:%d,steps:%d,i:%d,j:%d,k:%d, Te still <=Tmin, using MD-temp:%.4e\n",
                           myid, steps, i_global, j_global, k_global, l1[i][j][k].md_temp);
// #endif
//HOTIFX: still stiil< Tmin ?!?! (wegen pdecay,z.B.)
                    l2[i][j][k].temp = l1[i][j][k].temp = MAX(l1[i][j][k].md_temp, 0.003);
                  }

                }
              }
              else  // No neighbors? -> Get MD-temp
              {
                l1[i][j][k].temp = l1[i][j][k].md_temp;
                l2[i][j][k].temp = l1[i][j][k].temp;
// #if DEBUG_LEVEL>0
                printf("proc:%d,steps:%d,i:%d,j:%d,k:%d, Te is NaN or <=0, using md-temp=>Te=%f\n",
                       myid, steps, i_global, j_global, k_global, l1[i][j][k].temp);

// #endif
              }
            } //isnan ....
          } // endif isnan(temp) || temp<=0

        } // endif ..new cell activated...

        //l1[i][j][k].natoms_old=l2[i][j][k].natoms_old=l1[i][j][k].natoms;
      } //for k
    } //for j
  } //for i
}

// ****************************************************************************************************
// *   DO A SINGLE DIFFUSION STEP HERE
// *   tau gibt an, wie groß der ts sein soll
// *   in IMD-Einheiten (Möglichkeit in ein Strang-Operator Splitting Schema einzubauen)
// *****************************************************************************************************
void do_DIFF(double tau)
{
  int i, j, k;
  int xmin, xmax, /* these are neighboring indices     */      
      ymin, ymax, /* (to account for bc & deactivated cells) */
      zmin, zmax;

  double xi_fac=fd_vol/3.0/((double) diff_substeps);///((double) l1[i][j][k].natoms); //ORIGINAL
  //double xi_fac = 26.9815 * AMU / 3.0 * 1e30 / ((double) diff_substeps); //NEU

  //xi=1/fdsteps * sum_(n=1)^(fdsteps) {m*fd_g/(3*rho*k_b)*(Te-Ti)/Ti }
  // (k_B weglassen --> imd-units fuer temp., AMU muss rein, weil ich rho in kg/m^3 messe,
  // dazu muss "mal" 1e30, damit Volumen in Angstrom ist
  //ich ersetzte also das originale Vol/natoms durch mass/rho
  //Um zu vermeiden, dass bei sehr geringer Zahl an atomen, die Beschleunigung
  //Zu groß wird. Mit der lokalen "Umgebungsdichte" sollte das nicht passieren
  double Ce; //specific heat
  double invxsq = 1.0 / (fd_h.x * fd_h.x);
  double invysq = 1.0 / (fd_h.y * fd_h.y);
  double invzsq = 1.0 / (fd_h.z * fd_h.z);

  int i_global, j_global, k_global;
  double xmaxTe, xminTe, ymaxTe, yminTe, zmaxTe, zminTe; //temps
  double xmaxk, xmink, ymaxk, ymink, zmaxk, zmink; //kappas



#if DEBUG_LEVEL>1
//  if(DEBUG_LEVEL>0)
  printf("steps:%d,proc:%d,entered do_DIFF\n", steps, myid);
#endif

  ///////////////////////////
  // DIFFUSION             //
  ///////////////////////////
  for (i = 1; i < local_fd_dim.x - 1; i++)
  {
    i_global =  ((i - 1) + my_coord.x * (local_fd_dim.x - 2));
    for (j = 1; j < local_fd_dim.y - 1; j++)
    {
      j_global = ((j - 1) + my_coord.y * (local_fd_dim.y - 2));
      for (k = 1; k < local_fd_dim.z - 1; k++)
      {
        k_global =  ((k - 1) + my_coord.z * (local_fd_dim.z - 2));

        //compute absorbed laser-energy
        if (laser_active)
        {
          Eabs_local += l1[i][j][k].source * fd_vol * tau; //eV
        }

#ifdef FDTD
        SwapTTM(i, j, k);
#endif

        /* only do calculation if cell is not deactivated */
        if (l1[i][j][k].natoms < fd_min_atoms)   continue;
        if (l1[i - 1][j][k].natoms < fd_min_atoms)  xmin = i; else  xmin = i - 1;
        if (l1[i + 1][j][k].natoms < fd_min_atoms) xmax = i; else  xmax = i + 1;
        if (l1[i][j - 1][k].natoms < fd_min_atoms) ymin = j; else  ymin = j - 1;
        if (l1[i][j + 1][k].natoms < fd_min_atoms) ymax = j; else  ymax = j + 1;
        if (l1[i][j][k - 1].natoms < fd_min_atoms) zmin = k; else zmin = k - 1;
        if (l1[i][j][k + 1].natoms < fd_min_atoms)  zmax = k; else zmax = k + 1;

        xmaxTe = l1[xmax][j][k].temp;
        xminTe = l1[xmin][j][k].temp;
        ymaxTe = l1[i][ymax][k].temp;
        yminTe = l1[i][ymin][k].temp;
        zmaxTe = l1[i][j][zmax].temp;
        zminTe = l1[i][j][zmin].temp;

        xmaxk = l1[xmax][j][k].fd_k;
        xmink = l1[xmin][j][k].fd_k;
        ymaxk = l1[i][ymax][k].fd_k;
        ymink = l1[i][ymin][k].fd_k;
        zmaxk = l1[i][j][zmax].fd_k;
        zmink = l1[i][j][zmin].fd_k;

        /////////////////////////////////////////
        //dirichlet cases for outermost cells  //
        //ACHTUNG: Dirichlet Zellen werden     //
        //nicht in ttm-output geschrieben.     //
        //Sie sind "virtuell"      //
        /////////////////////////////////////////
#ifdef DIRICHLET
        if (i_global >= dirichlet_surfx_int) //don't cool the ablated material
        {
          if (dirichlet_maxy_global[i_global] == j_global)
          {
            ymaxTe = 0.025850926; //=RT
            ymaxk = 1.933442e+01; //getKappa(ymaxTe, l1[i][j][k].md_temp, l1[i][j][k].ne,double Z)
          }
          else if (dirichlet_miny_global[i_global] == j_global)
          {
            yminTe = 0.025850926;
            ymink = 1.933442e+01;
          }
          if (dirichlet_maxx_global[j_global] == i_global)
          {
            xmaxTe = 0.025850926;
            xmaxk = 1.933442e+01;
          }
        }
#endif

         Ce = l1[i][j][k].Ce; //eV/(eV*Angs^3)

        /***********************************************************************
              * Explicit diffusion with variable kappa   (Convervative formulation)  *
        ************************************************************************/    
        l2[i][j][k].temp=tau/Ce*
        //first diffusion terms
        ( 
          // dK/dx * d^2 T/dx^2
        ( ((l1[i][j][k].fd_k+xmaxk)/2 * (xmaxTe-l1[i][j][k].temp)*invxsq)
         -((l1[i][j][k].fd_k+xmink)/2 * (l1[i][j][k].temp-xminTe)*invxsq)
         // dK/dy * d^2 T/dy^2
         +((l1[i][j][k].fd_k+ymaxk)/2 * (ymaxTe-l1[i][j][k].temp)*invysq)
         -((l1[i][j][k].fd_k+ymink)/2 * (l1[i][j][k].temp-yminTe)*invysq)                                 
#ifndef FDTD   //weil bisher nur 1D oder 2D
         //
         +((l1[i][j][k].fd_k+zmaxk)/2 * (zmaxTe-l1[i][j][k].temp)*invzsq)
         -((l1[i][j][k].fd_k+zmink)/2 * (l1[i][j][k].temp-zminTe)*invzsq)
#endif
         )
              //now coupling+source term
          -l1[i][j][k].fd_g*(l1[i][j][k].temp-l1[i][j][k].md_temp)
          +l1[i][j][k].source
        ) +l1[i][j][k].temp;
        

        //via enerie statt temp und cv : Nicht vergessen W/K/m^3 in IMDU = 7.3739e-22 eV/IMDt/eV/Angstrom^3
        //                               --> teile durch density --> eV/IMDt/eV/u
        //                               --> mal 26.9815u ---> eV/IMDt/eV/Atom 
        //                               --->Bekomme interne energie in eV/Atom       
        // l2[i][j][k].U = l1[i][j][k].U + tau / (l1[i][j][k].dens*6.02214076e-4)*26.9815* // (l1[i][j][k].dens / AMU / 1e30) *
        //                 (
        //                   ((l1[i][j][k].fd_k + xmaxk) / 2 * (xmaxTe - l1[i][j][k].temp) * invxsq)
        //                 - ((l1[i][j][k].fd_k + xmink) / 2 * (l1[i][j][k].temp - xminTe) * invxsq)
                          
        //                   - l1[i][j][k].fd_g * (l1[i][j][k].temp - l1[i][j][k].md_temp)
        //                   + l1[i][j][k].source
        //                 );
                             
        // l2[i][j][k].temp = EOS_te_from_r_ee(l1[i][j][k].dens, l2[i][j][k].U / (26.9815 * AMU * J2eV)) / 11604.5; // eV/atom -> J/kg

        l2[i][j][k].U=l1[i][j][k].U + (l2[i][j][k].temp-l1[i][j][k].temp)*Ce*fd_vol/((double) l1[i][j][k].natoms); // eV

        l1[i][j][k].xi += (l2[i][j][k].temp-l1[i][j][k].md_temp)*xi_fac*l1[i][j][k].fd_g/l1[i][j][k].md_temp/((double) l1[i][j][k].natoms);//Original
        //l1[i][j][k].xi += (l2[i][j][k].temp - l1[i][j][k].md_temp) * xi_fac * l1[i][j][k].fd_g / l1[i][j][k].md_temp / l1[i][j][k].dens; // NE
        l2[i][j][k].xi = l1[i][j][k].xi;

        
        /*
        if(l2[i][j][k].temp>l1[i][j][k].temp*10)  //kann auch passieren wenn zelle aktiviert wird
        {
        //unrealist. große Heizrate. Evtl. Cv zu klein weil Te zu klein!
          printf("DIFFPROBLEM:Heizrate viel zu groß!!:%.4e,\n"
              "steps:%d,atoms:%d\n"
              "myid:%d,i:%d,j:%d,atoms:%d,dens:%.4e\n"
                                "Te:%.4e,Te_old:%.4e,fd_g:%e,FD_C:%e,fd_k:%e,T_i:%f\n"
                                "Txmax:%f,Txmin:%f,Tymax:%f,Tymin:%f,Tzmax:%f,Tzmin:%f\n"
                                "kxmax:%f,kxmin:%f,kymax:%f,kymin:%f,kzmax:%f,kzmin:%f\n"
                                "\n",
              l1[i][j][k].source,
              steps,
              l1[i][j][k].natoms,
              myid,i_global,j_global,l1[i][j][k].natoms,l1[i][j][k].dens,
                                l2[i][j][k].temp,l1[i][j][k].temp,
                                l1[i][j][k].fd_g,
                                l1[i][j][k].Ce,l1[i][j][k].fd_k,l1[i][j][k].md_temp,
                                l1[xmax][j][k].temp,l1[xmin][j][k].temp,
                                l1[i][ymax][k].temp,l1[i][ymin][k].temp,
                                l1[i][j][zmax].temp,l1[i][j][zmin].temp,
                                l1[xmax][j][k].fd_k,l1[xmin][j][k].fd_k,
                                l1[i][ymax][k].fd_k,l1[i][ymin][k].fd_k,
                                l1[i][j][zmax].fd_k,l1[i][j][zmin].fd_k
                                );

        }
        */

#if DEBUG_LEVEL>0
        // NaN Temp oder <0
        if (l2[i][j][k].temp <= 0 || isnan(l2[i][j][k].temp) != 0)
        {
          printf("TEMP IS NAN or <0 IN DIFFLOOP!!!! steps:%d,proc:%d,i:%d,j:%d,k:%d,T:%.4e\n"
                 "dens: %.15e,atoms:%d\n"
                 "Q:%e,fd_g:%e,FD_C:%e,fd_k:%e,Te_old:%.4f,T_i:%f\n"
                 "Txmax:%f,Txmin:%f,Tymax:%f,Tymin:%f,Tzmax:%f,Tzmin:%f\n"
                 "kxmax:%f,kxmin:%f,kymax:%f,kymin:%f,kzmax:%f,kzmin:%f\n"
                 "fdx:%f,fdy:%f,fdz:%f\n"
                 "coupling:%.4e\n"
                 "Diffx:%.4e, Diffy:%.4e, Diffz:%.4e\n", steps,
                 myid, i, j, k, l2[i][j][k].temp,
                 l1[i][j][k].dens, l1[i][j][k].natoms,
                 l1[i][j][k].source, l1[i][j][k].fd_g,
                 l1[i][j][k].Ce, l1[i][j][k].fd_k, l1[i][j][k].temp, l1[i][j][k].md_temp,
                 l1[xmax][j][k].temp, l1[xmin][j][k].temp,
                 l1[i][ymax][k].temp, l1[i][ymin][k].temp,
                 l1[i][j][zmax].temp, l1[i][j][zmin].temp,
                 l1[xmax][j][k].fd_k, l1[xmin][j][k].fd_k,
                 l1[i][ymax][k].fd_k, l1[i][ymin][k].fd_k,
                 l1[i][j][zmax].fd_k, l1[i][j][zmin].fd_k,
                 fd_h.x, fd_h.y, fd_h.z,
                 -tau / Ce * l1[i][j][k].fd_g * (l1[i][j][k].temp - l1[i][j][k].md_temp),
                 tau / Ce * (((l1[i][j][k].fd_k + xmaxk) / 2 * (xmaxTe - l1[i][j][k].temp)*invxsq)
                             - ((l1[i][j][k].fd_k + xmink) / 2 * (l1[i][j][k].temp - xminTe)*invxsq)),

                 tau / Ce * (((l1[i][j][k].fd_k + ymaxk) / 2 * (ymaxTe - l1[i][j][k].temp)*invysq)
                             - ((l1[i][j][k].fd_k + ymink) / 2 * (l1[i][j][k].temp - yminTe)*invysq)),

                 tau / Ce * (((l1[i][j][k].fd_k + zmaxk) / 2 * (zmaxTe - l1[i][j][k].temp)*invzsq)
                             - ((l1[i][j][k].fd_k + zmink) / 2 * (l1[i][j][k].temp - zminTe)*invzsq))

                );
          error("Temp got NaN or <=0 in do_DIFF during fd-loop");
        }
#endif

      } // for k
    } // for j
  } //for i

  /* take care - l1 must always be the updated lattice */
  l3 = l1;
  l1 = l2;
  l2 = l3;

#if DEBUG_LEVEL>1
// if(DEBUG_LEVEL>0)
  printf("steps:%d,proc:%d,do_DIFF complete\n", steps, myid);
#endif

}

// *****************************************
// * ROUTINE ZUM AUSSCHREIBEN DER TTM-DATEN
// ******************************************
void ttm_writeout(int number)
{
  int n, nlocal;
  int i, j, k;
  ttm_Element * lglobal;
  ttm_Element * llocal;
  n = global_fd_dim.x * global_fd_dim.y * global_fd_dim.z;
  nlocal = (local_fd_dim.x - 2) * (local_fd_dim.y - 2) * (local_fd_dim.z - 2);

#ifdef MPIIO
  MPI_File fh;
  MPI_Status iostatus;
  int i_global, j_global, k_global;
  char mpifname[255];
  char singleline[TTMOUTBUFLEN];
  // create buffer for MPI_File_write_all

  char * ttmoutbuf = NULL;
  ttmoutbuf = (char*) calloc((local_fd_dim.x - 2) * (local_fd_dim.y - 2) * (local_fd_dim.z - 2) * TTMOUTBUFLEN, sizeof(char));
  //ttmoutbuf=(char*) malloc((local_fd_dim.x-2)*(local_fd_dim.y-2)*(local_fd_dim.z-2)*300*sizeof(char));
  if (NULL == ttmoutbuf)
  {
    error("couldn't alloc mem for ttmoutbuf\n");
  }

  char outbufline[TTMOUTBUFLEN]; //temporary line-buffer
  int  len_of_outbufline = 0;
  int  padding_len = 0;
  long l = 0;
  //int l=0;



  for (i = 1; i < local_fd_dim.x - 1; i++)
  {
    i_global =  ((i - 1) + my_coord.x * (local_fd_dim.x - 2));
    for (j = 1; j < local_fd_dim.y - 1; j++)
    {
      j_global = ((j - 1) + my_coord.y * (local_fd_dim.y - 2));
      for (k = 1; k < local_fd_dim.z - 1; k++)
      {
        k_global =  ((k - 1) + my_coord.z * (local_fd_dim.z - 2));
        //write to buffer
#ifndef FDTD
        sprintf(outbufline, "%d %d %d %d %e %e %e %e %e %e %e %e %e %e %e %d %f",
#else
        sprintf(outbufline, "%d %d %d %d %e %e %e %e %e %e %e %e %e %e %e %d %f %e %e %e %e %e %e %e %e %e %e",
#endif
                i_global, j_global, k_global, l1[i][j][k].natoms, l1[i][j][k].temp,
                l1[i][j][k].md_temp, l1[i][j][k].xi,
                l1[i][j][k].source, l1[i][j][k].dens,
                l1[i][j][k].vcomx, l1[i][j][k].vcomy, l1[i][j][k].vcomz,
                l1[i][j][k].fd_k, l1[i][j][k].fd_g,
#ifndef FDTD
                l1[i][j][k].Z, l1[i][j][k].proc, l1[i][j][k].Ce
#else
                l1[i][j][k].Z, l1[i][j][k].proc, l1[i][j][k].Ce,
                l1[i][j][k].Ezx, l1[i][j][k].Ezy, l1[i][j][k].Hx, l1[i][j][k].Hy,
                l1[i][j][k].sigmax, l1[i][j][k].sigmay,
                l1[i][j][k].Hzx, l1[i][j][k].Hzy, l1[i][j][k].Ex, l1[i][j][k].Ey
#endif
               );

        len_of_outbufline = strlen(outbufline);
        len_of_outbufline++; //wegen linebreak
        if (len_of_outbufline > TTMOUTBUFLEN)
          error("ttm out-buf too small. Increase buffer!\n");
        /*
          padding_len=300-len_of_outbufline-2; //-1 wegen newline
          if(padding_len>0)
          {
           if(sprintf(&ttmoutbuf[l],"%s%*s\n",outbufline,padding_len,"")<0) //trailing padding
               error("err in sprintf in ttm_writeout (1)");
          }
          else if(padding_len==0)
          {
            if(sprintf(&ttmoutbuf[l],"%s\n",outbufline)<0)
              error("err in sprintf in ttm_writeout (2)");
          }
          else
             error("outbufline has more than 300 characters. Increase buffer!\n");
          l+=300;
        */
        if (sprintf(&ttmoutbuf[l], "%s\n", outbufline) < 0)
          error("err in sprintf in ttm_writeout (2)");
        l += len_of_outbufline;
      }
    }
  }

  MPI_Barrier(cpugrid);
  sprintf(mpifname, "%s.%d.ttm", outfilename, number);
  MPI_File_open(cpugrid, mpifname,
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  //Header schreiben
  int headofs;
  if (myid == 0)
  {
#ifndef FDTD
    sprintf(outbufline,
            "#x y z natoms temp md_temp xi source dens vx vy vz fd_k fd_g Z proc Ce");
#else
    sprintf(outbufline,
            "#x y z natoms temp md_temp xi source dens vx vy vz fd_k fd_g Z proc Ce Ezx Ezy Hx Hy sigmax sigmay Hzx Hzy Ex Ey");
#endif
//    padding_len=300-strlen(outbufline)-1;
//    sprintf(singleline,"%s%*s\n",outbufline,padding_len,"");
    sprintf(singleline, "%s\n", outbufline);
    MPI_File_write_at(fh, 0, singleline, strlen(singleline), MPI_CHAR, &iostatus);
    headofs = strlen(singleline);
  }
  MPI_Bcast(&headofs, 1, MPI_INT, 0, cpugrid);

  //Daten schreiben
//  long myofs=(long) 300+myid*300*(local_fd_dim.x-2)*(local_fd_dim.y-2)*(local_fd_dim.z-2);
  long myofs = 0;
  MPI_Exscan(&l, &myofs, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD); //undefined on proc 0

  myofs += (long) headofs;
  MPI_File_set_view(fh, myofs, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

//  MPI_File_write_all(fh, ttmoutbuf, 300*(local_fd_dim.x-2)*(local_fd_dim.y-2)*(local_fd_dim.z-2), MPI_CHAR, &iostatus);
  MPI_File_write_all(fh, ttmoutbuf, l, MPI_CHAR, &iostatus);

  MPI_File_close(&fh);
  //////////////////////
  //free buffer array
  /////////////////////
  if (ttmoutbuf == NULL)
  {
    error("can't free a NULL-pointer (ttm_writeout)");
  }
  else
    free(ttmoutbuf);
  ttmoutbuf = NULL;
  return;
#endif  //MPIIO

#ifdef DEBUG
  assert(nlocal == n / num_cpus);
#endif

#ifdef MPI2
  MPI_Alloc_mem ( nlocal * sizeof(ttm_Element), MPI_INFO_NULL, &llocal );
#else
  llocal = malloc(nlocal * sizeof(ttm_Element));
#endif //MPI2

  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {
    for (j = 1; j < local_fd_dim.y - 1; ++j)
    {
      for (k = 1; k < local_fd_dim.z - 1; ++k)
      { /* all the "-1" and "-2" because we don't store ghost layers
           and don't want to waste the space in llocal */
        llocal[ (i - 1) * (local_fd_dim.y - 2) * (local_fd_dim.z - 2)
                + (j - 1) * (local_fd_dim.z - 2)
                + k - 1 ] = l1[i][j][k];
      }
    }
  }

#ifdef MPI
  if (myid == 0)
  {
#ifdef MPI2
    MPI_Alloc_mem ( n * sizeof(ttm_Element), MPI_INFO_NULL, &lglobal );
#else
    lglobal = malloc (n * sizeof(ttm_Element));
#endif /*MPI2*/
  }

  /* Note: This only works because mpi_element2 includes the last element
   * of ttm_Element (v_com.z). If this is changed, you need to work with
   * displacements or something (e.g. MPI standard, Ex. 4.5) */
  MPI_Gather( llocal, nlocal, mpi_element2,
              lglobal, nlocal, mpi_element2, 0, cpugrid );

#else /* no MPI */
  lglobal = llocal;
#endif /* MPI */

  if (myid == 0)
  {
    FILE *outfile;
    char fname[255];
    sprintf(fname, "%s.%d.ttm", outfilename, number);
    outfile = fopen(fname, "w");
    if (NULL == outfile) error ("Cannot open ttm file for writing.\n");
#ifndef FDTD
    fprintf(outfile,
            "#x y z natoms temp md_temp xi source dens vx vy vz fd_k fd_g Z proc Ce\n");
#else
    fprintf(outfile,
            "#x y z natoms temp md_temp xi source dens vx vy vz fd_k fd_g Z proc Ce Ezx Ezy Hx Hy sigmax sigmay Hzx Hzy Ex Ey\n");
#endif
    for (i = 0; i < global_fd_dim.x; ++i)
    {
      for (j = 0; j < global_fd_dim.y; ++j)
      {
        for (k = 0; k < global_fd_dim.z; ++k)
        {
          int index;

#ifdef MPI
          ivektor from_process; /* we need to look in the data from
           the process with these grid coords */
          /* data from the processes is sorted by rank. */
          from_process.x = i / (local_fd_dim.x - 2);
          from_process.y = j / (local_fd_dim.y - 2);
          from_process.z = k / (local_fd_dim.z - 2);
          index = cpu_grid_coord(from_process) * nlocal +
                  (i % (local_fd_dim.x - 2)) * (local_fd_dim.y - 2) * (local_fd_dim.z - 2) +
                  (j % (local_fd_dim.y - 2)) * (local_fd_dim.z - 2) +
                  (k % (local_fd_dim.z - 2));
#else /* no MPI */
          index =
            i * (local_fd_dim.y - 2) * (local_fd_dim.z - 2) +
            j * (local_fd_dim.z - 2) +
            k;
#endif /* MPI*/

#ifndef FDTD
          fprintf(outfile, "%d %d %d %d %e %e %e %e %e %e %e %e %e %e %e %d %f\n",
#else
          fprintf(outfile, "%d %d %d %d %e %e %e %e %e %e %e %e %e %e %e %d %f %e %e %e %e %e %e %e %e %e %e\n",
#endif
                  i, j, k, lglobal[index].natoms, lglobal[index].temp,
                  lglobal[index].md_temp, lglobal[index].xi,
                  lglobal[index].source, lglobal[index].dens,
                  lglobal[index].vcomx, lglobal[index].vcomy, lglobal[index].vcomz,
                  lglobal[index].fd_k, lglobal[index].fd_g,
#ifndef FDTD
                  lglobal[index].Z, lglobal[index].proc, lglobal[index].Ce
#else
                  lglobal[index].Z, lglobal[index].proc, lglobal[index].Ce,
                  lglobal[index].Ezx, lglobal[index].Ezy, lglobal[index].Hx, lglobal[index].Hy,
                  lglobal[index].sigmax, lglobal[index].sigmay,
                  lglobal[index].Hzx, lglobal[index].Hzy, lglobal[index].Ex, lglobal[index].Ey
#endif
                 );
        }
      }
    }
    fclose(outfile);
  }

#ifdef MPI
  if (myid == 0)
  {
#ifdef MPI2
    MPI_Free_mem (lglobal);
#else
    free (lglobal);
#endif
  }
#endif

#ifdef MPI2
  MPI_Free_mem (llocal);
#else
  free(llocal);
#endif

}





/************************
 * MY MOD: TTM READ for *
*  restart              *
*  ACHTUNG: bisher nur  *
*  fuer MPI-nutzung     *
*  implementiert        *
*************************/
void ttm_read(int number)
{
  int i, j, k;
  int lines = global_fd_dim.x * global_fd_dim.y * global_fd_dim.z;
  ttm_Element *buf;

#ifdef MPI
#ifdef MPI2
  MPI_Alloc_mem(lines * sizeof(ttm_Element), MPI_INFO_NULL, &buf);
#else
  buf = malloc((lines * sizeof(ttm_Element)));
#endif
#endif

  if (myid == 0)
  {
    char fname[255];
    sprintf(fname, "%s.%d.ttm", outfilename, number);
    FILE *infile = fopen(fname, "r");
    if (infile == NULL)
    {
      error_str("Cannot open input ttm-file:%s\n", fname);
      //error("Cannot open input ttm-file");
    }
    fseek (infile, 0, SEEK_SET);
    printf("\nReading %d lines from ttm-file:%s\n", lines, fname);

//    char line[256];
    char line[400];
    char **tokens;
    size_t j, numtokens;

    int ig, jg, kg;

    for (i = 0; i < lines + 1; i++) //skip first line (comments)
    {
      //read data
      if (fgets (line, 400, infile) == NULL) {
        printf("error reading ttm-input-file %s in line %d\n", fname, i);
        error("error reading ttm input file...");
      }
      if (i > 0) //skip first line
      {
        tokens = strsplit(line, ", \t\n", &numtokens);

        sscanf(tokens[3], "%d",  &buf[i - 1].natoms);
        sscanf(tokens[4], "%lf", &buf[i - 1].temp);
        sscanf(tokens[5], "%lf", &buf[i - 1].md_temp);
        sscanf(tokens[6], "%lf", &buf[i - 1].xi);
        sscanf(tokens[7], "%lf", &buf[i - 1].source);
        sscanf(tokens[8], "%lf", &buf[i - 1].dens);
//  sscanf(tokens[9],"%lf", &buf[i-1].v_com.x);
//  sscanf(tokens[10],"%lf", &buf[i-1].v_com.y);
//  sscanf(tokens[11],"%lf", &buf[i-1].v_com.z);
        sscanf(tokens[12], "%lf", &buf[i - 1].fd_k);
        sscanf(tokens[13], "%lf", &buf[i - 1].fd_g);
        sscanf(tokens[14], "%lf", &buf[i - 1].Z);
        sscanf(tokens[15], "%d", &buf[i - 1].proc);
        sscanf(tokens[16], "%lf", &buf[i - 1].Ce);
#ifdef FDTD
        sscanf(tokens[17], "%lf", &buf[i - 1].Ezx);
        sscanf(tokens[18], "%lf", &buf[i - 1].Ezy);
        sscanf(tokens[19], "%lf", &buf[i - 1].Hx);
        sscanf(tokens[20], "%lf", &buf[i - 1].Hy);
        sscanf(tokens[21], "%lf", &buf[i - 1].sigmax);
        sscanf(tokens[22], "%lf", &buf[i - 1].sigmay);
        sscanf(tokens[23], "%lf", &buf[i - 1].Hzx);
        sscanf(tokens[24], "%lf", &buf[i - 1].Hzy);
        sscanf(tokens[25], "%lf", &buf[i - 1].Ex);
        sscanf(tokens[26], "%lf", &buf[i - 1].Ey);
#endif

        //printf("%d %d %d %d %d %d\n",i-1,ig,jg,kg,buf[i-1].natoms,buf[i-1].proc);

        for (j = 0; j < numtokens; j++) {
          free(tokens[j]);
        }
        if (tokens != NULL)
          free(tokens);
      }
    }
    fclose(infile);
  } //endif myid==0

  MPI_Bcast(buf, lines, mpi_element2, 0, MPI_COMM_WORLD);

  int l = 0;
  //int ig,jg,kg;
  i = j = k = 1;
  for (l = 0; l < lines; l++)
  {
    //ig =  ((i-1) + my_coord.x*(local_fd_dim.x-2));
    //jg =  ((j-1) + my_coord.y*(local_fd_dim.y-2));
    //kg =  ((k-1) + my_coord.z*(local_fd_dim.z-2));
    //printf("l:%d,proc:%d,myid:%d,i:%d,j:%d,k:%d\n",l,buf[l].proc,myid,i,j,k);
    if (buf[l].proc == myid)
    {
      l1[i][j][k].temp = buf[l].temp;
      l1[i][j][k].natoms = buf[l].natoms;
      l1[i][j][k].md_temp = buf[l].md_temp;
      l1[i][j][k].xi = buf[l].xi;
      l1[i][j][k].source = buf[l].source;
      l1[i][j][k].dens = buf[l].dens;
      l1[i][j][k].fd_k = buf[l].fd_k;
      l1[i][j][k].fd_g = buf[l].fd_g;
      l1[i][j][k].Z = buf[l].Z;
      l1[i][j][k].proc = buf[l].proc;
      l1[i][j][k].Ce = buf[l].Ce;
      l1[i][j][k].vcomx = buf[l].vcomx;
      l1[i][j][k].vcomy = buf[l].vcomy;
      l1[i][j][k].vcomz = buf[l].vcomz;
#ifdef FDTD
      l1[i][j][k].Ezx = buf[l].Ezx;
      l1[i][j][k].Ezy = buf[l].Ezy;
      l1[i][j][k].Hx = buf[l].Hx;
      l1[i][j][k].Hy = buf[l].Hy;
      l1[i][j][k].eps = eps0; //buf[l].eps;
      l1[i][j][k].mu = mu0; //buf[l].mu;
      l1[i][j][k].sigmax = buf[l].sigmax;
      l1[i][j][k].sigmay = buf[l].sigmay;
      l1[i][j][k].Hzx = buf[l].Hzx;
      l1[i][j][k].Hzy = buf[l].Hzy;
      l1[i][j][k].Ex = buf[l].Ex;
      l1[i][j][k].Ey = buf[l].Ey;
#endif

      k++;
      if (k > local_fd_dim.z - 2)
      {
        k = 1;
        j++;
      }
      if (j > local_fd_dim.y - 2)
      {
        j = 1;
        i++;
      }
      if (i > local_fd_dim.x - 2)
      {
        i = 1;
      }
    }
  }
#ifdef MPI
#ifdef MPI2
  MPI_Free_mem(buf);
#else
  free(buf);
#endif
#endif

//*****************CHECKEN OB KORREKT EINGELESEN WURDE (DEBUG PURPOSE)########
  /*
    int i_global,k_global,j_global;
    for (i=1; i<local_fd_dim.x-1; ++i)
    {
      i_global =  ((i-1) + my_coord.x*(local_fd_dim.x-2));
      for (j=1; j<local_fd_dim.y-1; ++j)
      {
        j_global =  ((j-1) + my_coord.y*(local_fd_dim.y-2));
        for (k=1; k<local_fd_dim.z-1; ++k)
        {
          k_global =  ((k-1) + my_coord.z*(local_fd_dim.z-2));
          if(i_global==110 && j_global==16 && k_global==21)
          {
            printf("DEBUG:natoms:%d\n",l1[i][j][k].natoms);
          }
        }
      }
    }
  */
}

// ******************************************************************************
// *         AUXILIARY FUNCTIONS FOR WIDE-RANGE PROPERTIES    *
// ******************************************************************************/

// double Cv(double Te, double ne) //specific heat of quasi-free electrons,Te in eV, ne in m^-3 //deprecated
// {
//   //approximatin according to mazhukin, S. 240 (Fermi-integral analyt. genähert)
//   //error less than 5%
//   double EF = fermi_E(ne);
//   //double Te_K=Te*11605; //in Kelvin
//   double Te_J = Te * 1.6021766e-19; //in Joule
//   //return 1.5*ne*BOLTZMAN*BOLTZMAN*Te_K/sqrt(Te_J*Te_J+pow(3*EF/M_PI/M_PI,2.0))*7.243227582e-8; //J/K/m^3 -> IMD-UNITS
//   return 2.401087548821963e-49 * Te * ne / sqrt(Te_J * Te_J + pow(EF * 0.303963550927013, 2.0)); //alle konstanten zus.gefasst

// }

/*
double nu_e_e(double Te, double EF, double Ne, double Na, double valence) //collision freq.Mazkhukin  //deprecated
{
  //chemical valence=nr of conduction electron per atom
  double l_e_e = 1.0 / Ne / sigma_e_e(Te, EF, Na, valence); //mean free path
  double collfreq = v_e(Te, EF) / l_e_e;
  return collfreq;
}

double v_e(double Te, double EF) //Thermal velocity of electrons, approx. by mazhukin //inactive
{
  double eta = BOLTZMAN * Te / EF; //reduced temp.
  return sqrt(3.0 * EF / EMASS * sqrt(eta * eta + 0.16));
}

// ********************************************* //
//    CROSS SECTIONS NACH MAZHUKIN
// ********************************************* //
// DEPRECATED ...
double sigma_e_ph(double Te, double Ti, double EF, double Na, double Z) //inactive
{
  double r = pow(3.0 / 4.0 / M_PI / Na, 1.0 / 3.0); //Wigner-Seitz radius
  double eta = BOLTZMAN * Te / EF;
  double eta_sq = eta * eta;
  double t = pow(9.0 * M_PI / 4.0, 4.0 / 3.0) * pow(Z, 1.0 / 3.0) * BOHR / r * sqrt((eta_sq + 0.16) * (eta_sq + 4.0 / 9.0));
  double mAl = 26.9815 * AMU;
  double sigma = 4.0 * M_PI * pow(4.0 / 9.0 / M_PI, 4.0 / 3.0) * mAl * BOLTZMAN * Ti / HBAR / HBAR * pow(r / BOHR, 2.0) * pow(r, 4.0) / pow(Z, 4.0 / 3.0);
  sigma = sigma * eta * (log(t + 1) - t / (t + 1)) / sqrt(eta_sq + 0.16) / sqrt(eta_sq + 4.0 / 9.0);
  return sigma; //~6e-22
}

double sigma_e_e(double Te, double EF, double Na, double valence) //Na=atomic density in 1/m^3 //inactive
{
  //electron-electron elastic scattering cross section
  //Mazhukin S.241
  double r = pow(3.0 / 4.0 / M_PI / Na, 1.0 / 3.0); //Wigner-Seitz radius
  double rsq = r * r;
  double rsqsq = rsq * rsq;
  double eta = BOLTZMAN * Te / EF; //reduced temperature
  double eta_sq = eta * eta;
  double BOHR_sq = BOHR * BOHR;
  double t = pow(9.0 * M_PI / 4.0, 4.0 / 3.0) * pow(valence, 1.0 / 3.0) * BOHR / r * sqrt((eta_sq + 0.16) * (eta_sq + 4.0 / 9.0));
  double sigma = 2.0 * M_PI / 9.0 * pow(4.0 / 9.0 / M_PI, 4.0 / 3.0) * pow(valence, -4.0 / 3.0) * rsqsq / BOHR_sq * eta_sq * (log(t + 1.0) - t / (t + 1.0)) / ((eta_sq + 4.0 / 9.0) * (eta_sq + 0.16));
  return sigma;
}

*/

double fermi_E(double ne) {                             //IN: [1/m^3] OUT: [J]

  return HBAR * HBAR * pow(3.0 * M_PI * M_PI * ne, 2.0 / 3.0) / 2.0 / EMASS;
}
double fermi_T(double ne) {                             //IN: [1/m^3] OUT: [K]
  return 2.0 * fermi_E(ne) / (3.0 * BOLTZMAN);
}

double omega_pl(double ne) {                            //IN: [1/m^3] OUT: [1/s]
  return sqrt(ne * ECHARGE * ECHARGE / EMASS / ECONST);
}
double bMin(double X, double omega_las, double Z, double ni, double Te) {       // OUT: [m]
  return MAX(Z * ECHARGE * ECHARGE / (BOLTZMAN * Te),
             HBAR / sqrt(2 * EMASS * BOLTZMAN * Te));
}
double bMax(double omega_las, double Z, double ni, double Te) {                 // OUT: [m]
  double teff = sqrt(Te * Te + fermi_T(ni * Z) * fermi_T(ni * Z));
  return MAX(r0(ni), sqrt(BOLTZMAN * teff / (EMASS)) / MAX(omega_las, omega_pl(ni * Z)));
  //ACHTUNG: Sobald laser aus --> omega_las=0
}
double coulomb_log(double X, double omega_las, double Z, double ni, double Te) { //OUTPUT in dimensionless^M
  return MAX(2.0, log(bMax(omega_las, Z, ni, Te) / bMin(X, omega_las, Z, ni, Te)));
}

double ne_cr(double omega_las) {                // OUT: [1/m3]
  return ECONST * EMASS * omega_las * omega_las / ECHARGE / ECHARGE;

  // aus polly-2t
  // omega_las * omega_las * ELECTRON_MASS / (ELECTRON_CHARGE * ELECTRON_CHARGE  * 4 * PI);
}


double fd_density(double natoms, double mass, double vol) { //OUTPUT in kg/m3
  return natoms    // number of atoms in the cell (can be gained by "ttmelement".natoms)
         *mass    // mass
         *AMU     // atomic mass unit
         / (vol * 1e-30);

}
double r0(double ni) {                                  //IN: [1/m3] OUT: [m]
  return pow(3.0 / (4 * M_PI * ni), 1.0 / 3.0);
}

/***********************************************************************
 * Mean-charge nach interpolationsformel wie in POLLY-2T (More-Fit)
 ***********************************************************************/
double MeanCharge(double temp, double rho, double z0, double am, int i, int j, int k) {

  if (isnan(temp) != 0)
  {
    printf("proc:%d,i:%d,j:%d,k:%d,rho:%f\n",
           myid, i, j, k, rho);
    error("temp is NaN (from MeanCharge)");
  }
  if (isnan(rho) != 0)
    error("rho is NaN (from MeanCharge)");
  if (isnan(z0) != 0)
    error("z0 is NaN (from MeanCharge)");

  temp *= 1e-3;           //Electron Temperature [K] -> [kK] //nicht nötig
  rho *= 1e-3;            //Density [kg/m3] -> [g/cm3]
  /*
          z0:             atom number=13 für Alu
          am:             atomic mass number=26.98 u für Al
          t0:             temperature(eV) /z0**(4/3)
          r0:             rho(g/ccm)/(z0*am)
  */
  double YA1 =  0.003323;
  double YA2 =  0.971832;
  double YA3 =  0.926148e-5;
  double YA4 =  3.10165;
  double YB0 =  -1.7630;
  double YB1 =  1.43175;
  double YB2 =  0.315463;
  double YC1 =  0.036666667;
  double YC2 =  0.983333333;
  double YAL =  0.143139e2;
  double YBE =  0.6624e0;

  double t0 = temp / pow(atomic_charge, 4. / 3.); //
  if (isnan(t0) != 0)
    error("t0 is NaN");
  double r0 = rho / (atomic_charge * am);
  if (isnan(r0) != 0)
    error("r0 is NaN");
  double ytf = t0 / (1.0 + t0);
  if (isnan(ytf) != 0)
  {
    printf("proc:%d,ytf:%e,t0:%e,atomic_charge:%f,atomic_weight:%f,temp:%f\n", myid, ytf, t0, atomic_charge, atomic_weight, temp);
    error("ytf is NaN");
  }
  double ya = YA1 * pow(t0, YA2) + YA3 * pow(t0, YA4);
  if (isnan(ya) != 0) {
    printf("YA1:%e,t0:%e,YA2:%e,YA3:%e,YA4:%e\n", YA1, t0, YA2, YA3, YA4);
    error("ya is NaN");
  }
  double yb = -exp(YB0 + YB1 * ytf + YB2 * ytf * ytf);
  if (isnan(yb) != 0)
    error("yb is Nan");
  double yc = YC1 * ytf + YC2;
  if (isnan(yc) != 0)
    error("yc is NaN");
  double yq1 = ya * pow(r0, yb);
  if (isnan(yq1) != 0)
    error("yq1 is Nan");
  double yq = pow(pow(r0, yc) + pow(yq1, yc), 1.0 / yc);
  if (isnan(yq) != 0)
    error("yq is NaN");
  double yx = YAL * pow(yq, YBE);
  if (isnan(yx) != 0)
    error("yx is NaN");
  double retval = z0 * yx / (1.0 + yx + sqrt(1.0 + 2.0 * yx));
  if (isnan(retval) != 0 && steps > 0)
  {
    //Kommt vor,wenn yx ist inf
    //retval=0.001;
    error("retval is NaN (MeanCharge)");
  }
  return retval;

}
/**************************************
*
* Effective collision frequencies
* in different regimes
*
***************************************/
/*
double nueff(double Te, double Ti, double ne, double ni, double ni0, double Ce, double lod) //Kombination aus Petrov & Mazhukin (nu_ee_sol) //inactive
{
  //all SI units
  double nu_ei_sol, nu_ei_liq, nu_ei_pl, nu_ei_wr;
  double nu_ee, nu_ee_pl, nu_ee_wr;
  double kappa_pl;
  double ve_sq;
  double Z;
  double prefac = 2.908744629161452e-32;
  double denom;
  double coulomb;
  double ratio;
  double nu_ei_sum;
  double EF;
  double eta;

  Z = ne / ni;
  coulomb = coulomb_log(1.5, omega_laser, Z, ni, Te);
  ratio = pow(ni / ni0, -1.3);

  nu_ei_sol = 3.2e11 * Ti; //@ rho0 bzw ni0
  nu_ei_liq = 0.72e14 / (130 + 0.0367 * Ti - 66700 / Ti);
  //ve_sq=2.0/EMASS*U; //U=internal eng.
  EF = fermi_E(ne);
  eta = BOLTZMAN * Te / EF; //reduced temp.
  ve_sq = 3.0 * EF / EMASS * sqrt(eta * eta + 0.16);

  //kappa_pl=prefac*0.7*BOLTZMAN*pow(BOLTZMAN*Te,2.5)/denom; <--bullshit
  kappa_pl = 0.7 * prefac * pow(BOLTZMAN * Te, 5.0 / 2.0) / (Z * pow(ECHARGE, 4.0) * coulomb);

  nu_ei_pl = 0.7 * Ce * ve_sq / 3 / kappa_pl; //Ce=specific heat (volumetric)
  nu_ei_sum = ratio * (lod * nu_ei_sol + (1 - lod) * nu_ei_liq);
  nu_ei_wr = pow(pow(nu_ei_sum, -2.0) + pow(nu_ei_pl, -2.0), -0.5);

  nu_ee = 1.0 / (0.2474e15 * Te * Te) + 1.0 / (1.15e15 * pow(Te, 0.28));
  nu_ee = 1.0 / nu_ee * ratio;
  nu_ee_pl = (1 - 0.7) * Ce * ve_sq / 3 / kappa_pl;
  nu_ee_wr = pow(pow(nu_ee, -2.0) + pow(nu_ee_pl, -2.0), -0.5);

  return nu_ee_wr;
}
*/

double numet(double A1, double A2, double Te, double Ti, double TF) //POVAR
{
  //IN: [1],[1],[K],[K],[K]
  //OUT: [1/s]
  //return (A1 * BOLTZMAN * Ti + A2 * BOLTZMAN * Te*Te / TF) / HBAR;
  //return BOLTZMAN*(A1*Ti+A2*Te*Te/TF)/HBAR;
  return 1.309202957843054e+11 * (A1 * Ti + A2 * Te * Te / TF);
}

double numax(double A3, double vf, double Te, double ni) //POVAR
{
  // INT: [1], [m/s],[K],[m^-3]
  //OUT: [1/s]
  //double r = r0(ni);
  return A3 / r0(ni) * sqrt(vf * vf + BOLTZMAN * Te / EMASS);
}
double nupl(double omega_las, double Z, double ni, double ne, double Te) //POVAR
{
  // IN: [1/s],[1],[m^-3],[m^-3],[K]
  //OUT: [1/s]
  double coulomb = coulomb_log(1.5, omega_las, Z, ni, Te);
  //return sqrt(2.0 / M_PI / EMASS) / 12.0 / M_PI / ECONST/ECONST * Z *ne * pow(ECHARGE, 4.0) * coulomb / pow(BOLTZMAN * Te, 3.0/2.0);
  return 1.863839656495274e-40 * Z * ne * coulomb / pow(BOLTZMAN * Te, 1.5);
}

/**************************************
* Kappa und Gamma nach Povar
***************************************/
double getGamma(double Te, double Ti, double Ne, double Z) //POVAR
{

//return 1.785961556211823e-04; //TESTCASE
  //ACHTUNG: gamma ist nach Povarnitsyn die Kopplungskonstante und NICHT der Wärmekap. Koeff

  Te *= 11604.5; //ECHARGE / BOLTZMAN;              // [eV] -> [K]
  Ti *= 11604.5; //ECHARGE / BOLTZMAN;              // [eV] -> [K]

  double A1g = 50;
  double A2g = 20;
  double A3g = 0.25;

  double  m_atom = 26.9815; //u

  double Ni = Ne / Z;                                     // [1/m3]
  double EF = fermi_E(Ne);                                // [J]
  //double TF = 2.0 * EF / (3.0 * BOLTZMAN);                // [K]
  double TF = 4.828648689433765e+22 * EF;
  //double VF = sqrt(2 * EF / EMASS);                       // [m/s]
  double VF = 1.481734876966785e+15 * sqrt(EF);
  double nu_eff = MIN(numet(A1g, A2g, Te, Ti, TF), MIN(numax(A3g, VF, Te, Ni), nupl(omega_laser, Z, Ni, Ne, Te)));
  //double gamma = 3.0 * BOLTZMAN * EMASS / (m_atom*AMU) * Ne * nu_eff;     // [W/(K*m^3)]
  double gamma = 2.783313120645740e-74 / (m_atom * AMU) * Ne * nu_eff; // IMDu  
  return gamma;
}

/*
double KappaInterpol(double rho, double Te, double Ti)
{
  double lgTe = log10(Te);
  double lgTi = log10(Ti);
  double Kappa = do_tricubinterp(&kappa_interp, rho, lgTe, lgTi);
  if (Kappa == -1) return -1;
  return Kappa;
}
double NueffInterpol(double rho, double Te, double Ti)
{
  double lgTe = log10(Te);
  double lgTi = log10(Ti);
  double nueff = do_tricubinterp(&nueff_interp, rho, lgTe, lgTi);
  if (nueff == -1) return -1;
  return nueff;
}
*/

double getKappa(double Te, double Ti, double Ne, double Z) //POVAR
{
//return 1.933442e+01;  //TESTCASE
  //check for NaNs

  //Te *= ECHARGE / BOLTZMAN;             // [eV] -> [K]
  //Ti *= ECHARGE / BOLTZMAN;             // [eV] -> [K]

  Te *= 11604.5; //ev->K
  Ti *= 11604.5; //ev->K

  double A1t = 2.95;
  double A2t = 0.5;
  double A3t = 0.16;
  double A4t = 1.2;

  //const double m_atom = 26.9815;

  //double rho = fd_density(natoms,m_atom,vol);       // [kg/m3]
  //double omega_las = (2.0 * M_PI * LIGHTSPEED / lambda);  // [1/s]
  //double Ne = Z / (1/rho * m_atom * AMU); // [1/m3]

  double Ni = Ne / Z;                     // [1/m3]
  double EF = fermi_E(Ne);                // [J]

  //double TF = 2.0 * EF / (3.0 * BOLTZMAN);  // [K]
  double TF = 4.828648689433765e+22 * EF;

  //double VF = sqrt(2 * EF / EMASS);         // [m/s]
  double VF = 1.481734876966785e+15 * sqrt(EF);

  double nu_eff = MIN(numet(A1t, A2t, Te, Ti, TF), numax(A3t, VF, Te, Ni));
  double coulomb = coulomb_log(1.5, omega_laser, Z, Ni, Te);
  //double kappa_met = M_PI*M_PI * BOLTZMAN*BOLTZMAN * Ne * Te/ (3.0 * EMASS * nu_eff);
  double kappa_met = 6.884236239621913e-16 * Ne * Te / nu_eff;

  //double kappa_pl = sqrt(2.0 / pow(M_PI, 7) / EMASS) * ECONST*ECONST * BOLTZMAN * pow(BOLTZMAN * Te, 5.0/2.0)
  //                /(Z * pow(ECHARGE, 4.0) * coulomb);

  double kappa_pl = 4.428788911416808e+43 * pow(BOLTZMAN * Te, 2.5) / Z / coulomb;
  double kappa_wr = kappa_pl + (kappa_met - kappa_pl) * exp(-A4t * Te / TF); // [W/(K*m)]
  //return kappa_wr * 10.18e-15 / (1e10 * BOLTZMAN);                // Convert from SI -> IMDu
  return kappa_wr * 0.073768115942029;
}

/*************************************************************
* MY MOD: Auxiliary function for string Parsing in TTM-Read  *
**************************************************************/
char **strsplit(const char* str, const char* delim, size_t* numtokens) {
  // copy the original string so that we don't overwrite parts of it
  // (don't do this if you don't need to keep the old line,
  // as this is less efficient)
  char *s = strdup(str);
  // these three variables are part of a very common idiom to
  // implement a dynamically-growing array
  size_t tokens_alloc = 1;
  size_t tokens_used = 0;
  char **tokens = calloc(tokens_alloc, sizeof(char*));
  char *token, *strtok_ctx;
  for (token = strtok_r(s, delim, &strtok_ctx);
       token != NULL;
       token = strtok_r(NULL, delim, &strtok_ctx)) {
    // check if we need to allocate more space for tokens
    if (tokens_used == tokens_alloc) {
      tokens_alloc *= 2;
      tokens = realloc(tokens, tokens_alloc * sizeof(char*));
    }
    tokens[tokens_used++] = strdup(token);
  }
  // cleanup
  if (tokens_used == 0) {
    free(tokens);
    tokens = NULL;
  } else {
    tokens = realloc(tokens, tokens_used * sizeof(char*));
  }
  *numtokens = tokens_used;
  free(s);
  return tokens;
}
// ************************************************
// * Determine max. time-step from CFL-criterion
// *********************************************
void CFL_maxdt(void)
{
  //vars for max diffusion timestep
  int imax, jmax, kmax;
  int i, j, k;
  double maxdttmp = 1e9; // temp.var for comp. of max. allowed ttm-timestep
  double dxsq = fd_h.x * fd_h.x;
  double dysq = fd_h.y * fd_h.y; // maxdt <= 0.5*dx^2/k_(i+1/2)*C (numerical recp. s. 1045)
  double dzsq = fd_h.z * fd_h.z;
  double khalf;

  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {
    for (j = 1; j < local_fd_dim.y - 1; ++j)
    {
      for (k = 1; k < local_fd_dim.z - 1; ++k)
      {
        if (l1[i][j][k].natoms > fd_min_atoms)
        {
          // ************************************************
          // * DETERMINE MAX. TIME-STEP FOR STABILITY   *
          // ************************************************
          //x-dir
          if (l1[i + 1][j][k].natoms >= fd_min_atoms) imax = i + 1; else imax = i;
          khalf = (l1[i][j][k].fd_k + l1[imax][j][k].fd_k); //k_(i+1/2)*2
          maxdttmp = MIN(maxdttmp, l1[i][j][k].Ce * dxsq / khalf);

          //y-dir
          if (l1[i][j + 1][k].natoms >= fd_min_atoms) jmax = j + 1; else jmax = j;
          khalf = (l1[i][j][k].fd_k + l1[i][jmax][k].fd_k);
          maxdttmp = MIN(maxdttmp, l1[i][j][k].Ce  * dysq / khalf);

          //z-dir
          if (l1[i][j][k + 1].natoms >= fd_min_atoms) kmax = k + 1; else kmax = k;
          khalf = (l1[i][j][k].fd_k + l1[i][j][kmax].fd_k);
          maxdttmp = MIN(maxdttmp, l1[i][j][k].Ce * dzsq / khalf);
        }
      }
    }
  }

  //COMMINICATE GLOBALLY TTM-TIMESTEP
  if (steps > 0)
  {
    MPI_Allreduce(&maxdttmp, &max_dt_ttm, 1, MPI_DOUBLE, MPI_MIN, cpugrid);
  }
  else
  {
    max_dt_ttm = timestep / ((double) fd_n_timesteps);
  }
  max_dt_ttm *= 0.9; //sicherheits-faktor

//  if(myid==0)
//     printf("mdstep:%d,maxdt_diffusion:%.4e,maxdt_advection:%.4e,steps:%f\n",steps,max_dt_ttm,maxdt_advect,timestep/max_dt_ttm);
}

double QfromT(double T, double rho)
{
//TESTCASE
//return 2.5;

  double lgT = log10(T);
  //double lgr=log10(rho);
  double Q = bcinterp(&QfromT_interp, rho, lgT); // erwartet T in eV und rho in kg/m^3
  return Q; //Wird -1 bei Fehler

}

double EOS_ee_from_r_te(double r, double t)
{
  //if (r < RHOMIN) return 0;
  r = MAX(r, RHOMIN);
  r = MIN(r, 3500); //CHEAT

  double tsqrt = sqrt(t);
  point pout;
  pout.x = r;
  pout.y = tsqrt;
  // nnhpi_interpolate(intp_e_from_r_tsqrt.interpolator, &pout); //naturla neigh, sibson-rule
  lpi_interpolate_point(intp_ee_from_r_tesqrt.interpolator, &pout); //linear
  //pout.z*=26.9815*AMU*6.2415091E18; // J/kg --> eV/Atom
  if(isnan(pout.z)!=0)
  { 
    printf("ee_from_r_te retunred NaN!.r:%.4e,t:%.4e",r,t);
    error("ee_from_r_te retunred NaN");
  }
  return pout.z;
}

double EOS_cve_from_r_te(double r, double t)
{

  r = MAX(r, RHOMIN);
  r = MIN(r, 3500); //CHEAT

  //double tsqrt = sqrt(t);
  
  point pout;
  pout.x = r;
  pout.y = t;

  // nnhpi_interpolate(intp_e_from_r_tsqrt.interpolator, &pout); //naturla neigh, sibson-rule
  lpi_interpolate_point(intp_cve_from_r_te.interpolator, &pout); //linear

  if(isnan(pout.z)!=0)
  { 
    printf("cve_from_r_te retunred NaN!.r:%.4e,t:%.4e",r,t);
    error("cve_from_r_te retunred NaN");
  }

  // double tupper=t+0.02*t;
  // double tlower=t-0.02*t;
  // double eupper=EOS_ee_from_r_te(r,tupper);
  // double elower=EOS_ee_from_r_te(r,tlower);

  // pout.z=(eupper-elower)/(tupper-tlower);

  pout.z *= r; // J/(K*kg) --> J/K/m^3
  pout.z *= 11604.5; // -->J/eV/m^3
  pout.z *= 1e-30; // --> J/eV/Angs^3
  pout.z *= J2eV; // --> eV/eV/A^3

  return pout.z;
}


double EOS_te_from_r_ee(double r, double e)
{
  //if (r < RHOMIN) return 0;
  r = MAX(r, RHOMIN);
  r=MIN(r,3500);
  //double eSI=e/(26.9815*AMU*6.2415091E18); // eV/Atom ---> J/kg
  double a = pow(intp_ee_from_r_tesqrt.ymin, 2.0);
  double b = pow(intp_ee_from_r_tesqrt.ymax, 2.0);
  double m = fminbnd(a, b, eeminfun, 1e-4, r, e);
  return m;
}
double eeminfun(double x, double r, double e)
{
  return ABS(EOS_ee_from_r_te(r, x) - e);
}





// Wenn FDTD und TTM gemeinsam aktiv, muss in der Diffusion-loop dafür gesorgt werden, dass
// die EM-Feldkomponenten des verangenen Zeitschritts nicht mit den aktuellen vertauscht werden
// -->kopieren und wiederherstellen
#ifdef FDTD
void SwapTTM(int i, int j, int k)
{
  double tmp;
  tmp = l2[i][j][k].Ezx;
  l2[i][j][k].Ezx = l1[i][j][k].Ezx;
  l1[i][j][k].Ezx = tmp;

  tmp = l2[i][j][k].Ezy;
  l2[i][j][k].Ezy = l1[i][j][k].Ezy;
  l1[i][j][k].Ezy = tmp;

  tmp = l2[i][j][k].Ex;
  l2[i][j][k].Ex = l1[i][j][k].Ex;
  l1[i][j][k].Ex = tmp;

  tmp = l2[i][j][k].Ey;
  l2[i][j][k].Ey = l1[i][j][k].Ey;
  l1[i][j][k].Ey = tmp;

  tmp = l2[i][j][k].Hx;
  l2[i][j][k].Hx = l1[i][j][k].Hx;
  l1[i][j][k].Hx = tmp;

  tmp = l2[i][j][k].Hy;
  l2[i][j][k].Hy = l1[i][j][k].Hy;
  l1[i][j][k].Hy = tmp;

  tmp = l2[i][j][k].Hzx;
  l2[i][j][k].Hzx = l1[i][j][k].Hzx;
  l1[i][j][k].Hzx = tmp;

  tmp = l2[i][j][k].Hzy;
  l2[i][j][k].Hzy = l1[i][j][k].Hzy;
  l1[i][j][k].Hzy = tmp;

  tmp = l2[i][j][k].Jzx;
  l2[i][j][k].Jzx = l1[i][j][k].Jzx;
  l1[i][j][k].Jzx = tmp;

  tmp = l2[i][j][k].Jzy;
  l2[i][j][k].Jzy = l1[i][j][k].Jzy;
  l1[i][j][k].Jzy = tmp;

  tmp = l2[i][j][k].Jx;
  l2[i][j][k].Jx = l1[i][j][k].Jx;
  l1[i][j][k].Jx = tmp;

  tmp = l2[i][j][k].Jy;
  l2[i][j][k].Jy = l1[i][j][k].Jy;
  l1[i][j][k].Jy = tmp;


  l2[i][j][k].sigmax = l1[i][j][k].sigmax; //noetig?
  l2[i][j][k].sigmay = l1[i][j][k].sigmay;
  l2[i][j][k].mu = l1[i][j][k].mu;
  l2[i][j][k].eps = l1[i][j][k].eps;


  //Lorentz
  tmp = l2[i][j][k].Pzx;
  l2[i][j][k].Pzx = l1[i][j][k].Pzx;
  l1[i][j][k].Pzx = tmp;

  tmp = l2[i][j][k].Pzy;
  l2[i][j][k].Pzy = l1[i][j][k].Pzy;
  l1[i][j][k].Pzy = tmp;

  tmp = l2[i][j][k].Px;
  l2[i][j][k].Px = l1[i][j][k].Px;
  l1[i][j][k].Px = tmp;

  tmp = l2[i][j][k].Py;
  l2[i][j][k].Py = l1[i][j][k].Py;
  l1[i][j][k].Py = tmp;

  tmp = l2[i][j][k].Pzx;
  l2[i][j][k].Jlzx = l1[i][j][k].Jlzx;
  l1[i][j][k].Jlzx = tmp;

  tmp = l2[i][j][k].Jlzy;
  l2[i][j][k].Jlzy = l1[i][j][k].Jlzy;
  l1[i][j][k].Jlzy = tmp;

  tmp = l2[i][j][k].Jlx;
  l2[i][j][k].Jlx = l1[i][j][k].Jlx;
  l1[i][j][k].Jlx = tmp;

  tmp = l2[i][j][k].Jly;
  l2[i][j][k].Jly = l1[i][j][k].Jly;
  l1[i][j][k].Jly = tmp;

//        int itmp=l2[i][j][k].natoms_old;
//        l2[i][j][k].natoms_old=l1[i][j][k].natoms_old;
//        l1[i][j][k].natoms_old=itmp;


}
#endif


#ifdef FDTD
int fitDL(int i, int j, int k)
{
//TESTCASE CONST
  /*
  epsinf=2.73;
  omega_plasma=2.2955e+16;
  gamma_p=1.1174e+15;
  omegapl_L=7.6595e+15;
  Omega0_L=2.4024e+15;
  Gamma_L=4.5199e+14;
  */

  /*
    l1[i][j][k].DL[0]=2.73;
    l1[i][j][k].DL[1]=1.1174e+15;
    l1[i][j][k].DL[2]=7.6595e+15;
    l1[i][j][k].DL[3]=2.4024e+15;
    l1[i][j][k].DL[4]=4.5199e+14;
    l1[i][j][k].DL[5]=2.2955e+16;
    return 0;
  */

  double Te = MAX(l1[i][j][1].temp, 0.0259);
  double Ti = MAX(l1[i][j][1].md_temp, 0.0259);

  double rho = l1[i][j][1].dens;
  //lg10Timin=-1.587505178144741
  double lgTe = log10(Te);
  double lgTi = log10(Ti);
  //HOTFIX
  lgTi = MAX(lgTi, -1.587504);
  lgTe = MAX(lgTe, -1.587504);

  lgTe = MIN(lgTe, Lop1i.ymax);
  lgTi = MIN(lgTi, Lop1i.zmax);

  l1[i][j][1].DL[0] = do_tricubinterp(&Lop1i, rho, lgTe, lgTi); //epsinf, dimlos
  l1[i][j][1].DL[1] = do_tricubinterp(&Lop2i, rho, lgTe, lgTi) / hbarev; //gamma_plasma
  l1[i][j][1].DL[2] = do_tricubinterp(&Lop3i, rho, lgTe, lgTi) / hbarev; //omegapl_L
  l1[i][j][1].DL[3] = do_tricubinterp(&Lop4i, rho, lgTe, lgTi) / hbarev; //Omega0_L
  l1[i][j][1].DL[4] = do_tricubinterp(&Lop5i, rho, lgTe, lgTi) / hbarev; //Gamma_L
  l1[i][j][1].DL[5] = sqrt(l1[i][j][1].ne * ECHARGE * ECHARGE / ECONST / EMASS); //Omega_plasma

  int foo;
  for (foo = 0; foo < 6; foo++)
    if (l1[i][j][1].DL[foo] == -1)
      return -1;

  return 0;
}
#endif

#ifdef MPI
void ttm_create_mpi_datatypes(void)
{
  { /* type for our basic struct */

    // ************************************************************************************************
    // * MY MOD: Only send what is needed during DIFFUSION calculation, i.e. electron-temp. & kappa only
    // * ACHTUNG: natoms wird auch kommuniziert, da jede Zelle wissen muss, ob die Nachbarzelle aus aktiv ist
    // *    Eigentlich würde es reichen diese Info nur jeden MD-steps zu kommunizieren.
    // ****************************************************************************************************

    //ACHTUNG: Hier wird mpi_element erzeigt.Wird benötigt für das ursprüngliche (Ullrich) Kommunikations-Schema via SendRecv
    ttm_Element tmpelement;
    ttm_Element * tmpelement_pointer;
    tmpelement_pointer = &tmpelement;
    MPI_Aint tmpaddr;
    int blockcounts[4] = {1, 1, 1, 1};
    MPI_Datatype types[4] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_UB};
    MPI_Aint displs[4];

    MPI_Address(tmpelement_pointer, &tmpaddr);
    MPI_Address(&tmpelement.natoms, &displs[0]);
    MPI_Address(&tmpelement.temp, &displs[1]);
    //MPI_Address(&tmpelement.u, &displs[1]);
    MPI_Address(&tmpelement.fd_k, &displs[2]);
    tmpelement_pointer++;
    MPI_Address(tmpelement_pointer, &displs[3]);

    displs[3] -= tmpaddr;
    displs[2] -= tmpaddr;
    displs[1] -= tmpaddr;
    displs[0] -= tmpaddr;

    MPI_Type_struct(4, blockcounts, displs, types, &mpi_element);
    if (MPI_Type_commit(&mpi_element) != MPI_SUCCESS)
      error("type mpi_element failed to commit");

  }
  {
    // ****************************************************************
    // * Jetzt der derived-datatype für kommunikation für ttm-output
    // ****************************************************************
    int i;
    ttm_Element tmpelement;
    ttm_Element * tmpelement_pointer;
    tmpelement_pointer = &tmpelement;
    MPI_Aint tmpaddr;

#ifndef FDTD
//natoms temp md_temp xi u source dens vcom.x vcom.y vcom.z fd_k fd_g Z proc Ce
    int blockcounts[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[15] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                              MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE,
                              MPI_UB
                             };



    MPI_Aint displs[15];
#else

// natoms temp md_temp U xi source dens vcom.x vcom.y vcom.z fd_k fd_g Z proc Ce Ezx Ezy Hx Hy sigmax sigmay Hzx Hzy Ex Ey
    int blockcounts[25] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, //natoms-fd_g
                           1, 1, 1, //Z-Ce
                           1, 1, 1, 1, 1, 1, // Ezx-sigmay
                           1, 1, 1, 1
                          }; //Hzx-Ey
    MPI_Datatype types[25] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                              MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE,
                              //ADDITIONAL FDTD-elements from here
                              MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                              MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                              MPI_UB
                             };
    MPI_Aint displs[25];

#endif
    MPI_Address(&tmpelement, &tmpaddr);
    MPI_Address(&tmpelement.natoms, &displs[0]);
    MPI_Address(&tmpelement.temp, &displs[1]);
    MPI_Address(&tmpelement.md_temp, &displs[2]);

    MPI_Address(&tmpelement.xi, &displs[3]);

    MPI_Address(&tmpelement.source, &displs[4]);
    MPI_Address(&tmpelement.dens, &displs[5]);

    MPI_Address(&tmpelement.vcomx, &displs[6]);
    MPI_Address(&tmpelement.vcomy, &displs[7]);
    MPI_Address(&tmpelement.vcomz, &displs[8]);

    MPI_Address(&tmpelement.fd_k, &displs[9]);
    MPI_Address(&tmpelement.fd_g, &displs[10]);
    MPI_Address(&tmpelement.Z, &displs[11]);
    MPI_Address(&tmpelement.proc, &displs[12]);
    MPI_Address(&tmpelement.Ce, &displs[13]);

#ifndef FDTD
    tmpelement_pointer++;
    MPI_Address(tmpelement_pointer, &displs[14]);

    for (i = 0; i < 15; ++i)
    {
      displs[i] -= tmpaddr;
    }

    MPI_Type_struct(15, blockcounts, displs, types, &mpi_element2);
#else
    MPI_Address(&tmpelement.Ezx, &displs[14]);
    MPI_Address(&tmpelement.Ezy, &displs[15]);
    MPI_Address(&tmpelement.Hx, &displs[16]);
    MPI_Address(&tmpelement.Hy, &displs[17]);
    MPI_Address(&tmpelement.sigmax, &displs[18]);
    MPI_Address(&tmpelement.sigmay, &displs[19]);

    MPI_Address(&tmpelement.Hzx, &displs[20]);
    MPI_Address(&tmpelement.Hzy, &displs[21]);
    MPI_Address(&tmpelement.Ex, &displs[22]);
    MPI_Address(&tmpelement.Ey, &displs[23]);

    tmpelement_pointer++;
    MPI_Address(tmpelement_pointer, &displs[24]);
    for (i = 0; i < 25; ++i)
    {
      displs[i] -= tmpaddr;
    }

    MPI_Type_struct(25, blockcounts, displs, types, &mpi_element2);
#endif
    if (MPI_Type_commit(&mpi_element2) != MPI_SUCCESS)
      error("type mpi_element2 failed to commit");
  }

  /* datatype for one string of elements along z (short of 2 lattice points) */
  MPI_Type_contiguous(local_fd_dim.z - 2, mpi_element, &mpi_zrow);
  //MPI_Type_commit(&mpi_zrow);
  if (MPI_Type_commit(&mpi_zrow) != MPI_SUCCESS)
    error("type mpi_zrow failed to commit");


  { /* add displacements to skip ghost layers (mpi_zrow_block) */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore = &tmpel;

    int blockcounts[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_LB, mpi_zrow, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore, &displs[0]);
    tmpstore += 1;
    MPI_Address(tmpstore, &displs[1]);
    tmpstore += 1 + (local_fd_dim.z - 2);
    MPI_Address(tmpstore, &displs[2]);

    displs[2] -= displs[0];
    displs[1] -= displs[0];
    displs[0] -= displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_zrow_block);
    //MPI_Type_commit(&mpi_zrow_block);
    if (MPI_Type_commit(&mpi_zrow_block) != MPI_SUCCESS)
      error("type mpi_zrow_block failed to commit");

  }


  /* datatype for one layer of elements perpendicular to x
   * (short of 2 strings along z-axis)                    */
  MPI_Type_contiguous(local_fd_dim.y - 2, mpi_zrow_block, &mpi_xplane);
  //MPI_Type_commit(&mpi_xplane);
  if (MPI_Type_commit(&mpi_xplane) != MPI_SUCCESS)
    error("type mpi_xplane failed to commit");


  { /* add displacements to skip ghost layers (mpi_xplane_block) */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore = &tmpel;

    int blockcounts[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_LB, mpi_xplane, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore, &displs[0]);
    tmpstore += 2 + (local_fd_dim.z - 2);
    MPI_Address(tmpstore, &displs[1]);
    tmpstore += (2 + (local_fd_dim.z - 2)) * ((local_fd_dim.y - 2) + 1);
    MPI_Address(tmpstore, &displs[2]);

    displs[2] -= displs[0];
    displs[1] -= displs[0];
    displs[0] -= displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_xplane_block);
    //MPI_Type_commit(&mpi_xplane_block);
    if (MPI_Type_commit(&mpi_xplane_block) != MPI_SUCCESS)
      error("type mpi_xplane_block failed to commit");
  }

  // **********************************************************************
  // * datatype for one layer of elements perpendicular to y (short of 2) *
  // **********************************************************************
  MPI_Type_vector((local_fd_dim.x - 2), 1, 2 + (local_fd_dim.y - 2),
                  mpi_zrow_block, &mpi_yplane);
  //MPI_Type_commit(&mpi_yplane);
  if (MPI_Type_commit(&mpi_yplane) != MPI_SUCCESS)
    error("type mpi_yplane failed to commit");


  { /* add displacements to skip ghost layers (mpi_yplane_block) */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore = &tmpel;

    int blockcounts[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_LB, mpi_yplane, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore, &displs[0]);
    tmpstore += (2 + (local_fd_dim.z - 2)) *
                (2 + (local_fd_dim.y - 2));
    MPI_Address(tmpstore, &displs[1]);
    tmpstore += (1 + (local_fd_dim.x - 2)) *
                (2 + (local_fd_dim.z - 2)) *
                (2 + (local_fd_dim.y - 2));
    MPI_Address(tmpstore, &displs[2]);

    displs[2] -= displs[0];
    displs[1] -= displs[0];
    displs[0] -= displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_yplane_block);
    //MPI_Type_commit(&mpi_yplane_block);
    if (MPI_Type_commit(&mpi_yplane_block) != MPI_SUCCESS)
      error("type mpi_yplane_block failed to commit");

  }

  // ************************************************************
  // * datatype for one string of elements along y (short of 2) *
  // ************************************************************
  MPI_Type_vector((local_fd_dim.y - 2), 1, 2 + (local_fd_dim.z - 2),
                  mpi_element, &mpi_yrow);
  //MPI_Type_commit(&mpi_yrow);
  if (MPI_Type_commit(&mpi_yrow) != MPI_SUCCESS)
    error("type mpi_yrow failed to commit");

  /* add displacements to create datatype from which mpi_zplane will be built
   * (again, skipping ghost layers) */
  { /* mpi_yrow_block */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore = &tmpel;

    int blockcounts[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_LB, mpi_yrow, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore, &displs[0]);
    tmpstore += 2 + (local_fd_dim.z - 2);
    MPI_Address(tmpstore, &displs[1]);
    tmpstore += (1 + (local_fd_dim.y - 2)) *
                (2 + (local_fd_dim.z - 2));
    MPI_Address(tmpstore, &displs[2]);

    displs[2] -= displs[0];
    displs[1] -= displs[0];
    displs[0] -= displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_yrow_block);
    //MPI_Type_commit(&mpi_yrow_block);
    if (MPI_Type_commit(&mpi_yrow_block) != MPI_SUCCESS)
      error("type mpi_yrow_block failed to commit");

  }

#ifdef DEBUG
  if (myid == 0)
  { /* output size/extent comparisons */
    int size;
    MPI_Aint extent;

    MPI_Type_size(mpi_zrow, &size);
    MPI_Type_extent(mpi_zrow, &extent);
    printf("Size / Extent of mpi_zrow: %d / %ld\n", size, (long)extent);

    MPI_Type_size(mpi_yrow, &size);
    MPI_Type_extent(mpi_yrow, &extent);
    printf("Size / Extent of mpi_yrow: %d / %ld\n", size, (long)extent);

    MPI_Type_size(mpi_yrow_block, &size);
    MPI_Type_extent(mpi_yrow_block, &extent);
    printf("Size / Extent of mpi_yrow_block: %d / %ld\n", size, (long)extent);

    MPI_Type_size(mpi_element, &size);
    MPI_Type_extent(mpi_element, &extent);
    printf("Size / Extent of mpi_element: %d / %ld\n", size, (long)extent);

  }
#endif /*DEBUG*/
  // ********************************************************************************
  // * datatype for one layer of elements perpendicular to z (short of 2 strings)
  // ********************************************************************************
  MPI_Type_contiguous((local_fd_dim.x - 2), mpi_yrow_block, &mpi_zplane);
  MPI_Type_commit(&mpi_zplane);

  { /* add displacements to skip ghost layers (mpi_zplane_block) */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore = &tmpel;

    int blockcounts[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_LB, mpi_zplane, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore, &displs[0]);
    tmpstore += (2 + (local_fd_dim.z - 2)) *
                (2 + (local_fd_dim.y - 2));
    MPI_Address(tmpstore, &displs[1]);
    tmpstore += (1 + (local_fd_dim.x - 2)) *
                (2 + (local_fd_dim.z - 2)) *
                (2 + (local_fd_dim.y - 2));
    MPI_Address(tmpstore, &displs[2]);

    displs[2] -= displs[0];
    displs[1] -= displs[0];
    displs[0] -= displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_zplane_block);
    MPI_Type_commit(&mpi_zplane_block);
  }
}

void ttm_fill_ghost_layers(void)
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
  int i, j, k;

  /////////////////
  // x direction //
  /////////////////
  if (pbc_dirs.x == 1 || (my_coord.x != 0 && my_coord.x != cpu_dim.x - 1) )
  {
    /* send left slice to left neighbor. */
    /* Simultaneously receive slice from right neighbor */
    MPI_Sendrecv(&l1[1][0][0], 1, mpi_xplane_block, nbeast, 7100,
                 &l1[(local_fd_dim.x - 2) + 1][0][0], 1, mpi_xplane_block, nbwest, 7100,
                 cpugrid, &stati[0]);


    /* send right slice to right neighbor. */
    /* Simultaneously receive slice from left neighbor */
    MPI_Sendrecv(&l1[(local_fd_dim.x - 2)][0][0], 1, mpi_xplane_block, nbwest, 7200,
                 &l1[0][0][0], 1, mpi_xplane_block, nbeast, 7200,
                 cpugrid, &stati[1]);
  }
  else /* no pbc and we are at the surface */
    if (my_coord.x == 0 && my_coord.x != cpu_dim.x - 1) /* left surface */
    {
      /* only receive from right */
      MPI_Recv(&l1[(local_fd_dim.x - 2) + 1][0][0], 1, mpi_xplane_block, nbwest, 7100,
               cpugrid, &stati[0]);

      /* only send to right */
      MPI_Send(&l1[(local_fd_dim.x - 2)][0][0], 1, mpi_xplane_block, nbwest, 7200,
               cpugrid);
// WOZU??? -->Macht schon Sinn: die Buffer-Zelle muss mit natoms=0 initialisiert werden
//        sonst steht hier irgendein Schrott drin,evtl. > fd_min_atoms und do_DIFF
//        behandelt die Zelle, als sei sie aktiv!
      // left ghost layer receives reflecting bc
      for (j = 1; j <= (local_fd_dim.y - 2); j++)
      {
        for (k = 1; k <= (local_fd_dim.z - 2); k++)
        {
          // no atoms -> no conduction.
          l1[0][j][k].natoms = 0;
        }
      }
    }
    else if (my_coord.x == cpu_dim.x - 1 && my_coord.x != 0) /* right surface */
    {

      /* only send to left */
      MPI_Send(&l1[1][0][0], 1, mpi_xplane_block, nbeast, 7100,
               cpugrid);

      /* only receive from left */
      MPI_Recv(&l1[0][0][0], 1, mpi_xplane_block, nbeast, 7200,
               cpugrid, &stati[1]);
      // right ghost layer receives reflecting bc
      for (j = 1; j <= (local_fd_dim.y - 2); j++)
      {
        for (k = 1; k <= (local_fd_dim.z - 2); k++)
        {
          // no atoms -> no conduction.
          l1[local_fd_dim.x - 1][j][k].natoms = 0;
        }
      }
    }
    else
    { // two surfaces, just apply reflecting bc
      for (j = 1; j <= (local_fd_dim.y - 2); j++)
      {
        for (k = 1; k <= (local_fd_dim.z - 2); k++)
        {
          l1[0][j][k].natoms = l1[local_fd_dim.x - 1][j][k].natoms = 0;
        }
      }
    }
  /////////////////
  // y direction //
  /////////////////
  if (pbc_dirs.y == 1 || (my_coord.y != 0 && my_coord.y != cpu_dim.y - 1) )
  {
    /* send left slice to left neighbor. */
    /* Simultaneously receive slice from right neighbor */
    MPI_Sendrecv(&l1[0][1][0], 1, mpi_yplane_block, nbnorth, 710,
                 &l1[0][(local_fd_dim.y - 2) + 1][0], 1, mpi_yplane_block, nbsouth, 710,
                 cpugrid, &stati[2]);
    /* send right slice to right neighbor. */
    /* Simultaneously receive slice from left neighbor */
    MPI_Sendrecv(&l1[0][(local_fd_dim.y - 2)][0], 1, mpi_yplane_block, nbsouth, 720,
                 &l1[0][0][0], 1, mpi_yplane_block, nbnorth, 720,
                 cpugrid, &stati[3]);
  }
  else /* no pbc and we are at the surface */
  {
    if (my_coord.y == 0  && my_coord.y != cpu_dim.y - 1) /* left surface */
    {
      /* only receive from right */
      MPI_Recv(&l1[0][(local_fd_dim.y - 2) + 1][0], 1, mpi_yplane_block, nbsouth, 710,
               cpugrid, &stati[2]);
      /* only send to right */
      MPI_Send(&l1[0][(local_fd_dim.y - 2)][0], 1, mpi_yplane_block, nbsouth, 720,
               cpugrid);
      for (i = 1; i <= (local_fd_dim.x - 2); i++)
      {
        for (k = 1; k <= (local_fd_dim.z - 2); k++)
        {
          /* no atoms -> no conduction. */
          l1[i][0][k].natoms = 0;
        }
      }
    }
    else if (my_coord.y == cpu_dim.y - 1 && my_coord.y != 0) /* right surface */
    {
      /* only send to left */
      MPI_Send(&l1[0][1][0], 1, mpi_yplane_block, nbnorth, 710,
               cpugrid);
      /* only receive from left */
      MPI_Recv(&l1[0][0][0], 1, mpi_yplane_block, nbnorth, 720,
               cpugrid, &stati[3]);

      /* right ghost layer receives reflecting bc */
      for (i = 1; i <= (local_fd_dim.x - 2); i++)
      {
        for (k = 1; k <= (local_fd_dim.z - 2); k++)
        {
          // no atoms -> no conduction.
          l1[i][local_fd_dim.y - 1][k].natoms = 0;
        }
      }
    }
    //else { error("This should be logically impossible.\n");} //<--doch klar!
  }
  /////////////////
  // z direction //
  /////////////////
#ifndef FDTD
  if (pbc_dirs.z == 1 || (my_coord.z != 0 && my_coord.z != cpu_dim.z - 1) ) //BULK
  {
    MPI_Sendrecv(&l1[0][0][1], 1, mpi_zplane_block, nbup, 71,
                 &l1[0][0][(local_fd_dim.z - 2) + 1], 1, mpi_zplane_block, nbdown, 71,
                 cpugrid, &stati[4]);
    MPI_Sendrecv(&l1[0][0][(local_fd_dim.z - 2)], 1, mpi_zplane_block, nbdown, 72,
                 &l1[0][0][0], 1, mpi_zplane_block, nbup, 72,
                 cpugrid, &stati[5]);
  }
  else //no pbc & we are at surface
  {
    if (my_coord.z == 0) // && my_coord.z!=cpu_dim.z-1) //lower surf
    {
      MPI_Recv(&l1[0][0][(local_fd_dim.z - 2) + 1], 1, mpi_zplane_block, nbdown, 71,
               cpugrid, &stati[4]);
      MPI_Send(&l1[0][0][(local_fd_dim.z - 2)], 1, mpi_zplane_block, nbdown, 72,
               cpugrid);
      for (i = 1; i <= (local_fd_dim.x - 2); i++)
      {
        for (j = 1; j <= (local_fd_dim.y - 2); j++)
        {
          /* no atoms -> no conduction. */
          l1[i][j][0].natoms = 0;
        }
      }
    }
    else if (my_coord.z == cpu_dim.z - 1) // && my_coord.z!=0) //upper surface //down=+z,up=-z
    {
      MPI_Send(&l1[0][0][1], 1, mpi_zplane_block, nbup, 71,
               cpugrid);
      MPI_Recv(&l1[0][0][0], 1, mpi_zplane_block, nbup, 72,
               cpugrid, &stati[5]);
      /* right ghost layer receives reflecting bc */
      for (i = 1; i <= (local_fd_dim.x - 2); i++)
      {
        for (j = 1; j <= (local_fd_dim.y - 2); j++)
        {
          /* no atoms -> no conduction. */
          l1[i][j][local_fd_dim.z - 1].natoms = 0;
        }
      }
    }
    else { error("This should be logically impossible.\n");}
  }

#endif
}
#endif /*MPI*/