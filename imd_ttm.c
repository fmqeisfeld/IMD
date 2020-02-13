#include "imd.h"
#ifdef DEBUG
#include <assert.h>
#endif
#ifdef BUFCELLS
#define NBUFFC 2
#else
#define NBUFFC 0
#endif /*BUFCELLS*/

// ****************************************************
// * ACHTUNG: Das ist eine modifizierte TTM-Version,
// *          die sich mit loadbalance verträgt
// *          Allerdings bisher nur für 1D-TTM
// ***************************************************

#include <sys/time.h> // Performance messen
#include <complex.h>
// #include "/usr/include/gsl/gsl_sf_fermi_dirac.h"  // Fuer fermi-dirac integral zur berechnung der FEG-internen energie
                                                   // als alternative zur interpol.Tabelle (falls gewünscht: EOSMODE = 0)
//gsl_sf_fermi_dirac deakt. (genau wie FEG_ee_from_r_ne_te), weil es auf hazel nicht zu finden ist


//TODO: 
// TTM-READ UPDATEN FUER FDTD-NUTZUNG

//#define DENSHOTFIX //erste Zelle soll für bestimmte Zeit festkörperdichte haben!
                   //Da meine umgebungsdichte-berechnung für die oberfläche eine geringere dichte liefert 
                   //als im Bulk. Im Kontinuum-Bild hat oberfläche aber auch bulk-dichte
                   //Fazit: --> kaum unterschied, nur geringfügig größere Refl.

#define DENSHOTFIXSTEP 10000 //bist zu diesem step bleibt surf-dens const.

#define EOSMODE 1  // 1=EOS-TABLE, 0 = Free Electron Gas
                   // ACHTUNG: FEG Momentan totaler BS --> t_from_e viel zu ungenau (+/- 10 %)
                   // Grund: weiss nicht genau aber vermute Fermi_dirac integral 
                   // Evtl. wäre "manuelles" integrieren zuverlässiger
#define ADVMODE 1  // 0=NO ADVECTION, 1=DISCRETE FLUX SOLVER (PREDICT ATOMIC FLUXES)
//#define ADVMODE2d  // FALLS y dim offen sein soll, müssen atomic-fluxes auch über kanten kommmuniziert werden

#define VLATTICE   //VIRTUAL LATTICE HINTER DER PROBE. ACHTUNG: NUR 1D !!!!
                    //Falls gewünscht kann mit vlatbuffer die zahl
                    //der zellen (vom ende der probe gezählt) angegeben werden, die NICHT im TTM berücksichtigt
                    //werden soll, da diese als Puffer dienen 
                    //Grund:md-temp (und damit auch elek.temp) nimmt wegen Schockabsorbierenden Randbed. 
                    //stark ab -> PROBLEM: Wide-range properties...
                    //zellen mit "-1" atomen sind echte, von VLATTICE deaktivierte Zellen
                    //Zellen mit "-2" atomen sind virtuelle Zellen
#define DEBUG_LEVEL 1




#define TTMOUTBUFLEN 400  //Wird bei ttm-output im MPIIO-Modus gebraucht. Entspricht max. Zeilenlänge
#define node  l1[i]       //Achtung: MPIIO-lange nicht mehr aktualisiert (z.b. fehlt colrad)
#define node2 l2[i]

#define RHOMIN       50


// ****************************************************
// *                MAIN FUNC
// *****************************************************
void calc_ttm()
{
  int i, j, k;
  double fd_subtime;
  double tau_DIFF; //integr. timestep-size in IMD-time units
  // double tau_FDTD; //integr. timestep-size in SI-units

  // diff_substeps=fd_n_timesteps;
  // tau_DIFF=(timestep)/((double) fd_n_timesteps);


  // int fdtd_substeps;
  //int diff_substeps;
  // int fdtd_div_diff; // fdtd_substeps/diff_substeps
  // ACHTUNG: integr.zeitschritt fuer Diffusion ist maximal so groß
  // wie integr.zeitschritt fuer maxwell-solver
  // Außerdem kann der benutzer via fd_n_timesteps eine untere grenze bestimmen

  // int telaps_md;     //fuer timing
  // int telaps_ttm;
  update_fd();
#if ADVMODE == 1  
  do_ADV(1.0);
#endif
  do_cell_activation();
#ifdef COLRAD
  do_colrad(timestep*10.18*1.0e-15);
#endif
  do_FILLMESH(); // calc. wide-range properties
  ttm_fill_ghost_layers();  

  CFL_maxdt();    
  diff_substeps = MAX(fd_n_timesteps, (int) (timestep / max_dt_ttm));
  tau_DIFF = timestep / ((double) diff_substeps);  


  for (i = 0; i < diff_substeps; i++)
  {
    do_tmm(tau_DIFF); //Helmholtz-solver, tau_DIFF spielt hier keine Rolle mehr
    tmm_time += tau_DIFF * 10.18 * 1.0e-15; //in sek
    do_DIFF(tau_DIFF);
    do_FILLMESH();
    ttm_fill_ghost_layers();            
  }
  //xi_arr global muss noch aus xiarrlocal reduziert werden
  //beginnend bei index=1 (brauche ghost-layer nicht)
  MPI_Allgather(&xiarr_local[1], local_fd_dim.x-2, MPI_DOUBLE, 
                &xiarr_global[0],local_fd_dim.x-2, MPI_DOUBLE, cpugrid);

{
    //Calc internal eng and updt new U after diff
    // tot_elec_energy_local =0;
    // for (i=1; i<local_fd_dim.x-1; ++i)
    // {
    // for (j=1; j<local_fd_dim.y-1; ++j)
    // {
    // for (k=1; k<local_fd_dim.z-1; ++k)
    // {
    //       if(node.natoms >= 1)
    //       {
    //         //node.U =EOS_ee_from_r_te(node.dens, node.temp * 11604.5) * 26.9815 * AMU * J2eV; //eV/Atom      
    //         tot_elec_energy_local += node.U*((double) node.natoms);
    //       }      
    //       else
    //       {
    //         node.U=0.0; //node2.U=0.0;
    //       }
    // }
    // }
    // }
}

MPI_Reduce(&Eabs_local, &Eabs_global, 1, MPI_DOUBLE, MPI_SUM, 0, cpugrid);

if(myid==0)
  printf("step:%d, It:%.2e, substeps:%d, Finc:%.4e, t-t0:%.4e, Refl:%.4e, laser_active:%d \n",
          steps,I_t, diff_substeps, Eabs_global * eV2J / laser_spot_area,(tmm_time - laser_t_0) * 1e15,tmm_refl, laser_active);

}

// *********************************************************************************************************************************************
// * UPDATE_FD COMPUTES MD-TEMPERATURES AND ATOMIC FLUXES ACROSS FD-CELLS  AND CLEARS OLD ARRAYS
// *********************************************************************************************************************************************
void update_fd()
{
// printf("myid:%d,entered update_fd\n",myid);  
  int i, j, k, i_global, j_global, k_global;
  //int fd_tags=(local_fd_dim.x-2)*(local_fd_dim.y-2)*(local_fd_dim.z-2)+1; //every MPI-message needs to be unique
  int loopvar, l;


  l = 0;
  cellptr p;

  //buffer arrays für allreduce
  int * natomslocal, *natomsglobal;
  int * totneighslocal, *totneighsglobal;

  //Folgende arrays zu globalen vars gemacht, da
  //imd integrate sie braucht

  // double *vcomxlocal,*vcomxglobal;
  // double *vcomylocal,*vcomyglobal;
  // double *vcomzlocal,*vcomzglobal;

  int *fluxfromrightlocal,*fluxfromleftlocal;
  int *fluxfromrightglobal,*fluxfromleftglobal;
  double * mdtemplocal,*mdtempglobal;

  alloc1darr(int, fluxfromrightlocal, global_fd_dim.x);
  alloc1darr(int, fluxfromleftlocal, global_fd_dim.x);

  alloc1darr(int, fluxfromrightglobal, global_fd_dim.x);
  alloc1darr(int, fluxfromleftglobal, global_fd_dim.x);

  alloc1darr(int, natomslocal, global_fd_dim.x);
  alloc1darr(int, natomsglobal, global_fd_dim.x);

  alloc1darr(int, totneighslocal, global_fd_dim.x);
  alloc1darr(int, totneighsglobal, global_fd_dim.x);

  //Nacht init_ttm verschoben
  // alloc1darr(double, vcomxlocal, global_fd_dim.x);
  // alloc1darr(double, vcomxglobal, global_fd_dim.x);

  // alloc1darr(double, vcomylocal, global_fd_dim.x);
  // alloc1darr(double, vcomyglobal, global_fd_dim.x);

  // alloc1darr(double, vcomzlocal, global_fd_dim.x);
  // alloc1darr(double, vcomzglobal, global_fd_dim.x);

  alloc1darr(double, mdtemplocal,  global_fd_dim.x);
  alloc1darr(double, mdtempglobal, global_fd_dim.x);

  //zero buffer arrays
  for(i=0;i<global_fd_dim.x;i++)
  {
    natomslocal[i]=totneighslocal[i]=fluxfromleftlocal[i]=fluxfromrightlocal[i]=0;
    vcomxlocal[i]=vcomylocal[i]=vcomzlocal[i]=mdtemplocal[i]=0.0;

    xiarr_global[i]=0.0;
  }


#if DEBUG_LEVEL>1
  printf("steps:%d,proc:%d,entered update_fd\n", steps, myid);
#endif


#ifdef VLATTICE
  last_active_cell_local=-5000;
#endif

#ifdef DENSHOTFIX
  first_active_cell_local=first_active_cell_global=999999;
#endif 

#if ADVMODE==1
// *********************************************************
// Clear edges & corners of fluxes & temp (erstmal nur 1D und 2D)  *
// *********************************************************
  for (k = 0; k < 2; k++)
  {
    l1[0].flux[k] = l1[local_fd_dim.x - 1].flux[k] = 0;
    l1[local_fd_dim.x - 1].flux[k] = l1[0].flux[k] = 0;
  }
  l1[0].temp = l1[local_fd_dim.x - 1].temp = 0;
  l1[local_fd_dim.x - 1].temp = l1[0].temp = 0;
#endif

  // **********************************************************************
  // * Clear arrays
  // ************************************************************************
  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {
    // i_global =  ((i - 1) + my_coord.x * (local_fd_dim.x - 2));    
    //CLEAR ARRAYS
    node.vcomx = 0.0;
    node.vcomy = 0.0;
    node.vcomz = 0.0;
    node.source = node2.source = 0.0;
#if ADVMODE==1
    node.flux[0] = node.flux[1] == 0;    
#endif
    node.natoms_old = node2.natoms_old = node.natoms;
    node.natoms = node2.natoms = 0; 
    node.xi = node2.xi = 0.0;
    node.source = node2.source = 0.0;

    xiarr_local[i]=0.0;
    node.dens=node2.dens=0;
  }
// **********************
//   1st LOOP OVER ATOMS
// **********************
// printf("myid:%d,1st loop start\n",myid);  
 for (k=0; k<NCELLS; ++k) 
 {
    p = CELLPTR(k);
    if(p->lb_cell_type != LB_REAL_CELL)
      continue;
    for (l=0; l<p->n; ++l) 
    {
      int i_global=(int) (ORT(p,l,X)/fd_h.x) ;
      // int i_local=i_global+1-my_coord.x*(local_fd_dim.x-2);

      natomslocal[i_global]++;
      totneighslocal[i_global]+= NUMNEIGHS(p, l);

      vcomxlocal[i_global]+=IMPULS(p,l,X);
      vcomylocal[i_global]+=IMPULS(p,l,Y);
      vcomzlocal[i_global]+=IMPULS(p,l,Z);

// ***********************************************
// ATOMIC FLUXES FOR ADVECTION, erstmal nur 1D/2D
// *************************************************
#if ADVMODE==1
      //REMINDER
      // flux[0] : teilchen erhalten von +x,y
      // flux[1] : teilchen erhalten von -x,y
      if (steps > 0)
      {
        if (p->fdi[l] != i_global)
        {
          // +x
          if (p->fdi[l] > i_global)
            fluxfromrightlocal[i_global]++;            
          //-x
          else if (p->fdi[l] < i_global)
            fluxfromleftlocal[i_global]++;
        }
      }      
      p->fdi[l] = i_global; //für advection
      p->fdj[l] = 1;
      p->fdk[l] = 1;
#endif
    }
  }
 // *************************
 //  REDUCE Buffer arrays
 // *************************
 // MPI_Barrier(cpugrid);
 // printf("myid:%d,before allred\n",myid);  

 MPI_Allreduce(natomslocal, natomsglobal, global_fd_dim.x, MPI_INT, MPI_SUM, cpugrid);

 MPI_Allreduce(fluxfromleftlocal,  fluxfromleftglobal,  global_fd_dim.x, MPI_INT, MPI_SUM, cpugrid); 
 MPI_Allreduce(fluxfromrightlocal, fluxfromrightglobal, global_fd_dim.x, MPI_INT, MPI_SUM, cpugrid); 

 MPI_Allreduce(totneighslocal, totneighsglobal, global_fd_dim.x, MPI_INT, MPI_SUM, cpugrid);

 MPI_Allreduce(vcomxlocal, vcomxglobal, global_fd_dim.x, MPI_DOUBLE, MPI_SUM, cpugrid);
 MPI_Allreduce(vcomylocal, vcomyglobal, global_fd_dim.x, MPI_DOUBLE, MPI_SUM, cpugrid);
 MPI_Allreduce(vcomzlocal, vcomzglobal, global_fd_dim.x, MPI_DOUBLE, MPI_SUM, cpugrid);
// printf("myid:%d,after allred\n",myid);  
 // MPI_Barrier(cpugrid);

// ****************************** 
//   vcom muss korrigiert werden
// *******************************
for(i_global=0;i_global<global_fd_dim.x;i_global++)
{
  // int i_local=i_global+1-my_coord.x*(local_fd_dim.x-2);
  // if(i_local < 0 || i_local > local_fd_dim.x-1)
  //   continue;

    double totmass=(double) natomsglobal[i_global]*atomic_weight;
    if(totmass>0)
    {
      vcomxglobal[i_global] /= totmass;
      vcomyglobal[i_global] /= totmass;
      vcomzglobal[i_global] /= totmass;
// printf("myid:%d,vcomx:%f\n",myid,vcomxglobal[i_global]);

    }
    else
    {
      vcomxglobal[i_global]=0.0;
      vcomyglobal[i_global]=0.0;
      vcomzglobal[i_global]=0.0;
    }  

}
// printf("myid:%d,2nd loop start\n",myid);  
// **********************
// 2nd loop for md-temp
// **********************
for (k=0; k<NCELLS; ++k) 
{
  p = CELLPTR(k);
  if(p->lb_cell_type != LB_REAL_CELL)
    continue;
  for (l=0; l<p->n; ++l) 
  {            
      int i_global=(int) (ORT(p,l,X)/fd_h.x);
      // int i_local=i_global+1-my_coord.x*(local_fd_dim.x-2);
      // if(i_local < 0 || i_local > local_fd_dim.x-1)
      //   continue;

      mdtemplocal[i_global] += MASSE(p, l) * SQR(IMPULS(p, l, X) / MASSE(p, l)
                                   - vcomxglobal[i_global]);

      mdtemplocal[i_global] += MASSE(p, l) * SQR(IMPULS(p, l, Y) / MASSE(p, l)
                                   - vcomyglobal[i_global]);

      mdtemplocal[i_global] += MASSE(p, l) * SQR(IMPULS(p, l, Z) / MASSE(p, l)
                                   - vcomzglobal[i_global]);
  }
}
MPI_Allreduce(mdtemplocal, mdtempglobal, global_fd_dim.x, MPI_DOUBLE, MPI_SUM, cpugrid);

// **********************************************************
// Reconstruct  local TTM-lattice from global buffer arrays
// *************************************************************
// printf("myid:%d,reconstruct\n",myid);  
for(i_global=0; i_global < global_fd_dim.x;i_global++)
{
  int i_local=i_global+1-myid*(local_fd_dim.x-2);

// printf("myid:%d,ig:%d,il:%d,max:%d\n",myid,i_global,i_local,local_fd_dim.x-1);
  if(i_local<1 || i_local > local_fd_dim.x-2)
    continue;


    i=i_local;
    node.natoms=natomsglobal[i_global];
    node2.natoms=node.natoms;
    node.flux[0]=node2.flux[0]=fluxfromrightglobal[i];
    node.flux[1]=node2.flux[1]=fluxfromleftglobal[i];

    if(node.natoms > 0)
    {
      node.dens=((double) totneighsglobal[i_global])/((double) node.natoms) * atomic_weight / neighvol * 1660.53907; //kg/m^3       
      if(node.dens==0)
        node.dens=(double) node.natoms * atomic_weight / fd_vol * 1660.53907;    

      node.md_temp=mdtempglobal[i_global];
      node.md_temp /= (3.0 * (double) node.natoms);

#ifdef VLATTICE
      if(node.natoms >= fd_min_atoms && node.dens > RHOMIN)
      {  
        last_active_cell_local=MAX(last_active_cell_local,i_global- vlatbuffer);
      }
#endif     

#ifdef DENSHOTFIX
      if(node.natoms >= fd_min_atoms && node.dens > RHOMIN)
      {
        first_active_cell_local=MIN(first_active_cell_local,i_global);
      }
#endif           
    }
    else 
    {
      node.dens=0.0;  
      node.md_temp=0.0;
      node.natoms=node2.natoms=0;
    }

    node2.md_temp=node.md_temp;      
    node2.dens=node.dens;    
    ////////
    // ghost layers können direkt hier gefüllt werden
    ///////////
    if(i_local==1)
    {
      l1[0].natoms=natomsglobal[i_global-1];
      if(l1[0].natoms > 0 )
        l1[0].dens=((double) totneighsglobal[i_global-1])/((double) l1[0].natoms) * atomic_weight / neighvol * 1660.53907; //kg/m^3       

      l2[0].natoms=l1[0].natoms;
      l2[0].dens=l1[0].dens;
    }
    else if (i_local==local_fd_dim.x-2)
    {
      l1[local_fd_dim.x-1].natoms=natomsglobal[i_global+1];

      if(l1[local_fd_dim.x-1].natoms > 0 )
        l1[local_fd_dim.x-1].dens=((double) totneighsglobal[i_global+1])/((double) l1[local_fd_dim.x-1].natoms) *
                                    atomic_weight / neighvol * 1660.53907; //kg/m^3  

      l2[local_fd_dim.x-1].natoms=l1[local_fd_dim.x-1].natoms;
      l2[local_fd_dim.x-1].dens  =l1[local_fd_dim.x-1].dens;

// if(myid==1)      
//   printf("\n\n\nFUCK:%d, you:%f\n\n\n\n", l2[local_fd_dim.x-2].natoms, l2[local_fd_dim.x-2].dens);
    }




  }

  free1darr(natomslocal);
  free1darr(natomsglobal);
  free1darr(totneighslocal);
  free1darr(totneighsglobal);
  // free1darr(vcomxlocal);
  // free1darr(vcomxglobal);
  // free1darr(vcomylocal);
  // free1darr(vcomyglobal);
  // free1darr(vcomzlocal);
  // free1darr(vcomzglobal);
  free1darr(fluxfromleftlocal);
  free1darr(fluxfromleftglobal);
  free1darr(fluxfromrightlocal);
  free1darr(fluxfromrightglobal);
  free1darr(mdtemplocal);
  free1darr(mdtempglobal);

  // ******************************************************************************
  // *  ONLY ONCE: INITIALIZE ELECTRON TEMPERATURE AND CHECK EOS PLAUSIBILITY
  // *****************************************************************************
  if(steps<1)
  {
    for(i=1;i<local_fd_dim.x-1;i++)
    {
      node.natoms_old = node2.natoms_old = node.natoms;
      if (node.natoms >= fd_min_atoms && node.dens > RHOMIN)
      {
         node.temp = node2.temp = node.md_temp;
      }      
    } //for i
  }

#if DEBUG_LEVEL>1
  printf("steps:%d,proc:%d,update_fd complete\n", steps, myid);
#endif
  

#ifdef VLATTICE
  MPI_Status vlatstatus;  
  MPI_Allreduce(&last_active_cell_local, &last_active_cell_global, 1, MPI_INT, MPI_MAX, cpugrid);
  //cur_vlattice_proc=last_active_cell_global.rank;)  
  cur_vlattice_proc=(int) (((double) last_active_cell_global*num_cpus) / ((double) global_fd_dim.x));

//printf("oldvlat:%d,nuvlat:%d, lastcell:%d\n",old_vlattice_proc, cur_vlattice_proc, last_active_cell_global.val);
  
  if(cur_vlattice_proc != old_vlattice_proc) //last-filled cell is now on other proc --> send to new proc
  {
    //MPI_Bcast(&vlattice1[0], vlatdim, mpi_element2, last_active_cell_global.rank, cpugrid);       
    if(steps>0) //Beim 0-ten step macht das keinen sinn!
       MPI_Sendrecv(&vlattice1[0], vlatdim, mpi_element2, cur_vlattice_proc, 101010, 
                    &vlattice1[0], vlatdim, mpi_element2, old_vlattice_proc, 101010, cpugrid, &vlatstatus);       
  }

  old_vlattice_proc=cur_vlattice_proc;

  //Jetzt müssen noch die zellen deaktiviert werden, die als "puffer" dienen
  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {
    i_global =  ((i - 1) + myid * (local_fd_dim.x - 2));
    if(i_global > last_active_cell_global)
      l1[i].natoms=l2[i].natoms=-1;  //So erkenne ich deaktivierte Zellen
  }
#endif //VLATTICE

#ifdef DENSHOTFIX
  MPI_Allreduce(&first_active_cell_local, &first_active_cell_global, 1, MPI_INT, MPI_MIN, cpugrid);
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
  j=k=1;
  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {
    i_global =  ((i - 1) + myid * (local_fd_dim.x - 2));

    if (node.natoms >= fd_min_atoms && node.dens > RHOMIN)
    {

#ifdef DENSHOTFIX
      if(i_global==first_active_cell_global)
      {
        if(steps< DENSHOTFIXSTEP)
          node.dens=node2.dens=MAX(node2.dens,2680);
      }
#endif          
          //////////////////////////////////////////////////////////////          
          //           IONISATIONSGRAD UND ELEK.DICHTE
          ////////////////////////////////////////////////////////////// 
#ifndef COLRAD          
          node.Z=MeanCharge(node.temp*11604.5, node.dens, atomic_charge, atomic_weight,i,j,k);
          // node.Z=QfromT(node.temp,node.dens); //Möglichkeit zur Interpolation aus Tabelle, falls More-Fit ungenügend
#if DEBUG_LEVEL > 0
          if (node.Z == -1.0)
          {
            printf("steps:%d,proc:%d,i:%d,j:%d,k:%d, ERROR during QfromT in FILLMESH, Te:%f (K), dens:%f (kg/m^3),atoms:%d\n",
                   steps, myid, i, j, k, node.temp * 11604.5, node.dens, node.natoms);
            error("ERROR during QfromT in FILLMESH");
          }
#endif
          node.ne = node.Z * node.dens / (atomic_weight * AMU);
          node2.ne = node.ne;
#endif
          //////////////////////////////////////////////////////////////          
          //                      Wärmekapazität
          //////////////////////////////////////////////////////////////      
#if EOSMODE == 1          
            node.Ce = EOS_cve_from_r_te(node.dens, node.temp * 11604.5); //Interpol.Tabelle                    
#else          
            node.Ce=  FEG_cve_from_ne_te(node.dens, node.ne,node.temp*11604.5);
            // pritf("rho:%.4e,temp:%.4e, ne:%.4e,Cv:%.4e,\n",node.dens,node.temp,node.ne,)
#endif
          node2.Ce=node.Ce;

#if DEBUG_LEVEL>0
          if (node.Ce < 0)
          {
            printf("steps:%d,proc:%d,i:%d,j:%d,k:%d, ERROR during CfromT in FILLMESH, Ce:%.4e, Te:%f (K), dens:%f (kg/m^3)\n",
                   steps, myid, i, j, k, node.Ce, node.temp * 11604.5, node.dens);
            error("ERROR during CfromT in FILLMESH");
          }
#endif
          //////////////////////////////////////////////////////////////          
          //                WÄRMELEITFÄHIGKEIT 
          //////////////////////////////////////////////////////////////          
          node.fd_k = getKappa(node.temp, node.md_temp, node.ne, node.Z); //Hardcoding ist faster
          node2.fd_k = node.fd_k;

          //Möglichkeit zur trikubischen Interpolation aus Tablle --> Beliebiges Material 
          //Dasselbe kann für Elek.phonon-Kopplung über einfaches copy/paste ergänzt werden
          /* //TRIKUBISCHE INTERPOLATION AUS TABELLE // trikub. interpol ist recht langsam 
           node.fd_k=KappaInterpol(node.dens,node.temp,node.md_temp);
           if(node.fd_k==-1.0)
           {
              printf("steps:%d,proc:%d,i:%d,j:%d,k:%d, ERROR during KappaInterpol in FILLMESH, Te:%f (eV), Ti:%f (eV) ,dens:%f (kg/m^3),atoms:%d\n",
                      steps,myid,i,j,k,node.temp,node.md_temp,node.dens,node.natoms);
              error("ERROR during KappaInterpol in FILLMESH");
           }
          */

#if DEBUG_LEVEL>0  
          //ACHTUNG: Das ist nur bei Hardcoding noetig. Bei Interpolation, sollte die Interpol-Funktion error-check machen
          //NaN kappa
          if (isnan(node.fd_k) != 0 || node.fd_k < 0) //&& steps>0 || node.fd_k<0 && steps>0)
          {
            char errstr[255];
            sprintf(errstr,"proc:%d,i:%d,iloc:%d,k:%d,steps:%d, atoms:%d,fd_k is NaN,Te=%.4e,dens=%.4e,Ti=%.4e,Z=%.4e,ne=%.4e\n", 
              myid, i_global, i, k, steps, node.natoms, node.temp, node.dens,node.md_temp,node.Z,node.ne);
            error(errstr);
          }
#endif
          //////////////////////////////////////////////////////////////   
          //                   GAMMA  (KOPPL.CONST)
          //////////////////////////////////////////////////////////////                    
          node.fd_g = getGamma(node.temp, node.md_temp, node.ne, node.Z);
          node2.fd_g = node.fd_g;
          //NaN gamma
#if DEBUG_LEVEL>0
          if (isnan(node.fd_g) != 0 || node.fd_g < 0)
          {
            char errstr[255];
            sprintf(errstr,"proc:%d,i:%d,j:%d,k:%d,steps:%d, atoms:%d,fd_g is NaN,Te=%.4e,dens=%.4e,Ti=%.4e\n", 
              myid, i, j, k, steps, node.natoms, node.temp, node.dens,node.md_temp);
            error(errstr);
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
                   steps, myid, i, j, k, node.temp * 11604.5, node.dens);
            error("ERROR during fitDL in FILLMESH");
          }
#endif //DEBUG LEVEL
#endif  //FDTD     
//

        // ******************************************************************************
        // *  ONLY ONCE: INITIALIZE ELECTRON TEMPERATURE AND CHECK EOS PLAUSIBILITY
        // *****************************************************************************
        if (steps < 1)
        {
          if (node.natoms >= fd_min_atoms && node.dens > RHOMIN)
          {        
#if EOSMODE==1               
             node.U =EOS_ee_from_r_te(node.dens, node.temp * 11604.5) * 26.9815 * AMU * J2eV; // eV/Atom
             node2.U=node.U;

             //PLAUSIBILITY EOS CHECK:
             double echeck= EOS_ee_from_r_te(node.dens, node.temp * 11604.5)*26.9815 * AMU * J2eV;
             double tcheck = EOS_te_from_r_ee(node.dens, echeck/26.9815/AMU * eV2J) / 11604.5;
             double tinit=node.temp;

             if(ABS(tcheck -tinit) > tinit*0.01) // 1% unterschied
             {
               char errstr[255];

               sprintf(errstr,"ERROR: EOS Plausibility check failed, TfromU != Tinit. Tinit:%.4e, TfromU:%.4e\n"
                               "Maybe Interpolation table too sparse or increase tolerance",tinit,tcheck); 
               error(errstr);
             }
#else          
//ACHTUNG: gsl_fermi_integral schwankt extrem  
            node.U =FEG_ee_from_r_ne_te(node.dens, node.ne, node.temp * 11604.5) * 26.9815 * AMU * J2eV; // eV/Atom
            node2.U=node.U;

             //PLAUSIBILITY EOS CHECK:
             double echeck= FEG_ee_from_r_ne_te(node.dens, node.ne, node.temp * 11604.5)*26.9815 * AMU * J2eV;
             double tcheck = FEG_te_from_r_ne_ee(node.dens, node.ne, echeck/26.9815/AMU * eV2J) / 11604.5;
             double tinit=node.temp;

             if(ABS(tcheck -tinit) > tinit*0.05) // 5% unterschied
             {
               char errstr[255];

               sprintf(errstr,"ERROR: EOS Plausibility check failed, TfromU != Tinit. Tinit:%.4e, TfromU:%.4e,ne:%.4e\n",
                               tinit*11604.5,tcheck*11604.5,node.ne); 
               error(errstr);
             }            

#endif     
          } 
        } // if steps < 1


        }// if >= min_atoms ....
        else
        {
          node.fd_k = 0.0;
          node.fd_g = 0.0;
          node.Z = 0.0;
          node.ne = 0.0;
          node.Ce = 0.0;
          node.temp=0.0;  //temp wird nun in do_DIFF genutzt zum zu prüfen ob
        }                 //Nachbarzellen aktiv oder nicht,sodass atomzahl nicht mehr kommuniziert werden muss

  }     // for i

#ifdef VLATTICE
        //if(last_active_cell_global.rank==myid) // I own last active cell        
        if(cur_vlattice_proc==myid)
        {
          
          for(i=0;i<vlatdim;i++)
          {
            vlattice1[i].Z= MeanCharge(vlattice1[i].temp*11604.5, vlattice1[i].dens, atomic_charge, atomic_weight,i,1,1);
            vlattice1[i].ne= vlattice1[i].Z * vlattice1[i].dens / (atomic_weight * AMU);
            vlattice1[i].fd_k = getKappa(vlattice1[i].temp,vlattice1[i].md_temp, vlattice1[i].ne, vlattice1[i].Z); 
#if EOSMODE==1
            vlattice1[i].Ce = EOS_cve_from_r_te(vlattice1[i].dens, vlattice1[i].temp * 11604.5);
#else
            vlattice1[i].Ce = FEG_cve_from_ne_te(vlattice1[1].dens, vlattice1[i].ne, vlattice1[i].temp * 11604.5);
#endif
            vlattice1[i].fd_g= getGamma(vlattice1[i].temp, vlattice1[i].md_temp, vlattice1[i].ne, vlattice1[i].Z); 

            vlattice2[i].Z=vlattice1[i].Z;
            vlattice2[i].ne=vlattice1[i].ne;
            vlattice2[i].fd_k=vlattice1[i].fd_k;
            vlattice2[i].fd_g=vlattice1[i].fd_g;
          }
        }
#endif    //VLATTICE


#if DEBUG_LEVEL>1
  printf("steps:%d,proc:%d,FILLMESH complete\n", steps, myid);
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////
// Kommunikation von diversen größen, die nur jeden MD-step gebraucht werden
// statt jeden fd-step (dafür gibts ja fill_ghost_cells)
// Bisher nur Größen die für Advection-step gebraucht werden.
///////////////////////////////////////////////////////////////////////////////////////////
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


  /** MPI communication (can occur before and/or during MD calculations?) */
  /* Remember:
   * east -> -x
   * west -> +x
   * north-> -y
   * south-> +y
   * up   -> -z
   * down -> +z
   * *************/

  //first comm. yz-plane (+x/-x)
  if (myid != 0 && myid != num_cpus - 1) //BULK
  {
    //flux
    MPI_Sendrecv(l1[local_fd_dim.x - 2].flux, 2, MPI_INT, myid+1, 7302, //+x
                 l1[0].flux, 2, MPI_INT, myid-1, 7302,                  //-x
                 cpugrid, &stati[0]);

    MPI_Sendrecv(l1[1].flux, 2, MPI_INT, myid-1, 7402,                  //-x
                 l1[local_fd_dim.x - 1].flux, 2, MPI_INT, myid+1, 7402, //+x
                 cpugrid, &stati[1]);

    //selbes fuer interne eng.
    MPI_Sendrecv(&l1[local_fd_dim.x - 2].U, 1, MPI_DOUBLE, myid+1, 5302, //+x
                 &l1[0].U, 1, MPI_DOUBLE, myid-1, 5302,                  //-x
                 cpugrid, &stati[0]);

    MPI_Sendrecv(&l1[1].U, 1, MPI_DOUBLE, myid-1, 5402,                  //-x
                 &l1[local_fd_dim.x - 1].U, 1, MPI_DOUBLE, myid+1, 5402, //+x
                 cpugrid, &stati[1]);
  #ifdef COLRAD
    //selbes fuer colrad
    MPI_Sendrecv(l1[local_fd_dim.x - 2].y, colrad_neq, MPI_DOUBLE, myid+1, 4302, //+x
                 l1[0].y, colrad_neq, MPI_DOUBLE, myid-1, 4302,                  //-x
                 cpugrid, &stati[0]);

    MPI_Sendrecv(l1[1].y, colrad_neq, MPI_DOUBLE, myid-1, 4402,                  //-x
                 l1[local_fd_dim.x - 1].y, colrad_neq, MPI_DOUBLE, myid+1, 4402, //+x
                 cpugrid, &stati[1]);     
  #endif
  }
  else if (myid == 0) //SURFACE:only comm. with nbwest
  {
    //flux
    MPI_Sendrecv(l1[local_fd_dim.x - 2].flux, 2, MPI_INT, myid+1, 7302, //+x
                 l1[local_fd_dim.x - 1].flux, 2, MPI_INT, myid+1, 7402, //+x
                 cpugrid, &stati[1]);
    //selbes fuer u
    MPI_Sendrecv(&l1[local_fd_dim.x - 2].U, 1, MPI_DOUBLE, myid+1, 5302, //+x
                 &l1[local_fd_dim.x - 1].U, 1, MPI_DOUBLE, myid+1, 5402, //+x
                 cpugrid, &stati[1]);
  #ifdef COLRAD
    MPI_Sendrecv(l1[local_fd_dim.x - 2].y, colrad_neq, MPI_DOUBLE, myid+1, 4302, //+x
                 l1[local_fd_dim.x - 1].y, colrad_neq, MPI_DOUBLE, myid+1, 4402, //+x
                 cpugrid, &stati[1]);        
  #endif        
  }
  else if (myid == num_cpus - 1) //SURFACE:only comm with east
  {
    //flux
    MPI_Sendrecv(l1[1].flux, 2, MPI_INT, myid-1, 7402,         //-x
                 l1[0].flux, 2, MPI_INT, myid-1, 7302,         //-x
                 cpugrid, &stati[0]);
    //selbes fuer u
    MPI_Sendrecv(&l1[1].U, 1, MPI_DOUBLE, myid-1, 5402,         //-x
                 &l1[0].U, 1, MPI_DOUBLE, myid-1, 5302,         //-x
                 cpugrid, &stati[0]);
  #ifdef COLRAD
    MPI_Sendrecv(l1[1].y, colrad_neq, MPI_DOUBLE, myid-1, 4402,         //-x
                 l1[0].y, colrad_neq, MPI_DOUBLE, myid-1, 4302,         //-x
                 cpugrid, &stati[0]);        
  #endif        
  }

#if DEBUG_LEVEL>1
//  if(DEBUG_LEVEL>0)
  printf("proc:%d,steps:%d,do_COMMFLUX complete\n", myid, steps);
#endif

}

/****************************************************************************
 *        TTM INIT                                  *
*****************************************************************************/
/* init_ttm(): initialize FD stuff, ttm_init */
void init_ttm()
{
  int i, j, k;
  cell *p;

  if (myid == 0)
    printf("imdrestart=%d\n", imdrestart);


  global_fd_dim.x=ttmdimx; //User supplied
  if( global_fd_dim.x % num_cpus != 0)
  {
    error("global_fd_dim.x mod num_cpus != 0.");
  }

  local_fd_dim.x=  global_fd_dim.x/num_cpus;
  local_fd_dim.x +=2; //Ghost layer

  local_fd_dim.y=3;
  local_fd_dim.z=3;

  fd_h.x=box_x.x/global_fd_dim.x;
  fd_h.y=box_y.y;  
  fd_h.z=box_z.z;

  fd_vol = fd_h.x * fd_h.y * fd_h.z;
  max_dt_ttm = timestep / ((double) fd_n_timesteps); //nur zu beginn...im weiteren verlauf adaptiv
  neighvol = pow(sqrt(pair_pot.end[0]), 3.0) * 4.0 / 3.0 * M_PI;


  /* Time to initialize our FD lattice... */
  l1 = (ttm_Element*) malloc( local_fd_dim.x * sizeof(ttm_Element) );
  l2 = (ttm_Element*) malloc( local_fd_dim.x * sizeof(ttm_Element) );


  // ******************************************
  // * AUXILIARAY ARRAYS FOR imd_integrate.c
  // ********************************************
  xiarr_global=(double*) malloc( global_fd_dim.x * sizeof(double)); //damit xi auf allen procs vorhanden ist!
  xiarr_local=(double*) malloc(  local_fd_dim.x * sizeof(double)); //damit xi auf allen procs vorhanden ist!

  alloc1darr(double, vcomxlocal, global_fd_dim.x);
  alloc1darr(double, vcomxglobal, global_fd_dim.x);

  alloc1darr(double, vcomylocal, global_fd_dim.x);
  alloc1darr(double, vcomyglobal, global_fd_dim.x);

  alloc1darr(double, vcomzlocal, global_fd_dim.x);
  alloc1darr(double, vcomzglobal, global_fd_dim.x);


  #ifdef VLATTICE //NUR 1D!!! (global_fd_dim.y und global_fd_dim.z=1)
  vlattice1= (ttm_Element*) malloc(sizeof(ttm_Element)* vlatdim);
  vlattice2= (ttm_Element*) malloc(sizeof(ttm_Element)* vlatdim);

  for(i=0;i<vlatdim;i++)
  {
    vlattice1[i].dens=vlattice2[i].dens=vlatdens;
    vlattice1[i].temp=vlattice2[i].temp=0.0264;
    vlattice1[i].md_temp=vlattice2[i].md_temp=0.0264;    
  }
  #endif

 for (k=0; k<NCELLS; ++k) 
 {
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) 
    {
            
      int i_global=(int) (ORT(p,i,X)/fd_h.x) ;
      int i_local=i_global+1-myid*(local_fd_dim.x-2);

      if(i_local> -1 && i_local< local_fd_dim.x+1 )
      {
        p->fd_cell_idx.x=i_local; 
        p->fd_cell_idx.y=1;
        p->fd_cell_idx.z=1;
      }
    }
 }

  for(i=0;i<local_fd_dim.x;i++)
  {
    node2.source = node.source = 0.0;
    node2.fd_k = node.fd_k = 0;
    node.fd_g = node.fd_g = 0;
    node2.proc = node.proc = myid;
    node2.xi = node.xi = 0.0;
    node2.ne = node.ne = 0.0;
    node2.temp = node.temp = 0.0;
    node2.dens = node.dens = 0.0;
    node2.natoms = node.natoms = 0;
    node2.natoms_old = node.natoms_old = 0;
    node2.U = node.U = 0.0;        
    xiarr_local[i]=0.0;
    #ifdef FDTD
    int bar;
    for (bar = 0; bar < 6; bar++)
      node2.DL[bar] = node.DL[bar] = 0.0;
    #endif
#if ADVMODE==1 
    int foo;
    for (foo = 0; foo < 8; foo++)
    {
      node.flux[foo] = node2.flux[foo] = 0;
    }
#endif
  }


#ifdef LASERYZ
  error("TTM with LB does not support LASERYZ!")
#endif

  //MPI_Barrier(cpugrid);
  if (myid == 0)
  {
    printf("***************************************************\n");
    printf("*            TWO-TEMPERATURE MODEL   	      *\n");
    printf("***************************************************\n");
    printf("Global FD cell array dimension: %d \n", global_fd_dim.x);
    printf("Local FD cell array dimensions: %d\n",local_fd_dim.x);    
    printf("fd_h.x:%f A, fd_h.y:%f A,fd_h.z:%f A\n", fd_h.x, fd_h.y, fd_h.z);
    printf("Volume of one FD cell: %e [cubic Angstroms]\n",
           fd_h.x * fd_h.y * fd_h.z );
    printf("Volume of whole sample: %e [cubic Angstroms]\n",
           fd_h.x * fd_h.y * fd_h.z * global_fd_dim.x );

#ifdef VLATTICE
    printf("Virtual lattice dim:%d\n",vlatdim);
    printf("Virtual lattice buffer-cells:%d\n", vlatbuffer);
#endif    

#if DEBUG_LEVEL>0
    printf("DEBUG_LEVEL>0\n");
#endif

#ifdef DIRICHLET
    printf("dirichlet_surfx=%.2f, dirichlet_surfx_int=%d\n", dirichlet_surfx, dirichlet_surfx_int);
#endif
  }

#ifdef MPI
  // create MPI datatypes
  ttm_create_mpi_datatypes();
#endif


  // *****************************************
  // * READ AND BCAST INTERPOLATION TABLES
  // ******************************************  
#if EOSMODE==1
  nn_read_table(&intp_cve_from_r_te, "EOS_cve_from_r_te.txt");
  nn_read_table(&intp_ee_from_r_tesqrt, "EOS_ee_from_r_tesqrt.txt");
#endif
  //read_tricub_interp(&kappa_interp,"kappa.txt"); //Hardcoding ist schneller
  //Lese Drude-Lorentz Interpolationstabellen


#ifdef FDTD
  //Read and bcast Drude-Lorentz param tables
  read_tricub_interp(&Lop1i, "DL1.txt");
  read_tricub_interp(&Lop2i, "DL2.txt");
  read_tricub_interp(&Lop3i, "DL3.txt");
  read_tricub_interp(&Lop4i, "DL4.txt");
  read_tricub_interp(&Lop5i, "DL5.txt");
  
  //Set Minimum allowed Te und dens according to interpolation tables
  Temin=MAX(Temin,pow(10.0,Lop1i.ymin));
  Temin=MAX(Temin,pow(10.0,Lop2i.ymin));
  Temin=MAX(Temin,pow(10.0,Lop3i.ymin));
  Temin=MAX(Temin,pow(10.0,Lop4i.ymin));
  Temin=MAX(Temin,pow(10.0,Lop5i.ymin));

  rhomin=MAX(rhomin,Lop1i.xmin);
  rhomin=MAX(rhomin,Lop2i.xmin);
  rhomin=MAX(rhomin,Lop3i.xmin);
  rhomin=MAX(rhomin,Lop4i.xmin);
  rhomin=MAX(rhomin,Lop5i.xmin);

#endif


//ttm_writeout(9999);
  /***********************
  *  MY MOD: restart ttm *
  ************************/
  if (imdrestart > 0)
  {

    int readstep = imdrestart * checkpt_int;
    int readttm = readstep / ttm_int;
    ttm_read(readttm);

#ifdef FDTD
    t_SI = (double) imdrestart * (double) checkpt_int * timestep * 10.18 / 1.0e15;
    if (myid == 0)
      printf("t_SI:%.4e s\n", t_SI);
#endif

#ifdef TMM
    tmm_time=(double) imdrestart * (double) checkpt_int * timestep * 10.18 / 1.0e15;
#endif    
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
// *************************************************
void do_ADV(double tau)
{
  if (steps < 1) return;
  int i, j, k;
  int i_global, j_global, k_global;

  k = 1; //erstmal nur 2D und 1D

  do_COMMFLUX();

  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {

#ifdef FDTD
      SwapTTM(i, j, k);
#endif

      i_global = ((i - 1) + myid * (local_fd_dim.x - 2));
      //j_global = ((j - 1) + my_coord.y * (local_fd_dim.y - 2));
      //k_global =  ((k-1) + my_coord.z*(local_fd_dim.z-2));

//REMINDER
// flux[0] : teilchen erhalten von +x,y
// flux[1] : teilchen erhalten von -x,y

      double Nold = (double) node.natoms_old;
      double Nnew = (double) node.natoms; //tmp;


      if (Nnew > 0)
      {
        node2.U= node.U * Nold / Nnew+ tau*(
                                  // +x/-x
        + (double) node.flux[0] * l1[i + 1].U //erhalten von +x,y
        - (double) l1[i + 1].flux[1] * node.U //nach +x,y abgeflossen

        + (double) node.flux[1] * l1[i - 1].U //erhalten von -x,y
        - (double) l1[i - 1].flux[0] * node.U //nach -x,y abgeflossen
      ) / Nnew;

#ifdef COLRAD
          int l=0;
          for(l=0;l<colrad_neq;l++)
          {
            Ith(node2.y,l) = Ith(node.y,l) * Nold / Nnew+ tau*(
                                      // +x/-x
            + (double) node.flux[0] * Ith(l1[i + 1].y,l) //erhalten von +x,y
            - (double) l1[i + 1].flux[1] * Ith(node.y,l) //nach +x,y abgeflossen

            + (double) node.flux[1] * Ith(l1[i - 1].y,l) //erhalten von -x,y
            - (double) l1[i - 1].flux[0] * Ith(node.y,l) //nach -x,y abgeflossen
          // +y/-y
          ) / Nnew;
        }
#endif        

          //Temp updaten wenn zelle aktiviert
          if(node.natoms >= fd_min_atoms && node.dens > RHOMIN)
          {
#if EOSMODE==1            
            node2.temp = EOS_te_from_r_ee(node.dens, node2.U / (26.9815 * AMU * J2eV)) / 11604.5;
#else
            node.Z=MeanCharge(node.temp*11604.5, node.dens, atomic_charge, atomic_weight,i,j,k);
            node.ne = node.Z * node.dens / (atomic_weight * AMU); //Assumption: Quasi-neutrality condition
            node2.ne = node.ne;            
            node2.temp=  FEG_te_from_r_ne_ee(node.dens,node.ne, node2.U / (26.9815 * AMU * J2eV)) / 11604.5;
#endif     
          }     
      }
      else if (Nnew < 1)
      {
        node2.U = 0.0;
        node2.temp = 0.0;
#ifdef COLRAD
        int l=0;
        for(l=0;l<colrad_neq;l++)
        {
          Ith(node2.y,l)=0.0;
        }
#endif      
      }
  } //for i

  l3 = l1;
  l1 = l2;
  l2 = l3;

}

// **********************************************
// * Falls zahl der atome seit letztem step auf >= fd_min_atoms gesprungen ist,
// * soll zelle aktiviert werden.
// * Dabei wird auch geprüft ob der advection-step erfolgreich war und die neue Zelle
// * eine gewisse mindest-temperatur hat (zu interolierende transport-eigenschaften)
// * Falls nicht, wird als fallback die alte Mittelugns-methode verwendet.
// * Solte auch dies scheitern wird die md-temp verwendet
// **********************************************
void do_cell_activation(void)
{
  int i,j,k;
  int i_global,j_global,k_global;
  j_global=k_global=1;

    //NOW CHECK IF NEW CELL ACTIVATED!!
  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {
    i_global = ((i - 1) + myid * (local_fd_dim.x - 2));
    if (node.natoms_old >= fd_min_atoms && node.natoms < fd_min_atoms && node.dens > RHOMIN)
    {
      // ZELLE DEAKTIVIERT
#if DEBUG_LEVEL > 0
      printf("Warning:FD cell deactivated on proc %d on step %d at i:%d,j:%d,k%d with %d atoms and temp:%.4e\n", myid, steps,
             i_global, j_global, k_global, node.natoms, node.temp);
#endif
      // Cell deactivated. Deduce its electronic energy from E_new_local
      node.xi = 0.0;
    }
    // ZELLE AKTIVIERT
    else if (node.natoms_old < fd_min_atoms && node.natoms >= fd_min_atoms && node.dens > RHOMIN)
    {
#if DEBUG_LEVEL > 0
      printf("Warning:New FD cell activated on proc %d at ig:%d,jg:%d,kg:%d with %d atoms on step:%d and T=%.4e, dens=%.4e, atoms_old:%d\n",
             myid, i, j, k, node.natoms, steps, node.temp, node.dens, node.natoms_old);
#endif
      // *****************************************************
      // * NEU AKTIVIERTE ZELLE MIT UNSINNIGER TEMPERATUR,   *
      // * d.h. ADVECTION HAT NICHT FUNKTIONIERT             *
      // * ---> WENDE ALTES SCHEMA AN UND BERECHNE MITTEL    *
      // * AUS NACHBARZELLEN               *
      // *****************************************************

      if (isnan(node.temp) != 0 || node.temp <= Temin) //Temp zu klein (etwa 35K) -> altes schema
      {
#if DEBUG_LEVEL>0
        printf("proc:%d,steps:%d,ig:%d,jg:%d,kg:%d WARNING: Freshly activated cell with Te=%.4e is NaN or < Tmin:%.4e,atoms:%d, dens: %.4e ,"
               "using neighbor cells or mdtemp\n",
               myid, steps, i_global, j_global, 0, node.temp, Temin,  node.natoms, node.dens);
#endif

        // Freshly activated cell. Gets avg. electron energy of active
        // neighbor cells, the created energy is added to E_new_local
        int n_neighbors = 0;
        double E_el_neighbors = 0.0;
        // 6 indices: -x,x,-y,y,-z,z //

        // Was folgt, ist das ursprüngliche schema, mit den neu aktivierten Zellen
        // eine elec-temp. aus dem mittel der nachbarzelen zugewiesen wird ---> verletzt Energieerhaltung
        if (l1[i + 1].natoms >= fd_min_atoms && l1[i+1].dens >= RHOMIN)  //eigentlich noch dense zu prüfen aber bisher nicht kommuniziert
        {
          E_el_neighbors += SQR(l1[i + 1].temp);
          n_neighbors++;
        }
        if (l1[i - 1].natoms >= fd_min_atoms && l1[i-1].dens >= RHOMIN)
        {
          E_el_neighbors += SQR(l1[i - 1].temp);
          n_neighbors++;
        }
        {
          if (n_neighbors != 0)
          {
            node.temp = sqrt(E_el_neighbors / ((double)n_neighbors));
            node2.temp = node.temp;
            
#if DEBUG_LEVEL>0
            printf("proc:%d,steps:%d,i:%d,j:%d,k:%d, Te is NaN or <=0, using neighbor cells=>Te=%f\n",
                   myid, steps, i_global, j_global, 0, node.temp);
#endif
//HOTFIX: still < Tmin? --> use md-temp
            if (node.temp < Temin)
            {
#if DEBUG_LEVEL > 0
              printf("proc:%d,steps:%d,i:%d,j:%d,k:%d, Te=%.4e still < Tmin=%.4e, using MD-temp:%.4e\n",
                     myid, steps, i_global, j_global, k_global, node.temp, Temin, node.md_temp);
#endif
              node2.temp = node.temp = node.md_temp;
            }

          }
          else  // No neighbors? -> Get MD-temp
          {
            node.temp = node.md_temp;
            node2.temp = node.temp;
#if DEBUG_LEVEL > 0
            printf("proc:%d,steps:%d,i:%d,j:%d,k:%d, WARNING: No neighbors in activated cell. Using md-temp=>Te=%f\n",
                   myid, steps, i_global, j_global, k_global, node.temp);

#endif
          }
        } //isnan ....
        //Interne Energie muss noch geupdatet werden (Falls fallback auf altes Schema)

#if EOSMODE==1            
        node.U = node2.U= EOS_ee_from_r_te(node.dens, node.temp * 11604.5) * 26.9815 * AMU * J2eV; // eV/Atom
#else
        node.U=node2.U=FEG_ee_from_r_ne_te(node.dens,node.ne,node.temp*11604.5) * 26.9815 * AMU * J2eV; // eV/Atom
#endif
      } // endif isnan(temp) || temp<=0

    } // endif ..new cell activated...

    //node.natoms_old=node2.natoms_old=node.natoms; //nach update_fd verschoben


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

  //habe xi-berechnung ersetzt: Da alle wide-range properties von der Umgebungsdichte abhängen und 
  //nicht von der zahl der atome pro Zelle, muss auch xi umgeschrieben werden da sonst
  //zellen mit wenigen Atomen viel zu stark erhitzt werden!    
  //double xi_fac=fd_vol/3.0/((double) diff_substeps);///((double) node.natoms); //ORIGINAL
  double xi_fac = 26.9815 * AMU / 3.0 * 1e30 / ((double) diff_substeps); //NEU


  double Ce; //specific heat
  double invxsq = 1.0 / (fd_h.x * fd_h.x);
  // double invysq = 1.0 / (fd_h.y * fd_h.y);
  // double invzsq = 1.0 / (fd_h.z * fd_h.z);

  int i_global, j_global, k_global;
  double xmaxTe, xminTe, ymaxTe, yminTe, zmaxTe, zminTe; //temps
  double xmaxk, xmink, ymaxk, ymink, zmaxk, zmink; //kappas


#if DEBUG_LEVEL>1
//  if(DEBUG_LEVEL>0)
  printf("steps:%d,proc:%d,entered do_DIFF\n", steps, myid);
#endif
  j_global=k_global=1;

  ///////////////////////////
  // DIFFUSION             //
  ///////////////////////////
  for (i = 1; i < local_fd_dim.x - 1; i++)
  {     
    i_global =  ((i - 1) + myid* (local_fd_dim.x - 2));
    //compute absorbed laser-energy
    if (laser_active)
    {
      Eabs_local += node.source * fd_vol * tau; //eV
    }
    
#ifdef COLRAD  //Weil y-vec wird nur 1 mal allokiiert. (Colrad braucht keine zusätzl. Kopie)
    node2.y=node.y;
#endif
#ifdef FDTD
    SwapTTM(i, j, k);
#endif

    // if(i==1 || i==local_fd_dim.x-2)
    // printf("myid:%d, i:%d, ig:%d, atoms:%d,at-1:%d,at+1:%d\n",myid,i,i_global, node.natoms, l1[i-1].natoms, l1[i+1].natoms);

    /* only do calculation if cell is not deactivated */
    //ACHTUNG: dens in buffer-zelle ist nicht bekannt, da nicht kommuniziert
    //         aber in fill_ghost_layers wird natoms genullt falls dens < RHOMIN    

    if (node.natoms < fd_min_atoms || node.dens < RHOMIN)   continue;

    if (l1[i - 1].natoms < fd_min_atoms || l1[i-1].dens < RHOMIN)  xmin = i; else  xmin = i - 1; //das geht nun, da ich update_fd
    if (l1[i + 1].natoms < fd_min_atoms || l1[i+1].dens < RHOMIN)  xmax = i; else  xmax = i + 1; //auch direkt ghost-layer befülle


    // if (l1[i - 1].temp <= 0 )  xmin = i; else  xmin = i - 1; //das geht nun, da ich update_fd
    // if (l1[i + 1].temp <= 0 )  xmax = i; else  xmax = i + 1; //auch direkt ghost-layer befülle

    xmaxTe = l1[xmax].temp;
    xminTe = l1[xmin].temp;

    xmaxk = l1[xmax].fd_k;
    xmink = l1[xmin].fd_k;

#ifdef VLATTICE
    //if(last_active_cell_global.rank==myid) // I own last active cell        
    if(cur_vlattice_proc==myid)
    {
      xmaxTe = vlattice1[0].temp;
      xmaxk  = vlattice1[0].fd_k;
    }
#endif          

    /////////////////////////////////////////
    //dirichlet cases for outermost cells  //
    //ACHTUNG: Dirichlet Zellen werden     //
    //nicht in ttm-output geschrieben.     //
    //Sie sind "virtuell"      //
    /////////////////////////////////////////
     Ce = node.Ce; //eV/(eV*Angs^3)
    /***********************************************************************
          * Explicit diffusion with variable kappa   (Convervative formulation)  *
    ************************************************************************/    
    node2.temp=tau/Ce*
    ( 
      // dK/dx * d^2 T/dx^2
      (  ((node.fd_k+xmaxk)/2 * (xmaxTe-node.temp)*invxsq)
        -((node.fd_k+xmink)/2 * (node.temp-xminTe)*invxsq)
      )
      //now coupling+source term
      -node.fd_g*(node.temp-node.md_temp)
      +node.source
    ) +node.temp;
    

// if(myid==3)
// {
//   printf("tnow:%f,told:%f,k:%f,kmax:%f,tmax:%f,kmin:%f,tmin:%f\n",
//           node2.temp,node.temp, node.fd_k, xmaxk, xmaxTe, xmink,xminTe);
// }
    //Folgende Zeile setzt vorraus, dass Cve und U kompatiblen Tabllen zugrunde liegen
    //Sonst macht T_from_E keinen sinn
    node2.U=node.U + (node2.temp-node.temp)*Ce*fd_vol/((double) node.natoms); // eV/atom

    //node.xi += (node2.temp-node.md_temp)*xi_fac*node.fd_g/node.md_temp/((double) node.natoms);//Original
    node.xi += (node2.temp - node.md_temp) * xi_fac * node.fd_g / node.md_temp / node.dens; // NEU
    node2.xi = node.xi;

    xiarr_local[i] += (node2.temp - node.md_temp) * xi_fac * node.fd_g / node.md_temp / node.dens;
    //xiarr_global soll auf allen procs (für alle ttm-zellen) vorhanden sein (für imd_integrate.c)


    if(isnan(node2.temp) != 0 || isinf(node2.temp) != 0)
    {
      char errstr[255];
      sprintf(errstr,"Te got NaN or Inf in diffloop.natoms:%d, Ce:%.4e, told:%f, fdk:%.2e, fdg:%.2e,mdtemp:%f, xmaxk:%.2e,xmink:%.2e,xmaxT:%.2e,xminT:%.2e\n",
                      node.natoms, Ce, node.temp, node.fd_k, node.fd_g, node.md_temp, xmaxk,xmink, xmaxTe,xminTe);
      error(errstr);
    }
  } //for i



#ifdef VLATTICE
  //Keine diffusion für md-temp, nur kopplung mit elek.temp
  //Weil gitter-temp kaum diffundiert im vgl. mit elek-temp
  //if(last_active_cell_global.rank==myid) // I own last active cell        
  if(cur_vlattice_proc==myid)
  {

    int ilocal= last_active_cell_global+1-myid*(local_fd_dim.x-2);
    xminTe = l1[ilocal].temp;
    xmink  = l1[ilocal].fd_k;
    // Ci AUS GEOS: rho: 2.665655433e+03 temp: 3.000000000e+02 cvi: 8.589449886e+02           
    double Ci=8.589449886e+02;
    Ci *= vlatdens; // J/(K*kg) --> J/K/m^3
    Ci *= 11604.5; // -->J/eV/m^3
    Ci *= 1e-30; // --> J/eV/Angs^3
    Ci *= J2eV; // --> eV/eV/A^3
    for(i=0;i<vlatdim;i++)
    {
      // Te-diffusion
      Ce = vlattice1[i].Ce;                        
      vlattice2[i].temp=tau/Ce*
      //first diffusion terms
      ( 
        // dK/dx * d^2 T/dx^2
        ((vlattice1[i].fd_k+xmaxk)/2 * (xmaxTe- vlattice1[i].temp)*invxsq)
       -((vlattice1[i].fd_k+xmink)/2 * (vlattice1[i].temp-xminTe)*invxsq)
      
        //now coupling+source term
        -vlattice1[i].fd_g*(vlattice1[i].temp- vlattice1[i].md_temp)
      ) +vlattice1[i].temp;            

      //lattice
      double dT=tau/Ci*vlattice1[i].fd_g*(vlattice1[i].temp- vlattice1[i].md_temp);
      vlattice2[i].md_temp=vlattice1[i].md_temp + dT;          
    }
  }
#endif   

  /* take care - l1 must always be the updated lattice */

  l3 = l1;
  l1 = l2;
  l2 = l3;


#ifdef VLATTICE
  vlattice3=vlattice1;
  vlattice1=vlattice2;
  vlattice2=vlattice3;
#endif

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
  n = global_fd_dim.x; // * global_fd_dim.y * global_fd_dim.z;
  nlocal = (local_fd_dim.x - 2);

//ACHTUNG: MPIIO-OUTPUT SCHON LANGE NICHT MEHR AKTUALISIESRT. COLRAD FEHLT 
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
                i_global, j_global, k_global, node.natoms, node.temp,
                node.md_temp, node.xi,
                node.source, node.dens,
                node.vcomx, node.vcomy, node.vcomz,
                node.fd_k, node.fd_g,
 #ifndef FDTD
                node.Z, node.proc, node.Ce
 #else
                node.Z, node.proc, node.Ce,
                node.Ezx, node.Ezy, node.Hx, node.Hy,
                node.sigmax, node.sigmay,
                node.Hzx, node.Hzy, node.Ex, node.Ey
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


#ifdef VLATTICE
    MPI_Status vlatstatus;
    if(steps>0)
    {
      if(myid==cur_vlattice_proc)
      {
        MPI_Send(&vlattice1[0], vlatdim, mpi_element2, 0, 101011, cpugrid);
      }

      if(myid==0)
      {
        MPI_Recv(&vlattice1[0], vlatdim, mpi_element2, cur_vlattice_proc, 101011, cpugrid, &vlatstatus);
      }        
    }
#endif 

#ifdef MPI2
  MPI_Alloc_mem ( nlocal * sizeof(ttm_Element), MPI_INFO_NULL, &llocal );
#else
  llocal = malloc(nlocal * sizeof(ttm_Element));
#endif //MPI2

  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {
   /* all the "-1" and "-2" because we don't store ghost layers
     and don't want to waste the space in llocal */
    llocal[ i-1 ] = node;    
  }

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

  //-->Das bedeutet: Letztes Element in mpi_elemnt2 muss vom type REAL sein!
  // printf("\n\n\n\nmyid:%d,gather start\n",myid);
  MPI_Gather( llocal, nlocal, mpi_element2,
              lglobal, nlocal, mpi_element2, 0, cpugrid );
// printf("\n\n\n\nmyid:%d,gather end\n",myid);


  if (myid == 0)
  {    
    FILE *outfile;
    char fname[255];
    sprintf(fname, "%s.%d.ttm", outfilename, number);

    outfile = fopen(fname, "w");
    if (NULL == outfile) error ("Cannot open ttm file for writing.\n");

    //Das wird immer rausgeschrieben
    fprintf(outfile,"#x y z natoms temp md_temp U xi source dens vx vy vz fd_k fd_g Z proc Ce");


#ifdef FDTD //Das kommt bei FDTD zusätzlich dazu
    fprintf(outfile, " Ezx Ezy Hx Hy sigmax sigmay Hzx Hzy Ex Ey");
#endif

#ifdef COLRAD
    fprintf(outfile, " P_EE P_EI P_MPI2 P_MPI3 P_RR");
#endif    

    //Linebreak immer am ende der ifdefs
    fprintf(outfile,"\n"); 
    j=k=0;
    for (i = 0; i < global_fd_dim.x; ++i)
    {
          int index=i; //In diesem Fall easy weil Gather nach cpu-nr sortiert

          fprintf(outfile, "%d %d %d %d %e %e %e %e %e %e %e %e %e %e %e %e %d %f",
                  i, j, k, lglobal[index].natoms, lglobal[index].temp,
                  lglobal[index].md_temp, lglobal[index].U, lglobal[index].xi,
                  lglobal[index].source, lglobal[index].dens,
                  lglobal[index].vcomx, lglobal[index].vcomy, lglobal[index].vcomz,
                  lglobal[index].fd_k, lglobal[index].fd_g,
                  lglobal[index].Z, lglobal[index].proc, lglobal[index].Ce);

#ifdef FDTD
          fprintf(outfile, " %e %e %e %e %e %e %e %e %e %e",
                  lglobal[index].Ezx, lglobal[index].Ezy, lglobal[index].Hx, lglobal[index].Hy,
                  lglobal[index].sigmax, lglobal[index].sigmay,
                  lglobal[index].Hzx, lglobal[index].Hzy, lglobal[index].Ex, lglobal[index].Ey
                 );
#endif          
#ifdef COLRAD
          fprintf(outfile," %e %e %e %e %e ",
                lglobal[index].P_EE,lglobal[index].P_EI,lglobal[index].P_MPI2,lglobal[index].P_MPI3, lglobal[index].P_RR);
#endif          
          fprintf(outfile,"\n");
    }

#ifdef VLATTICE
    for(i=0;i<vlatdim;i++)
    {

          fprintf(outfile, "%d %d %d %d %e %e %e %e %e %e %e %e %e %e %e %e %d %f",                  
                  i + last_active_cell_global+1, 0, 0, 
                  -2, //daran erkenne ich virtual-lattice im output file                  
                  vlattice1[i].temp,
                  vlattice1[i].md_temp,
                  0.0, //u
                  0.0, //xi
                  0.0, //source
                  vlatdens,
                  0.0,0.0,0.0, //vcom's
                  vlattice1[i].fd_k,
                  vlattice1[i].fd_g,
                  vlattice1[i].Z,
                  cur_vlattice_proc,
                  vlattice1[i].Ce);
          fprintf(outfile,"\n");
    }
#endif 


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
 * TTM READ for         *
*  restart              *
*  ACHTUNG: bisher nur  *
*  fuer MPI-nutzung     *
*  implementiert,d.h.   *
*  ohne mpi -> crash    *
*************************/
void ttm_read(int number)
{
  int i, j, k;
  int lines = global_fd_dim.x * global_fd_dim.y * global_fd_dim.z;
  ttm_Element *buf;

#ifdef VLATTICE
  lines+=vlatdim;
#endif


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
    char line[600];
    char **tokens;
    size_t j, numtokens;

    // int ig, jg, kg;

    for (i = 0; i < lines + 1; i++) //skip first line (comments)
    {
      //read data
      if (fgets (line, 400, infile) == NULL) {
        char errstr[255];
        sprintf(errstr,"Error Reading ttm-input-file: %s in line %d. Maybe shiftx not deactivated?\n", fname,i);
        error(errstr);
      }
      if (i > 0) //skip first line
      {
        tokens = strsplit(line, ", \t\n", &numtokens);

        sscanf(tokens[3], "%d",  &buf[i - 1].natoms);
        sscanf(tokens[4], "%lf", &buf[i - 1].temp);
        sscanf(tokens[5], "%lf", &buf[i - 1].md_temp);
        sscanf(tokens[6], "%lf", &buf[i - 1].U);
        sscanf(tokens[7], "%lf", &buf[i - 1].xi);
        sscanf(tokens[8], "%lf", &buf[i - 1].source);
        sscanf(tokens[9], "%lf", &buf[i - 1].dens);
        sscanf(tokens[10],"%lf", &buf[i-1].vcomx); 
        sscanf(tokens[11],"%lf", &buf[i-1].vcomy);
        sscanf(tokens[12],"%lf", &buf[i-1].vcomz);
        sscanf(tokens[13], "%lf", &buf[i - 1].fd_k);
        sscanf(tokens[14], "%lf", &buf[i - 1].fd_g);
        sscanf(tokens[15], "%lf", &buf[i - 1].Z);
        sscanf(tokens[16], "%d", &buf[i - 1].proc);
        sscanf(tokens[17], "%lf", &buf[i - 1].Ce);
#ifdef FDTD //TODO:  array-index dynamisch wie in ttm_write...
        sscanf(tokens[18], "%lf", &buf[i - 1].Ezx);
        sscanf(tokens[19], "%lf", &buf[i - 1].Ezy);
        sscanf(tokens[20], "%lf", &buf[i - 1].Hx);
        sscanf(tokens[21], "%lf", &buf[i - 1].Hy);
        sscanf(tokens[22], "%lf", &buf[i - 1].sigmax);
        sscanf(tokens[23], "%lf", &buf[i - 1].sigmay);
        sscanf(tokens[24], "%lf", &buf[i - 1].Hzx);
        sscanf(tokens[25], "%lf", &buf[i - 1].Hzy);
        sscanf(tokens[26], "%lf", &buf[i - 1].Ex);
        sscanf(tokens[27], "%lf", &buf[i - 1].Ey);
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

  MPI_Bcast(buf, lines, mpi_element2, 0, cpugrid);
  int l = 0;
  //int ig,jg,kg;
  i = j = k = 1;
#ifdef VLATTICE
  int v=0;
  last_active_cell_local=last_active_cell_global=99999;
#endif
  // ******************************************
  // *   NOW COPY BUFF TO LATTICE 
  // *****************************************
  for (l = 0; l < lines; l++)
  {

#ifdef VLATTICE
    if(buf[l].natoms==-2) //vlattice zellen
    {
      vlattice1[v].temp=buf[l].temp;
      vlattice1[v].natoms = buf[l].natoms;
      vlattice1[v].md_temp = buf[l].md_temp;
      vlattice1[v].xi = buf[l].xi;
      vlattice1[v].source = buf[l].source;
      vlattice1[v].dens = buf[l].dens;
      vlattice1[v].fd_k = buf[l].fd_k;
      vlattice1[v].fd_g = buf[l].fd_g;
      vlattice1[v].Z = buf[l].Z;
      vlattice1[v].proc = buf[l].proc;
      vlattice1[v].Ce = buf[l].Ce;
      vlattice1[v].vcomx = buf[l].vcomx;
      vlattice1[v].vcomy = buf[l].vcomy;
      vlattice1[v].vcomz = buf[l].vcomz;     
      v++;
      continue;
    }
#endif    
    //ig =  ((i-1) + my_coord.x*(local_fd_dim.x-2));
    //jg =  ((j-1) + my_coord.y*(local_fd_dim.y-2));
    //kg =  ((k-1) + my_coord.z*(local_fd_dim.z-2));
    //printf("l:%d,proc:%d,myid:%d,i:%d,j:%d,k:%d\n",l,buf[l].proc,myid,i,j,k);
#ifdef VLATTICE
    if(buf[l].natoms==-1)
      last_active_cell_local=MIN(last_active_cell_local,l-1);
      last_active_cell_global=last_active_cell_local;
#endif

    if (buf[l].proc == myid)
    {
      node.temp = buf[l].temp;
      node.natoms = buf[l].natoms;

      node.md_temp = buf[l].md_temp;
      node.U = buf[l].U;
      node.xi = buf[l].xi;
      node.source = buf[l].source;
      node.dens = buf[l].dens;
      node.vcomx= buf[l].vcomx;
      node.vcomy= buf[l].vcomy;
      node.vcomz= buf[l].vcomz;
      node.fd_k = buf[l].fd_k;
      node.fd_g = buf[l].fd_g;
      node.Z = buf[l].Z;
      node.proc = buf[l].proc;
      node.Ce = buf[l].Ce;
      node.vcomx = buf[l].vcomx;
      node.vcomy = buf[l].vcomy;
      node.vcomz = buf[l].vcomz;
#ifdef FDTD
      node.Ezx = buf[l].Ezx;
      node.Ezy = buf[l].Ezy;
      node.Hx = buf[l].Hx;
      node.Hy = buf[l].Hy;
      node.eps = eps0; //buf[l].eps;
      node.mu = mu0; //buf[l].mu;
      node.sigmax = buf[l].sigmax;
      node.sigmay = buf[l].sigmay;
      node.Hzx = buf[l].Hzx;
      node.Hzy = buf[l].Hzy;
      node.Ex = buf[l].Ex;
      node.Ey = buf[l].Ey;
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
#ifdef VLATTICE
  cur_vlattice_proc=vlattice1[0].proc;
  old_vlattice_proc=cur_vlattice_proc;
#endif


#ifdef MPI
#ifdef MPI2
  MPI_Free_mem(buf);
#else
  free(buf);
#endif
#endif

  MPI_Barrier(cpugrid);
  ttm_writeout(100000);  //only for debug

}


double fermi_E(double ne)
{                             //IN: [1/m^3] OUT: [J]

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
double MeanCharge(double temp, double rho, double z0, double am, int i, int j, int k) 
{
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
  double r0 = rho / (atomic_charge * am);
  double ytf = t0 / (1.0 + t0);
  double ya = YA1 * pow(t0, YA2) + YA3 * pow(t0, YA4);
  double yb = -exp(YB0 + YB1 * ytf + YB2 * ytf * ytf);

  double yc = YC1 * ytf + YC2;

  double yq1 = ya * pow(r0, yb);

  double yq = pow(pow(r0, yc) + pow(yq1, yc), 1.0 / yc);

  double yx = YAL * pow(yq, YBE);

  double retval = z0 * yx / (1.0 + yx + sqrt(1.0 + 2.0 * yx));
  return retval;

}

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
* MY MOD: Auxiliary function for string Parsing in TTM-Read  
* WIESO MACHE ICH MIR DEN AUFWAND? -->
* Diese routine ist extrem praktisch wenn man dataen einliest,
* von denen man vorher nicht weiss wie genau die struktur aussieht
* Siehe z.b. auch colrad_read, wo die Zahl der spalten davon abhängt,
* wieviele energie-zustände man haben möchte
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
  double khalf;

  for (i = 1; i < local_fd_dim.x - 1; ++i)
  {

        if (node.natoms >= fd_min_atoms && node.dens > RHOMIN)
        {
          // ************************************************
          // * DETERMINE MAX. TIME-STEP FOR STABILITY   *
          // ************************************************
          //x-dir
          if (l1[i + 1].natoms >= fd_min_atoms) imax = i + 1; else imax = i;
          khalf = (node.fd_k + l1[imax].fd_k); //k_(i+1/2)*2
          maxdttmp = MIN(maxdttmp, node.Ce * dxsq / khalf);        
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


double EOS_ee_from_r_te(double r, double t)
{
  //if (r < RHOMIN) return 0;

  double tsqrt = sqrt(t);

  point pout;
  pout.x = r;
  pout.y = tsqrt;

#if DEBUG_LEVEL > 0
if(tsqrt > intp_ee_from_r_tesqrt.ymax || tsqrt < intp_ee_from_r_tesqrt.ymin )
{
  char errstr[400];
  sprintf(errstr, "ERRROR in EOS_ee_from_r_te: tsqrt=%.4e exceeded interpolation range in EOS-table. ymin=%.4e,ymax=%.4e\n",
        tsqrt,intp_ee_from_r_tesqrt.ymin, intp_ee_from_r_tesqrt.ymax);
  printf("myid:%d, WARNING: %s. using bounds.\n",myid,errstr);
  tsqrt=MIN(intp_ee_from_r_tesqrt.ymax,tsqrt);
  tsqrt=MAX(intp_ee_from_r_tesqrt.ymin,tsqrt);
  //error(errstr);
}
if(r < intp_ee_from_r_tesqrt.xmin || r > intp_ee_from_r_tesqrt.xmax)
{
  char errstr[400];
  sprintf(errstr, "ERROR in EOS_ee_from_r_te: Density=%.4e exceeds interpolation range, rmin= %.4e, rmax=%.4e\n",
    r, intp_ee_from_r_tesqrt.xmin,intp_ee_from_r_tesqrt.xmax);
  printf("myid:%d, WARNING: %s. using bounds.\n",myid,errstr);
  r=MIN(intp_ee_from_r_tesqrt.xmax, r);
  r=MAX(r,intp_ee_from_r_tesqrt.xmin);
  //error(errstr);
}
#endif

  // nnhpi_interpolate(intp_e_from_r_tsqrt.interpolator, &pout); //naturla neigh, sibson-rule
  lpi_interpolate_point(intp_ee_from_r_tesqrt.interpolator, &pout); //linear
  //pout.z*=26.9815*AMU*6.2415091E18; // J/kg --> eV/Atom
#if DEBUG_LEVEL > 0
  if(isnan(pout.z)!=0)
  { 
    char errstr[400];
    sprintf(errstr, "ERROR in EOS_ee_from_r_te: ee_from_r_te retunred NaN!.r:%.4e,t:%.4e",r,t);
    error(errstr);
  }
#endif

  return pout.z;
}

double EOS_cve_from_r_te(double r, double t)
{

  //r = MAX(r, RHOMIN);
  //r = MIN(r, 3500); //CHEAT

  //double tsqrt = sqrt(t);
  
  point pout;
  pout.x = r;
  pout.y = t;
  
#if DEBUG_LEVEL > 0   
  if(r < intp_cve_from_r_te.xmin || r > intp_cve_from_r_te.xmax)
  {
    char errstr[400];
    sprintf(errstr, "ERROR in EOS_cve_from_r_te: Density=%.4e  exceeds interpolation range: xmin=%.4e, xmax=%.4e\n",
      r,intp_cve_from_r_te.xmin,intp_cve_from_r_te.xmax);

    printf("myid:%d, WARNING: %s. using bounds.\n",myid,errstr);
    r=MIN(intp_cve_from_r_te.xmax, r);
    r=MAX(r,intp_cve_from_r_te.xmin);    
    //error(errstr);
  }
  if(t < intp_cve_from_r_te.ymin || t > intp_cve_from_r_te.ymax)
  {
    char errstr[400];
    sprintf(errstr, "ERROR in EOS_cve_from_r_te: Te=%.4e  exceeds interpolation range: ymin=%.4e, ymax=%.4e\n",
      t,intp_cve_from_r_te.ymin,intp_cve_from_r_te.ymax);

    printf("myid:%d, WARNING: %s. using bounds.\n",myid,errstr);
    t=MIN(intp_cve_from_r_te.ymax, t);
    t=MAX(t,intp_cve_from_r_te.ymin);      
  }
#endif

  // nnhpi_interpolate(intp_e_from_r_tsqrt.interpolator, &pout); //naturla neigh, sibson-rule
  lpi_interpolate_point(intp_cve_from_r_te.interpolator, &pout); //linear
#if DEBUG_LEVEL > 0
  if(isnan(pout.z)!=0)
  { 
    char errstr[400];
    sprintf(errstr, "ERROR: cve_from_r_te retunred NaN!.r:%.4e,t:%.4e",r,t);
    error(errstr);
  }
  #endif

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
  //r = MAX(r, RHOMIN);
  //r=MIN(r,3500);
#if DEBUG_LEVEL > 0
  if(r < intp_ee_from_r_tesqrt.xmin || r > intp_ee_from_r_tesqrt.xmax)
  {
    char errstr[400];
    sprintf(errstr, "ERROR in EOS_te_from_r_ee: Density=%.4e exceeds interpolation range of table: xmin=%.4e, xmax=%.4e\n",
      r,intp_ee_from_r_tesqrt.xmin, intp_ee_from_r_tesqrt.xmax);

    printf("myid:%d, WARNING: %s. using bounds.\n",myid,errstr);
    r=MIN(r,intp_ee_from_r_tesqrt.xmax);
    r=MAX(r,intp_ee_from_r_tesqrt.xmin);  

    //error(errstr);
  }    
#endif

  //double eSI=e/(26.9815*AMU*6.2415091E18); // eV/Atom ---> J/kg
  double a = pow(intp_ee_from_r_tesqrt.ymin, 2.0);
  double b = pow(intp_ee_from_r_tesqrt.ymax, 2.0);
  double m = fminbnd(a, b, eeminfun, 1e-4, r, e);
  return m;
}

double eeminfun(double x, double r, double e) //Auxiliary function for Te_from_dens_U. This func. is minimized by brent's algo.
{
  return ABS(EOS_ee_from_r_te(r, x) - e);
}


//Selber für FEG-model
double FEG_eeminfun(double x, double r, double ne, double e) 
{
  return ABS(FEG_ee_from_r_ne_te(r,ne,x) - e);
}
double FEG_ee_from_r_ne_te(double r,double ne, double T) //free-electron internal energy from Fermi-integral, T in K
{
  return 0;
  /*
  double mu=chempot(ne,T);
  double F_3_half= gsl_sf_fermi_dirac_3half(mu/BOLTZMAN/T);
  double Gamma_5_half=3.0*sqrt(M_PI)/4.0;
  double C=EMASS/M_PI/M_PI/pow(HBAR,3.0)*sqrt(2*EMASS)*pow(BOLTZMAN*T,2.5)/Gamma_5_half;
  double result=C*F_3_half/r;
  return result;
  */
}
double FEG_te_from_r_ne_ee(double r,double ne,double e)
{
  //tmin=100 ?
  //tmax=1e6?
  double m=fminbnd2(100.0,1e5,FEG_eeminfun,1e-9,r,ne,e);
}





// Wenn FDTD und TTM gemeinsam aktiv, muss in der Diffusion-loop dafür gesorgt werden, dass
// die EM-Feldkomponenten des verangenen Zeitschritts nicht mit den aktuellen vertauscht werden
// -->kopieren und wiederherstellen
#ifdef FDTD
void SwapTTM(int i, int j, int k)
{
  double tmp;
  tmp = node2.Ezx;
  node2.Ezx = node.Ezx;
  node.Ezx = tmp;

  tmp = node2.Ezy;
  node2.Ezy = node.Ezy;
  node.Ezy = tmp;

  tmp = node2.Ex;
  node2.Ex = node.Ex;
  node.Ex = tmp;

  tmp = node2.Ey;
  node2.Ey = node.Ey;
  node.Ey = tmp;

  tmp = node2.Hx;
  node2.Hx = node.Hx;
  node.Hx = tmp;

  tmp = node2.Hy;
  node2.Hy = node.Hy;
  node.Hy = tmp;

  tmp = node2.Hzx;
  node2.Hzx = node.Hzx;
  node.Hzx = tmp;

  tmp = node2.Hzy;
  node2.Hzy = node.Hzy;
  node.Hzy = tmp;

  tmp = node2.Jzx;
  node2.Jzx = node.Jzx;
  node.Jzx = tmp;

  tmp = node2.Jzy;
  node2.Jzy = node.Jzy;
  node.Jzy = tmp;

  tmp = node2.Jx;
  node2.Jx = node.Jx;
  node.Jx = tmp;

  tmp = node2.Jy;
  node2.Jy = node.Jy;
  node.Jy = tmp;


  node2.sigmax = node.sigmax; //noetig?
  node2.sigmay = node.sigmay;
  node2.mu = node.mu;
  node2.eps = node.eps;


  //Lorentz
  tmp = node2.Pzx;
  node2.Pzx = node.Pzx;
  node.Pzx = tmp;

  tmp = node2.Pzy;
  node2.Pzy = node.Pzy;
  node.Pzy = tmp;

  tmp = node2.Px;
  node2.Px = node.Px;
  node.Px = tmp;

  tmp = node2.Py;
  node2.Py = node.Py;
  node.Py = tmp;

  tmp = node2.Pzx;
  node2.Jlzx = node.Jlzx;
  node.Jlzx = tmp;

  tmp = node2.Jlzy;
  node2.Jlzy = node.Jlzy;
  node.Jlzy = tmp;

  tmp = node2.Jlx;
  node2.Jlx = node.Jlx;
  node.Jlx = tmp;

  tmp = node2.Jly;
  node2.Jly = node.Jly;
  node.Jly = tmp;

//        int itmp=node2.natoms_old;
//        node2.natoms_old=node.natoms_old;
//        node.natoms_old=itmp;


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
    node.DL[0]=2.73;
    node.DL[1]=1.1174e+15;
    node.DL[2]=7.6595e+15;
    node.DL[3]=2.4024e+15;
    node.DL[4]=4.5199e+14;
    node.DL[5]=2.2955e+16;
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

  //mpi_element1 wird in der laodbalance variante nicht benötigt
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

    int  blockcountselements; //==alle meine elemnte + 1 für UB
    blockcountselements=15+1;
#ifdef FDTD    
    blockcountselements+=10,
#endif    
#ifdef COLRAD
    blockcountselements+=5; // P_EE P_EI P_MPI2 P_MPI3 P_RR
#endif    

    int *blockcounts;
    alloc1darr(int, blockcounts, blockcountselements);

    MPI_Datatype* types;
    alloc1darr(MPI_Datatype, types, blockcountselements);

    MPI_Aint* displs;
    alloc1darr(MPI_Aint, displs, blockcountselements);

    for(i=0;i<blockcountselements;i++)
      blockcounts[i]=1;    

    //erstmal die blockcount-elemente besetzen, die immer benutzt werden
    types[0]=MPI_INT;   //+1:   natoms
    for(i=1;i<=12;i++)  //+12:  temp, md_temp, xi, U, source, dens, vcom.x, vcom.y, vcom.z, fd_k, fd_g, Z
      types[i]=MPI_DOUBLE;

    types[13]=MPI_INT;    // +1: proc
    types[14]=MPI_DOUBLE; // +1: Ce
    types[blockcountselements-1]=MPI_UB;  //das letzte immer MPI_UB 
    int next_i=15;

#ifdef FDTD
    for(i=next_i; i < next_i + 10; i++) //+10: Ezx Ezy Hx Hy sigmax sigmay Hzx Hzy Ex Ey
      types[i]=MPI_DOUBLE;
    next_i+=10;
#endif    
#ifdef COLRAD
    for(i=next_i; i < next_i + 5; i++)
      types[i]=MPI_DOUBLE;    
#endif    
    MPI_Address(&tmpelement, &tmpaddr);
    MPI_Address(&tmpelement.natoms, &displs[0]);
    MPI_Address(&tmpelement.temp, &displs[1]);
    MPI_Address(&tmpelement.md_temp, &displs[2]);
    MPI_Address(&tmpelement.xi, &displs[3]);
    MPI_Address(&tmpelement.U, &displs[4]);
    MPI_Address(&tmpelement.source, &displs[5]);
    MPI_Address(&tmpelement.dens, &displs[6]);
    MPI_Address(&tmpelement.vcomx, &displs[7]);
    MPI_Address(&tmpelement.vcomy, &displs[8]);
    MPI_Address(&tmpelement.vcomz, &displs[9]);
    MPI_Address(&tmpelement.fd_k, &displs[10]);
    MPI_Address(&tmpelement.fd_g, &displs[11]);
    MPI_Address(&tmpelement.Z, &displs[12]);
    MPI_Address(&tmpelement.proc, &displs[13]);
    MPI_Address(&tmpelement.Ce, &displs[14]);
    next_i=15;
#ifdef FDTD
    MPI_Address(&tmpelement.Ezx, &displs[next_i++]);
    MPI_Address(&tmpelement.Ezy, &displs[next_i++]);
    MPI_Address(&tmpelement.Hx, &displs[next_i++]);
    MPI_Address(&tmpelement.Hy, &displs[next_i++]);
    MPI_Address(&tmpelement.sigmax, &displs[next_i++]);
    MPI_Address(&tmpelement.sigmay, &displs[next_i++]);

    MPI_Address(&tmpelement.Hzx, &displs[next_i++]);
    MPI_Address(&tmpelement.Hzy, &displs[next_i++]);
    MPI_Address(&tmpelement.Ex, &displs[next_i++]);
    MPI_Address(&tmpelement.Ey, &displs[next_i++]);    
#endif
#ifdef COLRAD
    MPI_Address(&tmpelement.P_EE, &displs[next_i++]);
    MPI_Address(&tmpelement.P_EI, &displs[next_i++]);
    MPI_Address(&tmpelement.P_MPI2, &displs[next_i++]);
    MPI_Address(&tmpelement.P_MPI3, &displs[next_i++]);    
    MPI_Address(&tmpelement.P_RR, &displs[next_i++]);    
#endif    


    tmpelement_pointer++;
    MPI_Address(tmpelement_pointer, &displs[blockcountselements-1]);
    for (i = 0; i < blockcountselements; ++i)
    {
      displs[i] -= tmpaddr;
    }

    MPI_Type_struct(blockcountselements, blockcounts, displs, types, &mpi_element2);

    if (MPI_Type_commit(&mpi_element2) != MPI_SUCCESS)
      error("type mpi_element2 failed to commit");
    MPI_Type_contiguous(local_fd_dim.z - 2, mpi_element, &mpi_zrow);
    //MPI_Type_commit(&mpi_zrow);
    if (MPI_Type_commit(&mpi_zrow) != MPI_SUCCESS)
      error("type mpi_zrow failed to commit");
      
  }

  /* datatype for one string of elements along z (short of 2 lattice points) */
  


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
  double buffs[2];
  double buffr[2];

  if(myid != 0 && myid !=num_cpus-1)
  {

    buffs[0]=l1[1].temp;
    buffs[1]=l1[1].fd_k;

    //links-recht
    MPI_Sendrecv(&buffs[0], 2, MPI_DOUBLE, myid-1, 7100,
                 &buffr[0], 2, MPI_DOUBLE, myid+1, 7100,
                 cpugrid, &stati[0]);    

    l1[(local_fd_dim.x - 2) + 1].temp=buffr[0];  
    l1[(local_fd_dim.x - 2) + 1].fd_k=buffr[1];  
    // l1[(local_fd_dim.x - 2) + 1].natoms=(int) buffr[2];  

    //rechts-links
    buffs[0]=l1[(local_fd_dim.x - 2)].temp;
    buffs[1]=l1[(local_fd_dim.x - 2)].fd_k;
    // buffs[2]=(double) l1[(local_fd_dim.x - 2)].natoms;
    // if(l1[(local_fd_dim.x - 2)].dens < RHOMIN) buffs[2]=0;

    MPI_Sendrecv(&buffs[0], 2, MPI_DOUBLE, myid+1, 7200,
                 &buffr[0], 2, MPI_DOUBLE, myid-1, 7200,
                 cpugrid, &stati[1]);      

    l1[0].temp=buffr[0];
    l1[0].fd_k=buffr[1];
    // l1[0].natoms=(int) buffr[2];
  }
  else
  {
    if(myid==0)
    {
      MPI_Recv(&buffr[0], 2, MPI_DOUBLE, myid+1, 7100,
               cpugrid, &stati[0]);   

      l1[(local_fd_dim.x - 2) + 1].temp=buffr[0];
      l1[(local_fd_dim.x - 2) + 1].fd_k=buffr[1];
      l1[(local_fd_dim.x - 2) + 1].natoms=(int) buffr[2];

      //MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
      buffs[0]=l1[(local_fd_dim.x - 2)].temp;
      buffs[1]=l1[(local_fd_dim.x - 2)].fd_k;
      // buffs[2]=(double) l1[(local_fd_dim.x - 2)].natoms;
      // if(l1[(local_fd_dim.x - 2)].dens < RHOMIN) buffs[2]=0;

      MPI_Send(&buffs[0], 2, MPI_DOUBLE, myid+1, 7200,
               cpugrid);      


      l1[0].natoms=0;
    }
    else if(myid==num_cpus-1)
    {
      
      buffs[0]=l1[1].temp;
      buffs[1]=l1[1].fd_k;
      // buffs[2]=(double) l1[1].natoms;
      // if(l1[1].dens < RHOMIN) buffs[2]=0.0;

      MPI_Send(&buffs, 2, MPI_DOUBLE, myid-1, 7100,
               cpugrid);

    
      MPI_Recv(&buffr[0], 2, MPI_DOUBLE, myid-1, 7200,
               cpugrid, &stati[1]);

      l1[0].temp=buffr[0];
      l1[0].fd_k=buffr[1];
      // l1[0].natoms=(int) buffr[2];

      l1[local_fd_dim.x - 1].natoms = 0;
    }
  }
}
#endif /*MPI*/

// **************************************************************************
// *  CHEMICAL POTENTIAL FROM SOMMERFELD EXPANSION
// *  TODO: BRENT ROOT-FINDING --> CHEMPOT FROM NORMALIZATION CONDITION
// **************************************************************************
double chempot(double ne,double Te)
{
  //double EF= HBAR * HBAR * pow(3.0 * M_PI * M_PI * ne, 2.0 / 3.0) / 2.0 / EMASS;
  double EF= 5.842256986370049e-38 *pow(ne,0.666666666666667);
  double mu=EF*(1.0-1.0/3.0*pow(M_PI*BOLTZMAN*Te/2/EF,2.0));
  return mu;

}


double FEG_cve_from_ne_te(double r,double ne,double T)//specific heat of quasi-free electrons
{
  //approximatin according to mazhukin, S. 240 (Fermi-integral analyt. genähert)
  //error less than 5%
  double EF = fermi_E(ne);
  // double Te_K=Te*11605; //in Kelvin
  // double Te_J = T * 1.6021766e-19; //in Joule
  // return 1.5*ne*BOLTZMAN*BOLTZMAN*Te_K/sqrt(Te_J*Te_J+pow(3*EF/M_PI/M_PI,2.0))*7.243227582e-8; //J/K/m^3 -> IMD-UNITS
  // return 2.401087548821963e-49 * T * ne / sqrt(Te_J * Te_J + pow(EF * 0.303963550927013, 2.0)); //alle konstanten zus.gefasst


  // Cv_class=@(T) ne(T).*1.5*boltzman;
  // Cv_deg=@(T) pi^2.*ne(T).*boltzman^2.*T./(2.*EF(T));
  // Cv_mix=@(T) 1./sqrt(1./Cv_deg(T).^2+1./Cv_class(T).^2);   

  double Cv_class=ne*1.5*BOLTZMAN;
  double Cv_deg=M_PI*M_PI*ne*BOLTZMAN*BOLTZMAN*T/2.0/EF;
  double Cv_mix=1.0/sqrt(1.0/Cv_deg/Cv_deg +1/Cv_class/Cv_class);

  double result=Cv_mix;
  result *= r; // J/(K*kg) --> J/K/m^3
  result *= 11604.5; // -->J/eV/m^3
  result *= 1e-30; // --> J/eV/Angs^3
  result *= J2eV; // --> eV/eV/A^3
  return result;


  // double tupper=T+T*0.02;
  // double tlower=T-T*0.02;

  // double eupper=FEG_ee_from_r_ne_te(r,ne,tupper);
  // double elower=FEG_ee_from_r_ne_te(r,ne,tlower);
  // double result=(eupper-elower)/(tupper-tlower);

  // result *= r; // J/(K*kg) --> J/K/m^3
  // result *= 11604.5; // -->J/eV/m^3
  // result *= 1e-30; // --> J/eV/Angs^3
  // result *= J2eV; // --> eV/eV/A^3


  // return result;

}


