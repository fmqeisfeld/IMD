
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_fix_cells_3d.c -- code for fixing cell distribution, three dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#define INDEXED_ACCESS
#include "imd.h"

/******************************************************************************
*
*  fix_cells
*
*  check if each atom is in the correct cell and on the correct CPU;
*  move atoms that have left their cell or CPU
*
*  this also uses Plimpton's communication scheme
*
******************************************************************************/

//ACHTUNG: Dieser Code braucht DRINGEND eine politur..
// Der jetzige Zustand ist katastrophal.
// Mit den ganzen wild in der gegend verstreutgen if..else ..#ifdefs etc. ist es unmöglich damit vernünftig zu arbeiten
// Selbst vim kann die Klammern-Paare nicht mehr zuordnen

//WARNUNG: multijump und loadbalance nicht kompatibel ausser wenn zusätzlich FILTER aktiviert
//         Warum? --> nicht die leiseste Ahnung...
void fix_cells(void)
{
  int i, j, k, l, clone, to_cpu, ii;
  minicell *p, *q;
  ivektor coord, lcoord;
  msgbuf *buf;

//MYMOD
#ifdef LOADBALANCE
#ifdef MULTIJUMP
  int lbmj = 0;
#endif
#endif
//ENDOF MYMID

#ifdef MPI
  empty_mpi_buffers();
#endif

  /* apply periodic boundary conditions */
  do_boundaries();

  /* for each cell in bulk */
  for (i = cellmin.x; i < cellmax.x; ++i)
    for (j = cellmin.y; j < cellmax.y; ++j)
      for (k = cellmin.z; k < cellmax.z; ++k) 
      {

        p = PTR_3D_V(cell_array, i, j, k, cell_dim);
#ifdef LOADBALANCE
        /* Only check content in real cells */
        if (p->lb_cell_type != LB_REAL_CELL) continue;
#endif
        /* loop over atoms in cell */
        l = 0;
        while ( l < p->n )
        {

          coord  = cell_coord( ORT(p, l, X), ORT(p, l, Y), ORT(p, l, Z) );
          lcoord = local_cell_coord( coord );

//MYMOD
#if defined(LOADBALANCE) && defined(MULTIJUMP)
          lbmj = 0;      //my flag
#endif
//ENDOF MYMOD          

#ifdef FILTER

          //if (ORT(p,l,X) < filter_min_x || ORT(p,l,X) > filter_max_x)
          if (FILTERME(p, l) == 1)
          {
            //printf("\n WARNING myid:%d,steps:%d,deleting nr:%d, x:%f, y:%f, z:%f \n",myid,steps,NUMMER(p,l),ORT(p,l,X),ORT(p,l,Y),ORT(p,l,Z));
            to_cpu = cpu_coord(coord);
            buf    = NULL;
            buf = &dump_buf;
            dump_buf.n = 0;
            natoms  -= nclones;
            nactive -= nclones * DIM;
            num_sort [ SORTE(p, l)] -= nclones;
            num_vsort[VSORTE(p, l)] -= nclones;

            if (buf != NULL)
            {
              copy_one_atom( buf, to_cpu, p, l, 1);
#ifdef CLONE
              if (l < p->n - nclones)
                for (clone = 1; clone < nclones; clone++)
                  copy_one_atom( buf, to_cpu, p, l + clone, 1);
              else
                for (clone = 1; clone < nclones; clone++)
                  copy_one_atom( buf, to_cpu, p, l, 1);
#endif  //CLONE
            }
          }
          else //if DELME...
          {
#endif  //FILTER
#ifdef LOADBALANCE
            /* Wrap around pbcs if necessary to get the correct index using the loadbalance scheme*/
            if (cpu_dim.x >= 2 && pbc_dirs.x == 1)
            {
              if (lcoord.x >= global_cell_dim.x) lcoord.x -= global_cell_dim.x;
              if (lcoord.x < 0) lcoord.x += global_cell_dim.x;
            }
            if (cpu_dim.y >= 2 && pbc_dirs.y == 1)
            {
              if (lcoord.y >= global_cell_dim.y) lcoord.y -= global_cell_dim.y;
              if (lcoord.y < 0) lcoord.y += global_cell_dim.y;
            }
            if (cpu_dim.z >= 2 && pbc_dirs.z == 1)
            {
              if (lcoord.z >= global_cell_dim.z) lcoord.z -= global_cell_dim.z;
              if (lcoord.z < 0) lcoord.z += global_cell_dim.z;
            }
#endif //LB

            /* see if atom is in wrong cell */
            if ((lcoord.x == i) && (lcoord.y == j) && (lcoord.z == k))
            {
// if(NUMMER(p,l)==35)
// printf("myid:%d, nr:%d in wrong coord\n",myid,NUMMER(p,l));
              l++;

            }
            else
            {
#ifdef LOADBALANCE
              if (lcoord.x < 0 || lcoord.x >= cell_dim.x ||
                  lcoord.y < 0 || lcoord.y >= cell_dim.y ||
                  lcoord.z < 0 || lcoord.z >= cell_dim.z ||
                  (PTR_VV(cell_array, lcoord, cell_dim))->lb_cell_type == LB_EMPTY_CELL)
              {
                //MYMOD :
#ifndef MULTIJUMP
                char errstr[256];
                sprintf(errstr, "LB: Illegal cell accessed, Atom %d jumped multiple CPUs\n", NUMMER(p, l));
                //error("LB: Illegal cell accessed, Atom jumped multiple CPUs");
                error(errstr);
#else
                lbmj = 1;                
                warning("Atom jumped multiple CPUs");
#endif
                //ENDOF MYMOD
              }
              to_cpu = (PTR_VV(cell_array, lcoord, cell_dim))->lb_cpu_affinity;

#else //LB
              to_cpu = cpu_coord(coord);
#endif //LB
              buf    = NULL;

              /* atom is on my cpu */
              if (to_cpu == myid)
              {
                q = PTR_VV(cell_array, lcoord, cell_dim);
                MOVE_ATOM(q, p, l);
#ifdef CLONE
                if (l < p->n - nclones)
                  for (clone = 1; clone < nclones; clone++)
                    MOVE_ATOM(q, p, l + clone);
                else /* we are dealing with the last in the stack */
                  for (clone = 1; clone < nclones; clone++)
                    MOVE_ATOM(q, p, l);
#endif //CLINE
              }
#ifdef MPI
#ifdef LOADBALANCE
              else
              {
                //MYMOD
#if defined(MULTIJUMP)
                //overwrite buf with dump_buf if multijump
                if (lbmj)
                {
                  if (l < p->n) copy_atom_cell_cell(p, l, p, p->n); //overwrite
                  p->n--;
                  natoms-=1;
                  nactive-=DIM;
                  buf=NULL;
                  lbmj=0;
                }
                else //normal..c&p die nächste zeile
                  buf = &lb_send_buf[(PTR_VV(cell_array, lcoord, cell_dim))->lb_neighbor_index];
#else //MJ
                buf = &lb_send_buf[(PTR_VV(cell_array, lcoord, cell_dim))->lb_neighbor_index]; //das stand vorher schon da.
#endif //MJ
                //ENDOF MYMOD
              }
#else //LB
              //folgende elseif's nur wenn nicht LB
              /* west */
              else if ((cpu_dim.x > 1) && //to cpu==myid
                       ((to_cpu == nbwest) || (to_cpu == nbnw)  || (to_cpu == nbws) ||
                        (to_cpu == nbuw  ) || (to_cpu == nbunw) || (to_cpu == nbuws) ||
                        (to_cpu == nbdw  ) || (to_cpu == nbdwn) || (to_cpu == nbdsw)))
              {
                buf = &send_buf_west;
              }

              /* east */
              else if ((cpu_dim.x > 1) &&
                       ((to_cpu == nbeast) || (to_cpu == nbse)  || (to_cpu == nben) ||
                        (to_cpu == nbue  ) || (to_cpu == nbuse) || (to_cpu == nbuen) ||
                        (to_cpu == nbde  ) || (to_cpu == nbdes) || (to_cpu == nbdne)))
              {
                buf = &send_buf_east;
              }

              /* south  */
              else if ((cpu_dim.y > 1) &&
                       ((to_cpu == nbsouth) || (to_cpu == nbus)  || (to_cpu == nbds)))
              {
                buf = &send_buf_south;
              }

              /* north  */
              else if ((cpu_dim.y > 1) &&
                       ((to_cpu == nbnorth) || (to_cpu == nbun)  || (to_cpu == nbdn)))
              {
                buf = &send_buf_north;
              }

              /* down  */
              else if ((cpu_dim.z > 1) && (to_cpu == nbdown))
              {
                buf = &send_buf_down;
              }

              /* up  */
              else if ((cpu_dim.z > 1) && (to_cpu == nbup))
              {
                buf = &send_buf_up;
              }

              else
              {
#if defined(SHOCK) || defined(MULTIJUMP)
                /* remove atom from simulation */
                buf = &dump_buf;
                dump_buf.n = 0;
                natoms  -= nclones;
                nactive -= nclones * DIM;
                num_sort [ SORTE(p, l)] -= nclones;
                num_vsort[VSORTE(p, l)] -= nclones;
                warning("Atom jumped multiple CPUs");
#else
                error("Atom jumped multiple CPUs");
#endif
              }
#endif /* not LOADBALANCE*/

              if (buf != NULL) 
              {
                copy_one_atom( buf, to_cpu, p, l, 1);
#ifdef CLONE
                if (l < p->n - nclones)
                  for (clone = 1; clone < nclones; clone++)
                    copy_one_atom( buf, to_cpu, p, l + clone, 1);
                else /* we are dealing with the last in the stack */
                  for (clone = 1; clone < nclones; clone++)
                    copy_one_atom( buf, to_cpu, p, l, 1);
#endif
              }
#ifdef FILTER
            } //if DELME ... else...?
#endif

#endif /* MPI */
          } //???
        } //while
      } // for k

#ifdef MPI
  /* send atoms to neighbbour CPUs */
  send_atoms();
#endif

#ifdef NBLIST
  /* tag neighbor list as outdated */
  have_valid_nbl = 0;
#endif

}

#ifdef MPI

#ifdef LOADBALANCE
void send_atoms()
{
  int i;
  MPI_Status stat;

  for (i = 0; i < lb_nTotalComms; ++i) {
    /*Send data away*/
    isend_buf(&lb_send_buf[i], lb_commIndexToCpu[i], &lb_req_send[i]);
    lb_requests[2 * i] = lb_req_send[i];
    lb_request_indices[2 * i] = -1; /* Indicates no processing required */
    /*Start receiving data*/
    irecv_buf(&lb_recv_buf[i], lb_commIndexToCpu[i], &lb_req_recv[i]);
    lb_requests[2 * i + 1] = lb_req_recv[i];
    lb_request_indices[2 * i + 1] = i;
  }

  /*Receive and process data as soon as something is available*/
  for (i = 2 * lb_nTotalComms; i > 0; i--) {
    int finished;
    MPI_Waitany(i, lb_requests, &finished, &stat);
    int ind = lb_request_indices[finished];
    if (ind != -1) {
      MPI_Get_count(&stat, REAL, &lb_recv_buf[ind].n);
      process_buffer( &lb_recv_buf[ind]);
    }
    lb_requests[finished] = lb_requests[i - 1];
    lb_request_indices[finished] = lb_request_indices[i - 1];
  }
}
#else

#ifdef SR

/******************************************************************************
*
* send_atoms  -  only used for fix_cells
*
******************************************************************************/

void send_atoms()
{
  MPI_Status  stat;

  if (cpu_dim.x > 1) {
    /* send east, receive west, move atoms from west to cells */
    sendrecv_buf( &send_buf_east, nbeast, &recv_buf_west, nbwest, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_west.n );
    process_buffer( &recv_buf_west );

    /* send west, receive east, move atoms from east to cells */
    sendrecv_buf( &send_buf_west, nbwest, &recv_buf_east, nbeast, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_east.n );
    process_buffer( &recv_buf_east );

    if (cpu_dim.y > 1) {
      /* append atoms from west and east to north send buffer */
      copy_atoms_buf( &send_buf_north, &recv_buf_west );
      copy_atoms_buf( &send_buf_north, &recv_buf_east );
      /* check special case cpu_dim.y==2 */
      if (nbsouth != nbnorth) {
        /* append atoms from west and east to south send buffer */
        copy_atoms_buf( &send_buf_south, &recv_buf_west );
        copy_atoms_buf( &send_buf_south, &recv_buf_east );
      }
    }

    if (cpu_dim.z > 1) {
      /* append atoms from west and east to up send buffer */
      copy_atoms_buf( &send_buf_up, &recv_buf_east );
      copy_atoms_buf( &send_buf_up, &recv_buf_west );
      /* check special case cpu_dim.z==2 */
      if (nbdown != nbup) {
        /* append atoms from west and east to down send buffer */
        copy_atoms_buf( &send_buf_down, &recv_buf_east );
        copy_atoms_buf( &send_buf_down, &recv_buf_west );
      }
    }
  }

  if (cpu_dim.y > 1) {
    /* send north, receive south, move atoms from south to cells */
    sendrecv_buf(  &send_buf_north, nbnorth, &recv_buf_south, nbsouth, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_south.n );
    process_buffer( &recv_buf_south );

    /* send south, receive north, move atoms from north to cells */
    sendrecv_buf( &send_buf_south, nbsouth, &recv_buf_north, nbnorth, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_north.n );
    process_buffer( &recv_buf_north );

    if (cpu_dim.z > 1) {
      /* append atoms from north and south to up send buffer */
      copy_atoms_buf( &send_buf_up, &recv_buf_north );
      copy_atoms_buf( &send_buf_up, &recv_buf_south );
      /* check special case cpu_dim.z==2 */
      if (nbdown != nbup) {
        /* append atoms from north and south to down send buffer */
        copy_atoms_buf( &send_buf_down, &recv_buf_north );
        copy_atoms_buf( &send_buf_down, &recv_buf_south );
      }
    }
  }

  if (cpu_dim.z > 1) {
    /* send up, receive down, move atoms from down to cells */
    sendrecv_buf( &send_buf_up, nbup, &recv_buf_down, nbdown, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_down.n );
    process_buffer( &recv_buf_down );

    /* send down, receive up, move atoms from up to cells */
    sendrecv_buf( &send_buf_down, nbdown, &recv_buf_up, nbup, &stat );
    MPI_Get_count( &stat, REAL, &recv_buf_up.n );
    process_buffer( &recv_buf_up );
  }

}

#else /* not SR */

/******************************************************************************
*
* send_atoms  -  only used for fix_cells
*
******************************************************************************/

void send_atoms()
{
  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  if (cpu_dim.x > 1) {
    /* send east */
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );

    /* send west */
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );

    /* Wait for atoms from west, move them to cells */
    MPI_Waitall(2, reqwest, statwest);
    MPI_Get_count( &statwest[1], REAL, &recv_buf_west.n );
    process_buffer( &recv_buf_west );

    /* Wait for atoms from east, move them to cells */
    MPI_Waitall(2, reqeast, stateast);
    MPI_Get_count( &stateast[1], REAL, &recv_buf_east.n );
    process_buffer( &recv_buf_east );

    if (cpu_dim.y > 1) {
      /* append atoms from west and east to north send buffer */
      copy_atoms_buf( &send_buf_north, &recv_buf_west );
      copy_atoms_buf( &send_buf_north, &recv_buf_east );
      /* check special case cpu_dim.y==2 */
      if (nbsouth != nbnorth) {
        /* append atoms from west and east to south send buffer */
        copy_atoms_buf( &send_buf_south, &recv_buf_east );
        copy_atoms_buf( &send_buf_south, &recv_buf_west );
      }
    }

    if (cpu_dim.z > 1) {
      /* append atoms from west and east to up send buffer */
      copy_atoms_buf( &send_buf_up, &recv_buf_east );
      copy_atoms_buf( &send_buf_up, &recv_buf_west );
      /* check special case cpu_dim.z==2 */
      if (nbdown != nbup) {
        /* append atoms from west and east to down send buffer */
        copy_atoms_buf( &send_buf_down, &recv_buf_east );
        copy_atoms_buf( &send_buf_down, &recv_buf_west );
      }
    }
  }

  if (cpu_dim.y > 1) {
    /* Send atoms north */
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );

    /* Send atoms south */
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );

    /* Wait for atoms from south, move them to cells */
    MPI_Waitall(2, reqsouth, statsouth);
    MPI_Get_count( &statsouth[1], REAL, &recv_buf_south.n );
    process_buffer( &recv_buf_south );

    /* Wait for atoms from north, move them to cells */
    MPI_Waitall(2, reqnorth, statnorth);
    MPI_Get_count( &statnorth[1], REAL, &recv_buf_north.n );
    process_buffer( &recv_buf_north );

    if (cpu_dim.z > 1) {
      /* append atoms from north and south to up send buffer */
      copy_atoms_buf( &send_buf_up, &recv_buf_north );
      copy_atoms_buf( &send_buf_up, &recv_buf_south );
      /* check special case cpu_dim.z==2 */
      if (nbdown != nbup) {
        /* append atoms from north and south to down send buffer */
        copy_atoms_buf( &send_buf_down, &recv_buf_north );
        copy_atoms_buf( &send_buf_down, &recv_buf_south );
      }
    }
  }

  if (cpu_dim.z > 1) {
    /* Send atoms up */
    irecv_buf( &recv_buf_down , nbdown, &reqdown[1]);
    isend_buf( &send_buf_up   , nbup  , &reqdown[0]);

    /* Send atoms down */
    irecv_buf( &recv_buf_up  , nbup  , &requp[1] );
    isend_buf( &send_buf_down, nbdown, &requp[0] );

    /* Wait for atoms from down, move them to cells */
    MPI_Waitall(2, reqdown, statdown);
    MPI_Get_count( &statdown[1], REAL, &recv_buf_down.n );
    process_buffer( &recv_buf_down );

    /* Wait for atoms from up, move them to cells */
    MPI_Waitall(2, requp, statup);
    MPI_Get_count( &statup[1], REAL, &recv_buf_up.n );
    process_buffer( &recv_buf_up );
  }

}

#endif /* not SR */

#endif /* not LOADBALANCE*/

#endif /* MPI */
