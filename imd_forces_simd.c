
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2011 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/*****************************************************************************
*
* imd_forces -- force loop
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#include "potaccess.h"

/******************************************************************************
*
*  do_forces, version for scalar processors
*
*  computes the forces between atoms in two given cells
*
******************************************************************************/

void do_forces(cell *p, cell *q, vektor pbc, real *Epot, real *Virial, 
               real *Vir_xx, real *Vir_yy, real *Vir_zz,
               real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  int i,j,k;
  vektor d;
  vektor tmp_d;
  vektor force;
  real r2, rho_h;
  real tmp_virial;
#ifdef P_AXIAL
  vektor tmp_vir_vect;
#endif
  real pot_zwi, pot_grad;
  int col, col2, is_short=0, inc = ntypes * ntypes;
  int jstart, q_typ, p_typ;
  real *qptr, *pfptr, *qfptr, *qpdptr, *ppdptr, *qpoptr, *ppoptr;

  tmp_virial     = 0.0;
#ifdef P_AXIAL
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
  tmp_vir_vect.z = 0.0;
#endif

  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = ORT(p,i,X) - pbc.x;
    tmp_d.y = ORT(p,i,Y) - pbc.y;
    tmp_d.z = ORT(p,i,Z) - pbc.z;

    p_typ  = SORTE(p,i);
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0))               ? i+1 : 0);
    qptr   = &ORT(q,jstart,X);
    /* for each atom in neighbouring cell */

    // *****************************
    // * TRY TO VECTORIZE THIS
    // ****************************
    const int nn=q->n;
    for (j = jstart; j < nn; ++j) 
    {
      /* calculate distance */

      d.x = *qptr - tmp_d.x; ++qptr;
      d.y = *qptr - tmp_d.y; ++qptr;
      d.z = *qptr - tmp_d.z; ++qptr;


      //q_typ = SORTE(q,j);
      q_typ=0;
/*
      col   = p_typ * ntypes + q_typ;
      col2  = q_typ * ntypes + p_typ;
      r2    = SPROD(d,d);
*/

#if defined(PAIR) || defined(KEATING)
      /* PAIR and KEATING are mutually exclusive */
#if defined(PAIR)
//      if (r2 <= pair_pot.end[col]) {
#ifdef LINPOT
//        PAIR_INT_LIN(pot_zwi, pot_grad, pair_pot_lin, col, inc, r2, is_short)
#else
//        PAIR_INT(pot_zwi, pot_grad, pair_pot, col, inc, r2, is_short)
#endif
#endif

        /* store force in temporary variable */
/*
        force.x = d.x * pot_grad;
        force.y = d.y * pot_grad;
        force.z = d.z * pot_grad;
*/

        /* accumulate forces */
        pfptr = &KRAFT(p,i,X);
        qfptr = &KRAFT(q,j,X);
/*
        *pfptr     += force.x; 
        *qfptr     -= force.x; 
        *(++pfptr) += force.y; 
        *(++qfptr) -= force.y; 
        *(++pfptr) += force.z; 
        *(++qfptr) -= force.z; 
*/
//        *Epot      += pot_zwi;

 //       pot_zwi *= 0.5;   /* avoid double counting */
 //       POTENG(p,i) += pot_zwi;
 //       POTENG(q,j) += pot_zwi;

/*
#ifdef P_AXIAL
        tmp_vir_vect.x -= d.x * force.x;
        tmp_vir_vect.y -= d.y * force.y;
        tmp_vir_vect.z -= d.z * force.z;
#else
        tmp_virial     -= r2 * pot_grad;
#endif
*/

//      }
#endif /* PAIR || KEATING */

#ifdef EAM2
      /* compute host electron density */
/*
      if (r2 < rho_h_tab.end[col])  {

        VAL_FUNC(rho_h, rho_h_tab, col,  inc, r2, is_short);
        EAM_RHO(p,i) += rho_h; 
#ifdef EEAM
        EAM_P(p,i) += rho_h*rho_h; 
#endif
      }
      if (p_typ==q_typ) {
        if (r2 < rho_h_tab.end[col]) 
          {EAM_RHO(q,j) += rho_h;
#ifdef EEAM
           EAM_P(q,j) += rho_h*rho_h;
#endif
          }
      } else {
        col2 = q_typ * ntypes + p_typ;
        if (r2 < rho_h_tab.end[col2]) {
          VAL_FUNC(rho_h, rho_h_tab, col2, inc, r2, is_short);
          EAM_RHO(q,j) += rho_h; 
#ifdef EEAM
          EAM_P(q,j) += rho_h*rho_h; 
#endif
        }
      }
*/
#endif
    } /* for j */
  } /* for i */

#ifdef DEBUG
  if (is_short==1) printf("\n Short distance!\n");
#endif
#ifdef P_AXIAL
  *Vir_xx += tmp_vir_vect.x;
  *Vir_yy += tmp_vir_vect.y;
  *Virial += tmp_vir_vect.x;
  *Virial += tmp_vir_vect.y;
  *Vir_zz += tmp_vir_vect.z;
  *Virial += tmp_vir_vect.z;
#else
  *Virial += tmp_virial;
#endif 

}
