
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
#ifndef SIMD
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
#ifndef TWOD
  tmp_vir_vect.z = 0.0;
#endif
#endif

//MYMOD
#ifdef LOD
  vektor qvec; //for summation over n.n.wave-vectors
  vektor d2;
  double qvec_fac=4.0*M_PI/alat;
#endif 
#ifdef NRB
  //displacement vectors for NRB
  vektor U_self; //self
  ///yx-plane
  vektor U_0; //neighbor at -1,+1,0
  vektor U_1; //-1,-1,0
  vektor U_2; //+1,-1,0
  vektor U_3; //+1,+1,0
  vektor V_0;
  vektor V_1;
  vektor V_2;
  vektor V_3;
  //yz-plane
  vektor U_4; //0,-1,1
  vektor U_5; //0,-1,-1
  vektor U_6; //0,+1,-1
  vektor U_7; //0,+1,+1
  vektor V_4;
  vektor V_5;
  vektor V_6;
  vektor V_7;
  //xz-plane
  vektor U_8; //-1,0,1
  vektor U_9; //-1,0,-1
  vektor U_10;//+1,0,-1
  vektor U_11;//+1,0,+1
  vektor V_8;
  vektor V_9;
  vektor V_10;
  vektor V_11;
 
  vektor U_dot;
  int	nrbactive;
#endif
//ENDOF MYMOD
    
  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = ORT(p,i,X) - pbc.x;
    tmp_d.y = ORT(p,i,Y) - pbc.y;
#ifndef TWOD
    tmp_d.z = ORT(p,i,Z) - pbc.z;
#endif

    p_typ  = SORTE(p,i);
#ifdef TWOD
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0))               ? i+1 : 0);
#else
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0)) ? i+1 : 0);
#endif
    qptr   = &ORT(q,jstart,X);
//MYMOD
#ifdef NRB
    if(NRBI(p,i,12)>0)
    {
      nrbactive=NRBI(p,i,12);
      U_self.x=REF_POS(p,i,X)-ORT(p,i,X);
      U_self.y=REF_POS(p,i,Y)-ORT(p,i,Y);
      U_self.z=REF_POS(p,i,Z)-ORT(p,i,Z);

      //MINIMUM IMAGE CONVENTION
      U_self.x=MINIMGX(U_self.x);
      U_self.y=MINIMGY(U_self.y);
      U_self.z=MINIMGZ(U_self.z);
  
 
      //Self-displacement beitrag
      //+x boundary (yz-plane, wave from -x)
      if(nrbactive==1)
      {
        U_dot.x=nrbk*4.0*U_self.x;
        U_dot.y=nrbk*dblsqrt2*U_self.y;
        U_dot.z=nrbk*dblsqrt2*U_self.z;
	IMPULS(p,i,X)+=U_dot.x*MASSE(p,i);
	IMPULS(p,i,Y)+=U_dot.y*MASSE(p,i);
	IMPULS(p,i,Z)+=U_dot.z*MASSE(p,i);
      }
    }
    else
      nrbactive=0;
#endif
    /* for each atom in neighbouring cell */
    for (j = jstart; j < q->n; ++j) 
    {
      /* calculate distance */
      d.x = *qptr - tmp_d.x; ++qptr;
      d.y = *qptr - tmp_d.y; ++qptr;
#ifndef TWOD
      d.z = *qptr - tmp_d.z; ++qptr;
#endif

      q_typ = SORTE(q,j);
      col   = p_typ * ntypes + q_typ;
      col2  = q_typ * ntypes + p_typ;
      r2    = SPROD(d,d);

#ifdef DEBUG
      if (0==r2) { char msgbuf[256];
        sprintf(msgbuf, "Distance is zero between particles %d and %d!\n",
                NUMMER(p,i), NUMMER(q,j));
        error(msgbuf);
      }
#endif

      /* compute pair interactions */
#if defined(PAIR) || defined(KEATING)
      /* PAIR and KEATING are mutually exclusive */
#if defined(PAIR)
      if (r2 <= pair_pot.end[col]) {

//MYMOD
/*
#ifdef NRB
      // ACHTUNG: die aller aussersten nrb-atome haben keine 12 neighs!
      // Hat z.B. das atom i keinen nachbarn auf position 3, so ist NRBI(p,i,3)==0 
      // In diesem fall wird die fallunterscheidung keinen Beitrag von diesem nachbar-atom liefern
      // frage: Muss der entsprechende self-beitrag dann auch entfernt werden?

      // nrbactive=1 beudeted es handelt sich um ein +x-bnd-atom, also schockwelle kommt von -x auf yz-Ebene
      // die noetigen nachbarn sind auf positionen: 0,1,8,9
      if(nrbactive==1) 
      {
	//Now check for neighs
	if(NUMMER(q,j)==NRBI(q,j,0))
	{
	  //calc displacement from refpos
	  U_0.x=REF_POS(q,j,X)-ORT(q,j,X);
          U_0.y=REF_POS(q,j,Y)-ORT(q,j,Y);
          U_0.z=REF_POS(q,j,Z)-ORT(q,j,Z);
	  U_0.x=MINIMGX(U_0.x);
	  U_0.y=MINIMGY(U_0.y);
          U_0.z=MINIMGZ(U_0.z);
	
	  V_0.x=IMPULS(q,j,X)/MASSE(q,j);
	  V_0.y=IMPULS(q,j,Y)/MASSE(q,j);
	  V_0.z=IMPULS(q,j,Z)/MASSE(q,j);
	  
	  U_dot.x=nrbk*U_0.x;
 	  U_dot.y=nrbk*sqrt2half*U_0.y;
	  U_dot.z=nrbk*sqrt2half*U_0.z;
          U_dot.x=-0.25*V_0.x;
          U_dot.y=-0.25*V_0.y;
          U_dot.z=-0.25*V_0.z;

          IMPULS(p,i,X)+=U_dot.x*MASSE(p,i);
          IMPULS(p,i,Y)+=U_dot.y*MASSE(p,i);
          IMPULS(p,i,Z)+=U_dot.z*MASSE(p,i);

	}
	else if(NUMMER(q,j)==NRBI(q,j,1))
	{
          U_1.x=REF_POS(q,j,X)-ORT(q,j,X);
          U_1.y=REF_POS(q,j,Y)-ORT(q,j,Y);
          U_1.z=REF_POS(q,j,Z)-ORT(q,j,Z);
          U_1.x=MINIMGX(U_1.x);
          U_1.y=MINIMGY(U_1.y);
          U_1.z=MINIMGZ(U_1.z);

          V_1.x=IMPULS(q,j,X)/MASSE(q,j);
          V_1.y=IMPULS(q,j,Y)/MASSE(q,j);
          V_1.z=IMPULS(q,j,Z)/MASSE(q,j);

          U_dot.x=nrbk*U_1.x;
          U_dot.y=nrbk*sqrt2half*U_1.y;
          U_dot.z=nrbk*sqrt2half*U_1.z;
          U_dot.x=-0.25*V_1.x;
          U_dot.y=-0.25*V_1.y;
          U_dot.z=-0.25*V_1.z;

          IMPULS(p,i,X)+=U_dot.x*MASSE(p,i);
          IMPULS(p,i,Y)+=U_dot.y*MASSE(p,i);
          IMPULS(p,i,Z)+=U_dot.z*MASSE(p,i);

	}
	else if(NUMMER(q,j)==NRBI(q,j,8))
	{
          U_8.x=REF_POS(q,j,X)-ORT(q,j,X);
          U_8.y=REF_POS(q,j,Y)-ORT(q,j,Y);
          U_8.z=REF_POS(q,j,Z)-ORT(q,j,Z);
          U_8.x=MINIMGX(U_8.x);
          U_8.y=MINIMGY(U_8.y);
          U_8.z=MINIMGZ(U_8.z);

          V_8.x=IMPULS(q,j,X)/MASSE(q,j);
          V_8.y=IMPULS(q,j,Y)/MASSE(q,j);
          V_8.z=IMPULS(q,j,Z)/MASSE(q,j);
          
          U_dot.x=nrbk*U_8.x;
          U_dot.y=nrbk*sqrt2half*U_8.y;
          U_dot.z=nrbk*sqrt2half*U_8.z;
          U_dot.x=-0.25*V_8.x;
          U_dot.y=-0.25*V_8.y;
          U_dot.z=-0.25*V_8.z;

          IMPULS(p,i,X)+=U_dot.x*MASSE(p,i);
          IMPULS(p,i,Y)+=U_dot.y*MASSE(p,i);
          IMPULS(p,i,Z)+=U_dot.z*MASSE(p,i);
	}
	else if(NUMMER(q,j)==NRBI(q,j,9))
	{
          U_9.x=REF_POS(q,j,X)-ORT(q,j,X);
          U_9.y=REF_POS(q,j,Y)-ORT(q,j,Y);
          U_9.z=REF_POS(q,j,Z)-ORT(q,j,Z);
          U_9.x=MINIMGX(U_9.x);
          U_9.y=MINIMGY(U_9.y);
          U_9.z=MINIMGZ(U_9.z);

          V_9.x=IMPULS(q,j,X)/MASSE(q,j);
          V_9.y=IMPULS(q,j,Y)/MASSE(q,j);
          V_9.z=IMPULS(q,j,Z)/MASSE(q,j);
          
          U_dot.x=nrbk*U_9.x;
          U_dot.y=nrbk*sqrt2half*U_9.y;
          U_dot.z=nrbk*sqrt2half*U_9.z;
          U_dot.x=-0.25*V_9.x;
          U_dot.y=-0.25*V_9.y;
          U_dot.z=-0.25*V_9.z;

          IMPULS(p,i,X)+=U_dot.x*MASSE(p,i);
          IMPULS(p,i,Y)+=U_dot.y*MASSE(p,i);
          IMPULS(p,i,Z)+=U_dot.z*MASSE(p,i);

	}
      }
      // *****************************************************
      // * -y boundary, d.h. welle kommt von +y auf xz-Ebene
      // neigh-positionen sind 0,3,7,6
      // ******************************************************   
      else if(nrbactive==2) 
      {
	if(NUMMER(q,j)==NRBI(q,j,0))
	{
	  //calc displacement from refpos
	  U_0.x=REF_POS(q,j,X)-ORT(q,j,X);
          U_0.y=REF_POS(q,j,Y)-ORT(q,j,Y);
          U_0.z=REF_POS(q,j,Z)-ORT(q,j,Z);
	  U_0.x=MINIMGX(U_0.x);
	  U_0.y=MINIMGY(U_0.y);
          U_0.z=MINIMGZ(U_0.z);

	  U_dotdot.x+=nrbk*(U_0.x-U_0.y);
 	  U_dotdot.y+=nrbk*(U_0.y-U_0.x);
	}
	else if(NUMMER(q,j)==NRBI(q,j,3))
	{
          U_3.x=REF_POS(q,j,X)-ORT(q,j,X);
          U_3.y=REF_POS(q,j,Y)-ORT(q,j,Y);
          U_3.z=REF_POS(q,j,Z)-ORT(q,j,Z);
          U_3.x=MINIMGX(U_3.x);
          U_3.y=MINIMGY(U_3.y);
          U_3.z=MINIMGZ(U_3.z);

          U_dotdot.x+=nrbk*(U_3.x+U_3.y);
          U_dotdot.y+=nrbk*(U_3.x+U_3.y);
	}
	else if(NUMMER(q,j)==NRBI(q,j,6))
	{
          U_6.x=REF_POS(q,j,X)-ORT(q,j,X);
          U_6.y=REF_POS(q,j,Y)-ORT(q,j,Y);
          U_6.z=REF_POS(q,j,Z)-ORT(q,j,Z);
          U_6.x=MINIMGX(U_6.x);
          U_6.y=MINIMGY(U_6.y);
          U_6.z=MINIMGZ(U_6.z);

          U_dotdot.y+=nrbk*(U_6.y-U_6.z);
          U_dotdot.z+=nrbk*(U_6.z-U_6.y);
	}
	else if(NUMMER(q,j)==NRBI(q,j,7))
	{
          U_7.x=REF_POS(q,j,X)-ORT(q,j,X);
          U_7.y=REF_POS(q,j,Y)-ORT(q,j,Y);
          U_7.z=REF_POS(q,j,Z)-ORT(q,j,Z);
          U_7.x=MINIMGX(U_7.x);
          U_7.y=MINIMGY(U_7.y);
          U_7.z=MINIMGZ(U_7.z);

          U_dotdot.y+=nrbk*(U_7.y+U_7.z);
          U_dotdot.z+=nrbk*(U_7.z+U_7.y);
	}
      }
      // *****************************************************
      // * +y boundary, d.h. welle kommt von -y auf xz-Ebene
      // neigh-positionen sind 1,2,4,5
      // ****************************************************** 
      else if(nrbactive==2)
      {
 	if(NUMMER(q,j)==NRBI(q,j,1))
	{
          U_1.x=REF_POS(q,j,X)-ORT(q,j,X);
          U_1.y=REF_POS(q,j,Y)-ORT(q,j,Y);
          U_1.z=REF_POS(q,j,Z)-ORT(q,j,Z);
          U_1.x=MINIMGX(U_1.x);
          U_1.y=MINIMGY(U_1.y);
          U_1.z=MINIMGZ(U_1.z);

          U_dotdot.x+=nrbk*(U_1.x+U_1.y);
          U_dotdot.y+=nrbk*(U_1.x+U_1.y);
	}

	else if(NUMMER(q,j)==NRBI(q,j,2))
	{
          U_2.x=REF_POS(q,j,X)-ORT(q,j,X);
          U_2.y=REF_POS(q,j,Y)-ORT(q,j,Y);
          U_2.z=REF_POS(q,j,Z)-ORT(q,j,Z);
          U_2.x=MINIMGX(U_2.x);
          U_2.y=MINIMGY(U_2.y);
          U_2.z=MINIMGZ(U_2.z);

          U_dotdot.x+=nrbk*(U_2.x-U_2.y);
          U_dotdot.y+=nrbk*(U_2.y-U_2.x);
	}
	else if(NUMMER(q,j)==NRBI(q,j,4))
	{
          U_4.x=REF_POS(q,j,X)-ORT(q,j,X);
          U_4.y=REF_POS(q,j,Y)-ORT(q,j,Y);
          U_4.z=REF_POS(q,j,Z)-ORT(q,j,Z);
          U_4.x=MINIMGX(U_4.x);
          U_4.y=MINIMGY(U_4.y);
          U_4.z=MINIMGZ(U_4.z);

          U_dotdot.y+=nrbk*(U_4.y-U_4.z);
          U_dotdot.z+=nrbk*(U_4.z-U_4.y);
	}
	else if(NUMMER(q,j)==NRBI(q,j,5))
	{
          U_5.x=REF_POS(q,j,X)-ORT(q,j,X);
          U_5.y=REF_POS(q,j,Y)-ORT(q,j,Y);
          U_5.z=REF_POS(q,j,Z)-ORT(q,j,Z);
          U_5.x=MINIMGX(U_5.x);
          U_5.y=MINIMGY(U_5.y);
          U_5.z=MINIMGZ(U_5.z);

          U_dotdot.y+=nrbk*(U_5.y+U_5.z);
          U_dotdot.z+=nrbk*(U_5.z+U_5.y);
	}
      }
#endif // NRB
*/

#ifdef TTM
        NUMNEIGHS(p,i)++;
        NUMNEIGHS(q,j)++;
#endif

#ifdef LOD	
	//erstmal nur fcc
	//Wichtig: sowohl für atom i als auch für j berechnen
	//damit zählt man zwar doppelt aber das verhindert probleme
	//in ghost-layern, da hier sonst lod größer ausfallen würde, weil
	//beide procs von "beiden" seiten ins selbe atom schreiben!
	d2.x=-d.x;
	d2.y=-d.y;
	d2.z=-d.z;
	//100
	qvec.x=qvec_fac;qvec.y=0;qvec.z=0;	
 	LODP(p,i)+=cexp(I*SPROD(qvec,d));
        LODP(q,j)+=cexp(I*SPROD(qvec,d2));
	
	//010
        qvec.x=0;qvec.y=qvec_fac;qvec.z=0;
        LODP(p,i)+=cexp(I*SPROD(qvec,d));
	LODP(q,j)+=cexp(I*SPROD(qvec,d2));

	//001
        qvec.x=0;qvec.y=0;qvec.z=1;
        LODP(p,i)+=cexp(I*SPROD(qvec,d));
	LODP(q,j)+=cexp(I*SPROD(qvec,d2));

	//110
        qvec.x=qvec_fac;qvec.y=qvec_fac;qvec.z=0;
        LODP(p,i)+=cexp(I*SPROD(qvec,d));
	LODP(q,j)+=cexp(I*SPROD(qvec,d2));

	//011
        qvec.x=0;qvec.y=qvec_fac;qvec.z=qvec_fac;
        LODP(p,i)+=cexp(I*SPROD(qvec,d));
	LODP(q,j)+=cexp(I*SPROD(qvec,d2));

	//101
        qvec.x=qvec_fac;qvec.y=0;qvec.z=qvec_fac;
        LODP(p,i)+=cexp(I*SPROD(qvec,d));
	LODP(q,j)+=cexp(I*SPROD(qvec,d2));

#endif
//ENDOF MYMOD

#ifdef LINPOT
        PAIR_INT_LIN(pot_zwi, pot_grad, pair_pot_lin, col, inc, r2, is_short)
#else
        PAIR_INT(pot_zwi, pot_grad, pair_pot, col, inc, r2, is_short)
#endif
#elif defined(KEATING)
      if (r2 < keat_r2_cut[p_typ][q_typ]) {
	PAIR_INT_KEATING(pot_zwi, pot_grad, p_typ, q_typ, r2)
#endif

        /* store force in temporary variable */
        force.x = d.x * pot_grad;
        force.y = d.y * pot_grad;
#ifndef TWOD
        force.z = d.z * pot_grad;
#endif

        /* accumulate forces */
        pfptr = &KRAFT(p,i,X);
        qfptr = &KRAFT(q,j,X);
        *pfptr     += force.x; 
        *qfptr     -= force.x; 
        *(++pfptr) += force.y; 
        *(++qfptr) -= force.y; 
#ifndef TWOD
        *(++pfptr) += force.z; 
        *(++qfptr) -= force.z; 
#endif

        *Epot      += pot_zwi;

#ifndef MONOLJ
        pot_zwi *= 0.5;   /* avoid double counting */
#ifdef NNBR
        if (r2 < nb_r2_cut[col ]) NBANZ(p,i)++;
        if (r2 < nb_r2_cut[col2]) NBANZ(q,j)++;
#endif
#ifdef ORDPAR
        if (r2 < op_r2_cut[col ]) POTENG(p,i) += op_weight[col ] * pot_zwi;
        if (r2 < op_r2_cut[col2]) POTENG(q,j) += op_weight[col2] * pot_zwi;
#else
        POTENG(p,i) += pot_zwi;
        POTENG(q,j) += pot_zwi;
#endif
#endif

#ifdef P_AXIAL
        tmp_vir_vect.x -= d.x * force.x;
        tmp_vir_vect.y -= d.y * force.y;
#ifndef TWOD
        tmp_vir_vect.z -= d.z * force.z;
#endif
#else
        tmp_virial     -= r2 * pot_grad;
#endif

#ifdef STRESS_TENS
        if (do_press_calc) {
          /* avoid double counting of the virial */
          force.x *= 0.5;
          force.y *= 0.5;
#ifndef TWOD
          force.z *= 0.5;
#endif
          PRESSTENS(p,i,xx) -= d.x * force.x;
          PRESSTENS(q,j,xx) -= d.x * force.x;
          PRESSTENS(p,i,yy) -= d.y * force.y;
          PRESSTENS(q,j,yy) -= d.y * force.y;
          PRESSTENS(p,i,xy) -= d.x * force.y;
          PRESSTENS(q,j,xy) -= d.x * force.y;
#ifndef TWOD
          PRESSTENS(p,i,zz) -= d.z * force.z;
          PRESSTENS(q,j,zz) -= d.z * force.z;
          PRESSTENS(p,i,yz) -= d.y * force.z;
          PRESSTENS(q,j,yz) -= d.y * force.z;
          PRESSTENS(p,i,zx) -= d.z * force.x;
          PRESSTENS(q,j,zx) -= d.z * force.x;
#endif
	}
#endif
      }
#endif /* PAIR || KEATING */

#ifdef EAM2
      /* compute host electron density */
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
#endif

#if defined(COVALENT)
      /* make neighbor tables for covalent systems */
      if (r2 <= neightab_r2cut[col]) {

        neightab *neigh;
        real  *tmp_ptr;

        /* update neighbor table of particle i */
        neigh = NEIGH(p,i);
        if (neigh->n_max <= neigh->n) {
          increase_neightab( neigh, neigh->n_max + NEIGH_LEN_INC );
        }
        neigh->typ[neigh->n] = q_typ;
        neigh->cl [neigh->n] = q;
        neigh->num[neigh->n] = j;
        tmp_ptr  = &neigh->dist[3*neigh->n];
        *tmp_ptr = d.x; ++tmp_ptr; 
        *tmp_ptr = d.y; ++tmp_ptr; 
        *tmp_ptr = d.z;
        neigh->n++;

        /* update neighbor table of particle j */
        neigh = NEIGH(q,j);
        if (neigh->n_max <= neigh->n) {
          increase_neightab( neigh, neigh->n_max + NEIGH_LEN_INC );
        }
        neigh->typ[neigh->n] = p_typ;
        neigh->cl [neigh->n] = p;
        neigh->num[neigh->n] = i;
        tmp_ptr  = &neigh->dist[3*neigh->n];
        *tmp_ptr = -d.x; ++tmp_ptr; 
        *tmp_ptr = -d.y; ++tmp_ptr; 
        *tmp_ptr = -d.z;
        neigh->n++;
      }
#endif  /* COVALENT */

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
#ifndef TWOD
  *Vir_zz += tmp_vir_vect.z;
  *Virial += tmp_vir_vect.z;
#endif
#else
  *Virial += tmp_virial;
#endif 

}

#endif //SIMD
