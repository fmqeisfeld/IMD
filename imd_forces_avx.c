#include "imd.h"
#include "potaccess.h"
#include <immintrin.h>

//MAKRO dass die vecs befuellen soll aber darauf achtet dass "fake" atome ignoriert werden,
//d.h. dass z.b. ORT(p,i,X) nicht geht wenn i>maxatoms. In diesem Fall wird "0" zurueck
//gegeben
// ireal, d.h. die tatsächliche zahl von atomen muss vorher definiert worden sein
// macfun2=makro mit 2 argumenten, cptr=cell-ptr (p oder q), ii=index
//
//VECBOUND MAKROS
#define VBND2(macfun2,cptr,ii) (ii)>(nreal)-1 ? 0 : macfun2(cptr,(ii))  
#define VBND3(macfun3,cptr,ii,arg3) (ii)>(nreal)-1 ? 0: macfun3(cptr,(ii),arg3)
	

void do_forces_avx(cell *p, cell *q, vektor pbc, real *Epot, real *Virial, 
               real *Vir_xx, real *Vir_yy, real *Vir_zz,
               real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  int i,j,l;
  vektor d;
  vektor tmp_d;
  vektor force;
  //real r2, rho_h;
  real tmp_virial;
#ifdef P_AXIAL
  vektor tmp_vir_vect;
#endif
  //real pot_zwi, pot_grad;
  int col, col2, is_short=0, inc = ntypes * ntypes;
  int jstart, q_typ, p_typ;
  real  *qptr, *pfptr, *qfptr, *qpdptr, *ppdptr, *qpoptr, *ppoptr;

  tmp_virial     = 0.0;
#ifdef P_AXIAL
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
  tmp_vir_vect.z = 0.0;
#endif

  //const vecs
    __m256 onevec= _mm256_set_ps(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0);
    __m256 halfvec=_mm256_set_ps(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5);
    __m256 zerovec = _mm256_setzero_ps();
    __m256i zeroveci= _mm256_set_epi32(0,0,0,0,0,0,0,0);
    __m256 neg2 = _mm256_set_ps(-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0);
    __m256 twovec=_mm256_set_ps(2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0);
   
    float* ffx,*ffy,*ffz; //for force-accumulation

  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) {   //theoretich auch vektorisierbar

    tmp_d.x = ORT(p,i,X) - pbc.x;
    tmp_d.y = ORT(p,i,Y) - pbc.y;
    tmp_d.z = ORT(p,i,Z) - pbc.z;

    p_typ  = SORTE(p,i);
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0))               ? i+1 : 0);
    qptr   = &ORT(q,jstart,X);
    /* for each atom in neighbouring cell */
    // *****************************
    // * POS VECTORS 
    // ****************************   
    __m256 ix = _mm256_set_ps(tmp_d.x,tmp_d.x,tmp_d.x,tmp_d.x,
                               tmp_d.x,tmp_d.x,tmp_d.x,tmp_d.x);

    __m256 iy = _mm256_set_ps(tmp_d.y,tmp_d.y,tmp_d.y,tmp_d.y,
                               tmp_d.y,tmp_d.y,tmp_d.y,tmp_d.y);

    __m256 iz = _mm256_set_ps(tmp_d.z,tmp_d.z,tmp_d.z,tmp_d.z,
                               tmp_d.z,tmp_d.z,tmp_d.z,tmp_d.z);


    const int nn=(int) (ceil(((float) q->n)/8-0)*8); //nächst höherer, durch 8 teilbarer integer
    const int nreal=q->n;                            //tatsächliche zahl an atomen
    const int resid=nn-nreal;                        //die letzten "resid" atome bzw. vektor-elemente sind fake
    // ************************************
    // * INNER FORCE-LOOP VECTORIZATION
    // *************************************
    for (j=jstart;j<nn;j+=8)
    {
    
      //pos
      __m256 jx = _mm256_set_ps(VBND3(ORT,q,j,X),VBND3(ORT,q,j+1,X),VBND3(ORT,q,j+2,X),VBND3(ORT,q,j+3,X),
				VBND3(ORT,q,j+4,X),VBND3(ORT,q,j+5,X),VBND3(ORT,q,j+6,X),VBND3(ORT,q,j+7,X));

      __m256 jy = _mm256_set_ps(VBND3(ORT,q,j,Y),VBND3(ORT,q,j+1,Y),VBND3(ORT,q,j+2,Y),VBND3(ORT,q,j+3,Y),
                                VBND3(ORT,q,j+4,Y),VBND3(ORT,q,j+5,Y),VBND3(ORT,q,j+6,Y),VBND3(ORT,q,j+7,Y));

      __m256 jz = _mm256_set_ps(VBND3(ORT,q,j,Z),VBND3(ORT,q,j+1,Z),VBND3(ORT,q,j+2,Z),VBND3(ORT,q,j+3,Z),
                                VBND3(ORT,q,j+4,Z),VBND3(ORT,q,j+5,Z),VBND3(ORT,q,j+6,Z),VBND3(ORT,q,j+7,Z));

      //distance
      __m256  dx = _mm256_sub_ps(jx, ix);
      __m256  dy = _mm256_sub_ps(jy, iy);
      __m256  dz = _mm256_sub_ps(jz, iz);

      //type
//      __m256i jt = _mm256_set_epi32(SORTE(q,j),SORTE(q,j),SORTE(q,j),SORTE(q,j),
//				   SORTE(q,j),SORTE(q,j),SORTE(q,j),SORTE(q,j));
      __m256i jt = zeroveci;
                          
     // 8 dot products at a time (instead of SPROD(d,d)
     __m256 dx2 = _mm256_mul_ps(dx,dx);
     __m256 dy2 = _mm256_mul_ps(dy,dy);
     __m256 dz2 = _mm256_mul_ps(dz,dz);
     __m256 r2_tmp  = _mm256_add_ps(dx2,dy2);
     __m256 r2 = _mm256_add_ps(r2_tmp,dz2);

     //now we need poteng and potgrad, as vecs
     //so far only quadratic interpol, c.f. PAIT_INT2
     //start by building tmp-vars for interpol. coeffs
     float  ptends=pair_pot.end[0]; //bisher nur 1 atomtyp
     __m256 ptend   = _mm256_set_ps(ptends,ptends,ptends,ptends,
				   ptends,ptends,ptends,ptends);
     __m256 r2a     = _mm256_min_ps(ptend,r2);
     r2a	    = _mm256_max_ps(r2a,zerovec);
     float invstep=pair_pot.invstep[col];
     __m256 istep= _mm256_set_ps(invstep,invstep,invstep,invstep,
                               invstep,invstep,invstep,invstep);
     r2a=_mm256_mul_ps(r2a,istep); // r2a=r2a*istep

     __m256i k	    = _mm256_cvtps_epi32(r2a); //float to int cast
     __m256 kf	    = _mm256_cvtepi32_ps(k);   //kf ist float-variante von k
     __m256  chi    = _mm256_sub_ps(r2a,kf);

     int* indx = (int*)&k;
     int inc=1;  
     int col=0;
     __m256 p0	    = _mm256_set_ps(pair_pot.table[indx[0]*inc+col],pair_pot.table[indx[1]*inc+col],
				    pair_pot.table[indx[2]*inc+col],pair_pot.table[indx[3]*inc+col],
				    pair_pot.table[indx[4]*inc+col],pair_pot.table[indx[5]*inc+col],
				    pair_pot.table[indx[6]*inc+col],pair_pot.table[indx[7]*inc+col]);	

     __m256 p1      = _mm256_set_ps(pair_pot.table[indx[0]*inc+col+inc],pair_pot.table[indx[1]*inc+col+inc],
                                    pair_pot.table[indx[2]*inc+col+inc],pair_pot.table[indx[3]*inc+col+inc],
                                    pair_pot.table[indx[4]*inc+col+inc],pair_pot.table[indx[5]*inc+col+inc],
                                    pair_pot.table[indx[6]*inc+col+inc],pair_pot.table[indx[7]*inc+col+inc]);

     __m256 p2      = _mm256_set_ps(pair_pot.table[indx[0]*inc+col+2*inc],pair_pot.table[indx[1]*inc+col+2*inc],
                                    pair_pot.table[indx[2]*inc+col+2*inc],pair_pot.table[indx[3]*inc+col+2*inc],
                                    pair_pot.table[indx[4]*inc+col+2*inc],pair_pot.table[indx[5]*inc+col+2*inc],
                                    pair_pot.table[indx[6]*inc+col+2*inc],pair_pot.table[indx[7]*inc+col+2*inc]);
    __m256 dv   = _mm256_sub_ps(p1,p0);
    __m256 dv2  = _mm256_mul_ps(neg2,p1);
    dv2		= _mm256_add_ps(dv2,p0);
    dv2		= _mm256_add_ps(dv2,p2);
  
    //Now calc pot: pot=p0+chi*dv + 0.5*chi*(chi-1)*d2v
    __m256 pot  =  _mm256_fmadd_ps(chi,dv,p0); //a*b+c bzw. chi*dv+p0
    __m256 chi_mone=_mm256_sub_ps(chi,onevec); //(chi - 1)
    __m256 chi_half=_mm256_mul_ps(chi,halfvec); //chi*0.5
    __m256 tmpvec  =_mm256_mul_ps(chi_mone,dv2); // (chi-1)*dv2
    pot= _mm256_fmadd_ps(tmpvec,chi_half,pot); // ((chi-1)*dv2*0.5*chi) +(chi*dv+p0)

   //now calc grad
   __m256 chi_mhalf=_mm256_sub_ps(chi,halfvec); //chi-0.5
   tmpvec=_mm256_fmadd_ps(chi_mhalf,dv2,dv); // (chi-0.5)*dv2+dv
   __m256 grad=_mm256_mul_ps(twovec,istep); // 2.0*istep
   grad=_mm256_mul_ps(grad,tmpvec); // 2*istep*((chi-0.5)*dv2+dv)
  
   //now calc forces
   __m256 fx=_mm256_mul_ps(grad,dx);
   __m256 fy=_mm256_mul_ps(grad,dy);
   __m256 fz=_mm256_mul_ps(grad,dz);


  //now accumulate forces (scalar step)
  //maybe extract is faster??

  ffx = (float*)&fx;
  ffy = (float*)&fy;
  ffz = (float*)&fz;
  int lmax=MIN(8,nreal-i); 
  for(l=0;l<lmax;l++)
  {
    KRAFT(p,i,X)+=ffx[l];
    KRAFT(p,i,Y)+=ffy[l];
    KRAFT(p,i,Z)+=ffz[l];

    KRAFT(q,j+l,X)-=ffx[l];
    KRAFT(q,j+l,Y)-=ffy[l];
    KRAFT(q,j+l,Z)-=ffz[l];
  }
  // ****************************************************************** 
  // * Now EAM-part: ALLES exakt dasselbe AUSSER table --> rho_h_tab
  // * ****************************************************************
#ifdef EAM2
  //compute host electron density: Eigentlich nur fuer r2<rho_h_tab.end
  //Spaeter mit einem mask-vektor resetten

  ptends=rho_h_tab.end[0];
  ptend   = _mm256_set_ps(ptends,ptends,ptends,ptends,
			  ptends,ptends,ptends,ptends);  
  invstep=rho_h_tab.invstep[col];
  istep= _mm256_set_ps(invstep,invstep,invstep,invstep,
                               invstep,invstep,invstep,invstep);
  r2a     = _mm256_min_ps(ptend,r2);
  r2a     = _mm256_max_ps(r2a,zerovec); 
  r2a=_mm256_mul_ps(r2a,istep); // r2a=r2a*istep

  k      = _mm256_cvtps_epi32(r2a); //float to int cast
  kf      = _mm256_cvtepi32_ps(k);   //kf ist float-variante von k
  chi    = _mm256_sub_ps(r2a,kf);

  indx = (int*)&k;
  p0      = _mm256_set_ps(rho_h_tab.table[indx[0]*inc+col],rho_h_tab.table[indx[1]*inc+col],
                                    rho_h_tab.table[indx[2]*inc+col],rho_h_tab.table[indx[3]*inc+col],
                                    rho_h_tab.table[indx[4]*inc+col],rho_h_tab.table[indx[5]*inc+col],
                                    rho_h_tab.table[indx[6]*inc+col],rho_h_tab.table[indx[7]*inc+col]);

  p1      = _mm256_set_ps(rho_h_tab.table[indx[0]*inc+col+inc],rho_h_tab.table[indx[1]*inc+col+inc],
                                    rho_h_tab.table[indx[2]*inc+col+inc],rho_h_tab.table[indx[3]*inc+col+inc],
                                    rho_h_tab.table[indx[4]*inc+col+inc],rho_h_tab.table[indx[5]*inc+col+inc],
                                    rho_h_tab.table[indx[6]*inc+col+inc],rho_h_tab.table[indx[7]*inc+col+inc]);

  p2      = _mm256_set_ps(rho_h_tab.table[indx[0]*inc+col+2*inc],rho_h_tab.table[indx[1]*inc+col+2*inc],
                                    rho_h_tab.table[indx[2]*inc+col+2*inc],rho_h_tab.table[indx[3]*inc+col+2*inc],
                                    rho_h_tab.table[indx[4]*inc+col+2*inc],rho_h_tab.table[indx[5]*inc+col+2*inc],
                                    rho_h_tab.table[indx[6]*inc+col+2*inc],rho_h_tab.table[indx[7]*inc+col+2*inc]);
  dv   = _mm256_sub_ps(p1,p0);
  dv2  = _mm256_mul_ps(neg2,p1);
  dv2  = _mm256_add_ps(dv2,p0);
  dv2  = _mm256_add_ps(dv2,p2);

  //now calc val. (bzw. ich nenns weiterhin pot)
  pot  =  _mm256_fmadd_ps(chi,dv,p0); //a*b+c bzw. chi*dv+p0
  chi_mone=_mm256_sub_ps(chi,onevec); //(chi - 1)
  chi_half=_mm256_mul_ps(chi,halfvec); //chi*0.5
  tmpvec  =_mm256_mul_ps(chi_mone,dv2); // (chi-1)*dv2
  pot= _mm256_fmadd_ps(tmpvec,chi_half,pot); // ((chi-1)*dv2*0.5*chi) +(chi*dv+p0)
 
  //if(p_typ==t_typ){
  //Estmal nur fuer p_typ==q_typ
  float* eamrho = (float*)&pot;
  for(l=0;l<lmax;l++)
  {
    EAM_RHO(q,j+l)+=eamrho[l];
    EAM_RHO(p,i)+=eamrho[l];
  }
  //}
  //
  //Jetzt muessten die Atome eigentlich not mit poteng versehen werden, aber fuer die
  //MD spielt das keine rolle, nur fuer auswertung, visualisierung, etc.

#endif

    } /* for j */
  } /* for i */
}

// ********************************************************************
// * 2nd FORCE-LOOP FOR EAM-POT IN VECTORIZED AVX2 FORMAT USING INTRINSICS
// * CALCULATES FORCE AND ENERGY CAUSED BY THE EMBEDDING ELECTRON DENSITY
// * USES PHI(r2), Rho(r2), F(rho) AND ITS DERIVATIVES
// *******************************************************************
void do_forces_eam2_avx(cell *p, cell *q, vektor pbc, real *Virial,
                    real *Vir_xx, real *Vir_yy, real *Vir_zz,
                    real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  int i,j,l;
  vektor d;
  vektor tmp_d;
  vektor force;
  real tmp_virial;
  //real pot_zwi, pot_grad;
  int col, col2, is_short=0, inc = ntypes * ntypes;
  int jstart, q_typ, p_typ;
  real  *qptr, *pfptr, *qfptr, *qpdptr, *ppdptr, *qpoptr, *ppoptr;

  tmp_virial     = 0.0;

  //const vecs
  __m256 zerovec = _mm256_setzero_ps();
  __m256 neg2 = _mm256_set_ps(-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0);
  __m256 halfvec=_mm256_set_ps(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5);
  __m256 twovec=_mm256_set_ps(2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0);
  __m256 onevec= _mm256_set_ps(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0);

  float* ffx,*ffy,*ffz;

  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) 
  {
    tmp_d.x = ORT(p,i,X) - pbc.x;
    tmp_d.y = ORT(p,i,Y) - pbc.y;
    tmp_d.z = ORT(p,i,Z) - pbc.z;

    p_typ  = SORTE(p,i);
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0))               ? i+1 : 0);
    qptr   = &ORT(q,jstart,X);
    /* for each atom in neighbouring cell */
    // *****************************
    // * POS VECTORS 
    // ****************************   
    __m256 ix = _mm256_set_ps(tmp_d.x,tmp_d.x,tmp_d.x,tmp_d.x,
                               tmp_d.x,tmp_d.x,tmp_d.x,tmp_d.x);

    __m256 iy = _mm256_set_ps(tmp_d.y,tmp_d.y,tmp_d.y,tmp_d.y,
                               tmp_d.y,tmp_d.y,tmp_d.y,tmp_d.y);

    __m256 iz = _mm256_set_ps(tmp_d.z,tmp_d.z,tmp_d.z,tmp_d.z,
                               tmp_d.z,tmp_d.z,tmp_d.z,tmp_d.z);
   
    real rho_i_strich; 
    const int nn=(int) (q->n/8.0);    
    //for (j = jstart; j < nn; ++j) 
    // ************************************
    // * INNER FORCE-LOOP VECTORIZATION
    // *************************************
    for (j=jstart;j<nn;j+=8)
    {
      //pos
      __m256 jx = _mm256_set_ps(ORT(q,j,X),ORT(q,j+1,X),ORT(q,j+2,X),ORT(q,j+3,X),
				ORT(q,j+4,X),ORT(q,j+5,X),ORT(q,j+6,X),ORT(q,j+7,X));

      __m256 jy = _mm256_set_ps(ORT(q,j,Y),ORT(q,j+1,Y),ORT(q,j+2,Y),ORT(q,j+3,Y),
                                ORT(q,j+4,Y),ORT(q,j+5,Y),ORT(q,j+6,Y),ORT(q,j+7,Y));

      __m256 jz = _mm256_set_ps(ORT(q,j,Z),ORT(q,j+1,Z),ORT(q,j+2,Z),ORT(q,j+3,Z),
                                ORT(q,j+4,Z),ORT(q,j+5,Z),ORT(q,j+6,Z),ORT(q,j+7,Z));

      //distance
      __m256  dx = _mm256_sub_ps(jx, ix);
      __m256  dy = _mm256_sub_ps(jy, iy);
      __m256  dz = _mm256_sub_ps(jz, iz);

      //type
      __m256i jt = _mm256_set_epi32(SORTE(q,j),SORTE(q,j),SORTE(q,j),SORTE(q,j),
				   SORTE(q,j),SORTE(q,j),SORTE(q,j),SORTE(q,j));

     // 8 dot products at a time (instead of SPROD(d,d)
     __m256 dx2 = _mm256_mul_ps(dx,dx);
     __m256 dy2 = _mm256_mul_ps(dy,dy);
     __m256 dz2 = _mm256_mul_ps(dz,dz);
     __m256 r2_tmp  = _mm256_add_ps(dx2,dy2);
     __m256 r2 = _mm256_add_ps(r2_tmp,dz2);

     // rho_strich_i(r_ij). c.f. DERIV_FUNC2(...)
     // -->Genau wie grad in PAIR_INT2, aber mit rho_h_tab
     float  ptends=rho_h_tab.end[0]; //bisher nur 1 atomtyp
     __m256 ptend   = _mm256_set_ps(ptends,ptends,ptends,ptends,
                                   ptends,ptends,ptends,ptends);
     __m256 r2a     = _mm256_min_ps(ptend,r2);
     r2a            = _mm256_max_ps(r2a,zerovec);
     float invstep=rho_h_tab.invstep[col];
     __m256 istep= _mm256_set_ps(invstep,invstep,invstep,invstep,
                               invstep,invstep,invstep,invstep);
     r2a=_mm256_mul_ps(r2a,istep); // r2a=r2a*istep

     __m256i k      = _mm256_cvtps_epi32(r2a); //float to int cast
     __m256 kf      = _mm256_cvtepi32_ps(k);   //kf ist float-variante von k
     __m256  chi    = _mm256_sub_ps(r2a,kf);

     int* indx = (int*)&k;
     int inc=1;
     int col=0;
     __m256 p0      = _mm256_set_ps(rho_h_tab.table[indx[0]*inc+col],rho_h_tab.table[indx[1]*inc+col],
                                    rho_h_tab.table[indx[2]*inc+col],rho_h_tab.table[indx[3]*inc+col],
                                    rho_h_tab.table[indx[4]*inc+col],rho_h_tab.table[indx[5]*inc+col],
                                    rho_h_tab.table[indx[6]*inc+col],rho_h_tab.table[indx[7]*inc+col]);

     __m256 p1      = _mm256_set_ps(rho_h_tab.table[indx[0]*inc+col+inc],rho_h_tab.table[indx[1]*inc+col+inc],
                                    rho_h_tab.table[indx[2]*inc+col+inc],rho_h_tab.table[indx[3]*inc+col+inc],
                                    rho_h_tab.table[indx[4]*inc+col+inc],rho_h_tab.table[indx[5]*inc+col+inc],
                                    rho_h_tab.table[indx[6]*inc+col+inc],rho_h_tab.table[indx[7]*inc+col+inc]);

     __m256 p2      = _mm256_set_ps(rho_h_tab.table[indx[0]*inc+col+2*inc],rho_h_tab.table[indx[1]*inc+col+2*inc],
                                    rho_h_tab.table[indx[2]*inc+col+2*inc],rho_h_tab.table[indx[3]*inc+col+2*inc],
                                    rho_h_tab.table[indx[4]*inc+col+2*inc],rho_h_tab.table[indx[5]*inc+col+2*inc],
				    rho_h_tab.table[indx[6]*inc+col+2*inc],rho_h_tab.table[indx[7]*inc+col+2*inc]);


    __m256 dv   = _mm256_sub_ps(p1,p0);
    __m256 dv2  = _mm256_mul_ps(neg2,p1);
    dv2         = _mm256_add_ps(dv2,p0);
    dv2         = _mm256_add_ps(dv2,p2);

   //calc. grad
   __m256 chi_mhalf=_mm256_sub_ps(chi,halfvec); //chi-0.5
   __m256 tmpvec=_mm256_fmadd_ps(chi_mhalf,dv2,dv); // (chi-0.5)*dv2+dv
   //__m256 grad=_mm256_mul_ps(twovec,istep); // 2.0*istep
   __m256 grad=_mm256_mul_ps(istep,tmpvec); // 1*istep*((chi-0.5)*dv2+dv) // ACHTUNG: statt 2*istep ... weil nachher eh halbiert wird!

   //grad=rho_i_strich==rho_j_strich weil type identisch ist
   
   //now comp. eam2_force = 0.5 * (EAM_DF(p,i)*rho_j_strich+EAM_DF(q,j)*rho_i_strich); 
   //EAM_DF(p,i) is computed in do_embedding_energy
   // 
   // ich berechne eam2force=(EAM_DF(p,i)+EAM_DF(q,j))*rho_i_strich
   float dfi_f=EAM_DF(p,i);   
   __m256 dfi=_mm256_set_ps(dfi_f,dfi_f,dfi_f,dfi_f,
			   dfi_f,dfi_f,dfi_f,dfi_f);

   __m256 dfj=_mm256_set_ps(EAM_DF(q,j),EAM_DF(q,j+1),EAM_DF(q,j+2),EAM_DF(p,j+3),
			  EAM_DF(q,j+4),EAM_DF(q,j+5),EAM_DF(q,j+6),EAM_DF(p,j+7));
   
   __m256 eam2f=_mm256_mul_ps(_mm256_add_ps(dfi,dfj),grad);

   //now calc forces
   __m256 fx=_mm256_mul_ps(eam2f,dx);
   __m256 fy=_mm256_mul_ps(eam2f,dy);
   __m256 fz=_mm256_mul_ps(eam2f,dz);

   //accumulate forces
   ffx = (float*)&fx;
   ffy = (float*)&fy;
   ffz = (float*)&fz;

   for(l=0;l<8;l++)
   {
     KRAFT(p,i,X)+=ffx[l];
     KRAFT(p,i,Y)+=ffy[l];
     KRAFT(p,i,Z)+=ffz[l];

     KRAFT(q,j+l,X)-=ffx[l];
     KRAFT(q,j+l,Y)-=ffy[l];
     KRAFT(q,j+l,Z)-=ffz[l];
   }
   
   } // loop 2nd atom
 } //loop 1st atom
}


// ******************************************************************************
// *
// * compute embedding energy and its derivative for all atoms
// *
// *******************************************************************************/
				
void do_embedding_energy_avx(void)
{
  int k2;
  int l;
  //const vecs
  __m256 onevec= _mm256_set_ps(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0);
  __m256 halfvec=_mm256_set_ps(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5);
  __m256 zerovec = _mm256_setzero_ps();
  __m256 neg2 = _mm256_set_ps(-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0);
  __m256 twovec=_mm256_set_ps(2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0);
  
  for (k2=0; k2<NCELLS; k2++) {
    int  i, idummy=0;
    real pot;
    cell *p;
    p = CELLPTR(k2);

    //const int nn=(int) (p->n/8.0);
    const int nn=(int) (ceil(((float) p->n)/8-0)*8); //nächst höherer, durch 8 teilbarer integer
    const int nreal=p->n;			     //tatsächliche zahl an atomen
    const int resid=nn-nreal;			     //die letzten "resid" atome bzw. vektor-elemente sind fake

    for (i = 0; i < nn; i+=8)    
    // ************************************
    // * INNER FORCE-LOOP VECTORIZATION
    // *************************************
    {
      //Hier passiert dasselbe wie in force-loop 1, nur mit embed_pot und EAM_RHO(p,i) steht fuer r2
      __m256 eamrho = _mm256_set_ps(VBND2(EAM_RHO,p,i), VBND2(EAM_RHO,p,i+1), VBND2(EAM_RHO,p,i+2), VBND2(EAM_RHO,p,i+3),
				  VBND2(EAM_RHO,p,i+4), VBND2(EAM_RHO,p,i+5), VBND2(EAM_RHO,p,i+6), VBND2(EAM_RHO,p,i+7));

     float  ptends=embed_pot.end[0]; //bisher nur 1 atomtyp
     __m256 ptend   = _mm256_set_ps(ptends,ptends,ptends,ptends,
                                   ptends,ptends,ptends,ptends);

     __m256 r2a     = _mm256_min_ps(ptend,eamrho);
     r2a            = _mm256_max_ps(r2a,zerovec);

     int inc=1;  //weil nur 1 atomtyp bisher
     int col=0;

     float invstep=embed_pot.invstep[col];
     __m256 istep= _mm256_set_ps(invstep,invstep,invstep,invstep,
                               invstep,invstep,invstep,invstep);
     r2a=_mm256_mul_ps(r2a,istep); // r2a=r2a*istep

     __m256i k      = _mm256_cvtps_epi32(r2a); //float to int cast
     __m256 kf      = _mm256_cvtepi32_ps(k);   //kf ist float-variante von k
     __m256  chi    = _mm256_sub_ps(r2a,kf);

     int* indx = (int*)&k;
     __m256 p0      = _mm256_set_ps(embed_pot.table[indx[0]*inc+col],embed_pot.table[indx[1]*inc+col],
                                    embed_pot.table[indx[2]*inc+col],embed_pot.table[indx[3]*inc+col],
                                    embed_pot.table[indx[4]*inc+col],embed_pot.table[indx[5]*inc+col],
                                    embed_pot.table[indx[6]*inc+col],embed_pot.table[indx[7]*inc+col]);

     __m256 p1      = _mm256_set_ps(embed_pot.table[indx[0]*inc+col+inc],embed_pot.table[indx[1]*inc+col+inc],
                                    embed_pot.table[indx[2]*inc+col+inc],embed_pot.table[indx[3]*inc+col+inc],
                                    embed_pot.table[indx[4]*inc+col+inc],embed_pot.table[indx[5]*inc+col+inc],
                                    embed_pot.table[indx[6]*inc+col+inc],embed_pot.table[indx[7]*inc+col+inc]);

     __m256 p2      = _mm256_set_ps(embed_pot.table[indx[0]*inc+col+2*inc],embed_pot.table[indx[1]*inc+col+2*inc],
                                    embed_pot.table[indx[2]*inc+col+2*inc],embed_pot.table[indx[3]*inc+col+2*inc],
                                    embed_pot.table[indx[4]*inc+col+2*inc],embed_pot.table[indx[5]*inc+col+2*inc],
                                    embed_pot.table[indx[6]*inc+col+2*inc],embed_pot.table[indx[7]*inc+col+2*inc]);
    __m256 dv   = _mm256_sub_ps(p1,p0);
    __m256 dv2  = _mm256_mul_ps(neg2,p1);
    dv2         = _mm256_add_ps(dv2,p0);
    dv2         = _mm256_add_ps(dv2,p2);

     //calc pot
     __m256 pot  =  _mm256_fmadd_ps(chi,dv,p0); //a*b+c bzw. chi*dv+p0
     __m256 chi_mone=_mm256_sub_ps(chi,onevec); //(chi - 1)
     __m256 chi_half=_mm256_mul_ps(chi,halfvec); //chi*0.5
     __m256 tmpvec  =_mm256_mul_ps(chi_mone,dv2); // (chi-1)*dv2
     pot= _mm256_fmadd_ps(tmpvec,chi_half,pot); // ((chi-1)*dv2*0.5*chi) +(chi*dv+p0)

     //calc grad
     __m256 chi_mhalf=_mm256_sub_ps(chi,halfvec); //chi-0.5
     tmpvec=_mm256_fmadd_ps(chi_mhalf,dv2,dv); // (chi-0.5)*dv2+dv
     __m256 grad=_mm256_mul_ps(twovec,istep); // 2.0*istep
     grad=_mm256_mul_ps(grad,tmpvec); // 2*istep*((chi-0.5)*dv2+dv)

     // ORIGINAL 
     /*
      PAIR_INT( pot, EAM_DF(p,i), embed_pot, SORTE(p,i),
                ntypes, EAM_RHO(p,i), idummy);
      POTENG(p,i)    += pot;
      tot_pot_energy += pot;
     */
     float* potf =(float*)&pot;
     float* gradf=(float*)&grad;

     int lmax=MIN(8,nreal-i); //"virtuelle" atome auslassen
     for(l=0;l<lmax;l++)
     {
       POTENG(p,i+l) += potf[l];
       EAM_DF(p,i+l)  = gradf[l];
       tot_pot_energy+=potf[l];
     }
    } // loop over i
  } // loop k
}

