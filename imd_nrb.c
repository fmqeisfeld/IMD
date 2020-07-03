#include "imd.h"
// ***********************************************
// *  NON-REFLECTING BOUNDARY CONDITIONS         *
// *  Nach Comput. Mech. 50:645-655 (2012)       *
// *  Bisher nur für fcc (100) und inplane-atome *
// ***********************************************
// DEV-LOG:
//*****************************************************************************************************************************************************************************
// 
// 1.1.20: imd_nrb geht nun auch mit LB, aber nicht im init step wegen inverse_send_cells
//         Das macht aber nichts, denn man kann eifach 1-step simulation anlaufen lassen,
//         die dann das nrb-outfile schreibt.
//         anschließend kann die simulation mit LB fortgesetzt werden
//
// Wichtig: Es muss darauf geachtet werden, dass die Probe in z-Richtung nicht zu dünn ist.
// 	    Es muss verhindert werden dass es in Z-Richtung nur 1 MD-Zelle gibt.
// 	    Trifft dies zu, dann hat diese Zelle 2 Nachbar-Bufferzellen (in z-Richtung) die genau dieselben Atome enthalten.
// 	    Dadurch kann es passieren, dass ein bnd-neigh-paar doppelt auftritt und der Kraftbeitrag doppelt addiert wird.
// 	    Z.B.: das bnd-atom und sein nachbar-atom befinden sich am Rand von der MD-Zelle mit local.dim.z=1.
// 	    	  in den Nachbarzellem local.dim.z=2 und local.dim.z.=0 befinden sich kopien dieser atome.
// 	    	  Hat nun das bnd-atom eine solche z-Koordinate, dass seine kopie in local.dim.z=2 auch in der Nachbarliste vom
// 	    	  neigh-atom in local.dim.z=1 ist, dann kriegt das bnd-atom und seine kopie in der buffer-zelle den kraftbeitrag ab.
// 	    	  Sendforces sorgt dann dafür, dass das originale bnd-atom die doppelte kraft bekommt.
// 	    Kurzgesagt: Dasselbe bnd-atom darf nicht doppelt in der nachbarliste seines neigh-atoms auftreten
// 	    		Genauso darf dasselbe neigh-atom nicht doppelt in der Nachbarliste seines bnd-atoms auftreten
//
// 	    
// Ausserdem: NRB wurde bisher nur für den NVE-integrator implementiert (siehe imd_integrate.c)
//
// Update:  in nrb_init werden die Randzonen nun automatisch bestimmt. Die Richtigkeit kann allerdings nur für das perfekte Gitter garantiert werden.
//	    Deswegen empfiehlt es sich das ganze bereits vor equilibrierung kurz anlaufen zu lassen, damit ein nrb-restart file geschrieben werden kann, 
//	    welches nach der equilibrierung (ohne nrb) eingelesen wird. 
//	    Dabei sollten dann die GGew. Positionen (von vor der equilibierung) mit den aktuellen (von nach der equilib.) ersetzt werden.
//	    Natürlich entsprechen diese Positionen nicht wirklich den GGew. Positionen. (dafür müsste man z.B. ein Ensemblemittel oder Zeitmittel berechnen)
//	    Aber für den Anfang soll das erstmal genügen
//
//
// UPDATE: Ausschreibe- und Einleseroutinen ergänzt, damit Simulation fortgesetzt werden kann. Details hierzu weiter unten bei den Routinen
//	   Es kann entweder "automatisch" aus dem restart-file (selbes invervall wie chkpt-file) eingelesen werden oder aber ein anderes input-file "nrb_infile" erzwungen werden
//	   Ferner gibt es die Option "nrb_overwrite", mit der die REF_POS aus dem eingelesenen file mit den aktuellen überschrieben werden
//
// TODO: Automatisch nrb-layers identifizieren 
// UPDATE:erledigt
//	   z.B. xhi-region geht von MAX(atom.x)-alat/4 bis MAX(atom.x).  
// 	   der faktor (1/4) ist dabei willkürlich gewählt. es würde auch 1/3 fuktionieren. Wichtig ist nur, dass nicht mehr als 1 layer in der bnd-region landet
//
// TIPP: nbl_margin und evtl. nbl_size vergrößern, damit selterner rebuild
//
// UPDATE: nrb_forces komplett neu konzipiert: es wird nun keine nrb_ifromid-liste mehr benötigt.
//	   Zusammengehörige bnd/neigh-paare werden nun während der kraftschleife identifiziert.
//	   Außerdem wurde das Kommunikationsschema direkt in das Standardschema (send_cells/send_forces) in imd_comm_force_3d.c implementiert
//	   Problematisch sind weiterhin Eck- bzw. Randatome, da deren Nachbarn teilweise selber bnd-atome sind.
//         Bisher wird das so geregelt, dass deren Kraftbeitrag dann ignoriert wird ---> nicht optimal. Hier bilden sich kleinere "Stoßwellen"-Quellen
//	   VORTEIL:  Es funktioniert endlich
//	   NACHTEIL: Wahrscheinlich nicht sonderlich performant + Kommunikationsmonster (alle NRB-infos werden JEDEN step kommuniziert. send_cells vor der kraftschleife und send_forces danach)
//********************************************************************************************************************************************************************************
#define DEBUG_LEVEL 0
#define WATCHME 72259
//#define WATCH ( (NUMMER(bndcell,bndi)==38400) || (NUMMER(bndcell,bndi)==38084) )
//#define WATCH ( (NUMMER(bndcell,bndi)==9764) )
#define PRINTSTEP 10

//13189
//#define nrb_xhi  140  // alles rechts davon sind +x-bnd-atome
/*
#define nrb_xhi  261  // alles rechts davon sind +x-bnd-atome
#define nrb_ylo  22   // alle atome mit y<y_lo sind -y-bnd atome
#define nrb_yhi  261
*/

int approx(double x,double x0,double deps)
{
  if( ABS(x-x0)<=deps) return 1;
  return 0;
}

//Comparator func for qsort for 2d arr 
int compare( const void* pa, const void* pb)
{
    const int *a = *(const int **)pa;
    const int *b = *(const int **)pb;
    if(a[0] == b[0])
        return a[1] - b[1];
    else
        return a[0] - b[0];
}


//recursive binary search (wird seit 02.08.19 nicht mehr benötigt)
int nrb_binarySearch( int lo, int hi, int x,int n,int** arr) 
{ 
//first occurence search
    if(hi >= lo) 
    { 
        int mid = lo + (hi - lo)/2; 
        //if( ( mid == 0 || x > nrb_ifromid[mid-1][0]) && nrb_ifromid[mid][0] == x) 
        //if( ( mid == 0 || x > nrb_ifromid2[mid-1][0]) && nrb_ifromid2[mid][0] == x) 
        if( ( mid == 0 || x > arr[mid-1][0]) && arr[mid][0] == x) 
            return mid; 
        //else if(x > nrb_ifromid2[mid][0]) 
        else if(x > arr[mid][0]) 
            return nrb_binarySearch((mid + 1), hi, x, n,arr); 
        else
            return nrb_binarySearch(lo, (mid -1), x, n,arr); 
    } 

/*
    if(hi >= lo) 
    { 
        int mid = lo + (hi - lo)/2; 
        if( ( mid == 0 || x > nrb_ifromid[mid-1][0]) && nrb_ifromid[mid][0] == x) 
            return mid; 
        else if(x > nrb_ifromid[mid][0]) 
            return first(arr, (mid + 1), high, x, n); 
        else
            return first(arr, low, (mid -1), x, n); 
    } 
*/

    // We reach here when element is not 
    // present in array 
    return -1; 
} 
// *********************************************************************************************************************+
int init_nrb() //nrb_eps=Abstands-toleranz in Angstrom,alat=Lattice const.
{
#ifndef NBL
#ifdef PAIR
    #ifndef AR
    error("If nrb is used without NBL it must be used with AR")
    #endif
#endif
#endif
  int i,j,k,n,m;
  int same_cell,jstart;
  int cmax=0;
  int c=0;
  invbox_x=1.0/box_x.x;
  invbox_y=1.0/box_y.y;
  invbox_z=1.0/box_z.z;

  double xmaxlocal,xmaxglobal,ymaxlocal,ymaxglobal,yminlocal,yminglobal;
  xmaxlocal=-9e9;
  ymaxlocal=-9e9;
  yminlocal=9e9;
  
  int  nrb_neighs; //zum checken ob alle neighs gefunden
  real *qptr;
  real nd=nrb_alat/2.0; //component-wise neighbor-distance
  
//  nrbk=0.3; //0.82905; // u/imdtime^2 ? u^2/imdtime^2
//  nrbk=0.9;
//  nrbk=1.2;
  //der coeff fuer die matrizen ist sqrt(k/m)
  nrbk=sqrt(nrbk/26.9815);

  if(imdrestart>0)
  {
    int readstep=imdrestart;
    sprintf(nrb_restart_file,"%s.%05d.nrb", outfilename,readstep); 
  }

  double r2;
  if(myid==0)
  {
    printf("*************************************************\n");
    printf("* NON-REFLECTING  BOUNDARY  CONDITIONS\n");
#ifdef NBL
    printf("* USING NBL\n");    
#else    
    printf("* USING PAIR\n");    
#endif    

    printf("* nrb_alat:%.4e\n",nrb_alat);
    printf("* nrb_eps: %.4e\n",nrb_eps);
    printf("* nrb_k: %.4e\n",pow(nrbk,2.0)*26.9815);
    if(nrb_readfile==1)
    {
      printf("* reading non-restart input-file: %s\n",nrb_input_file);
      //printf("* nrb_overwrite:%d\n",nrb_overwrite);
    }
    if(imdrestart>0 && nrb_readfile==0)
    {
      printf("* reading restart-file:%s\n",nrb_restart_file);
      printf("* nrb_overwrite:%d\n",nrb_overwrite); //REFPOS mit aktuellen überschreiben (z.b. nach equilib.)?    
      //printf("*************************************************\n");
    }
  }


  int r;
  if( (imdrestart==0 && nrb_readfile==0))
  {
#ifndef NBL
    error("Initial NRB setup currently only works with NBL.");

#endif    

    vektor d;
    //SET REFPOS
    // Hier auch direkt xmax,ymax,ymin suchen
    for(k=0;k<nallcells;k++)
    {
      cell*p =cell_array+k;
      for(i=0;i<p->n;i++)
      {	
	REF_POS(p,i,X)=ORT(p,i,X);
	REF_POS(p,i,Y)=ORT(p,i,Y);
	REF_POS(p,i,Z)=ORT(p,i,Z);

	if(ORT(p,i,X) > xmaxlocal) xmaxlocal=ORT(p,i,X);
	if(ORT(p,i,Y) > ymaxlocal) ymaxlocal=ORT(p,i,Y); //ACHTUNG: das ist natuerlich im 1D-Fall nicht noetig
	if(ORT(p,i,Y) < yminlocal) yminlocal=ORT(p,i,Y);
      }
    } 
    MPI_Allreduce(&xmaxlocal,  &xmaxglobal,  1, REAL, MPI_MAX, cpugrid);
    MPI_Allreduce(&ymaxlocal,  &ymaxglobal,  1, REAL, MPI_MAX, cpugrid);
    MPI_Allreduce(&yminlocal,  &yminglobal,  1, REAL, MPI_MIN, cpugrid);
   
    nrb_xhi=xmaxglobal-nrb_alat/4; //-nrb_alat/4.0; // -nrb_alat/2.0 würde evtl. noch 2-te Lage mit einschließen
    nrb_yhi=ymaxglobal-nrb_alat/4; //-nrb_alat/4.0; 
    nrb_ylo=yminglobal+nrb_alat/4; //+nrb_alat/4.0;

//HOTFIX fuer 1D    
if(pbc_dirs.y==1)
{
  nrb_yhi=10e9;
  nrb_ylo=-10e9;
}
  
    if(myid==0) 
    {
     printf("* nrb_xhi:%.2f, nrb_yhi:%.2f, nrb_ylo:%.2f\n", nrb_xhi,nrb_yhi,nrb_ylo);
    }

#ifdef NBL 

    //geschieht nur 1 mal zu Beginn-->muss nicht unbedingt vektorisiert werden
    for (k=0; k<ncells2;k++) 
    {
      //cell *p = CELLPTR(k); //cell_array + cnbrs[k].np;
      int c1=cnbrs[k].np;
      cell *p = cell_array+c1;
      for (i=0; i<p->n; i++) 
      {
	vektor d1;
	int    m, it, nb = 0;
	// positon of 1st atom
	d1.x = ORT(p,i,X);
	d1.y = ORT(p,i,Y);
	d1.z = ORT(p,i,Z);

	//TO DO: BND-ATOME ÜBER "REGION" (wie in lammps) festlegen
	//PROBEWEISE ERST EINMAL HARDCODING
	if(d1.x >nrb_xhi)
	{
	  NRBBND(p,i)=1;  //x bnd-hat vorrang vor y-bnd
	}     

	else if(d1.y < nrb_ylo)
	{
	  NRBBND(p,i)=2;
	}
	else if(d1.y>nrb_yhi)
	{
	  NRBBND(p,i)=3;
	}
	else
	  NRBBND(p,i)=0;

  #if DEBUG_LEVEL>4
    if(NRBBND(p,i)>0) printf("myid:%d, index:%d, p,i, IDENTIFIED AS BND\n",myid,NUMMER(p,i));
  #endif

	//REFPOS needed to calculate displacement from "equilib." positions
  /*
	REF_POS(p,i,X)=ORT(p,i,X); 
	REF_POS(p,i,Y)=ORT(p,i,Y);
	REF_POS(p,i,Z)=ORT(p,i,Z);
  */
	// **********************
	// loop over neighbors //
	// **********************
	//for (m=tl[n]; m<tl[n+1]; m++) 
	for(m=0;m<NNBCELL;m++) //die "Nachbarzelle" kann auch dieselbe Zelle sein!
	{
	  int c2=cnbrs[k].nq[m];
//	  if(c2<0) {printf("myid:%d, c2<0..skip\n",myid);continue;}
	  if(c2<0) continue;
	  cell* q=cell_array+c2;
	  int jstart;
	  if (c2==c1) jstart = i+1;
	  else        jstart = 0;

	for(j=jstart;j<q->n;j++)             
	{
	  vektor d2,d;

	  // position of 2nd atom
	  d2.x = ORT(q,j,X);
	  d2.y = ORT(q,j,Y);
	  d2.z = ORT(q,j,Z);
	  
	  if(d2.x>nrb_xhi)
	    NRBBND(q,j)=1;
	  else if(d2.y<nrb_ylo)	
	    NRBBND(q,j)=2;
	  else if(d2.y>nrb_yhi)
	    NRBBND(q,j)=3;
	  else
	    NRBBND(q,j)=0;
  #if DEBUG_LEVEL>4
    if(NRBBND(q,j)>0) printf("myid:%d, index:%d, q,j, IDENTIFIED AS BND\n",myid,NUMMER(q,j));
  #endif

	  d.x=d2.x-d1.x;
	  d.y=d2.y-d1.y;
	  d.z=d2.z-d1.z;

	  r2    = SPROD(d,d); //dist

	  //REFPOS needed to calculate displacement from "equilib." positions
  /*
	  REF_POS(q,j,X)=ORT(q,j,X);    
	  REF_POS(q,j,Y)=ORT(q,j,Y);
	  REF_POS(q,j,Z)=ORT(q,j,Z);
  */

	  /////////////
	  //XY-PLANE // 
	  /////////////
  //SITE 0
  if(NRBBND(p,i)==1 || NRBBND(p,i)==2)
  {
	  if(approx(d.x,-nd,nrb_eps)==1 && approx(d.y,nd,nrb_eps)==1 && approx(d.z,0,nrb_eps)==1) //pos 0
	  {
	    NRBI(p,i,0)=NUMMER(q,j);
	    NRBNEIGH(q,j)=1;
	    nrb_neighs++;

	    //if(NRBBND(q,j)==0) //nur übernehmen, wenn das atom KEIN bnd-atom ist (bnd-atome haben vorrang)
	    //NRBI(q,j,0)=NUMMER(p,i); //neigh-atom kennt sein bnd-atom
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 0 for nri:%d,ci:%d,p:%p -> nrj:%d,cj:%d,q:%p,ptype:%d,qtype:%d\n", myid,NUMMER(p,i),c1,p,NUMMER(q,j),c2,q,p->celltype,q->celltype);
  #endif
	  }
  }

  //SITE 1
  if(NRBBND(p,i)==1 || NRBBND(p,i)==3)
  {
	  if(approx(d.x,-nd,nrb_eps)==1 && approx(d.y,-nd,nrb_eps)==1 && approx(d.z,0,nrb_eps)==1) //pos 1
	  {
	    NRBI(p,i,1)=NUMMER(q,j);
	    NRBNEIGH(q,j)=1;
	    nrb_neighs++;

	    //if(NRBBND(q,j)==0)
	    //NRBI(q,j,1)=NUMMER(p,i);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 1 for nri:%d,ci:%d,p:%p -> nrj:%d,cj:%d,q:%p,ptype:%d,qtype:%d\n", myid,NUMMER(p,i),c1,p,NUMMER(q,j),c2,q,p->celltype,q->celltype);
  #endif

	  }
  }
  //SITE 2
  if(NRBBND(p,i)==3)
  {
	  if(approx(d.x,nd,nrb_eps)==1 && approx(d.y,-nd,nrb_eps)==1 && approx(d.z,0,nrb_eps)==1) //pos 2
	  {
	    NRBI(p,i,2)=NUMMER(q,j);
	    NRBNEIGH(q,j)=1;
	    nrb_neighs++;
		
	    //if(NRBBND(q,j)==0)
	    //NRBI(q,j,2)=NUMMER(p,i);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 2 for nri:%d,ci:%d,p:%p -> nrj:%d,cj:%d,q:%p,ptype:%d,qtype:%d\n", myid,NUMMER(p,i),c1,p,NUMMER(q,j),c2,q,p->celltype,q->celltype);
  #endif
	  }
  }

  //SITE 3
  if(NRBBND(p,i)==2)
  {
	  if(approx(d.x,nd,nrb_eps)==1 && approx(d.y,nd,nrb_eps)==1 && approx(d.z,0,nrb_eps)==1) //pos 3
	  {
	    NRBI(p,i,3)=NUMMER(q,j);
	    NRBNEIGH(q,j)=1;
	    nrb_neighs++;

	    //if(NRBBND(q,j)==0)
	    //NRBI(q,j,3)=NUMMER(p,i);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 3 for nri:%d,ci:%d,p:%p -> nrj:%d,cj:%d,q:%p,ptype:%d,qtype:%d\n", myid,NUMMER(p,i),c1,p,NUMMER(q,j),c2,q,p->celltype,q->celltype);
  #endif
	  }
  }
	  ///////////////
	  // XZ-PLANE  //
	  ///////////////
  //SITE 9
  if(NRBBND(p,i)==1)
  {
	  if(approx(d.x,-nd,nrb_eps)==1 && approx(d.y,0,nrb_eps)==1 && approx(-d.z,nd,nrb_eps)==1) //pos 9
	  {
	    NRBI(p,i,9)=NUMMER(q,j);          
	    NRBNEIGH(q,j)=1;
	    nrb_neighs++;

	    //if(NRBBND(q,j)==0)
	    //NRBI(q,j,9)=NUMMER(p,i);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 9 for nri:%d,ci:%d,p:%p -> nrj:%d,cj:%d,q:%p,ptype:%d,qtype:%d\n", myid,NUMMER(p,i),c1,p,NUMMER(q,j),c2,q,p->celltype,q->celltype);
  #endif
	  }
  }

  //SITE 8
  if(NRBBND(p,i)==1)
  {
	  if(approx(d.x,-nd,nrb_eps)==1 && approx(d.y,0,nrb_eps)==1 && approx(d.z,nd,nrb_eps)==1) //pos 8
	  {
	    NRBI(p,i,8)=NUMMER(q,j);
	    NRBNEIGH(q,j)=1;
	    nrb_neighs++;

	    //if(NRBBND(q,j)==0)
	    //NRBI(q,j,8)=NUMMER(p,i);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 8 for nri:%d,ci:%d,p:%p -> nrj:%d,cj:%d,q:%p,ptype:%d,qtype:%d\n", myid,NUMMER(p,i),c1,p,NUMMER(q,j),c2,q,p->celltype,q->celltype);
  #endif
	  }
  }
	  /////////////
	  //YZ-PLANE // 
	  /////////////
  //SITE 4
  if(NRBBND(p,i)==3)
  {
	  if(approx(d.x,0,nrb_eps)==1 && approx(d.y,-nd,nrb_eps)==1 && approx(d.z,nd,nrb_eps)==1) //pos 4
	  {
	    NRBI(p,i,4)=NUMMER(q,j);
	    NRBNEIGH(q,j)=1;
	    nrb_neighs++;

	    //if(NRBBND(q,j)==0)
	    //NRBI(q,j,4)=NUMMER(p,i);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 4 for nri:%d,ci:%d,p:%p -> nrj:%d,cj:%d,q:%p,ptype:%d,qtype:%d\n", myid,NUMMER(p,i),c1,p,NUMMER(q,j),c2,q,p->celltype,q->celltype);
  #endif
	  }
  }
  //SITE 5
  if(NRBBND(p,i)==3)
  {
	  if(approx(d.x,0,nrb_eps)==1 && approx(d.y,-nd,nrb_eps)==1 && approx(d.z,-nd,nrb_eps)==1) //pos 5
	  {
	    NRBI(p,i,5)=NUMMER(q,j);
	    NRBNEIGH(q,j)=1;
	    nrb_neighs++;

	    //if(NRBBND(q,j)==0)
	    //NRBI(q,j,5)=NUMMER(p,i);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 5 for nri:%d,ci:%d,p:%p -> nrj:%d,cj:%d,q:%p,ptype:%d,qtype:%d\n", myid,NUMMER(p,i),c1,p,NUMMER(q,j),c2,q,p->celltype,q->celltype);
  #endif
	  }
  }
  //SITE 7
  if(NRBBND(p,i)==2)
  {
	  if(approx(d.x,0,nrb_eps)==1 && approx(d.y,nd,nrb_eps)==1 && approx(d.z,nd,nrb_eps)==1) //pos 7
	  {
	    NRBI(p,i,7)=NUMMER(q,j);
	    NRBNEIGH(q,j)=1;
	    nrb_neighs++;

	    //if(NRBBND(q,j)==0)
	    //NRBI(q,j,7)=NUMMER(p,i);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 7 for nri:%d,ci:%d,p:%p -> nrj:%d,cj:%d,q:%p,ptype:%d,qtype:%d\n", myid,NUMMER(p,i),c1,p,NUMMER(q,j),c2,q,p->celltype,q->celltype);
  #endif
	  }
  }
  //SITE 6
  if(NRBBND(p,i)==2)
  {
	  if(approx(d.x,0,nrb_eps)==1 && approx(d.y,nd,nrb_eps)==1 && approx(d.z,-nd,nrb_eps)==1) //pos 6
	  {
	    NRBI(p,i,6)=NUMMER(q,j);
	    NRBNEIGH(q,j)=1;
	    nrb_neighs++;

	    //if(NRBBND(q,j)==0)
	    //NRBI(q,j,6)=NUMMER(p,i);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 6 for nri:%d,ci:%d,p:%p -> nrj:%d,cj:%d,q:%p,ptype:%d,qtype:%d\n", myid,NUMMER(p,i),c1,p,NUMMER(q,j),c2,q,p->celltype,q->celltype);
  #endif
	  }
  }
  // ****************************
  // * SAME FOR q,j
  // ****************************

     /////////////////
     //   XY-PLANE  //
     /////////////////
  //SITE 0 q,j
  if(NRBBND(q,j)==1 || NRBBND(q,j)==2)
  {
	  if(approx(d.x,nd,nrb_eps)==1 && approx(d.y,-nd,nrb_eps)==1 && approx(d.z,0,nrb_eps)==1) //pos 0
	  {
	    NRBI(q,j,0)=NUMMER(p,i);
	    NRBNEIGH(p,i)=1;
	    nrb_neighs++;

	    //if(NRBBND(p,i)==0)
	    //NRBI(p,i,0)=NUMMER(q,j);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 0 for nrj:%d,cj:%d,q:%p -> nri:%d,ci:%d,p:%p,ptype:%d,qtype:%d\n", myid,NUMMER(q,j),c2,q, NUMMER(p,i),c1,p,p->celltype,q->celltype);
  #endif
	  }
  }
  //SITE 1 q,j
  if(NRBBND(q,j)==1 || NRBBND(q,j)==3)
  {
	  if(approx(d.x,nd,nrb_eps)==1 && approx(d.y,nd,nrb_eps)==1 && approx(d.z,0,nrb_eps)==1) //pos 1
	  {
	    NRBI(q,j,1)=NUMMER(p,i);
	    NRBNEIGH(p,i)=1;
	    nrb_neighs++;

	    //if(NRBBND(p,i)==0)
	    //NRBI(p,i,1)=NUMMER(q,j);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 1 for nrj:%d,cj:%d,q:%p -> nri:%d,ci:%d,p:%p,ptype:%d,qtype:%d\n", myid,NUMMER(q,j),c2,q, NUMMER(p,i),c1,p,p->celltype,q->celltype);
  #endif
	  }
  }
  //SITE 2 q,j
  if(NRBBND(q,j)==3)
  {
	  if(approx(d.x,-nd,nrb_eps)==1 && approx(d.y,nd,nrb_eps)==1 && approx(d.z,0,nrb_eps)==1) //pos 2
	  {
	    NRBI(q,j,2)=NUMMER(p,i);
	    NRBNEIGH(p,i)=1;
	    nrb_neighs++;

	    //if(NRBBND(p,i)==0)
	    //NRBI(p,i,2)=NUMMER(q,j);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 2 for nrj:%d,cj:%d,q:%p -> nri:%d,ci:%d,p:%p,ptype:%d,qtype:%d\n", myid,NUMMER(q,j),c2,q, NUMMER(p,i),c1,p,p->celltype,q->celltype);
  #endif
	  }
  }


  //SITE 3 q,j
  if(NRBBND(q,j)==2)
  {
	  if(approx(d.x,-nd,nrb_eps)==1 && approx(d.y,-nd,nrb_eps)==1 && approx(d.z,0,nrb_eps)==1) //pos 3
	  {
	    NRBI(q,j,3)=NUMMER(p,i);
	    NRBNEIGH(p,i)=1;
	    nrb_neighs++;

	    //if(NRBBND(p,i)==0)
	    //NRBI(p,i,3)=NUMMER(q,j);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 3 for nrj:%d,cj:%d,q:%p -> nri:%d,ci:%d,p:%p,ptype:%d,qtype:%d\n", myid,NUMMER(q,j),c2,q, NUMMER(p,i),c1,p,p->celltype,q->celltype);
  #endif
	  }
  }
  ////////////////
  // XZ-PLANE   //
  ////////////////
  //SITE 9 q,j
  if(NRBBND(q,j)==1)
  {
	  if(approx(d.x,nd,nrb_eps)==1 && approx(d.y,0,nrb_eps)==1 && approx(d.z,nd,nrb_eps)==1) //pos 9
	  {
	    NRBI(q,j,9)=NUMMER(p,i);
	    NRBNEIGH(p,i)=1;
	    nrb_neighs++;

	    //if(NRBBND(p,i)==0)
	    //NRBI(p,i,9)=NUMMER(q,j);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 9 for nrj:%d,cj:%d,q:%p -> nri:%d,ci:%d,p:%p,ptype:%d,qtype:%d\n", myid,NUMMER(q,j),c2,q, NUMMER(p,i),c1,p,p->celltype,q->celltype);
  #endif
	  }
  }
  //SITE 8 q,j
  if(NRBBND(q,j)==1)
  {
	  if(approx(d.x,nd,nrb_eps)==1 && approx(d.y,0,nrb_eps)==1 && approx(d.z,-nd,nrb_eps)==1) //pos 8
	  {
	    NRBI(q,j,8)=NUMMER(p,i);
	    NRBNEIGH(p,i)=1;
	    nrb_neighs++;

	    //if(NRBBND(p,i)==0)
	    //NRBI(p,i,8)=NUMMER(q,j);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 8 for nrj:%d,cj:%d,q:%p -> nri:%d,ci:%d,p:%p,ptype:%d,qtype:%d\n", myid,NUMMER(q,j),c2,q, NUMMER(p,i),c1,p,p->celltype,q->celltype);
  #endif

	  }
  }
  //////////////////
  //  YZ-PLANE
  /////////////////
  //SITE 4 q,j
  if(NRBBND(q,j)==3)
  {
	  if(approx(d.x,0,nrb_eps)==1 && approx(d.y,nd,nrb_eps)==1 && approx(d.z,-nd,nrb_eps)==1) //pos 4
	  {
	    NRBI(q,j,4)=NUMMER(p,i);
	    NRBNEIGH(p,i)=1;
	    nrb_neighs++;

	    //if(NRBBND(p,i)==0)
	    //NRBI(p,i,4)=NUMMER(q,j);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 4 for nrj:%d,cj:%d,q:%p -> nri:%d,ci:%d,p:%p,ptype:%d,qtype:%d\n", myid,NUMMER(q,j),c2,q, NUMMER(p,i),c1,p,p->celltype,q->celltype);
  #endif
	  }
  }
  //SITE 5 q,j
  if(NRBBND(q,j)==3)
  {
	  if(approx(d.x,0,nrb_eps)==1 && approx(d.y,nd,nrb_eps)==1 && approx(d.z,nd,nrb_eps)==1) //pos 5
	  {
	    NRBI(q,j,5)=NUMMER(p,i);
	    NRBNEIGH(p,i)=1;
	    nrb_neighs++;

	    //if(NRBBND(p,i)==0)
	    //NRBI(p,i,5)=NUMMER(q,j);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 5 for nrj:%d,cj:%d,q:%p -> nri:%d,ci:%d,p:%p,ptype:%d,qtype:%d\n", myid,NUMMER(q,j),c2,q, NUMMER(p,i),c1,p,p->celltype,q->celltype);
  #endif
	  }
  }


  //SITE 7 q,j
  if(NRBBND(q,j)==2)
  {
	  if(approx(d.x,0,nrb_eps)==1 && approx(d.y,-nd,nrb_eps)==1 && approx(d.z,-nd,nrb_eps)==1) //pos 7
	  {
	    NRBI(q,j,7)=NUMMER(p,i);
	    NRBNEIGH(p,i)=1;
	    nrb_neighs++;

	    //if(NRBBND(p,i)==0)
	    //NRBI(p,i,7)=NUMMER(q,j);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 7 for nrj:%d,cj:%d,q:%p -> nri:%d,ci:%d,p:%p,ptype:%d,qtype:%d\n", myid,NUMMER(q,j),c2,q, NUMMER(p,i),c1,p,p->celltype,q->celltype);
  #endif

	  }
  }
  //SITE 6 q,j
  if(NRBBND(q,j)==2)
  {
	  if(approx(d.x,0,nrb_eps)==1 && approx(d.y,-nd,nrb_eps)==1 && approx(d.z,nd,nrb_eps)==1) //pos 6
	  {
	    NRBI(q,j,6)=NUMMER(p,i);
	    NRBNEIGH(p,i)=1;
	    nrb_neighs++;
	    //if(NRBBND(p,i)==0)
	    //NRBI(p,i,2)=NUMMER(q,j);
  #if DEBUG_LEVEL>3
  printf("myid:%d,FOUND 6 for nrj:%d,cj:%d,q:%p -> nri:%d,ci:%d,p:%p,ptype:%d,qtype:%d\n", myid,NUMMER(q,j),c2,q, NUMMER(p,i),c1,p,p->celltype,q->celltype);
  #endif

	  }
  }
	 } //for j
	 } //for m
     } //for i
    } //for k
#endif //NBL


    // Folgende 2 send_cells aufrufe nur 1 mal zu beginn nötig
    // um nrb-info auf allen beteiligten procs zu vereinheitlichen
    // **********************************************************************
    have_valid_nrb=1;
#ifdef LOADBALANCE
sync_cells_direct(copy_nrb_max,pack_nrb,unpack_nrb_max,0); //acumm. results
#else
nrb_send_cells(copy_nrb_max,pack_nrb,unpack_nrb_max); //acumm. results 
#endif

    // nrb_send_cells(copy_nrb_max,pack_nrb,unpack_nrb_max); //acumm. results 
nrb_inverse_send_cells(copy_nrb_max,pack_nrb,unpack_nrb_max);
    // **********************************************************************

    //korrektur-loop: atom darf nicht bnd-atom und neigh-atom gleichzeitig sein! 
    //WARUM? --> weil bei nrb_forces zunächst die impulse aller bnd-atome genullt werden muessen.
    //Da bnd-atome aber die geschwindigkeiten ihrer neighatome brauchen (die in diesem fall selber bnd-atome sein können),
    //führt dies zu fehlenden kraftbeiträgen, wenn die Geschwindigkeiten Null sind!

/*
    for(k=0;k<nallcells;k++)
    {
      cell* p=cell_array+k;
      for(i=0;i<p->n;++i)
	if(NRBBND(p,i)>0) NRBNEIGH(p,i)=0;
    }
*/

    nrb_writerestart(0); 
  } // if(imdrestart==0 && nrb_readfile==0)
  else 
  {
    have_valid_nrb=1;    
    if(nrb_readfile==1) nrb_readinputfile(nrb_input_file);    //anderes such-schema als bei restart (aufwändiger)
    else if(imdrestart>0 && nrb_readfile==0) nrb_readrestart(nrb_restart_file);  
    //d.h. nrb_input_file hat vorrang über restart-file
		

#ifdef LOADBALANCE
sync_cells_direct(copy_nrb_max,pack_nrb,unpack_nrb_max,0); //acumm. results //ACHTUNG: sync_cells macht alles kaputt --> NRB und LB nicht kompatibel
#else
nrb_send_cells(copy_nrb_max,pack_nrb,unpack_nrb_max); //acumm. results 
#endif

#ifdef NBL  //ERROR WENN NUR PAIR OHNE NBL
    //nrb_inverse_send_cells(copy_nrb_max,pack_nrb,unpack_nrb_max); // brauche ich nur im init step
#endif    

    nrb_build_ifromid();
  }


/*
    //int m; //bereits definiert
    int l;//,n; 
  // in imd_geom_3d.c
    for(l=0;l<cell_dim.x;l++)
      for(m=0;m<cell_dim.y;m++)
        for(n=0;n<cell_dim.z;n++)
        {
          cell* p=PTR_3D_V(cell_array, l, m, n, cell_dim);
          if(l==0 || l==cell_dim.x-1 || m==0 || m==cell_dim.y-1 || n==0 || n==cell_dim.z-1)
            p->celltype=2;
          else
            p->celltype=1;
          
        }
*/
    
  //printf("* identified %d bnd-atoms and %d neigh-atoms on proc:%d\n",mybndatoms,myneighatoms,myid);
  MPI_Barrier(cpugrid);
  if(myid==0)  printf("*************************************************\n");

  nrb_writerestart(0); //In jedem Fall ausschreiben zum checken!
  return 0;
}

// ********************************************************************************************************************************
// ROUTINE zum basteln eines 2d-arrays um aus dem atom-index die jeweilige atom-nummer zu in der zellen-liste zu finden
// Aufzurufen nach jedem make_nblist
// ********************************************************************************************************************************
int nrb_build_ifromid(void)
{
  int flag=0; 
  int flag_all=0;

  int k;
  int mybndatomsnew=0;
  int myneighatomsnew=0;
  int i;
  int m;
  int n;

  //Die komplette ROUTINE wird seit 02.08.19 nicht mehr benötigt


//   nrb_send_cells(copy_nrb,pack_nrb,unpack_nrb); //<-- muss. obwohl im normalen send_cells (siehe imd_comm_force_3d.c und imd_mpi_util.c) bereits kommuniziert wird
						//warum? --> k.P., aber führt zu fehlern, wenn hier oder dort auskommentiert.... 
// ********************************************************
  return 0;
}
// **********************************************************************************************************************************
// *  COMP NEW FORCES FOR BND-ATOMS
// **********************************************************************************************************************************
int nrb_forces(void)
{

// if(myid==0) printf("\nsteps:%d, entered nrb_forces\n",steps);

  vektor U_self; //self
  ///yx-plane
  vektor U[12];
  vektor V[12];

  double neighfac=1.0;

  int k;

  double mass=26.9815; // TODO: weiter unten einfach mit MASSE(p,i) tauschen.Caution: checken ob auch sinnvolle werte geliefert werden
  int i;               //       denn manche zellen sind nur buffer-zellen
  for(i=0;i<12;i++)
  {
    U[i].x=V[i].x=0.0;
    U[i].y=V[i].y=0.0;
    U[i].z=V[i].z=0.0;
  }

  vektor U_dot;
  int n,m;
  //CLEAR MOMENTA and SELF-CONTRIBUTIONS
  for (k=0;k<nallcells;k++)
  {
    cell *p=cell_array+k;
    for (i=0; i<p->n; i++)
    {
      IMPULS_ALT(p,i,X)=IMPULS(p,i,X);
      IMPULS_ALT(p,i,Y)=IMPULS(p,i,Y);
      IMPULS_ALT(p,i,Z)=IMPULS(p,i,Z);

      if(NRBBND(p,i)>0)
      {
        IMPULS(p,i,X)=0.0;
        IMPULS(p,i,Z)=0.0;
        IMPULS(p,i,Y)=0.0;

// if(NUMMER(p,i)==WATCHME)
//   printf("myid:%d, in:%d, bnd:%d, neigh:%d, celltype:%d,lbtype:%d\n",
//       myid,NUMMER(p,i),NRBBND(p,i), NRBNEIGH(p,i), p->celltype,p->lb_cell_type);
        // *************************************
        // * SELF CONTRIB only for real atoms
        // **************************************

#ifdef LOADBALANCE
        if(p->lb_cell_type==LB_REAL_CELL) //==1
#else
        if(p->celltype==1)
#endif      
        {
          //Comp. U-self contrib.
          U_self.x=ORT(p,i,X)-REF_POS(p,i,X);
          U_self.y=ORT(p,i,Y)-REF_POS(p,i,Y);
          U_self.z=ORT(p,i,Z)-REF_POS(p,i,Z);
          U_self.z=MINIMGZ(U_self.z);
          if(pbc_dirs.y==1)          
            U_self.y=MINIMGY(U_self.y);

          if(NRBBND(p,i)==1)
          {
            U_dot.x=-nrbk*4.0*U_self.x; //muss 
            U_dot.y=-nrbk*dblsqrt2*U_self.y;
            U_dot.z=-nrbk*dblsqrt2*U_self.z;
          }
          else if(NRBBND(p,i)==2 || NRBBND(p,i)==3)
          {
            U_dot.x=-nrbk*dblsqrt2*U_self.x; //muss 
            U_dot.y=-nrbk*U_self.y*4;
            U_dot.z=-nrbk*dblsqrt2*U_self.z;
          }

          IMPULS(p,i,X)+=U_dot.x*MASSE(p,i);
          IMPULS(p,i,Y)+=U_dot.y*MASSE(p,i);
          IMPULS(p,i,Z)+=U_dot.z*MASSE(p,i);

#if DEBUG_LEVEL>0
	 cell*bndcell=p;
	 int bndi=i;
//if(WATCH)
if(NUMMER(p,i)==WATCHME && (steps % PRINTSTEP ==0) )
printf("myid:%d,steps:%d, DOFORCE2, ind:%d, bnd:%d, dpx:%f dpy:%f dpz:%f px:%f py:%f pz:%f x:%f y:%f z:%f rx:%f,ry:%f,rz:%f\n",
        myid,steps,NUMMER(p,i),NRBBND(p,i),

        U_dot.x*MASSE(p,i),
        U_dot.y*MASSE(p,i),
        U_dot.z*MASSE(p,i),

        IMPULS(p,i,X),IMPULS(p,i,Y),IMPULS(p,i,Z),
	      ORT(p,i,X),ORT(p,i,Y),ORT(p,i,Z),
	      REF_POS(p,i,X),REF_POS(p,i,Y),REF_POS(p,i,Z));
//        KRAFT(p,i,X),KRAFT(p,i,Y),KRAFT(p,i,Z));
#endif

        }
      }
    }
  }
  // *************************
  // NOW NEIGH-CONTRIBUTIONS
  // *************************
  n=0;
  cell* bndcell;
  cell* neighcell;
  int bndi,neighi;

#ifdef NBL  
  for (k=0; k<ncells; k++) 
  {
    cell *p = cell_array + cnbrs[k].np;

    for (i=0; i<p->n; i++) 
    {
      if(NRBBND(p,i)==0 && NRBNEIGH(p,i)==0) //weder bnd noch neigh --> dont waste time
      {
	n++;continue;
      }
      // loop over neighbors
      for (m=tl[n]; m<tl[n+1]; m++) 
      {
        cell   *q;
        int    c, j;

        c = cl_num[ tb[m] ];
        j = tb[m] - cl_off[c];

        q = cell_array + c;
        int r=0;
	for(r=0;r<12;r++)
	{	  
	  if(NRBBND(p,i)>0 && NUMMER(q,j)==NRBI(p,i,r)) 
	  {
	    bndcell=p;
	    bndi=i;
	    neighcell=q;
	    neighi=j; 
	  }
	  else if(NRBBND(q,j)>0 && NUMMER(p,i)==NRBI(q,j,r))
	  {
	    bndcell=q;
	    bndi=j;
	    neighcell=p;
	    neighi=i;
	  }
	  else continue;
/*
	  if(NRBBND(neighcell,neighi)>0) 
	  {
	    continue;
	  } //quick'n dirty weil impulse gelöscht bzw. evtl bereits self-term (sehr groß) erhalalten
*/						   //lässt sich nur sauber lösen wenn atome gleichzeitig bnd- und neighatome sein können
          U_dot.x=U_dot.y=U_dot.z=0.0;

          //displacement from refpos
          U[r].x=ORT(neighcell,neighi,X)-REF_POS(neighcell,neighi,X);
          U[r].y=ORT(neighcell,neighi,Y)-REF_POS(neighcell,neighi,Y);
          U[r].z=ORT(neighcell,neighi,Z)-REF_POS(neighcell,neighi,Z);
 
          U[r].z=MINIMGZ(U[r].z); //ACHTUNG: Bei 1D-Simulation, muss U[r].y auch mittels MINIMGY berechnet werden!
          if(pbc_dirs.y==1)          
            U[r].y=MINIMGY(U[r].y);

          //velocities
//          V[r].x=IMPULS(neighcell,neighi,X)/mass;
//          V[r].y=IMPULS(neighcell,neighi,Y)/mass;
//          V[r].z=IMPULS(neighcell,neighi,Z)/mass;

          V[r].x=IMPULS_ALT(neighcell,neighi,X)/mass;
          V[r].y=IMPULS_ALT(neighcell,neighi,Y)/mass;
          V[r].z=IMPULS_ALT(neighcell,neighi,Z)/mass;


 	  if(NRBBND(bndcell,bndi)==1)
   	  {
            U_dot.x=nrbk*U[r].x;
            U_dot.y=nrbk*sqrt2half*U[r].y;
            U_dot.z=nrbk*sqrt2half*U[r].z;
	  }
	  else if(NRBBND(bndcell,bndi)==2 || NRBBND(bndcell,bndi)==3) //BND-typ 1 hat vorrang
	  {
            U_dot.x=nrbk*sqrt2half*U[r].x;
            U_dot.y=nrbk*U[r].y;
            U_dot.z=nrbk*sqrt2half*U[r].z;
	  }


          U_dot.x-=0.25*V[r].x*neighfac; //push from neigh
          U_dot.y-=0.25*V[r].y*neighfac;
          U_dot.z-=0.25*V[r].z*neighfac;

          IMPULS(bndcell,bndi,X)+=U_dot.x*mass;
          IMPULS(bndcell,bndi,Y)+=U_dot.y*mass;
          IMPULS(bndcell,bndi,Z)+=U_dot.z*mass;


#if DEBUG_LEVEL>0
if(NUMMER(bndcell,bndi)==WATCHME && (steps % PRINTSTEP ==0 ))
//if(WATCH) //NUMMER(bndcell,bndi)==WATCHME && steps>PRINTSTEP)
printf("myid:%d,steps:%d, DOFORCE, ind:%d, bnd:%d,nei:%d,r:%d, dpx:%f dpy:%f dpz:%f px:%f py:%f pz:%f btyp:%d,ntyp:%d,x:%f,y:%f,z:%f,rx:%f,ry:%f,rz:%f\n",
        myid,steps,NUMMER(bndcell,bndi),NRBBND(bndcell,bndi),NUMMER(neighcell,neighi),r,
        U_dot.x*mass,
        U_dot.y*mass,
        U_dot.z*mass,

        IMPULS(bndcell,bndi,X),IMPULS(bndcell,bndi,Y),IMPULS(bndcell,bndi,Z),
        bndcell->celltype,
        neighcell->celltype,
        ORT(neighcell,neighi,X),ORT(neighcell,neighi,Y),ORT(neighcell,neighi,Z),
        REF_POS(neighcell,neighi,X),REF_POS(neighcell,neighi,Y),REF_POS(neighcell,neighi,Z)); //vom nachbar-atom
#endif

	}//for r
      } //for m
      n++;
    }//for i
  } //for k


  #endif //NBL

#ifndef NBL
#ifdef PAIR
  for (n=0; n<nlists; ++n) 
  {
    for (k=0; k<npairs[n]; ++k) 
    {
      pair *P;
      P = pairs[n] + k;
      cell*p=cell_array+P->np;
      cell*q=cell_array+P->nq;

      vektor pbc;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;


      for (i=0; i<p->n; ++i) 
      {
         int j,jstart;
         if(NRBBND(p,i)==0 && NRBNEIGH(p,i)==0) //weder bnd noch neigh --> dont waste time
          continue;

         jstart = (((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0)) ? i+1 : 0);

         //Loop over neighs
         for (j = jstart; j < q->n; ++j)
         {

            int r=0;
            for(r=0;r<12;r++)
            {   
              if(NRBBND(p,i)>0 && NUMMER(q,j)==NRBI(p,i,r)) 
              {
                bndcell=p;
                bndi=i;
                neighcell=q;
                neighi=j; 
              }
              else if(NRBBND(q,j)>0 && NUMMER(p,i)==NRBI(q,j,r))
              {
                bndcell=q;
                bndi=j;
                neighcell=p;
                neighi=i;
              }
              else continue;

              U_dot.x=U_dot.y=U_dot.z=0.0;

              //displacement from refpos
              U[r].x=ORT(neighcell,neighi,X)-REF_POS(neighcell,neighi,X);
              U[r].y=ORT(neighcell,neighi,Y)-REF_POS(neighcell,neighi,Y);
              U[r].z=ORT(neighcell,neighi,Z)-REF_POS(neighcell,neighi,Z);
     
              U[r].z=MINIMGZ(U[r].z); //ACHTUNG: Bei 1D-Simulation, muss U[r].y auch mittels MINIMGY berechnet werden!
              if(pbc_dirs.y==1)          
                U[r].y=MINIMGY(U[r].y);


              V[r].x=IMPULS_ALT(neighcell,neighi,X)/mass;
              V[r].y=IMPULS_ALT(neighcell,neighi,Y)/mass;
              V[r].z=IMPULS_ALT(neighcell,neighi,Z)/mass;


              if(NRBBND(bndcell,bndi)==1)
              {
                      U_dot.x=nrbk*U[r].x;
                      U_dot.y=nrbk*sqrt2half*U[r].y;
                      U_dot.z=nrbk*sqrt2half*U[r].z;
              }
              else if(NRBBND(bndcell,bndi)==2 || NRBBND(bndcell,bndi)==3) //BND-typ 1 hat vorrang
              {
                      U_dot.x=nrbk*sqrt2half*U[r].x;
                      U_dot.y=nrbk*U[r].y;
                      U_dot.z=nrbk*sqrt2half*U[r].z;
              }


              U_dot.x-=0.25*V[r].x*neighfac; //push from neigh
              U_dot.y-=0.25*V[r].y*neighfac;
              U_dot.z-=0.25*V[r].z*neighfac;

              IMPULS(bndcell,bndi,X)+=U_dot.x*mass;
              IMPULS(bndcell,bndi,Y)+=U_dot.y*mass;
              IMPULS(bndcell,bndi,Z)+=U_dot.z*mass;

#if DEBUG_LEVEL>0
if(NUMMER(bndcell,bndi)==WATCHME && (steps % PRINTSTEP ==0 ))
//if(WATCH) //NUMMER(bndcell,bndi)==WATCHME && steps>PRINTSTEP)
printf("myid:%d,steps:%d, DOFORCE, ind:%d, bnd:%d,nei:%d,r:%d, dpx:%f dpy:%f dpz:%f px:%f py:%f pz:%f btyp:%d,ntyp:%d,x:%f,y:%f,z:%f,rx:%f,ry:%f,rz:%f\n",
        myid,steps,NUMMER(bndcell,bndi),NRBBND(bndcell,bndi),NUMMER(neighcell,neighi),r,
        U_dot.x*mass,
        U_dot.y*mass,
        U_dot.z*mass,

        IMPULS(bndcell,bndi,X),IMPULS(bndcell,bndi,Y),IMPULS(bndcell,bndi,Z),
        bndcell->celltype,
        neighcell->celltype,
        ORT(neighcell,neighi,X),ORT(neighcell,neighi,Y),ORT(neighcell,neighi,Z),
        REF_POS(neighcell,neighi,X),REF_POS(neighcell,neighi,Y),REF_POS(neighcell,neighi,Z)); //vom nachbar-atom
#endif

            }//for r
          } //for j
        }//for i
      }//for k
    }//for n
  
#endif
#endif


  return 0;
}



// *************************************************************************************************************
// copy func wird aufgerufen bei kommunikation mit sich selbst (einfach kopieren)
// copy_nrb_max wird nur 1 mal zu beginn von nrb_inverse_send_cells aufgerufen 
// hier wird geprüft ob ergebnisse aus buf-cells in real-cells übernommen werden sollen oder nicht
// ***********************************************************************************************************
void copy_nrb_max( int k, int l, int m, int r, int s, int t, vektor v )
{
  //Schickt immer von realcell nach buffercell
  int i;
  minicell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);


  int r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11;
  int isbnd,isneigh;

  for (i=0; i<to->n; ++i) {
    r0=MAX(NRBI(from,i,0),NRBI(to,i,0));
    r1=MAX(NRBI(from,i,1),NRBI(to,i,1));
    r2=MAX(NRBI(from,i,2),NRBI(to,i,2));
    r3=MAX(NRBI(from,i,3),NRBI(to,i,3));
    r4=MAX(NRBI(from,i,4),NRBI(to,i,4));
    r5=MAX(NRBI(from,i,5),NRBI(to,i,5));
    r6=MAX(NRBI(from,i,6),NRBI(to,i,6));
    r7=MAX(NRBI(from,i,7),NRBI(to,i,7));
    r8=MAX(NRBI(from,i,8),NRBI(to,i,8));
    r9=MAX(NRBI(from,i,9),NRBI(to,i,9));
    r10=MAX(NRBI(from,i,10),NRBI(to,i,10));
    r11=MAX(NRBI(from,i,11),NRBI(to,i,11));
    isbnd=MAX(NRBBND(from,i),NRBBND(to,i));
    isneigh=MAX(NRBNEIGH(from,i),NRBNEIGH(to,i));

     NRBI(to,i,0)  =r0; 
     NRBI(to,i,1)  =r1; 
     NRBI(to,i,2)  =r2; 
     NRBI(to,i,3)  =r3; 
     NRBI(to,i,4)  =r4; 
     NRBI(to,i,5)  =r5; 
     NRBI(to,i,6)  =r6; 
     NRBI(to,i,7)  =r7; 
     NRBI(to,i,8)  =r8; 
     NRBI(to,i,9)  =r9; 
     NRBI(to,i,10) =r10;
     NRBI(to,i,11) =r11; 
     NRBBND(to,i)  =isbnd; 
     NRBNEIGH(to,i)=isneigh; 

     REF_POS(to,i,X)=REF_POS(from,i,X);
     REF_POS(to,i,Y)=REF_POS(from,i,Y);
     REF_POS(to,i,Z)=REF_POS(from,i,Z);  
  }
}
// ************************************************* 
//    normales copy (nicht max)
//  wird aufgerufen von sendcells 
// ****************************************************
void copy_nrb( int k, int l, int m, int r, int s, int t, vektor v )
{
  //Schickt immer von realcell nach buffercell
  int i;
  minicell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

  for (i=0; i<to->n; ++i) {
     NRBI(to,i,0)  =NRBI(from,i,0);
     NRBI(to,i,1)  =NRBI(from,i,1);//=r1;
     NRBI(to,i,2)  =NRBI(from,i,2);//=r2;
     NRBI(to,i,3)  =NRBI(from,i,3);//=r3;
     NRBI(to,i,4)  =NRBI(from,i,4);//=r4;
     NRBI(to,i,5)  =NRBI(from,i,5);//=r5;
     NRBI(to,i,6)  =NRBI(from,i,6);//=r6;
     NRBI(to,i,7)  =NRBI(from,i,7);//=r7;
     NRBI(to,i,8)  =NRBI(from,i,8);//=r8;
     NRBI(to,i,9)  =NRBI(from,i,9);//=r9;
     NRBI(to,i,10) =NRBI(from,i,10);//=r10;
     NRBI(to,i,11) =NRBI(from,i,11);//=r11;
     NRBBND(to,i)  =NRBBND(from,i);
     NRBNEIGH(to,i)=NRBNEIGH(from,i);

     REF_POS(to,i,X)=REF_POS(from,i,X);
     REF_POS(to,i,Y)=REF_POS(from,i,Y);
     REF_POS(to,i,Z)=REF_POS(from,i,Z);  
  }
}


// ***************************************************************
// * void pack_nrb( msgbuf *b, int k, int l, int m) //forces
// * aufgerufen von send_cells und nrb_inverse_send_cells
// ****************************************************************
void pack_nrb( msgbuf *b, int k, int l, int m, vektor v )
{
  int i, j = b->n;
  minicell *from;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  b->data[ j++ ] = (real) from->n;

  for (i=0; i<from->n; ++i) {
    b->data[ j++]  = NRBI(from,i,0);
    b->data[ j++]  = NRBI(from,i,1);
    b->data[ j++]  = NRBI(from,i,2);
    b->data[ j++]  = NRBI(from,i,3);
    b->data[ j++]  = NRBI(from,i,4);
    b->data[ j++]  = NRBI(from,i,5);
    b->data[ j++]  = NRBI(from,i,6);
    b->data[ j++]  = NRBI(from,i,7);
    b->data[ j++]  = NRBI(from,i,8);
    b->data[ j++]  = NRBI(from,i,9);
    b->data[ j++]  = NRBI(from,i,10);
    b->data[ j++]  = NRBI(from,i,11);
    b->data[j++]   = NRBBND(from,i); 
    b->data[j++]   = NRBNEIGH(from,i);

    b->data[j++]=REF_POS(from,i,X); 
    b->data[j++]=REF_POS(from,i,Y);
    b->data[j++]=REF_POS(from,i,Z);
  }

  b->n = j;
  if (b->n_max < b->n)
  {    
    char errstr[255];
    sprintf(errstr,"Buffer overflow in pack_nrb - increase msgbuf_size. b->n_max:%d, b->n:%d",b->n_max,b->n);
    error(errstr);
  }
  
}



// *****************************************************************************
// *  unpack cell from MPI buffer to buffer cell
// *  wird aufgerufen von nrb_inverse_send_cells
// *  genau wird nrb_copy_max wird hier geprüft ob daten aus buf-cells in real-cells 
// *  übernommen werden sollen.
// *  geschieht nur 1 mal zu beginn
// ******************************************************************************/
void unpack_nrb_max( msgbuf *b, int k, int l, int m )
{
  int i, j = b->n,tmp_n;;
  minicell *to;

  to = PTR_3D_V(cell_array, k, l, m, cell_dim);
  
  tmp_n = (int) b->data[ j++ ];

  // increase minicell size if necessary  
  if (tmp_n > to->n_max) {  
    to->n = 0;
    ALLOC_MINICELL(to, tmp_n);
  }

  /* copy indices and atoms */
  to->n = tmp_n;


  int r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11;
  int isbnd,isneigh;

  for (i=0; i<to->n; ++i) {

    r0=b->data[ j++];
    r1=b->data[ j++];
    r2=b->data[ j++];
    r3=b->data[ j++];
    r4=b->data[ j++];
    r5=b->data[ j++];
    r6=b->data[ j++];
    r7=b->data[ j++];
    r8=b->data[ j++];
    r9=b->data[ j++];
    r10=b->data[ j++];
    r11=b->data[ j++];
    isbnd=b->data[ j++];
    isneigh=b->data[j++];

    r0=MAX(r0,NRBI(to,i,0));
    r1=MAX(r1,NRBI(to,i,1));
    r2=MAX(r2,NRBI(to,i,2));
    r3=MAX(r3,NRBI(to,i,3));
    r4=MAX(r4,NRBI(to,i,4));
    r5=MAX(r5,NRBI(to,i,5));
    r6=MAX(r6,NRBI(to,i,6));
    r7=MAX(r7,NRBI(to,i,7));
    r8=MAX(r8,NRBI(to,i,8));
    r9=MAX(r9,NRBI(to,i,9));
    r10=MAX(r10,NRBI(to,i,10));
    r11=MAX(r11,NRBI(to,i,11));
    isbnd=MAX(isbnd,NRBBND(to,i));
    isneigh=MAX(isneigh,NRBNEIGH(to,i));

    NRBI(to,i,0)= r0;
    NRBI(to,i,1)= r1;
    NRBI(to,i,2)= r2;
    NRBI(to,i,3)= r3;
    NRBI(to,i,4)= r4;
    NRBI(to,i,5)= r5;
    NRBI(to,i,6)= r6;
    NRBI(to,i,7)= r7;
    NRBI(to,i,8)= r8; 
    NRBI(to,i,9)= r9;
    NRBI(to,i,10)= r10;
    NRBI(to,i,11)= r11;
    NRBBND(to,i) = isbnd;
    NRBNEIGH(to,i)=isneigh;
    REF_POS(to,i,X)=b->data[j++];
    REF_POS(to,i,Y)=b->data[j++];
    REF_POS(to,i,Z)=b->data[j++];

  }

  b->n = j;
  if (b->n_max < b->n)
    error("Buffer overflow in unpack_nrb_max - increase msgbuf_size");
}


// Normales unpack
void unpack_nrb( msgbuf *b, int k, int l, int m )
{
  int i, j = b->n,tmp_n;;
  minicell *to;

  to = PTR_3D_V(cell_array, k, l, m, cell_dim);
  
  tmp_n = (int) b->data[ j++ ];

  /* increase minicell size if necessary */

  if (tmp_n > to->n_max) { //ACHTUNG:Sollte nicht passieren. Die Zelle (und das atom) ist ja schon da
    to->n = 0;
    ALLOC_MINICELL(to, tmp_n);
  }

  /* copy indices and atoms */
  to->n = tmp_n;

  for (i=0; i<to->n; ++i) {
    NRBI(to,i,0)= b->data[j++];
    NRBI(to,i,1)= b->data[j++];
    NRBI(to,i,2)= b->data[j++];
    NRBI(to,i,3)= b->data[j++];
    NRBI(to,i,4)= b->data[j++];
    NRBI(to,i,5)= b->data[j++];
    NRBI(to,i,6)= b->data[j++];
    NRBI(to,i,7)= b->data[j++];
    NRBI(to,i,8)= b->data[j++];
    NRBI(to,i,9)= b->data[j++];
    NRBI(to,i,10)= b->data[j++];
    NRBI(to,i,11)= b->data[j++];
    NRBBND(to,i) = b->data[j++];
    NRBNEIGH(to,i)=b->data[j++];
    REF_POS(to,i,X)=b->data[j++];
    REF_POS(to,i,Y)=b->data[j++];
    REF_POS(to,i,Z)=b->data[j++];
  }

  b->n = j;
  if (b->n_max < b->n)
  {
    char errstr[255];
    sprintf(errstr,"Buffer overflow in unpack_nrb - increase msgbuf_size. b->n_max:%d, b->n:%d",b->n_max,b->n);
    error(errstr);
  }
}

// ******************************************************
// * SENDFUNC zum mergen mit buffcell-nrb-data
// * bedingtes kopieren von buffer-cells nach real-cells um 
// * nrb-info zu vereinheitlichen 
// *******************************************************
void nrb_inverse_send_cells(void (*copy_func)  (int, int, int, int, int, int, vektor),
                	    void (*pack_func)  (msgbuf*, int, int, int, vektor),
	                    void (*unpack_func)(msgbuf*, int, int, int))
{

  /* Remember:
   * east -> -x
   * west -> +x
   * north-> -y
   * south-> +y
   * up   -> -z
   * down -> +z
   * *************/


  int i,j;

  vektor uvec={0,0,0}, dvec={0,0,0};
  vektor nvec={0,0,0}, svec={0,0,0};
  vektor evec={0,0,0}, wvec={0,0,0};

#ifdef MPI
  MPI_Status  stat;
  empty_mpi_buffers();
#endif

#ifdef NBLIST
  if (pbc_dirs.x==1) {
    if (my_coord.x==0) evec = box_x;
    if (my_coord.x==cpu_dim.x-1) {
      wvec.x = -box_x.x; wvec.y = -box_x.y; wvec.z = -box_x.z;
    }
  }
  if (pbc_dirs.y==1) {
    if (my_coord.y==0) nvec = box_y;
    if (my_coord.y==cpu_dim.y-1) {
      svec.x = -box_y.x; svec.y = -box_y.y; svec.z = -box_y.z;
    }
  }
  if (pbc_dirs.z==1) {
    if (my_coord.z==0) uvec = box_z;
    if (my_coord.z==cpu_dim.z-1) {
      dvec.x = -box_z.x; dvec.y = -box_z.y; dvec.z = -box_z.z;
    }
  }
#endif
  /* exchange up/down */
  // up: -z, down: +z
  if (cpu_dim.z==1) {
    // simply copy up/down atoms to buffer cells
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j) {
        //(*copy_func)( i, j, 1, i, j, cell_dim.z-1, uvec );
        //(*copy_func)( i, j, cell_dim.z-2, i, j, 0, dvec );
	(*copy_func) (i,j,0,i,j,cell_dim.z-2,uvec);
	(*copy_func) (i,j,cell_dim.z-1,i,j,1,dvec);
      }
  }
#ifdef MPI
  else {
    /* copy up atoms into send buffer */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_up, i, j, 0, uvec );

    /* send up, receive down */
    sendrecv_buf(&send_buf_up, nbup, &recv_buf_down, nbdown, &stat);

    /* upack atoms from down */
    recv_buf_down.n = 0;
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_down, i, j, cell_dim.z-2 );

    /* copy down atoms into send buffer */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_down, i, j, cell_dim.z-1, dvec );

    /* send down, receive up */
    sendrecv_buf(&send_buf_down, nbdown, &recv_buf_up, nbup, &stat);

    /* unpack atoms from up */
    recv_buf_up.n = 0;
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_up, i, j, 1 );
  }
#endif
  /* exchange north/south */
  // north: -y
  // south: +y
  if (cpu_dim.y==1) {
    /* simply copy north/south atoms to buffer cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*copy_func)( i, 0, j, i, cell_dim.y-2, j, nvec );
        (*copy_func)( i, cell_dim.y-1, j, i, 1, j, svec );
      }
  }
#ifdef MPI
  else {
    /* copy north atoms into send buffer */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_north, i, 0, j, nvec );

    /* send north, receive south */
    sendrecv_buf(&send_buf_north, nbnorth, &recv_buf_south, nbsouth, &stat);

    /* unpack atoms from south */
    recv_buf_south.n = 0;
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_south, i, cell_dim.y-2, j );

    /* copy south atoms into send buffer */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_south, i, cell_dim.y-1, j, svec );

    /* send south, receive north */
    sendrecv_buf(&send_buf_south, nbsouth, &recv_buf_north, nbnorth, &stat);
    /* unpack atoms from north */
    recv_buf_north.n = 0;
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_north, i, 1, j );
  }
#endif

  /* exchange east/west */
  // east : -x
  // west:  +x
  if (cpu_dim.x==1) {
    /* simply copy east/west atoms to buffer cells */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*copy_func)( 0, i, j, cell_dim.x-2, i, j, evec );
//#if !defined(AR) || defined(COVALENT)
        (*copy_func)( cell_dim.x-1, i, j, 1, i, j, wvec );
//#endif
      }
  }
#ifdef MPI
  else {
    //  copy east atoms into send buffer
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_east, 0, i, j, evec );

    // send east, receive west 
    sendrecv_buf(&send_buf_east, nbeast, &recv_buf_west, nbwest, &stat);

    // unpack atoms from west
    recv_buf_west.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_west, cell_dim.x-2, i, j );
//#if !defined(AR) || defined(COVALENT)
    // copy west atoms into send buffer 
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_west, cell_dim.x-1, i, j, wvec );

    // send west, receive east 
    sendrecv_buf(&send_buf_west, nbwest, &recv_buf_east, nbeast, &stat);

    // unpack atoms from east
    recv_buf_east.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_east, 1, i, j );
//#endif
  }
#endif
}
// ***********************************
// * NORMALES SEND CELL, OHNE DEN QUATSCH MIT IFDEF AR ...
// ******************************************************
void nrb_send_cells(void (*copy_func)  (int, int, int, int, int, int, vektor),
                void (*pack_func)  (msgbuf*, int, int, int, vektor),
                void (*unpack_func)(msgbuf*, int, int, int))
{
  int i,j;

  vektor uvec={0,0,0}, dvec={0,0,0};
  vektor nvec={0,0,0}, svec={0,0,0};
  vektor evec={0,0,0}, wvec={0,0,0};

#ifdef MPI
  MPI_Status  stat;
  empty_mpi_buffers();
#endif

#ifdef VEC
  atoms.n_buf = atoms.n;
#endif
#ifdef NBLIST
  if (pbc_dirs.x==1) {
    if (my_coord.x==0) evec = box_x;
    if (my_coord.x==cpu_dim.x-1) {
      wvec.x = -box_x.x; wvec.y = -box_x.y; wvec.z = -box_x.z;
    }
  }
  if (pbc_dirs.y==1) {
    if (my_coord.y==0) nvec = box_y;
    if (my_coord.y==cpu_dim.y-1) {
      svec.x = -box_y.x; svec.y = -box_y.y; svec.z = -box_y.z;
    }
  }
  if (pbc_dirs.z==1) {
    if (my_coord.z==0) uvec = box_z;
    if (my_coord.z==cpu_dim.z-1) {
      dvec.x = -box_z.x; dvec.y = -box_z.y; dvec.z = -box_z.z;
    }
  }
#endif
  /* exchange up/down */
  if (cpu_dim.z==1) {
    /* simply copy up/down atoms to buffer cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j) {
        (*copy_func)( i, j, 1, i, j, cell_dim.z-1, uvec );
        (*copy_func)( i, j, cell_dim.z-2, i, j, 0, dvec );
      }
  }
#ifdef MPI
  else {
    /* copy up atoms into send buffer */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_up, i, j, 1, uvec );

    /* send up, receive down */
    sendrecv_buf(&send_buf_up, nbup, &recv_buf_down, nbdown, &stat);

    /* upack atoms from down */
    recv_buf_down.n = 0;
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_down, i, j, cell_dim.z-1 );

    /* copy down atoms into send buffer */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_down, i, j, cell_dim.z-2, dvec );

    /* send down, receive up */
    sendrecv_buf(&send_buf_down, nbdown, &recv_buf_up, nbup, &stat);

    /* unpack atoms from up */
    recv_buf_up.n = 0;
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_up, i, j, 0 );
  }
#endif
  /* exchange north/south */
  if (cpu_dim.y==1) {
    /* simply copy north/south atoms to buffer cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*copy_func)( i, 1, j, i, cell_dim.y-1, j, nvec );
        (*copy_func)( i, cell_dim.y-2, j, i, 0, j, svec );
      }
  }
#ifdef MPI
  else {
    /* copy north atoms into send buffer */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_north, i, 1, j, nvec );

    /* send north, receive south */
    sendrecv_buf(&send_buf_north, nbnorth, &recv_buf_south, nbsouth, &stat);

    /* unpack atoms from south */
    recv_buf_south.n = 0;
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_south, i, cell_dim.y-1, j );

    /* copy south atoms into send buffer */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_south, i, cell_dim.y-2, j, svec );

    /* send south, receive north */
    sendrecv_buf(&send_buf_south, nbsouth, &recv_buf_north, nbnorth, &stat);

    /* unpack atoms from north */
    recv_buf_north.n = 0;
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_north, i, 0, j );
  }
#endif
  /* exchange east/west */
  if (cpu_dim.x==1) {
    /* simply copy east/west atoms to buffer cells */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*copy_func)( 1, i, j, cell_dim.x-1, i, j, evec );
//#if !defined(AR) || defined(COVALENT)
        (*copy_func)( cell_dim.x-2, i, j, 0, i, j, wvec );
//#endif
      }
  }
#ifdef MPI
  else {
    /* copy east atoms into send buffer */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_east, 1, i, j, evec );

    /* send east, receive west */
    sendrecv_buf(&send_buf_east, nbeast, &recv_buf_west, nbwest, &stat);

    /* unpack atoms from west */
    recv_buf_west.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_west, cell_dim.x-1, i, j );

//#if !defined(AR) || defined(COVALENT)
    /* copy west atoms into send buffer */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_west, cell_dim.x-2, i, j, wvec );

    /* send west, receive east */
    sendrecv_buf(&send_buf_west, nbwest, &recv_buf_east, nbeast, &stat);

    /* unpack atoms from east */
    recv_buf_east.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_east, 0, i, j );
//#endif
  }
#endif
}





// *******************************************************
// *  ANALOG zu send_forces (mit SR statt isend/irecv)
// ******************************************************
void nrb_sendmomenta(void (*add_func)   (int, int, int, int, int, int),
              void (*pack_func)  (msgbuf*, int, int, int),
              void (*unpack_func)(msgbuf*, int, int, int))
{
  int i,j;

#ifdef MPI
  MPI_Status  stat;
  empty_mpi_buffers();
#endif

  /* send forces east/west */
  if (cpu_dim.x==1) {
    /* simply add east/west forces to original cells */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) {
//#if defined(COVALENT)
        (*add_func)( 0, i, j, cell_dim.x-2, i, j );
//#endif
        (*add_func)( cell_dim.x-1, i, j, 1, i, j );
      }
  }
#ifdef MPI
  else {
//#if defined(COVALENT)
    /* copy east forces into send buffer */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_east, 0, i, j );

    /* send east, receive west */
    sendrecv_buf(&send_buf_east, nbeast, &recv_buf_west, nbwest, &stat);

    /* add forces from west to original cells */
    recv_buf_west.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_west, cell_dim.x-2, i, j );
//#endif
    /* copy west forces into send buffer */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_west, cell_dim.x-1, i, j );

    /* send west, receive east */
    sendrecv_buf(&send_buf_west, nbwest, &recv_buf_east, nbeast, &stat);

    /* add forces from east to original cells */
    recv_buf_east.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_east, 1, i, j );
  }
#endif
  /* send forces north/south */
  if (cpu_dim.y==1) {
    /* simply add north/south forces to original cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*add_func)( i, 0, j, i, cell_dim.y-2, j );
        (*add_func)( i, cell_dim.y-1, j, i, 1, j );
      }
  }
#ifdef MPI
  else {
    /* copy north forces into send buffer */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_north, i, 0, j );

    /* send north, receive south */
    sendrecv_buf(&send_buf_north, nbnorth, &recv_buf_south, nbsouth, &stat);

    /* add forces from south to original cells */
    recv_buf_south.n = 0;
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_south, i, cell_dim.y-2, j );

    /* copy south forces into send buffer */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_south, i, cell_dim.y-1, j );

    /* send south, receive north */
    sendrecv_buf(&send_buf_south, nbsouth, &recv_buf_north, nbnorth, &stat);

    /* add forces from north to original cells */
    recv_buf_north.n = 0;
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_north, i, 1, j );
  }
#endif
  /* send forces up/down */
  if (cpu_dim.z==1) {
    /* simply add up/down forces to original cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j) {
        (*add_func)( i, j, 0, i, j, cell_dim.z-2 );
        (*add_func)( i, j, cell_dim.z-1, i, j, 1 );
      }
  }
#ifdef MPI
  else {
    /* copy up forces into send buffer */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_up, i, j, 0 );

    /* send up, receive down */
    sendrecv_buf(&send_buf_up, nbup, &recv_buf_down, nbdown, &stat);

    /* add forces from down to original cells */
    recv_buf_down.n = 0;
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_down, i, j, cell_dim.z-2 );

    /* copy down forces into send buffer */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_down, i, j, cell_dim.z-1 );

    /* send down, receive up */
    sendrecv_buf(&send_buf_down, nbdown, &recv_buf_up, nbup, &stat);

    /* add forces from up to original cells */
    recv_buf_up.n = 0;
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_up, i, j, 1 );
  }
#endif
}

// ************************************************************************
// *      IMPULSE --> ANALOG zu SENDFORCES 
// ***********************************************************************
// ****************************
// * ANALOG ZU add_forces
// *****************************
void nrb_addmomenta( int k, int l, int m, int r, int s, int t )
{
  int i;
  minicell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

  real fpx,fpy,fpz;
  real tpx,tpy,tpz;

  for (i=0; i<to->n; ++i) 
  {
    if(NRBBND(to,i)>0)
    {
      fpx=IMPULS(from,i,X);
      fpy=IMPULS(from,i,Y);
      fpz=IMPULS(from,i,Z);
      	
      IMPULS(to,i,X)  += fpx;
      IMPULS(to,i,Y)  += fpy;
      IMPULS(to,i,Z)  += fpz;
    }
  }
}

/*
void nrb_addmomentav( int k, int l, int m, int r, int s, int t, vektor v )
{
  int i;
  minicell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

  real fpx,fpy,fpz;
  real tpx,tpy,tpz;

  for (i=0; i<to->n; ++i)
  {
    if(NRBBND(to,i)>0)
    {
      fpx=IMPULS(from,i,X);fpy=IMPULS(from,i,Y);fpz=IMPULS(from,i,Z);
        
      IMPULS(to,i,X)  += fpx;
      IMPULS(to,i,Y)  += fpy;
      IMPULS(to,i,Z)  += fpz;
    }
  }
}
*/



void nrb_copymomenta( int k, int l, int m, int r, int s, int t, vektor v )
{
  //Schickt immer von realcell nach buffercell
  int i;
  minicell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

  if (from->n > to->n_max) {
    to->n = 0;
    ALLOC_MINICELL(to, from->n);
  }
    
  to->n = from->n;

  for (i=0; i<to->n; ++i) 
  {
    if(NRBBND(to,i)>0)
    {
      IMPULS(to,i,X)=IMPULS(from,i,X);
      IMPULS(to,i,Y)=IMPULS(from,i,Y);
      IMPULS(to,i,Z)=IMPULS(from,i,Z);
    }
  }
}



void nrb_packmomentav( msgbuf *b, int k, int l, int m, vektor v )
{
  int i, j = b->n;
  minicell *from;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  b->data[ j++ ] = (real) from->n;

  for (i=0; i<from->n; ++i) {
    b->data[j++]=IMPULS(from,i,X);
    b->data[j++]=IMPULS(from,i,Y);
    b->data[j++]=IMPULS(from,i,Z);
  }

  b->n = j;
  if (b->n_max < b->n)
    error("Buffer overflow in nrb_packmomentav - increase msgbuf_size");
}

void nrb_packmomenta( msgbuf *b, int k, int l, int m)
{
  int i, j = b->n;
  minicell *from;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);

  for (i=0; i<from->n; ++i) {

    b->data[j++]=IMPULS(from,i,X);
    b->data[j++]=IMPULS(from,i,Y);
    b->data[j++]=IMPULS(from,i,Z);

//if(NUMMER(from,i)==WATCHME && steps>PRINTSTEP)
//printf("PACKING,myid:%d,steps:%d,index:%d,i:%d, px:%f,py:%f,pz:%f,fromtype:%d,j:%d,n:%d\n",
//         myid,steps,NUMMER(from,i),i,IMPULS(from,i,X),IMPULS(from,i,Y),IMPULS(from,i,Z),from->celltype,j,b->n);


  }

  b->n = j;
  if (b->n_max < b->n)
    error("Buffer overflow in nrb_packmomenta - increase msgbuf_size");
}

void nrb_unpackmomentaadd( msgbuf *b, int k, int l, int m )
{
  int i, j = b->n;
  minicell *to;
  to = PTR_3D_V(cell_array, k, l, m, cell_dim);

  double px,py,pz;

  for (i=0; i<to->n; ++i) {
    px=b->data[j++];
    py=b->data[j++];
    pz=b->data[j++];

//if(NUMMER(to,i)==WATCHME && steps>PRINTSTEP)
//printf("UNPACKING,myid:%d,steps:%d, index:%d,nr:%d, px:%f,py:%f,pz:%f,dpx:%f,dpy:%f,dpz:%f,totype:%d,j:%d,n:%d\n",
//	 myid,steps,NUMMER(to,i),i, IMPULS(to,i,X),IMPULS(to,i,Y),IMPULS(to,i,Z), px,py,pz,to->celltype,j,b->n);

    if(NRBBND(to,i)>0)
    {
      IMPULS(to,i,X)  += px;
      IMPULS(to,i,Y)  += py;
      IMPULS(to,i,Z)  += pz;
    }
  }
  b->n=j;
  if (b->n_max < b->n)
    error("Buffer overflow in nrb_unpackmomentaadd - increase msgbuf_size");

}


void nrb_unpackmomentacopy( msgbuf *b, int k, int l, int m )
{
  int i, j = b->n,tmp_n;;
  minicell *to;
  to = PTR_3D_V(cell_array, k, l, m, cell_dim);
  tmp_n = (int) b->data[ j++ ];

  if (tmp_n > to->n_max) {  //real-cell wird gelöscht???
    to->n = 0;
    ALLOC_MINICELL(to, tmp_n);
  }
  to->n = tmp_n;

  real px,py,pz;

  for (i=0; i<to->n; ++i) {
    px=b->data[j++];
    py=b->data[j++];
    pz=b->data[j++];

//Muss für alle atome geschehen, damit auch impulse der neighatome in bufcells bekannt sind
//    if(NRBBND(to,i)>0)
//    {
      IMPULS(to,i,X)  = px;
      IMPULS(to,i,Y)  = py;
      IMPULS(to,i,Z)  = pz;
//    }
  }
  b->n = j;
  if (b->n_max < b->n)
    error("Buffer overflow in nrb_packmomentacopy - increase msgbuf_size");

}

int sendrecv_buf(msgbuf *send, int to_cpu,
                 msgbuf *recv, int from_cpu, MPI_Status *status)
{
  return MPI_Sendrecv(send->data, send->n, REAL, to_cpu, BUFFER_TAG,
                      recv->data, recv->n_max, REAL, from_cpu, BUFFER_TAG,
                      cpugrid, status );
}

// *********************************************************************************************************************************************
// * Zum Fortsetzen von Simulationen ist ein restart-file mit allen bnd-und neigh-infos nötig
// *
// * Achtung:nur realtypes rausschreiben danach nrb_sendcells 
// *
// * Frage:Wie übernehme ich in der Einleseroutine die eingelesenen per-atom größen?
// * 	   1. für jedes atom in der liste über allcells loopen und  NRBBND(..), NRBNEIGH(...) setzen <--langsam
// *	   2. zusätzlich zellen-index sowie internen atom-index (d.h. "i") alle chkpt_int steps rausschreiben <-- fortgesetzte Sim. muss selbe cpu-dims aufweisen
// *	(*)3. analog zu writechkpt und readchkpt vorgehen, aber auf bnd-atome beschränken. Die Listen nrb_ifromid und nrb_ifromid2 in einem 2ten schritt bauen 
// *		
// * MPI/IO oder proc0 read & send? --> proc0 read & send ausreichend da A~V^(2/3), d.h. selbst bei 1e9 atomen, 
// * weniger als 1e6 oberflächenatome (weil -x-bnd offen). Dazu kommen allerdings nochmals genauso viele neigh-atome (1 atomlage tiefer)
// * --> keine performance/speicher-probleme zu erwarten
// *
// * Format?
// * index x y z bndtype neightype r0 r1 r2 r3 r4 r6 r7 r8 r9 r10 r11 refpos.x refpos.y refpos.z
// *
// * Evtl. zusätzliche optione "refoverwrite" : refpos wird beim einlesen mit aktuelle pos überschrieben
// * Vorteil: bnd-und neigh atome können vor equilib.simul. leicht zugewiesen werden (da perfektes gitter). Nach equilib. können dann die refpos. einfach überschrieben werden
// *
// * vorteil bei angabe von x,y,z (im gegensatz zu proc-nr): fortgesetzte sim. darf andere cpu_dims haben
// *
// * Anmerkung: Hier wurde zwecks Einlesen das normale Einlesen von chkpt-files adaptiert, was recht aufwändig ist.
// *		Unter der Voraussetzung dass es sich um relativ wenige bnd/neigh-atome handelt hätte das ganze viel einfacher gelöst werden können:
// *		Dazu ließt jeder proc das gesamte restart-file ein und bestimmt mittels der x,y,z-Werte ob es sich um "sein" atom handelt.
// *		Trifft dies zu, wird über local_cell_coord die entsprechende MD-zelle identifiziert, die anschließend nach der atom-nummer (index) durchsucht wird.
// *		Bei einem treffer werden danach einfach die eingelesenen werte übernommen
// *
// * Problem:   Dieses Vorgehen kann nur für restarts funktionieren, denn hier entsprechen die x,y,z-koords im nrb-restart-file auch den tatsächlichen coords der atome
// *	 	Sollen die nrb-infos aus einem anderen nrb-file erzwungen werden (das vor der equilib. erstellt wurde) ist dies nicht mehr möglich, da sich die atome bewegt haben
// *		und nun in anderen Zellen oder auch auf anderen procs gelandet sind. 
// *		Hier ist also eine gründlichere Suche nötig...jeder proc braucht die komplette nrb-info und loopt dann für jedes atom in dieser liste über alle seine zellen und atome!
// *		-->deutlich langsamer!
// *		
// * Was nun folgt, sind all die Routinen die für's Einlesen und Ausschreiben nötig sind
// * Update: Das Einlesen erfolgt analog zum Einlesen des chkpt-files, nur dass bei copy_atom_cell_cell und copy_atom_buf_cell keine neuen atome/zellen erzeugt werden.
// *	     Stattdessen wird über x,y,z-position und cell_coord die entsprechende, bereits vorhandenen Zelle nach dem eingelesenen Index durchsucht und die nrb-infos
// *	     werden übernommen
// * TODO:   Fehler-Handling? 
// *		
// ***************************************************************************************************************************************************************
void nrb_writerestart(int fzhlr)
{
  MPI_Barrier(cpugrid);
  outbuf=(char*) calloc(outbuf_size,sizeof(char));

  FILE *out=NULL;
  str255 fname;
  sprintf(fname,"%s.%05d.nrb", outfilename,fzhlr);
 
  if(myid==0)
  {
    out = fopen(fname,"w");
    if (NULL == out) error_str("Cannot open output for nrb %s",fname);
  }
  
  //(*write_atoms_fun)(out);
  // **********************
  write_nrb_fun(out);
  // **********************  

  if ((myid == my_out_id) && (out_grp_size > 1)) 
  {
    MPI_Status status;
    int m=1, len, source;
    while (m < out_grp_size)  //Recv-loop
    {
      MPI_Recv(outbuf, outbuf_size, MPI_CHAR, MPI_ANY_SOURCE,
             MPI_ANY_TAG, cpugrid, &status);

      MPI_Get_count(&status, MPI_CHAR, &len);
      if (status.MPI_TAG==OUTBUF_TAG+1) m++;
      if (len>1) fwrite(outbuf, 1, len-1, out);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (out) fclose(out);
  free(outbuf); 
}

void write_nrb_fun(FILE* out)
{
  int i, k, n, len=0;
  int l;

  int to_cpu;
  ivektor cellc;

//Schreibe bnd-atoms
  for(k=0;k<ncells;k++)
  {
    cell*p=CELLPTR(k);
    for(i=0;i<p->n;i++)
    {
      if(NRBBND(p,i)==0 && NRBNEIGH(p,i)==0) continue;

      cellc=cell_coord(ORT(p,i,X),ORT(p,i,Y),ORT(p,i,Z));
      to_cpu=cpu_coord(cellc);


      len += sprintf(outbuf+len, "%d %d %d %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf",
	   NUMMER(p,i), NRBBND(p,i), NRBNEIGH(p,i),
	   ORT(p,i,X), ORT(p,i,Y), ORT(p,i,Z),
	   NRBI(p,i,0),NRBI(p,i,1),NRBI(p,i,2),NRBI(p,i,3),NRBI(p,i,4),NRBI(p,i,5),
	   NRBI(p,i,6),NRBI(p,i,7),NRBI(p,i,8),NRBI(p,i,9),NRBI(p,i,10),NRBI(p,i,11),
	   REF_POS(p,i,X),REF_POS(p,i,Y),REF_POS(p,i,Z));

      len += sprintf(outbuf+len,"\n");
	  if (len > outbuf_size - 256) flush_outbuf(out,&len,OUTBUF_TAG);
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1); //imd_io.c, zeile 39
}
// ************************************************************************************
void nrb_readrestart(str255 infilename)
{

  MPI_Barrier(cpugrid);

  header_info_t info;
  cell *input;
  FILE *infile=NULL;
  char buf[1024];
  long addnumber = 0;
  int  p=0, k; 
  int  i, n[20], to_cpu, count,s;
  vektor   pos;
  ivektor  cellc;
  real     m, d[30];
  minicell *to;
  msgbuf   *input_buf=NULL, *b;

  int count_nrb=0;
  int count_neigh=0;


   if ((0 == myid) && (0==myrank)) {
//    printf("Reading atoms from %s.\n", infilename);
//    fflush(stdout);
    infile = fopen(infilename,"r");
    if (NULL == infile) error_str("Cannot open input-file for nrb %s",infilename);
  }

  /* size of temporary input buffers */
  if (myid != 0) //my_inp_id) 
  {
    nrb_recv_atoms();  //recv-schleife
    //read_atoms_cleanup();
    return;
  } 

  /* allocate temporary input buffer */
 // if (inp_grp_size > 1) 
  {
    input_buf = (msgbuf *) malloc( num_cpus * sizeof(msgbuf) );
    if (NULL==input_buf) error("cannot allocate input buffers");
    input_buf[0].data  = NULL;
    input_buf[0].n_max = 0;
    input_buf[0].n     = 0;
    for (i=0; i<num_cpus; i++) 
    {
      input_buf[i].data  = NULL;
      input_buf[i].n     = 0;
      input_buf[i].n_max = 0;
      if (i != myid)
        alloc_msgbuf(input_buf+i, inbuf_size);
    }
  }
 
  /* Set up 1 atom input cell */
  //if (inp_grp_size > 1) inbuf_size /= (sizeof(real) * (inp_grp_size-1));
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Cannot allocate input cell.");
  input->n_max=0;
  alloc_cell(input, 1);
  natoms  = 0;
  nactive = 0;

//printf("myid:%d, reading infile,inp_grp_size:%d\n",myid,inp_grp_size);
  while(!feof(infile))
  {
      p = 1;
      if (NULL==fgets(buf,sizeof(buf),infile)) p=0;
      if(p==0) break;

      p = sscanf( buf,
// index, nrbbnd, nrbneigh, x,y,z, r0,r1, .., r11, refx,refy,refz
	"%d %d %d %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf",
        &n[0], &n[1], &n[2],
	&d[0],&d[1],&d[2],
	&n[3],&n[4],&n[5],&n[6],&n[7],&n[8],
	&n[9],&n[10],&n[11],&n[12],&n[13],&n[14],
	&d[4],&d[5],&d[6]);

      if(p>0)
      {
	input->n=1;
	count=0;
	pos.x=d[0]+shiftx_front;
	pos.y=d[1]+shifty_front;
	pos.z=d[2]+shiftz_front;
//ver-realen?
	if(NRBBND(input,0)>0) count_nrb++;
 	if(NRBNEIGH(input,0)==1) count_neigh++;

	NUMMER(input,0)=n[0];
	NRBBND(input,0)=n[1];
	NRBNEIGH(input,0)=n[2];
	//n[2] ist nun die cpu-nr (visualisierung / debug) //<--nicht mehr
        NRBI(input,0,0)=n[3];
        NRBI(input,0,1)=n[4];
        NRBI(input,0,2)=n[5];
        NRBI(input,0,3)=n[6];
        NRBI(input,0,4)=n[7];
        NRBI(input,0,5)=n[8];
        NRBI(input,0,6)=n[9];
        NRBI(input,0,7)=n[10];
        NRBI(input,0,8)=n[11];
        NRBI(input,0,9)=n[12];
        NRBI(input,0,10)=n[13];
        NRBI(input,0,11)=n[14];
	
	REF_POS(input,0,X)=d[4]+shiftx_front;
	REF_POS(input,0,Y)=d[5]+shifty_front;
	REF_POS(input,0,Z)=d[6]+shiftz_front;

        ORT(input,0,X)=pos.x;
	ORT(input,0,Y)=pos.y;
	ORT(input,0,Z)=pos.z;

	cellc = cell_coord(pos.x,pos.y,pos.z);
        to_cpu = cpu_coord(cellc);

//printf("read atom nr:%d, x:%f,y:%f,z:%f, tocpu:%d\n",NUMMER(input,0),pos.x,pos.y,pos.z,to_cpu);

        /* to_cpu is in my input group, but not myself */
        //if ((inp_grp_size > 1) && (myid != to_cpu)) 
        if (myid != to_cpu) 
        {
          b = input_buf + to_cpu;
          if (b->data != NULL) 
	  {
            nrb_copy_atom_cell_buf(b, to_cpu, input, 0);  // message-buffer befüllen
            if (b->n_max - b->n < 22) //atom_size) 
	     {

//printf("sending to cpu:%d\n",to_cpu);

              MPI_Send(b->data, b->n, REAL, to_cpu, INBUF_TAG, cpugrid);
              b->n = 0;
            }
          }
        } 
	else
	{
          /* with per CPU input files, make sure we keep all atoms */
          //if ((to_cpu==myid) || ((1==parallel_input) && (1==addnumber))) 
          if (to_cpu==myid) 
          {
            cellc = local_cell_coord(cellc);
            to = PTR_VV(cell_array,cellc,cell_dim);
//          INSERT_ATOM(to, input, 0); //==move_atom  
	    nrb_move_atom(to,input,0);  // intra-cpu kopiervorgang
	  }
        }
      } // p>0
  }// !feof
  if(infile != NULL)
   fclose(infile);

  //if (inp_grp_size > 1) 
  {
    /* The last buffer is sent with a different tag, which tells the
       target CPU that reading is finished; we increase the size by
       one, so that the buffer is sent even if it is empty */
    for (s=0; s<num_cpus; s++) {
      b = input_buf + s;
      if (b->data) {
        b->data[b->n++] = 0.0;
        MPI_Send(b->data, b->n, REAL, s, INBUF_TAG+1, cpugrid);
        free_msgbuf(b);
      }
    }
    free(input_buf);
  }
  printf("* read %d nrb-atoms and %d neigh-atoms\n",count_nrb,count_neigh);
  //MPI_Barrier(cpugrid); //<-- hier nicht, weil ausser proc 0 sonst keiner hier ankommt (siehe "return" weiter oben)
} //nrb_readrestart
// ******************************************************
void nrb_move_atom(cell *to, cell *from, int index)
{
  /* check the parameters */
  if ((0 > index) || (index >= from->n))
    error("nrb_move_atom: index argument out of range.");
  /* see if we need more space */
  //if (to->n >= to->n_max) alloc_cell(to,to->n_max+incrsz);
  nrb_copy_atom_cell_cell(to, to->n, from, index);
  //++to->n;
  /* delete atom in original cell, by moving the last one to the empty slot */
  --from->n;
  if (index < from->n) nrb_copy_atom_cell_cell(from, index, from, from->n);
} // nrb_move_atom
// *******************************************************
void nrb_copy_atom_cell_cell(cell *to, int i, cell *from, int j)
{
  int l;
  int r;

  for(l=0;l<to->n;l++)
  {
    if(NUMMER(to,l)==from->nummer[j])
    {

      to->isnrbbnd[l]=from->isnrbbnd[j];
      to->isnrbneigh[l]=from->isnrbneigh[j];

      for(r=0;r<12;r++)
      {
        to->nrbid[l*12+r]=from->nrbid[j*12+r];
      }

      if(nrb_overwrite==0) // aus file übernehmen
      {
        to->refpos X(l) = from->refpos X(j);
        to->refpos Y(l) = from->refpos Y(j);
        to->refpos Z (l) = from->refpos Z(j);  
      }
      else //überschreiben
      {
        to->refpos X(l) = to->ort X(l);
        to->refpos Y(l) = to->ort Y(l);
        to->refpos Z (l) = to->ort Z(l);	
      }

    }
  }
} // nrb_copy_atom_cell_cell
// *******************************************************
void nrb_copy_atom_cell_buf(msgbuf *to, int to_cpu, cell *p, int ind )  // cell ---> msgbuf
{
//  if (to->n + atom_size > to->n_max) {
//    error("buffer overflow in nrb_copy_atom_cell_buf");
//  }

  ///data sollte alles real sein
  to->data[ to->n++ ] = to_cpu;
  to->data[ to->n++ ] = NUMMER(p,ind);
  to->data[ to->n++ ] = ORT(p,ind,X); //process_buffer braucht ort wegen cell_coord_local etc.
  to->data[ to->n++ ] = ORT(p,ind,Y);
  to->data[ to->n++ ] = ORT(p,ind,Z);
  to->data[ to->n++ ] = NRBBND(p,ind);
  to->data[ to->n++ ] = NRBNEIGH(p,ind);
  to->data[ to->n++ ] = NRBI(p,ind,0);
  to->data[ to->n++ ] = NRBI(p,ind,1);
  to->data[ to->n++ ] = NRBI(p,ind,2);
  to->data[ to->n++ ] = NRBI(p,ind,3);
  to->data[ to->n++ ] = NRBI(p,ind,4);
  to->data[ to->n++ ] = NRBI(p,ind,5);
  to->data[ to->n++ ] = NRBI(p,ind,6);
  to->data[ to->n++ ] = NRBI(p,ind,7);
  to->data[ to->n++ ] = NRBI(p,ind,8);
  to->data[ to->n++ ] = NRBI(p,ind,9);
  to->data[ to->n++ ] = NRBI(p,ind,10);
  to->data[ to->n++ ] = NRBI(p,ind,11);
  to->data[ to->n++ ] = REF_POS(p,ind,X); 
  to->data[ to->n++ ] = REF_POS(p,ind,Y); 
  to->data[ to->n++ ] = REF_POS(p,ind,Z); 

} // nrb_copy_atom_cell_buf
// ******************************************************
void nrb_recv_atoms(void)
{

 //printf("myid:%d, entered recv-loop\n",myid);

  MPI_Status status;
  int finished = 0, src = 0; //my_inp_id;
  msgbuf b = {NULL,0,0};

  alloc_msgbuf(&b, inbuf_size);
  do  //Recv-Schleife 
  {
    MPI_Recv(b.data, inbuf_size, REAL, src, MPI_ANY_TAG, cpugrid, &status);
//printf("myid:%d, received data, tag:%d\n",myid, status.MPI_TAG);

    MPI_Get_count(&status, REAL, &b.n);
    if (status.MPI_TAG==INBUF_TAG+1) { b.n--; finished=1; } /* last buffer */
    //process_buffer( &b );
    nrb_recv_process_buffer(&b);
  } while (0==finished);
  free_msgbuf(&b);
//printf("myid:%d,fin recv_loop\n",myid);
} // nrb_recv_atoms
// *****************************************
void nrb_recv_process_buffer(msgbuf *b)
{
  int i;
  minicell *to;
  ivektor coord;

  for (i=0; i<b->n; i+=22) //atom_size) 
  {
    if (myid != (int) (b->data[i] + 0.1)) continue; // ???
 
    coord = cell_coord( b->data[i+2], b->data[i+3], b->data[i+4] ); // data[i]= cpunr, data[i+1]=atom-nr
    coord = local_cell_coord( coord );

    to = PTR_VV(cell_array, coord, cell_dim);

//printf("myid:%d, process buffer, x:%f,y:%f,z:%f\n", myid, b->data[i+2], b->data[i+3], b->data[i+4]);
    nrb_copy_atom_buf_cell(to, b, i);
  }
} // nrb_recv_process_buffer
// **********************************************
void nrb_copy_atom_buf_cell(minicell *p, msgbuf *b, int start)  //msgbuf --> cell
{
  int  ind;
  int  j = start + 1;  // the first entry is the CPU number 
  int nr=b->data[j++]; // 2nd entry is atom-index (NUMMER)

  cell *to;

  to = p;
  //if (to->n >= to->n_max) alloc_cell(to,to->n_max+incrsz);
  //ind = to->n++;
  int i;
  j++;j++;j++; //x,y,z,brauche ich hier nicht
  double refx,refy,refz;	

  for(i=0;i<to->n;i++)
  {
    if(NUMMER(to,i)==nr)
    {
      NRBBND(to,i)=b->data[j++];
      NRBNEIGH(to,i)=b->data[j++];
      NRBI(to,i,0)=b->data[j++];
      NRBI(to,i,1)=b->data[j++];
      NRBI(to,i,2)=b->data[j++];
      NRBI(to,i,3)=b->data[j++];
      NRBI(to,i,4)=b->data[j++];
      NRBI(to,i,5)=b->data[j++];
      NRBI(to,i,6)=b->data[j++];
      NRBI(to,i,7)=b->data[j++];
      NRBI(to,i,8)=b->data[j++];
      NRBI(to,i,9)=b->data[j++];
      NRBI(to,i,10)=b->data[j++];
      NRBI(to,i,11)=b->data[j++];

      refx=b->data[j++];
      refy=b->data[j++];
      refz=b->data[j++];

      if(nrb_overwrite==0) //aus file übernehmen
      {
        REF_POS(to,i,X)=refx;
        REF_POS(to,i,Y)=refy;
        REF_POS(to,i,Z)=refz;
      }
      else  //überschreiben
      {
	REF_POS(to,i,X)=ORT(to,i,X);
	REF_POS(to,i,Y)=ORT(to,i,Y);
	REF_POS(to,i,Z)=ORT(to,i,Z);
      }
    }
  }
} // nrb_copy_atom_buf_cell
// **************************************************************
void nrb_readinputfile(str255 fname)
{
  FILE* infile;
  int nrblines=0;
  int i,k,l;

  if(myid==0)
  {
    infile=fopen(fname,"r");
    if(infile==NULL)
    {
        error("Error:infile in nrb_inputfile not found.");
    }
    char buf[300];
    while(fgets(buf,sizeof(buf),infile) != NULL)
     nrblines++;      
  }
  MPI_Barrier(cpugrid);
  MPI_Bcast(&nrblines,1,MPI_INT,0,cpugrid);

//printf("bcast lines %d\n",nrblines);

  int*  buf1d=(int*) malloc(nrblines*15*sizeof(int)); //index bnd neigh r0 r1 ...r11 (refpos wird überschrieben)
  if(myid==0)
  {
    i=0;    
//    fclose(infile);
//    infile=fopen(fname,"r");
    fseek (infile, 0, SEEK_SET);
    double tmpx,tmpy,tmpz,rx,ry,rz;
    for(l=0;l<nrblines;l++)
    {
      int tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15;
      fscanf(infile,"%d %d %d %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf",	     
             &tmp1,&tmp2,&tmp3,
             &tmpx,&tmpy,&tmpz,
	     &tmp4,&tmp5,&tmp6,&tmp7,&tmp8,&tmp9,
	     &tmp10,&tmp11,&tmp12,&tmp13,&tmp14,&tmp15,
             &rx,&ry,&rz);
/*
	     &buf1d[i++],&buf1d[i++],&buf1d[i++],
	     &tmpx,&tmpy,&tmpz,
	     &buf1d[i++],&buf1d[i++],&buf1d[i++],&buf1d[i++],&buf1d[i++],&buf1d[i++],
	     &buf1d[i++],&buf1d[i++],&buf1d[i++],&buf1d[i++],&buf1d[i++],&buf1d[i++],
	     &rx,&ry,&rz);
*/
	     //NOTE : k.P. warum ich diesen Umweg gehen muss....
	     //Der kommentierte block führt zu größtenteils leerem buf1d
	     buf1d[i]=tmp1;
	     buf1d[i+1]=tmp2;
	     buf1d[i+2]=tmp3; //das ist ab jetzt cpu-nr (nur visual/debug-purpose) statt nrbneigh (was sich auch so bestimmten lässt)
	     buf1d[i+3]=tmp4;
	     buf1d[i+4]=tmp5;
	     buf1d[i+5]=tmp6;
	     buf1d[i+6]=tmp7;
	     buf1d[i+7]=tmp8;
	     buf1d[i+8]=tmp9;
	     buf1d[i+9]=tmp10;
	     buf1d[i+10]=tmp11;
	     buf1d[i+11]=tmp12;
	     buf1d[i+12]=tmp13;
	     buf1d[i+13]=tmp14;
	     buf1d[i+14]=tmp15;
	
//printf("l:%d,tmp1:%d,tmp2:%d,tmp3:%d,tmp4:%d,tmp5:%d,tmp6:%d,tmp7:%d,tmp8:%d,tmp9:%d,tmp11:%d,tmp12:%d,tmp13:%d,tmp14:%d,tmp15:%d\n",
//	l,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15);
	i+=15;
    }
    fclose(infile);
  }
  MPI_Bcast(&buf1d[0],nrblines*15,MPI_INT,0,cpugrid);

  //jetzt nrb atome suchen und info übernehmen
  int ifoundbnd=0;
  int ifoundneigh=0;

  for(l=0;l<nrblines*15;l+=15)
  {
    int search=buf1d[l];
//if(myid==0) printf("myid:%d, searching %d,bnd:%d,neigh:%d\n",myid,buf1d[l],buf1d[l+1],buf1d[l+2]);
    for(k=0;k<ncells;k++)
    {
      cell*p=CELLPTR(k);
      for(i=0;i<p->n;i++)
      {
	if(NUMMER(p,i)==search)
	{
	  NRBBND(p,i)=buf1d[l+1];
	  if(NRBBND(p,i)>0) 
	    NRBNEIGH(p,i)=0;
	  else
	    NRBNEIGH(p,i)=1;

	  //NRBNEIGH(p,i)=buf1d[l+2]; //das ist nun cpu-nr (vom output-file)
	  if(NRBBND(p,i)>0) ifoundbnd++;	
 	  if(NRBNEIGH(p,i)>0) ifoundneigh++;
//printf("myid:%d found index:%d,bnd:%d,neigh:%d,bufneigh:%d\n", myid,NUMMER(p,i),NRBBND(p,i),NRBNEIGH(p,i),buf1d[l+2]);
	  NRBI(p,i,0)=buf1d[l+3];
	  NRBI(p,i,1)=buf1d[l+4];
	  NRBI(p,i,2)=buf1d[l+5];
	  NRBI(p,i,3)=buf1d[l+6];
	  NRBI(p,i,4)=buf1d[l+7];
	  NRBI(p,i,5)=buf1d[l+8];
	  NRBI(p,i,6)=buf1d[l+9];
	  NRBI(p,i,7)=buf1d[l+10];
	  NRBI(p,i,8)=buf1d[l+11];
	  NRBI(p,i,9)=buf1d[l+12];
	  NRBI(p,i,10)=buf1d[l+13];
	  NRBI(p,i,11)=buf1d[l+14];
	  REF_POS(p,i,X)=ORT(p,i,X);
	  REF_POS(p,i,Y)=ORT(p,i,Y);
	  REF_POS(p,i,Z)=ORT(p,i,Z);
	}
      }
    }
  }
  printf("* proc %d found %d bnd-atoms and %d neigh-atoms\n",myid,ifoundbnd,ifoundneigh);
  free(buf1d);  
}



