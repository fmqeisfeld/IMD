#include "imd.h"
//TODO: nachbarsuche mittels geometrischer loop über nachbarzellen für +y/-y-filteratome so wie für x-filter anpassen,
//	d.h. prüfen ob atom näher an +y oder -y grenze 
//	bei +y: for(j=jmin;j<=jmax;j++)
//	bei -y: for(j=jmax;j>=jmin;j--)
//
//	Weiterhin Problematisch: atome werden im imd_fix_cells_3d.c gelöscht. 
//				 Das Löschen geschieht nur wenn ein Atom in eine andere Zelle übergeht.
//				 Wenn z.B. die -x Filter-Grenze bei x=0 liegt, befindet sie sich auf dem 0-ten Proc (in 1D) in der 0-ten Zelle
//				 Atome, die erst hinter der Filter-Grenze (z.b. bei -100) zum Löschen freigegeben werden,
//				 werden von fix_cells ignoriert und nicht gelöscht!  
//
// UPDATE 22.08.19:  Der rekursive algorithmus wurde nun durch einen iterativen algorithmus mit explizitem "stack" (eigentlich heap) ersetzt, 
//		     da kaum vorherzusagen ist, wieviel platz auf dem stack noch frei ist
//		     Es wird zunächst ein pseudo-stack auf dem heap allokiert. 
//		     Für jedes Atom, dass geprüft wird, werden anschließend alle nachbaratome auf diesen "stack" geschoben und nacheinander abgearbeitet.
//		     Läuft der algorithmus in eine Sackgasse, wird das zuletzt auf den "stack" geschobene atom entfernt und der algorithmus springt zum vorherigen nachbarn etc.
//		     
//		     Der algo wurde so gestaltet, dass die komplette Nachbarkette, bis zum Abbruchkriterium im Speicher bleibt. 
//		     Alternativ könnte man auch das zuletzs auf den stack geschobene Atom mit dem jeweils neu gefundenen Nachbarn ersetzen, womit man 
//	             sehr viel speicher sparen würde. 
//		     Der Vorteil der ersten Vorgehensweise jedoch ist, dass, nachdem in der Nachbarkette ein nicht-Filter-Atom gefunden wurde, die
//		     komplette Nachbarkette als "nicht löschen" getaggt werden kann.
//		     Dies beschleunigt die folgenden Iterationen erheblich.
//		    
//		     Der Parameter NEIGHDEPTH gibt an, wie groß der Zwischenspeicher auf dem heap sein soll.
//		     Wird nach NEIGHDEPTH iterationen kein Nicht-Löschen-Kriterium gefunden, wird das Atom gelöscht. 
//		     Je größer NEIGHDEPTH umso langsamer der Algorithmus 		    
// *************************************
// * 1st LOOP TO INIT FILTER VARIABLES
// *************************************
//
#define NEIGHDEPTH 1000

void filter_init(void)
{
   int i,k;
   for(k=0;k<nallcells;k++)
   {
     cell*p=cell_array+k; //CELLPTR(k);
     for(i=0;i<p->n;i++)
     {
       p->filterflag[i]=0;

       if(ORT(p,i,X)<filter_min_x)
       {
	 FILTERME(p,i)=1;
	 KEEPME(p,i)=0;
	 continue;      
       }
       if(ORT(p,i,Y)<filter_min_y)
       {
	 FILTERME(p,i)=1;
	 KEEPME(p,i)=0;
	 continue;
       } 
       if(ORT(p,i,Y)>filter_max_y)
       {
	 FILTERME(p,i)=1;
	 KEEPME(p,i)=0;
	 continue;
       }
       //sonst
       FILTERME(p,i)=0;
       KEEPME(p,i)=1;            
     }   
   }

   //NOW FILL BUFFER-CELLS AND EXCHANGE KEEPME-VARIABLE
   //AMONG NEIGHBOURING CPUS
   //send_cells,....   
}
// *****************************
// * ACCUMULATE KEEPME-variable
// ****************************
void filter_forces(void)
{
   int i,k,m,n;
   n=0;

#ifdef NBL
   for(k=0;k<ncells;k++)
   {
     cell *p = cell_array + cnbrs[k].np;
     for (i=0; i<p->n; i++)
     {  
       for (m=tl[n]; m<tl[n+1]; m++) 
       {

         cell   *q;
         int    c, j;

         c = cl_num[ tb[m] ];
         j = tb[m] - cl_off[c];
         q = cell_array + c;

         real keepi=KEEPME(p,i);
         real keepj=KEEPME(q,j);

         KEEPME(p,i)+=keepj;
         KEEPME(q,j)+=keepi;
       }
       n++;
     }
   }
#else

#ifdef AR 
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
         jstart = (((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0)) ? i+1 : 0);
         for (j = jstart; j < q->n; ++j)
         {
           real keepi=KEEPME(p,i);
           real keepj=KEEPME(q,j);

           KEEPME(p,i)+=keepj;
           KEEPME(q,j)+=keepi;
         }
      }
    }
  }
#else
error("FILTER only works with NBL or PAIR LISTS WITH AR");
#endif

#endif


}
// **********************************
// * MAIN FILTER ROUTINE TO DETERMINE
// * IF ATOMS IS TO BE DELETED OR NOT
// **********************************
void filter_atoms(void)
{

  MPI_Barrier(cpugrid);

  //printf("myid:%d,steps:%d,entered filter_atoms,%d\n",myid,steps,have_valid_nbl);

  int i,j,k,l;
  filter_init();
  filter_forces();
  send_forces(add_filter,pack_filter,unpack_filter);
  //send_cells(...) buffer-zellen brauchen die infos auch!
  int delcnt=0;
  int fltcnt=0;

  filter_stacki=(int*)  malloc((NEIGHDEPTH)*sizeof(int));
  filter_stackq=(cell**) malloc((NEIGHDEPTH)*sizeof(cell*));
  for(k=0;k<ncells;k++)
  {
    //cell*p = PTR_3D_V(cell_array, i, j, k, cell_dim);
    cell *p=CELLPTR(k);
    for(l=0;l<p->n;l++)
    {
      if(FILTERME(p,l)==1 && KEEPME(p,l)==0)
      {

	//reset filterflag
	int c,r;
	for(c=0;c<nallcells;c++)
	{
	cell *foo=cell_array+c;
	for(r=0;r<foo->n;r++)
	  foo->filterflag[r]=0;
	}


	fltcnt++;
	p->filterflag[l]=1;
	DELME(p,l)=filter_check_neighs(p,l);

	if(DELME(p,l)==1)
	{
	  delcnt++;
//DEBUG zu Visualisierungszwecken
//MASSE(p,l)=26.91;
//DELME(p,l)=0;
	}
      }
    }
  }
  
  if(filter_stacki!= NULL) free(filter_stacki);
  if(filter_stackq!=NULL) free(filter_stackq);
  filter_stacki=NULL;
  filter_stackq=NULL;


//  printf("myid:%d,delcnt:%d,fltcnt:%d\n",myid,delcnt,fltcnt);
}
// ****************************************
// RECURISVELY CALLED ROUTINE TO CHECK
// NEIGHBOURS, NEIGHS OF NEIGHS, and so on
// ******************************************
int filter_check_neighs(cell* p,int l)
{
  int filter_stackindex=0;
  int filter_stackindex_before=0; //zum checken ob algo in dead-end läuft

  ivektor coord,lcoord;
  int i,j,k;
 
  //ausgangsatom in den stack schieben 
  filter_stacki[0]=l;
  filter_stackq[0]=p;
  filter_stackindex=0;

//printf("myid:%d, entered filter_check_neighs\n",myid);

  while(filter_stackindex >=0 ) //Stack abarbeiten
  {

//printf("myid:%d,index:%d\n",myid,filter_stackindex) ;

    filter_stackindex_before=filter_stackindex;   //Fuer dead-end check (nach neighloop)

    int   l2=filter_stacki[filter_stackindex]; //latest i
    cell*  q=filter_stackq[filter_stackindex];  //latest q

    q->filterflag[l2]=1; //als bereits besucht markieren
   
    if(KEEPME(q,l2)>0 || FILTERME(q,l2)==0) //Non-filter-atom found
    {
      KEEPME(p,l)=1;
      //Komplette Nachbarkette für weitere suchen taggen -> performance
      int ltmp=0;
      cell* qtmp;
      int loopvar;
      for(loopvar=0;loopvar<=filter_stackindex;loopvar++)
      {
	qtmp=filter_stackq[loopvar];
	ltmp=filter_stacki[loopvar];
	KEEPME(qtmp,ltmp)=1;
      }
      return 0; // nicht löschen
    }

    coord  = cell_coord( ORT(q,l2,X), ORT(q,l2,Y), ORT(q,l2,Z) );
    lcoord = local_cell_coord( coord );

    int imin=MAX(0,lcoord.x-1);
    int imax=MIN(lcoord.x+1,cellmax.x);

    int jmin=MAX(0,lcoord.y-1);
    int jmax=MIN(lcoord.y+1,cellmax.y);

    int kmin=MAX(0,lcoord.z-1);
    int kmax=MIN(lcoord.z+1,cellmax.z);

    vektor d1;
    d1.x=ORT(q,l2,X);
    d1.y=ORT(q,l2,Y);
    d1.z=ORT(q,l2,Z);
    
    int it=SORTE(q,l2);
    real r2,rcut2;

    //loop neighs and add to list
    for(i=imin;i<=imax;i++) 
    {		          
      for(j=jmin;j<=jmax;j++)
      {
	for(k=kmin;k<=kmax;k++)
	{
	  cell* q2 = PTR_3D_V(cell_array, i, j, k, cell_dim);
	  int l3;
	  for(l3=0;l3<q2->n;l3++)
	  {

	    if(q2->filterflag[l3]==1) continue; //bereits besuchte ignorieren
	    
	    vektor d;
	    d.x=ORT(q2,l3,X)-d1.x;
	    d.y=ORT(q2,l3,Y)-d1.y;
	    d.z=ORT(q2,l3,Z)-d1.z;

	    r2  = SPROD(d,d);
	    int jt=SORTE(q2,l3);
	    int col = it * ntypes + jt;
//HOTFIX    
col=0;
	    if(r2<=pair_pot.end[col]) //NEIGH gefunden
	    {	    
	      filter_stackindex++; //atom an nächste position in der liste schieben
	      if(filter_stackindex >NEIGHDEPTH-1) //stack voll
	        return 1;  // <--nach NEIGHDEPTH iterationen wurde kein atom in der nachbarkette gefunden, dass 
			   //    nicht gelöscht werden soll --> atom zum löschen freigeben


	      filter_stacki[filter_stackindex]=l3;
	      filter_stackq[filter_stackindex]=q2; // zur liste hinzufuegen
	    }
	  } //for l
	}//for k
      }//for j
    }//for i

    //Auf Sackgasse prüfen und evtl. zum vorherigen Atom zurückspringen
    if(filter_stackindex==filter_stackindex_before) //d.h. keine neuen neighs dazu gekommen
      filter_stackindex--; //Sackgasse --> 1 knoten zurück

  }//while
  return 1;
}
// **************************************
// * COMM
// **************************************
//
void add_filter( int k, int l, int m, int r, int s, int t )
{
  int i;
  minicell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

  for (i=0; i<to->n; ++i) 
  {
    KEEPME(to,i)  += KEEPME(from,i);
  }

}

void pack_filter( msgbuf *b, int k, int l, int m)
{
  int i, j = b->n;
  minicell *from;
  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  for (i=0; i<from->n; ++i) {
    b->data[ j++ ] = KEEPME(from,i); 
  }
  b->n = j;
  if (b->n_max < b->n)
    error("Buffer overflow in pack_filter - increase msgbuf_size");
}


void unpack_filter( msgbuf *b, int k, int l, int m )
 {
  int i, j = b->n;
  minicell *to;

  to = PTR_3D_V(cell_array, k, l, m, cell_dim);
  for (i=0; i<to->n; ++i) {
    KEEPME(to,i)  += b->data[ j++ ];
  }
  b->n = j;
  if (b->n_max < b->n)
    error("Buffer overflow in unpack_filter - increase msgbuf_size");
}


