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
// *************************************
// * 1st LOOP TO INIT FILTER VARIABLES
// *************************************
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
/*  
  int c;
  for (c=0; c<ncells2; c++) 
  {

    int i, c1 = cnbrs[c].np;
    cell *p   = cell_array + c1;

    for (i=0; i<p->n; i++) 
    {
      int m;
      for (m=0; m<NNBCELL; m++) 
      {  

        int    c2, jstart, j;
        real   r2;
        cell   *q;

        c2 = cnbrs[c].nq[m];
        if (c2<0) continue;
        if (c2==c1) jstart = i+1;
        else        jstart = 0;
        q = cell_array + c2;
	for (j=jstart; j<q->n; j++) 
	{
         real keepi=KEEPME(p,i);
         real keepj=KEEPME(q,j);

         KEEPME(p,i)+=keepj;
         KEEPME(q,j)+=keepi;

	}
      }      
    }
  }
*/

   int i,k,m,n;
   n=0;
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


  for (i=cellmin.x; i<cellmax.x; ++i)
  {
    for (j=cellmin.y; j<cellmax.y; ++j)
    {
      for (k=cellmin.z; k<cellmax.z; ++k)
      {
	cell*p = PTR_3D_V(cell_array, i, j, k, cell_dim);
	for(l=0;l<p->n;l++)
	{
	  filter_rd=0;

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


            //filter_rd=0; //reset recursion depth
	    fltcnt++;
	    p->filterflag[l]=1;
	    DELME(p,l)=filter_check_neighs(p,l,p,l);

	    if(DELME(p,l)==1)
	    {
	      delcnt++;
	  //    MASSE(p,l)=26.91;
	    }
	  }
	}
      }
    }
  }
  


//  printf("myid:%d,delcnt:%d,fltcnt:%d\n",myid,delcnt,fltcnt);
}
// ****************************************
// RECURISVELY CALLED ROUTINE TO CHECK
// NEIGHBOURS, NEIGHS OF NEIGHS, and so on
// ******************************************
int filter_check_neighs(cell* p,int l,cell* pstart,int lstart)
{

  ivektor coord,lcoord;
  int i,j,k,l2;
  coord  = cell_coord( ORT(p,l,X), ORT(p,l,Y), ORT(p,l,Z) );
  lcoord = local_cell_coord( coord );

  int imin=MAX(0,lcoord.x-1);
  int imax=MIN(lcoord.x+1,cellmax.x);

  int jmin=MAX(0,lcoord.y-1);
  int jmax=MIN(lcoord.y+1,cellmax.y);

  int kmin=MAX(0,lcoord.z-1);
  int kmax=MIN(lcoord.z+1,cellmax.z);

  vektor d1;
  d1.x=ORT(p,l,X);
  d1.y=ORT(p,l,Y);
  d1.z=ORT(p,l,Z);
  
  cell*q;  
  int it=SORTE(p,l);

  real r2,rcut2;
  
  filter_rd++;

  int watchme=1410;

  if(filter_rd>1000) // i.wann ist das sowieso erreicht...Achtung:höchst kritischer parameter
  {
//    if(NUMMER(pstart,lstart)==watchme)
//      printf("\n\n WAAA myid:%d, del %d due to rd\n\n\n",myid,NUMMER(pstart,lstart));   
    return 1;
  }

  for(i=imax;i>=imin;i--) //viel höhere chance auf ein nicht-filter atom in der nachbarkette zu stoßen, bevor
  {		          //max. recursion-depth erreicht ist
    for(j=jmin;j<=jmax;j++)
    {
      for(k=kmin;k<=kmax;k++)
      {
	q = PTR_3D_V(cell_array, i, j, k, cell_dim);
	for(l2=0;l2<q->n;l2++)
	{

	  if(q->filterflag[l2]==1) continue;
          if(NUMMER(q,l2)==NUMMER(pstart,lstart)) continue;
          if(NUMMER(q,l2)==NUMMER(p,l)) continue;       
	  
	  vektor d;
	  d.x=ORT(q,l2,X)-d1.x;
	  d.y=ORT(q,l2,Y)-d1.y;
	  d.z=ORT(q,l2,Z)-d1.z;

	  r2  = SPROD(d,d);
	  int jt=SORTE(q,l2);
	  int col = it * ntypes + jt;
	  if(r2<=pair_pot.end[col]) //NEIGH gefunden
	  {	    
 	    q->filterflag[l2]=1; //bereits besucht

	    /*
            if(DELME(q,l2)==1)
            {
              if(ORT(p,l,X)>70)
                printf("myid:%d, del %d due to neigh:%d\n",myid,NUMMER(p,l),NUMMER(q,l2));
	      
              return 1;
            }
	    */

	    if(FILTERME(q,l2)==0 || KEEPME(q,l2)>1)
	    {
	      //if(ORT(p,l,X)<70)
	      //if(ORT(p,l,X)<filter_min_x)
//	      if(NUMMER(pstart,lstart)==watchme)
//	      printf("\n WEEEE myid:%d, undel %d due to neigh:%d\n\n",myid,watchme,NUMMER(q,l2));
	      KEEPME(p,l)=1;	      
      	      return 0;
	    }

	    if(FILTERME(q,l2)==1)
	    {
 	      int checked=filter_check_neighs(q,l2,pstart,lstart);
	      if(checked==2) 
		return filter_check_neighs(p,l,pstart,lstart); //sackgasse->andere neighs anschauen
	      else
	      {
		//filter_rd--;  //MOD: 21.08.19, aufm hlrs mit mehr procs scheint rek.depth von 1000 immernoch zu viel?!?
		return checked;
	      }
	    }
	  }
	} //for l
      }//for k
    }//for j
  }//for i

  //keine ungetagged-ten neighs mehr gefunden?
  //zurück zum vorherigen atom  und neuen pfad nehmen
//  if(NUMMER(pstart,lstart)==watchme)
//  printf("SACKGASSE für atom %d bei %d...rd:%d\n",NUMMER(pstart,lstart),NUMMER(p,l),filter_rd);
//
  if(NUMMER(pstart,lstart)==NUMMER(p,l)) //Am Anfang der Kette angelangt--> löschen
    return 1;
  else
    return 2; //neuer versuch
  //return 1;

  /*
  for (m=tl[n]; m<tl[n+1]; m++) 
  {

    cell   *q;
    int    c, j;

    c = cl_num[ tb[m] ];
    j = tb[m] - cl_off[c];
    q = cell_array + c;   

    if(FILTERME(q,j) ==0  || KEEPME(q,j)>1)
    {
	return 0;
    }
    if(FILTERME(q,j)==1)
      DELME(q,j)=filter_check_neighs(c,j);
  }
  
  return 1;
  */
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


