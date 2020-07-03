#include "imd.h"
#include "tricub_coeffmat.h"  //A[64][64] matrix fuer tricubic interpolation

// **************************************************
// * AUXILIARY FUNCTIONS FOR BICUBIC INTERPOLATION  *
// * ACHTUNG: das funzt nur für regular grids!
// **************************************************
int read_bc_interp(struct bicubinterp *bi,const char *fname)
{
    // READ TABLE AND SETUP COEFFICIENTS
    int i,j,k,l;
    int is,js; //nr of rhos, and us
    double eng,rho,temp;
    double xmin,xmax;
    double ymin,ymax;
    FILE* myfile;

    double *buf;

if(myid==0)
{
    printf("Reading interpolation table %s\n",fname);
    myfile=fopen(fname,"r");
    if(myfile==NULL)
    {
        error("Error: Input Interpol-table  not found.");
    }

    fscanf(myfile,"%d %d",&is,&js);
//printf("rhos:%d,us:%d\n", rhos,us);
    fscanf(myfile,"%lf %lf %lf %lf",&xmin,&xmax,&ymin,&ymax);
}
    //ALLOCATE BUFFER MEM
    MPI_Bcast(&is,1,MPI_INT,0,cpugrid);
    MPI_Bcast(&js,1,MPI_INT,0,cpugrid);

if(myid==0)
  printf("Memory for table:%f MB\n", (double) is*js*4*8/1e6);

    //allocate buffer 1D array for bcast
    buf= malloc(sizeof(double)*is*js);


if(myid==0)
{
    for(i=0;i<is;i++)
    {
      for(j=0;j<js;j++)
      {
        fscanf(myfile,"%lf %lf %lf",&rho,&eng,&temp);
        //bi->arr[i][j]=temp;
	buf[i*js+j]=temp;
      }
    }
    fclose(myfile);
}
    MPI_Bcast(&xmax,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&xmin,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&ymax,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&ymin,1,MPI_DOUBLE,0,cpugrid);

    bi->is=is;
    bi->js=js;
    bi->xmin=xmin;
    bi->xmax=xmax;
    bi->ymin=ymin;
    bi->ymax=ymax;

    bi->arr= (double**) malloc(bi->is*sizeof(double*));
    for(i=0;i<is;i++){
      bi->arr[i]=(double*) malloc(bi->js*sizeof(double));
    }
    //NOW RECONSTRUCT arr from buf
    MPI_Bcast(buf,is*js,MPI_DOUBLE,0,cpugrid);
    for(i=0;i<is;i++)
      for(j=0;j<js;j++)
	bi->arr[i][j]=buf[i*js+j];
    free(buf);

    bi->dzdx=(double**) malloc(bi->is*sizeof(double*));
    bi->dzdy=(double**) malloc(bi->is*sizeof(double*));
    bi->dzdxdy=(double**) malloc(bi->is*sizeof(double*));
  for(i=0;i<bi->is;i++)
  {
    bi->dzdx[i]=(double*) malloc(bi->js*sizeof(double));
    bi->dzdy[i]=(double*) malloc(bi->js*sizeof(double));
    bi->dzdxdy[i]=(double*) malloc(bi->js*sizeof(double));
  }
  int ic,jc; //icell, jcell
  bi->cmat=(double****) malloc((bi->is-1)*sizeof(double***)); // nr of cells in x-dir=xodes-1
  for(ic=0;ic<bi->is-1;ic++)
  {
    bi->cmat[ic]=(double***) malloc((bi->js-1)*sizeof(double**)); //nr of cells in y-dir=ynodes-1    
    for(jc=0;jc<bi->js-1;jc++)
    {
        bi->cmat[ic][jc]=(double**) malloc(4*sizeof(double*));
        for(k=0;k<4;k++)
        {
            bi->cmat[ic][jc][k]=(double*) malloc(4*sizeof(double));
        }
    }
  }
  //NOW SETUP COEFFS (wird bei scattered data nicht benutzt!)
  bi->dx=(bi->xmax- bi->xmin)/((double) bi->is-1);
  bi->dy=(bi->ymax- bi->ymin)/((double) bi->js-1);

  double dyt=bi->dy;
  double dxt=bi->dx;

  //comp. gradients
  int imin,imax,jmin,jmax;
  for(i=0;i<bi->is;i++)
  {
    for(j=0;j<bi->js;j++)
    {
      imin=MAX(i-1,0);
      imax=MIN(i+1,bi->is-1);
      jmin=MAX(j-1,0);
      jmax=MIN(j+1,bi->js-1);
//fprintf(screen,"me:%d,i:%d,j:%d,imax:%d,jmax:%d,imin:%d,jmin:%d\n",me,i,j,imax,jmax,imin,jmin);

      bi->dzdx[i][j]=(bi->arr[imax][j] - bi->arr[imin][j])/dxt/2;
      bi->dzdy[i][j]=(bi->arr[i][jmax] - bi->arr[i][jmin])/dyt/2;
    }
  }
  //mixed deriv
  for(i=0;i<is;i++)
  {
    for(j=0;j<js;j++)
    {
      jmin=MAX(j-1,0);
      jmax=MIN(j+1,js-1);
      bi->dzdxdy[i][j]=(bi->dzdx[i][jmax] - bi->dzdx[i][jmin])/dyt/2;
    }
  }
  // Now matr. multipl. for each cell (ic,jc)
  double wt[256]=
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        -3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
        2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
        0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
        -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
        9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
        -6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
        2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
        -6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
        4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1};

  //temp/aux vars
  double Cdzdx[4], Cdzdy[4], Cdzdxdy[4], Cz[4];  //order of neigh-points:
                                          //lower left, lower right, upper right, upper left
  for(ic=0;ic<bi->is-1;ic++)
          {
    for(jc=0;jc<bi->js-1;jc++)
    {
        Cz[0]=bi->arr[ic][jc];
        Cz[1]=bi->arr[ic+1][jc];
        Cz[2]=bi->arr[ic+1][jc+1];
        Cz[3]=bi->arr[ic][jc+1];

        Cdzdx[0]=bi->dzdx[ic][jc];
        Cdzdx[1]=bi->dzdx[ic+1][jc];
        Cdzdx[2]=bi->dzdx[ic+1][jc+1];
        Cdzdx[3]=bi->dzdx[ic][jc+1];

        Cdzdy[0]=bi->dzdy[ic][jc];
        Cdzdy[1]=bi->dzdy[ic+1][jc];
        Cdzdy[2]=bi->dzdy[ic+1][jc+1];
        Cdzdy[3]=bi->dzdy[ic][jc+1];
        Cdzdxdy[0]=bi->dzdxdy[ic][jc];
        Cdzdxdy[1]=bi->dzdxdy[ic+1][jc];
        Cdzdxdy[2]=bi->dzdxdy[ic+1][jc+1];
        Cdzdxdy[3]=bi->dzdxdy[ic][jc+1];

        double xtmp[16]={
                          Cz[0],Cz[1],Cz[2],Cz[3],
                          Cdzdx[0]*dxt,Cdzdx[1]*dxt,Cdzdx[2]*dxt,Cdzdx[3]*dxt,
                          Cdzdy[0]*dyt,Cdzdy[1]*dyt,Cdzdy[2]*dyt,Cdzdy[3]*dyt,
                          Cdzdxdy[0]*dxt*dyt,Cdzdxdy[1]*dxt*dyt,Cdzdxdy[2]*dxt*dyt,Cdzdxdy[3]*dxt*dyt
                        };
        double c1[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        //matr.multpl.
        int i2=0;
        double xx;
        for(i2=0;i2<16;i2++)
        {
          xx=0.0;
          for(k=0;k<16;k++)
          {
            xx+=wt[16*(i2)+k]*xtmp[k];
          }
          c1[i2]=xx;
        }
        //unpack coeffs into coeff-matrices
        l=0;
        for(i2=0;i2<4;i2++)
        {
          for(k=0;k<4;k++)
          {
            bi->cmat[ic][jc][i2][k]=c1[l];
            l++;
          }
        }
    }   //j-cell
  }     //i-cell


  //free mem, which is not needed anymore (alles auser coeff-matrix)
  for(i = 0; i < is;i++)
  {
     free(bi->arr[i]);
     free(bi->dzdx[i]);
     free(bi->dzdy[i]);
     free(bi->dzdxdy[i]);
   }
   free(bi->arr);
   free(bi->dzdx);
   free(bi->dzdy);
   free(bi->dzdxdy);

  return 0;
}

void free_bc_interp(struct bicubinterp *bi) //macht das OS automatisch am ende
{
    //free array
    int i,j,k,l;
    int is=bi->is;
    int js=bi->js;

    for(i = 0; i < is;i++)
    {
     free(bi->arr[i]);
     free(bi->dzdx[i]);
     free(bi->dzdy[i]);
     free(bi->dzdxdy[i]);
   }
   free(bi->arr);
   free(bi->dzdx);
   free(bi->dzdy);
   free(bi->dzdxdy);

  for(i=0;i<is-1;i++)
   {
        for(j=0;j<js-1;j++)
        {
          for(k=0;k<4;k++)
            free(bi->cmat[i][j][k]);
          free(bi->cmat[i][j]);
        }
        free(bi->cmat[i]);
   }
   free(bi->cmat);
}

double bcinterp(struct bicubinterp *bi,double xp,double yp)
{
  double xmin=bi->xmin;
  double xmax=bi->xmax;
  double ymin=bi->ymin;
  double ymax=bi->ymax;
  int exc=0;

  if(xp > xmax || xp<xmin || yp > ymax || yp<ymin)
  {
    printf("proc:%d, INTERPOLATION TABLE EXCEEDED:xp:%.4e,yp:%.4e,xmin:%.4e,xmax:%.4e,ymin:%.4e,ymax:%.4e,using closest match.\n",myid,
            xp,yp,xmin,xmax,ymin,ymax);
    return -1;
  }

  int ic,jc; //which cell?
  int is,js;
  is=bi->is;
  js=bi->js;

  //HOTFIX fuer t<tmin
  ic=floor((xp-(xmin+bi->dx/2))/(xmax-xmin)*(is)); //is ist auch zahl der zellen
  jc=floor((yp-(ymin+bi->dy/2))/(ymax-ymin)*(js)); //(0-te zelle zählt auch)

//HOTFIX
if(xp<xmin)
{
  ic=1; //CLOSEMATCH
  xp=xmin+bi->dx;
}
if(yp<ymin)
{
  jc=1;
  yp=ymin+bi->dy;
}
if(yp>xmax)
{
  jc=js-1;
  yp=ymax-bi->dy;
}

  double xl=xmin+((double) ic)*(bi->dx);      //lower, upper etc.
  double xu=xl+bi->dx;
  double yl=ymin+((double) jc)*(bi->dy);
  double yu=yl+bi->dy;

  double t=(xp-xl)/(bi->dx);
  double u=(yp-yl)/(bi->dy);

  double s=0;
  int i;

  for(i=3;i>-1;i--)
  {
    s=t*s+((bi->cmat[ic][jc][i][3]*u+bi->cmat[ic][jc][i][2])*u+bi->cmat[ic][jc][i][1])*u+bi->cmat[ic][jc][i][0];
  }

  if(s!=s)
  {
    printf("proc:%d, s in bcinterp is NaN,xp:%.4e,yp:%.4e,xmax:%.4e,xmin:%.4e,ymax:%.4e,ymin:%.4e\n",
                    myid,xp,yp,xmax,xmin,ymax,ymin);
    //return -1.0;
    error("bcinterp returned NaN");
  }
  return s;
}
// ************************* BILIN INTERPOL ****************************
double lininterp(struct lninterp *bi,double xp,double yp,int option) //option=1 --> special case for TfromE ("trapez"-interpolation)
								     //option=0 --> general case (standard bilinear)
{
  // *************************************************************************
  // * Diese Routine hat 2 Optionen: 0 entspricht  der Standart-Bilinearen Interpolation
  // * Option 1 ist fuer den Fall das x-linear verteilt sind und 
  // * "scatterad data" in y-Richtung (d.h. die Werte sitzen nicht auf einem
  // * "strukturierten Gitter, sondern sind zerstreut)
  // * Die Einlesefunktion hat aus diesem grund auch 2 optionen,da im ersten fall
  // * für y ein 1D-array ausreicht, während im anderen fall y 2d sein muss
  // *
  // * Deswegen sind hier insg. 3 Suchen notwendig, da jedes "x" sein eigenes
  // * Paar von y_upper und y_lower hat (bzgl. dem geforderten yp)
  // * Ebenso müssen y1 und y2 interpoliert werden. Erst dann kann im 3ten
  // * Interpolationsschritt (entlang y-Achse) z(xp,yp) berechnet werden
  // * Im Prinzip kann hier mit wenigen Modifikationen mit höherer Ordnung (z.B. Bikubisch)
  // * interpoliert werden. Dazu sind allerdings mehr Nachbarpunkte nötig und
  // * dementsprechend zusätzliches Suchen in den arrays
  // * Mit wenigen Modifikationen kann dies ebenfalls auf beliebiges "quadrilaterales"
  // * Interpolieren erweitert werden, d.h. x und y sind beide "scattered"
  // * Hierzu muessen insgesamt 6 binäre suchen ausgeführt werden,also
  // * (xlower ->ylowe(xlower) & yupper(xlower) sowie 
  // * xupper -> ylower(xupper) & yupper(xupper)
  // *
  // * Normalerweise wird "scattered data" mit einem ganz anderen Schema 
  // * interpoliert (-->Nearest-Neighbor interpolation), wozu allerdings
  // * Eine voronoi-tesselation nötig wäre (nicht trivial)
  // * Vorerst, soll jedoch mit diesem simplem schema vorlieb genommen werden
  // *
  // **************************************************************************

  double xmin=bi->xmin;
  double xmax=bi->xmax;
  double ymin=bi->ymin;
  double ymax=bi->ymax;

  //bool z3_extr=false; //z3=f(x2,y22). wird aktiviert falls yp<ymin(rho_upper)
  //bool z2_extr=false; //z2=f(x2,y21). extrapolation fuer yp>ymax(rho_upper)
  //UPDATE: ---> Keine lin. extrapol. Funktioniert leider überhaupt nicht gut 

  if(xp>xmax || xp<xmin || yp>ymax || yp<ymin)
  {
    printf("ERROR: Interpolation table exceeded, xp:%.6e, xmin:%.6e,xmax:%.6e, yp:%.6e,ymin:%.6e,ymax:%.6e\n",
            xp,xmin,xmax,yp,ymin,ymax);
    return -1;
  }

  int ic,jc; //closest cell?
  int is,js;

  int i1,i2; //lower,upper
  i1=0;i2=0;
  int j11,j12;
  int j21,j22; //lower,upper 

  is=bi->is;
  js=bi->js;


  int L,U; //lower upper, temporary indices for algo
  int m; //middle index
  L=0;
  U=is-1;
  int it=0;
  ///////////////////////////
  // binary search x-array
  ///////////////////////////
  m=floor(((float)(L+U))/2.0);
  while(L<=U)
  {
//printf("it:%d, m:%d,xp:%.8f, x[m]:%.8f\n",it,m,xp,bi->x[m]);
    if(xp > bi->x[m])
        L=m+1;
    else
        U=m-1;
    it++;
    m=floor(((float)(L+U))/2.0);
  }
  i1=MAX(0,m);
  i2=MIN(is-1,i1+1);

  //printf("xp:%f,xlower:%f,xupper:%f,i1:%d,i2:%d\n",xp,bi->x[i1],bi->x[i2],i1,i2);

  ////////////////////////////////////////////////
  // binary search y-array
  // 2.Searches: for lower and upper rho's
  // Numerical Complexity = log_2(n)
  // TODO: Suchalgo. bekommt "Gedächtnis" verpasst
  //       (errinnert sich ans suchergebnis der letzten Suche)
  //       ---> Hunt-Phase & Bisection-Phase
  ////////////////////////////////////////////////
  L=0;
  U=js-1;
  it=0;
  double *ytmp;

if(option==1)
  ytmp=bi->y[i1];
else
  ytmp=bi->ys;
  //bin.search, leftmost-element
  m=floor(((float)(L+U))/2.0);
  while(L<=U)
  {
    if(yp > ytmp[m])
        L=m+1;
    else
        U=m-1;
    m=floor(((float)(L+U))/2.0);
  }
  j11=m;

  j11=MAX(0,j11);
  j12=MIN(js-1,j11+1);
//j12=0;
//printf("j11:%d,j12:%d\n",j11,j12);

  ///////////////////////
  // Same for i2
  ///////////////////////
if(option==1)
{
  ///////////////////////////////////////////////////////////////////////////
  // es kann passieren dass yp > max(y[i2][..]) bzw. yp < min(y[i2][...])
  // In diesem Fall extrapolieren und y[i2][...] gleich yp setzten 
  ///////////////////////////////////////////////////////////////////////////
  L=0;
  U=js-1;
  it=0;
  if(yp>bi->y[i2][js-1]) //kommt sowieso nicht vor
  { 
    printf("Warning:yp = %.15e > bi->y[i2][js-1]=%.15e, i.e. exceeding y-table for x=%.6f.\n"
	   "xl=%.6f, xu=%.6f\nUsing closest match.\n",yp,bi->y[i2][js-1],xp,bi->x[i1],bi->x[i2]);
    j22=js-1;
    j21=js-2;
  }
  else if(yp<bi->y[i2][0])
  {
    printf("Warning:yp = %.15e  < bi->y[i2][0]=%.15e, i.e. exceeding y-table for x=%.6f , xl=%.6f, xu=%.6f\n",yp,bi->y[i2][0],xp,bi->x[i1],bi->x[i2]);
    j21=0;
    j22=1;
return -1;
  }
  else
  {
    m=floor(((float)(L+U))/2.0);
    while(L<=U)
    {
      if(yp > bi->y[i2][m])
        L=m+1;
      else
        U=m-1;
      it++;
      //printf("it:%d,y[i2][%d]:%.15e,L:%d,U:%d\n",it,m,bi->y[i2][m],L,U);
      m=floor(((float)(L+U))/2.0);
    }
    j21=MAX(0,m);
    j22=MIN(js-1,j21+1);
  }
}

  // ***********************
  // * NOW INTERPOLATE     *
  // ***********************
if(option==1)
{
  //Quadrilateral interpolation
  double z11,z21,z22,z12; //lower-lower, lower-upper,upper-upper,upper-lower neigh. points
  double y11,y21,y22,y12;
  double x1,x2,y1,y2;
  double zxy1,zxy2;
  double zinterp;

  x1=bi->x[i1];
  x2=bi->x[i2];
  z11=bi->z[i1][j11];  //z11=z1
  z21=bi->z[i2][j21];  //z21=z2
  z12=bi->z[i1][j12];  //z12=z4
  z22=bi->z[i2][j22];  //z22=z3;

  y11=bi->y[i1][j11];  //y1  
  y21=bi->y[i2][j21];  //y2
  y12=bi->y[i1][j12];  //y4
  y22=bi->y[i2][j22];  //y3


//EXTRAPOL: Punkt 1,2 wird zu 1,1
//          Punkt 1,1 wird extrapoliert damit y11==yp --> draus wird z11 linear extrapoliert
//          Achtung:funzt garnicht mal so gut --> kann zu negativem z11 führen --> HOTFIX: MAX(z11,0)
//              |\
//              | \
//              |  \
//              |   \
//              \    \
//               \    \
//                \   |
//               * \  |
//                  \ |
//                   \|
if(yp < y11)
{
  double s=(z12-z11)/(y12-y11);
  double y11nu=yp;
  double z11nu=z11-s*(y11-yp);

  y12=y11;
  z12=z11;

  y11=yp;
  z11=MAX(z11nu,0); //negative quadratische temperaturen gibts net!

}
  //biline. interpol. on trapez: Mapping function (x,y)-->(l,m)
  //
  //ansatz: x=a_1+a_2*a*l
  //        y=b_1+b_2*l+b_3*m+b_4*l*m

  double a1,a2,b1,b2,b3,b4;
  double l,m;
  a1=x1;
  a2=-x1+x2;
  b1=y11;
  b2=-y11+y21;
  b3=-y11+y12;
  b4=y11-y21+y22-y12;

  l=(xp-a1)/a2;
  if(b4*l+b3==0)
  {
    error("ERROR: b4*l+b3==0!\n");
  }

  m=(-b2*l-b1+yp)/(b4*l+b3);
  //now weights = A_i/A, wobei A_i=Fläche des, dem Knoten gegen-
  //              uerberliegenden Rechtecks und A=Gesamtfläche(=l*m=1^2)
  double w11,w12,w21,w22;
  w11=(1-l)*(1-m);
  w12=(1-l)*m;
  w22=l*m;
  w21=l*(1-m);
  zinterp=w11*z11+w12*z12+w22*z22+w21*z21;

//zinterp=(xp-x1)/(x2-x1) * z11 + (x2-xp)/(x2-x1) * ( (yp-y12)/(y22-y21)*z12 + (y22-yp)/(y22-y21)*z22);

//  printf("x1:%f,y11:%f,y12:%f,z11:%.15e,w11:%f, z12:%.15e, w12:%f\n"
//         "x2:%f y21:%f,y22:%f,z22:%f,w22:%f, z21:%f, w21:%f\n",
//          x1,y11,y12,z11,w11,z12,w12,
//          x2,y21,y22,z22,w22,z21,w21);
//  printf("xp:%.15f,yp:%.15f,l:%f,m:%f,zinterp:%.4e\n", xp,yp,l,m,zinterp);
  return zinterp;
}
else //standard scheme
{
  double zxy1,zxy2;
  double x1,x2,y1,y2;
  double zinterp;
  //interpol in x-dir
  x1=bi->x[i1]; //lower 
  x2=bi->x[i2]; //upper
  y1=bi->ys[j11]; //j11=j21 for non-scattered y
  y2=bi->ys[j12]; //j22=j21 "..."

  zxy1=(x2-xp)/(x2-x1)*bi->z[i1][j11] + (xp-x1)/(x2-x1)*bi->z[i2][j11];
  zxy2=(x2-xp)/(x2-x1)*bi->z[i1][j12] + (xp-x1)/(x2-x1)*bi->z[i2][j12];

//printf("xp:%.f,x1:%f,x2:%f,yp:%f,y1:%.4e,y2:%.4e,zxy1:%.15e,zxy2:%.15e\n",
//        xp, x1,x2,yp,y1,y2,zxy1,zxy2);


  zinterp=(y2-yp)/(y2-y1)*zxy1+(yp-y1)/(y2-y1)*zxy2;
  return zinterp;
}
}




int read_lin_interp(struct lninterp *bi,const char *fname,int option) //option 0 --> allg., option 1-->trapez (y scatterd)
{
    // READ TABLE AND SETUP COEFFICIENTS
    int i,j,k,l;
    int is,js; //nr of rhos, and us

    double x,y,z;
    double xmin,xmax;
    double ymin,ymax;
    FILE* myfile=NULL;

if(myid==0)
{
    myfile=fopen(fname,"r");
    printf("Reading interpolation table %s\n",fname);
    myfile=fopen(fname,"r");
    if(myfile==NULL)
    {
        error("Error: Input Interpol-table  not found.");
    }
    fscanf(myfile,"%d %d",&is,&js);
    fscanf(myfile,"%lf %lf %lf %lf",&xmin,&xmax,&ymin,&ymax);
}
    MPI_Bcast(&is,1,MPI_INT,0,cpugrid);
    MPI_Bcast(&js,1,MPI_INT,0,cpugrid);
    MPI_Bcast(&xmax,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&xmin,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&ymax,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&ymin,1,MPI_DOUBLE,0,cpugrid);

    bi->is=is;
    bi->js=js;
    bi->xmin=xmin;
    bi->xmax=xmax;
    bi->ymin=ymin;
    bi->ymax=ymax;

//DEBUG
//printf("proc:%d,is:%d,js:%d,xmin:%f,xmax:%f,ymin:%.5e,ymax:%.5e\n",myid,bi->is,bi->js,bi->xmin,bi->xmax,bi->ymin,bi->ymax);

    bi->z= (double**) malloc(is*sizeof(double*));
    bi->x=(double*) malloc(is*sizeof(double));

if(option==1)
    bi->y=(double**) malloc(is*sizeof(double*));
else
    bi->ys=(double*) malloc(js*sizeof(double));

    for(i=0;i<is;i++)
    {    
      bi->z[i]=(double*) malloc(js*sizeof(double));
if(option==1)
      bi->y[i]=(double*) malloc(js*sizeof(double));
    }

    //allocate buffer 1D arrays for bcast of y and z (x can be bcastet directly)
    double *bufz,*bufy;
    bufz= malloc(sizeof(double)*is*js);

if(option==1) //bufy only needed for y-scattered data (=option 1)
    bufy= malloc(sizeof(double)*is*js);

    //read data on proc 0
if(myid==0)
{
    for(i=0;i<is;i++)
    {
      for(j=0;j<js;j++)
      {
        fscanf(myfile,"%lf %lf %lf",&x,&y,&z);
        //bi->z[i][j]=z;
        //bi->y[i][j]=y;
	bufz[i*js+j]=z;

if(option==1)
        bufy[i*js+j]=y;  //special case (TfromE with scattered y-data)
else
	bi->ys[j]=y; 	//general case
      }
      bi->x[i]=x;
    }
    fclose(myfile);

  double reqmem=(is*8+2*is*js*8)/1e6;
  printf("Interpolation table %s read. Total memory:%f MB\n",fname,reqmem);
}
    //NOW RECONSTRUCT arr from buf
    MPI_Bcast(bi->x,is,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(bufz,is*js,MPI_DOUBLE,0,cpugrid);
if(option==1)
    MPI_Bcast(bufy,is*js,MPI_DOUBLE,0,cpugrid);
else
   MPI_Bcast(bi->ys,js,MPI_DOUBLE,0,cpugrid);

    for(i=0;i<is;i++)
    {
      for(j=0;j<js;j++)
      {
        bi->z[i][j]=bufz[i*js+j];
if(option==1)
	bi->y[i][j]=bufy[i*js+j];
      }
    }
    free(bufz);
if(option==1)
    free(bufy);
  
  return 0;
}

int free_lin_interp(struct lninterp *bi) //Macht das OS automatisch
{
    //free array
    int i,j,k,l;
    int is=bi->is;
    int js=bi->js;

    for(i = 0; i < is;i++)
    {
     free(bi->z[i]);
     free(bi->y[i]);
   }
   free(bi->z);
   free(bi->y);

   return 1;
}

// ***********************************************************************
// * ROUTINEN FUER TRIKUBISCHE INTERPOLATION
// *
// * Nach F. Lekien, Int. J. Numer. Meth. Engng 2005; 63:455–471
// * Auch hier nur für regular grids möglich!
// ******************************************************************
double do_tricubinterp(struct tricubinterp *tri,double xp,double yp,double zp)
{
  double xmin=tri->xmin;
  double xmax=tri->xmax;

  double ymin=tri->ymin;
  double ymax=tri->ymax;

  double zmin=tri->zmin;
  double zmax=tri->zmax;

  double dx=tri->dx;
  double dy=tri->dy;
  double dz=tri->dz;


  if(xp>xmax || xp<xmin || yp>ymax || yp<ymin || zp<zmin || zp>zmax)
  {
    printf("ERROR: Interpolation table exceeded,\n xp:%.6e, xmin:%.6e,xmax:%.6e\n yp:%.6e,ymin:%.6e,ymax:%.6e\nzp:%.6e,zmin:%.6e,zmax:%.6e\n",
	    xp,xmin,xmax,yp,ymin,ymax,zp,zmin,zmax);
    return -1;
  }

  int i1,i2; //lower,upper cell
  int j1,j2;
  int k1,k2;
  int is,js,ks;

  is=tri->is;
  js=tri->js;
  ks=tri->ks; 
  double isd=(double) is;
  double jsd=(double) js;
  double ksd=(double) ks;

  i1=floor((xp-(xmin+dx/2))/(xmax-xmin)*isd); 
  j1=floor((yp-(ymin+dy/2))/(ymax-ymin)*jsd); 
  k1=floor((zp-(zmin+dz/2))/(zmax-zmin)*ksd);

  //ACHTUNG: i1,j1,k1 an den rändern checken!
  i1=MAX(i1,0);
  j1=MAX(j1,0);
  k1=MAX(k1,0);
  
  i2=i1+1;
  j2=j1+1;
  k2=k1+1;
  //limit of array-indices
  i2=MIN(i2,is-1);
  i2=MAX(i2,0);

  j2=MIN(j2,js-1);
  j2=MAX(j2,0);

  k2=MIN(k2,ks-1);
  k2=MAX(k2,0);

  int i1min=MAX(i1-1,0);
  int i1max=MIN(i1+1,is-1);
  int i2min=MAX(i2-1,0);
  int i2max=MIN(i2+1,is-1);

  int j1min=MAX(j1-1,0);
  int j1max=MIN(j1+1,js-1);
  int j2min=MAX(j2-1,0);
  int j2max=MIN(j2+1,js-1);
  
  int k1min=MAX(k1-1,0);
  int k1max=MIN(k1+1,ks-1);
  int k2min=MAX(k2-1,0);
  int k2max=MIN(k2+1,ks-1);
  

  //printf("dx:%.2e,dy:%.2e,dz:%.2e\n",dx,dy,dz); 
  //printf("\n\n");
  double x1,x2,y1,y2,z1,z2;
  x1=tri->x[i1];
  x2=tri->x[i2];
  y1=tri->y[j1];
  y2=tri->y[j2];
  z1=tri->z[k1];
  z2=tri->z[k2];

//DEBUG
/*
  printf("i1:%d,j1:%d k1:%d\nx1:%f,y1:%f,z1:%f\nx2:%f,y2:%f,z2:%f\n",
	  i1,j1,k1,
	  tri->x[i1],tri->y[j1],tri->z[k1],
	  tri->x[i2],tri->y[j2],tri->z[k2]);
*/

  //Now comp. derivatives
///  0: x=0; y=0; z=0;
///  1: x=1; y=0; z=0;
///  2: x=0; y=1; z=0;
///  3: x=1; y=1; z=0;
///  4: x=0; y=0; z=1;
///  5: x=1; y=0; z=1;
///  6: x=0; y=1; z=1;
///  7: x=1; y=1; z=1;
  double ***f=tri->f;
  double fs[8]= {f[i1][j1][k1],   // 0 0 0
  		 f[i2][j1][k1],   // 1 0 0
		 f[i1][j2][k1],   // 0 1 0
		 f[i2][j2][k1],   // 1 1 0
		 f[i1][j1][k2],   // 0 0 1
		 f[i2][j1][k2],   // 1 0 1
		 f[i1][j2][k2],   // 0 1 1
		 f[i2][j2][k2]};  // 1 1 1 
 
  double dfdx[8]={ (f[i1max][j1][k1] - f[i1min][j1][k1]), // 0 0 0
		   (f[i2max][j1][k1] - f[i2min][j1][k1]), // 1 0 0
		   (f[i1max][j2][k1] - f[i1min][j2][k1]), // 0 1 0
		   (f[i2max][j2][k1] - f[i2min][j2][k1]), // 1 1 0
		   (f[i1max][j1][k2] - f[i1min][j1][k2]), // 0 0 1
		   (f[i2max][j1][k2] - f[i2min][j1][k2]), // 1 0 1
		   (f[i1max][j2][k2] - f[i1min][j2][k2]), // 0 1 1
		   (f[i2max][j2][k2] - f[i2min][j2][k2])  // 1 1 1
	         };

  double dfdy[8]={ (f[i1][j1max][k1] - f[i1][j1min][k1]), // 0 0 0
                   (f[i2][j1max][k1] - f[i2][j1min][k1]), // 1 0 0
                   (f[i1][j2max][k1] - f[i1][j2min][k1]), // 0 1 0
                   (f[i2][j2max][k1] - f[i2][j2min][k1]), // 1 1 0
                   (f[i1][j1max][k2] - f[i1][j1min][k2]), // 0 0 1
                   (f[i2][j1max][k2] - f[i2][j1min][k2]), // 1 0 1
                   (f[i1][j2max][k2] - f[i1][j2min][k2]), // 0 1 1
                   (f[i2][j2max][k2] - f[i2][j2min][k2])  // 1 1 1
                 };

  double dfdz[8]={ (f[i1][j1][k1max] - f[i1][j1][k1min]), // 0 0 0
                   (f[i2][j1][k1max] - f[i2][j1][k1min]), // 1 0 0
                   (f[i1][j2][k1max] - f[i1][j2][k1min]), // 0 1 0
                   (f[i2][j2][k1max] - f[i2][j2][k1min]), // 1 1 0
                   (f[i1][j1][k2max] - f[i1][j1][k2min]), // 0 0 1
                   (f[i2][j1][k2max] - f[i2][j1][k2min]), // 1 0 1
                   (f[i1][j2][k2max] - f[i1][j2][k2min]), // 0 1 1
                   (f[i2][j2][k2max] - f[i2][j2][k2min])  // 1 1 1
                 };

  double d3fdxdydz[8]={ (((f[i1max][j1max][k1max]-f[i1min][j1max][k1max]) - (f[i1max][j1min][k1max]-f[i1min][j1min][k1max]))  // 0 0 0
                        -((f[i1max][j1max][k1min]-f[i1min][j1max][k1min]) - (f[i1max][j1min][k1min]-f[i1min][j1min][k1min]))),// 0 0 0
//
			(((f[i2max][j1max][k1max]-f[i2min][j1max][k1max]) - (f[i2max][j1min][k1max]-f[i2min][j1min][k1max]))  // 1 0 0
			-((f[i2max][j1max][k1min]-f[i2min][j1max][k1min]) - (f[i2max][j1min][k1min]-f[i2min][j1min][k1min]))),// 1 0 0
//
			(((f[i1max][j2max][k1max]-f[i1min][j2max][k1max]) - (f[i1max][j2min][k1max]-f[i1min][j2min][k1max]))  // 0 1 0
			-((f[i1max][j2max][k1min]-f[i1min][j2max][k1min]) - (f[i1max][j2min][k1min]-f[i1min][j2min][k1min]))),// 0 1 0
//
			(((f[i2max][j2max][k1max]-f[i2min][j2max][k1max]) - (f[i2max][j2min][k1max]-f[i2min][j2min][k1max]))  // 1 1 0
			-((f[i2max][j2max][k1min]-f[i2min][j2max][k1min]) - (f[i2max][j2min][k1min]-f[i2min][j2min][k1min]))),// 1 1 0
//
			(((f[i1max][j1max][k2max]-f[i1min][j1max][k2max]) - (f[i1max][j1min][k2max]-f[i1min][j1min][k2max]))  // 0 0 1
			-((f[i1max][j1max][k2min]-f[i1min][j1max][k2min]) - (f[i1max][j1min][k2min]-f[i1min][j1min][k2min]))),// 0 0 1
//
			(((f[i2max][j1max][k2max]-f[i2min][j1max][k2max]) - (f[i2max][j1min][k2max]-f[i2min][j1min][k2max]))  // 1 0 1
			-((f[i2max][j1max][k2min]-f[i2min][j1max][k2min]) - (f[i2max][j1min][k2min]-f[i2min][j1min][k2min]))),// 1 0 1
//
			(((f[i1max][j2max][k2max]-f[i1min][j2max][k2max]) - (f[i1max][j2min][k2max]-f[i1min][j2min][k2max]))  // 0 1 1
			-((f[i1max][j2max][k2min]-f[i1min][j2max][k2min]) - (f[i1max][j2min][k2min]-f[i1min][j2min][k2min]))),// 0 1 1
//
			(((f[i2max][j2max][k2max]-f[i2min][j2max][k2max]) - (f[i2max][j2min][k2max]-f[i2min][j2min][k2max]))  // 1 1 1
			-((f[i2max][j2max][k2min]-f[i2min][j2max][k2min]) - (f[i2max][j2min][k2min]-f[i2min][j2min][k2min]))) // 1 1 1	 
		      };
  
  double d2fdydz[8]= {  ((f[i1][j1max][k1max]-f[i1][j1min][k1max]) - (f[i1][j1max][k1min]-f[i1][j1min][k1min])), // 0 0 0 
                        ((f[i2][j1max][k1max]-f[i2][j1min][k1max]) - (f[i2][j1max][k1min]-f[i2][j1min][k1min])), // 1 0 0
                        ((f[i1][j2max][k1max]-f[i1][j2min][k1max]) - (f[i1][j2max][k1min]-f[i1][j2min][k1min])), // 0 1 0
			((f[i2][j2max][k1max]-f[i2][j2min][k1max]) - (f[i2][j2max][k1min]-f[i2][j2min][k1min])), // 1 1 0
                        ((f[i1][j1max][k2max]-f[i1][j1min][k2max]) - (f[i1][j1max][k2min]-f[i1][j1min][k2min])), // 0 0 1
                        ((f[i2][j1max][k2max]-f[i2][j1min][k2max]) - (f[i2][j1max][k2min]-f[i2][j1min][k2min])), // 1 0 1
                        ((f[i1][j2max][k2max]-f[i1][j2min][k2max]) - (f[i1][j2max][k2min]-f[i1][j2min][k2min])), // 0 1 1
                        ((f[i2][j2max][k2max]-f[i2][j2min][k2max]) - (f[i2][j2max][k2min]-f[i2][j2min][k2min]))  // 1 1 1
		     };

  double d2fdxdy[8]= {  ((f[i1max][j1max][k1]-f[i1min][j1max][k1]) - (f[i1max][j1min][k1]-f[i1min][j1min][k1])), // 0 0 0
		        ((f[i2max][j1max][k1]-f[i2min][j1max][k1]) - (f[i2max][j1min][k1]-f[i2min][j1min][k1])), // 1 0 0
		        ((f[i1max][j2max][k1]-f[i1min][j2max][k1]) - (f[i1max][j2min][k1]-f[i1min][j2min][k1])), // 0 1 0
			((f[i2max][j2max][k1]-f[i2min][j2max][k1]) - (f[i2max][j2min][k1]-f[i2min][j2min][k1])), // 1 1 0
		        ((f[i1max][j1max][k2]-f[i1min][j1max][k2]) - (f[i1max][j1min][k2]-f[i1min][j1min][k2])), // 0 0 1
			((f[i2max][j1max][k2]-f[i2min][j1max][k2]) - (f[i2max][j1min][k2]-f[i2min][j1min][k2])), // 1 0 1
			((f[i1max][j2max][k2]-f[i1min][j2max][k2]) - (f[i1max][j2min][k2]-f[i1min][j2min][k2])), // 0 1 1
			((f[i2max][j2max][k2]-f[i2min][j2max][k2]) - (f[i2max][j2min][k2]-f[i2min][j2min][k2]))  // 1 1 1		       
		     }; 

  double d2fdxdz[8]= {  ((f[i1max][j1][k1max]-f[i1min][j1][k1max]) - (f[i1max][j1][k1min]-f[i1min][j1][k1min])), // 0 0 0
                        ((f[i2max][j1][k1max]-f[i2min][j1][k1max]) - (f[i2max][j1][k1min]-f[i2min][j1][k1min])), // 1 0 0
                        ((f[i1max][j2][k1max]-f[i1min][j2][k1max]) - (f[i1max][j2][k1min]-f[i1min][j2][k1min])), // 0 1 0
			((f[i2max][j2][k1max]-f[i2min][j2][k1max]) - (f[i2max][j2][k1min]-f[i2min][j2][k1min])), // 1 1 0
                        ((f[i1max][j1][k2max]-f[i1min][j1][k2max]) - (f[i1max][j1][k2min]-f[i1min][j1][k2min])), // 0 0 1
                        ((f[i2max][j1][k2max]-f[i2min][j1][k2max]) - (f[i2max][j1][k2min]-f[i2min][j1][k2min])), // 1 0 1
                        ((f[i1max][j2][k2max]-f[i1min][j2][k2max]) - (f[i1max][j2][k2min]-f[i1min][j2][k2min])), // 0 1 1
                        ((f[i2max][j2][k2max]-f[i2min][j2][k2max]) - (f[i2max][j2][k2min]-f[i2min][j2][k2min]))  // 1 1 1
		     };
  // Die Re-Skalierung kann ich mir sparen, wenn ich bei der Berechnung der gradienten 
  // dx,dy,dz weglasse
 
  double a[64];
  tricub_get_coeff(a,fs,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz); 
  double f2;
  double xs,ys,zs; //scaled coords.
  if(x2-x1==0)
    xs=1.0;
  else
    xs=(xp-x1)/(x2-x1);
  if(y2-y1==0)
    ys=1.0;
  else
    ys=(yp-y1)/(y2-y1);
  if(z2-z1==0)
    zs=1.0;
  else
    zs=(zp-z1)/(z2-z1);
  f2=tricub_eval(a,xs,ys,zs);

//printf("xs:%.2e, ys:%.2e, zs:%.2e\n",xs,ys,zs);
  return f2;
}

int read_tricub_interp(struct tricubinterp *tri,const char *fname)
{
    // READ TABLE AND SETUP COEFFICIENTS
    int i,j,k,l;

    FILE* myfile;
if(myid==0)
{
    myfile=fopen(fname,"r");
    printf("Reading %s...\n",fname);
    if(myfile==NULL)
    {
        //error
        error("Error:table not found.");
    }
}

    int is,js,ks;; //nr of rhos, and us

    double x,y,z,f;
    double xmin,xmax;
    double ymin,ymax;
    double zmin,zmax;
if(myid==0)
{
    fscanf(myfile,"%d %d %d",&is,&js,&ks);
}
    MPI_Bcast(&is,1,MPI_INT,0,cpugrid);
    MPI_Bcast(&js,1,MPI_INT,0,cpugrid);
    MPI_Bcast(&ks,1,MPI_INT,0,cpugrid);
//printf("rhos:%d,us:%d\n", rhos,us);
    tri->is=is;
    tri->js=js;
    tri->ks=ks;
if(myid==0)
    fscanf(myfile,"%lf %lf %lf %lf %lf %lf",&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);

//printf("rhomin:%.4e,rhomax:%.4e,Umin:%.4e,Umax:%.4e\n",
//        rhomin,rhomax,Ushiftlogmin,Ushiftlogmax);
    MPI_Bcast(&xmax,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&xmin,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&ymax,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&ymin,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&zmax,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&zmin,1,MPI_DOUBLE,0,cpugrid);

    tri->xmin=xmin;
    tri->xmax=xmax;
    tri->ymin=ymin;
    tri->ymax=ymax;
    tri->zmin=zmin;
    tri->zmax=zmax;

    tri->dx=(xmax-xmin)/((double) is);
    tri->dy=(ymax-ymin)/((double) js);
    tri->dz=(zmax-zmin)/((double) ks);

    tri->z=(double*) malloc(ks*sizeof(double));
    tri->x=(double*) malloc(is*sizeof(double));
    tri->y=(double*) malloc(js*sizeof(double));
    tri->f=(double***) malloc(is*sizeof(double**));

    //alloc buffer mem for Bcast
    double* buf= malloc(sizeof(double)*is*js*ks);

    for(i=0;i<is;i++){
      tri->f[i]=(double**) malloc(js*sizeof(double*));
      for(j=0;j<js;j++)
        tri->f[i][j]=(double*) malloc(ks*sizeof(double));
    }

if(myid==0)
{
    for(i=0;i<is;i++)
    {
      for(j=0;j<js;j++)
      {
	for(k=0;k<ks;k++)
	{
          fscanf(myfile,"%lf %lf %lf %lf",&x,&y,&z,&f);
          tri->z[k]=z;
	  //tri->f[i][j][k]=f; 
	  buf[i + j * is + k * is * js]=f;	
        }
        tri->y[j]=y;
      }
      tri->x[i]=x;
    }
    fclose(myfile);
}

  MPI_Bcast(tri->x,is,MPI_DOUBLE,0,cpugrid);
  MPI_Bcast(tri->y,js,MPI_DOUBLE,0,cpugrid);
  MPI_Bcast(tri->z,ks,MPI_DOUBLE,0,cpugrid);

  double reqmem=(is*js*ks*8)/1e6;
if(myid==0)
  printf("Tricub-Interpolation table %s read. Total memory:%f MB\n",fname,reqmem);
    
  //Bcast buff-array and reconstruct 3d-arrays on every proc
  MPI_Bcast(buf,is*js*ks,MPI_DOUBLE,0,cpugrid);
/*
  for(l=0;k<is*js*ks;l++)
  {
    i= l % is;
    j=(l/ is) % js;
    k=l/(is*js);
printf("i:%d,j:%d,k:%d\n",i,j,k);
    tri->f[i][j][k]=buf[l];  
  }
*/

    for(i=0;i<is;i++)
    {
      for(j=0;j<js;j++)
      {
        for(k=0;k<ks;k++)
        {
	  tri->f[i][j][k]=buf[i + j * is + k * is * js];
	}
      }
    }

  free(buf);

  return 0;
}


void tricub_get_coeff(double a[64], double f[8], double dfdx[8], double dfdy[8], double dfdz[8],
                        double d2fdxdy[8], double d2fdxdz[8], double d2fdydz[8], double d3fdxdydz[8])
{
  int i;
  double x[64];
  #pragma ivdep
  for (i=0;i<8;i++) {
    x[0+i]=f[i];
    x[8+i]=dfdx[i];
    x[16+i]=dfdy[i];
    x[24+i]=dfdz[i];
    x[32+i]=d2fdxdy[i];
    x[40+i]=d2fdxdz[i];
    x[48+i]=d2fdydz[i];
    x[56+i]=d3fdxdydz[i];
  }
  tricub_get_coeff_stacked(a,x);
}


double tricub_eval(double a[64], double x, double y, double z)
{
  int i,j,k;
  int t1,t2,t3,t4;
  int t5,t6,t7;
  double s1=0;
  double s2=0;
  double s3=0;
  double s4=0;

  double ret=0.0;
  for (i=0;i<4;i++) 
  {    
    for (j=0;j<4;j++) 
    {
      t1=4*j;
      #pragma ivdep
      for (k=0;k<4;k++) 
      {
        //ret+=a[tricub_ijk2n(i,j,k)]*pow(x,i)*pow(y,j)*pow(z,k);
        //ret+=a[i+4*j+16*k]*pow(x,i)*pow(y,j)*pow(z,k);
        t2=16*k;
        t3=i+t1+t2;
        s1=a[t3];

        t4=pow(x,i);
        t5=pow(y,j);
        t6=pow(z,k);
        t7=t4*t5*t6;
        s2=t4*t5*t6;

        s3=s1*s2;
        s4=s4+s3;
      }
    }
  }
  ret=ret+s4;
  return(ret);
}



void tricub_get_coeff_stacked(double a[64], double x[64]) 
{
  int i,j;
  for(i=0;i<64;i++)
    a[i]=0.0;

  for (i=0;i<64;i++) {
    //a[i]=(double)(0.0);
    #pragma ivdep
    for (j=0;j<64;j++) 
    {
      a[i]+=tricub_A[i][j]*x[j];
    }
  }
}

int tricub_ijk2n(int i, int j, int k) 
{
  return(i+4*j+16*k);
}



// ******************  ROUTINEN FUER NATURAL NEIGHBOR INTERPOLATION ************************
int nn_read_table(nn_interp *nn, const char *fname)
{
    // READ TABLE AND SETUP COEFFICIENTS
    int i,j,k,l;
    int is,js; //nr of rhos, and us
    double x,y,z;
    double xmin,xmax;
    double ymin,ymax;
    FILE* myfile;

    double *xbuf,*ybuf,*zbuf;
if(myid==0)
{
    printf("Reading interpolation table %s\n",fname);
    myfile=fopen(fname,"r");
    if(myfile==NULL)
    {
        error("Error: Input Interpol-table for nn not found.");
	      MPI_Abort(cpugrid,0);
    }

    fscanf(myfile,"%d %d",&is,&js);
    fscanf(myfile,"%lf %lf %lf %lf",&xmin,&xmax,&ymin,&ymax);
}
    MPI_Bcast(&is,1,MPI_INT,0,cpugrid);
    MPI_Bcast(&js,1,MPI_INT,0,cpugrid);
    MPI_Bcast(&xmax,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&xmin,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&ymax,1,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(&ymin,1,MPI_DOUBLE,0,cpugrid);

    nn->xmin=xmin;
    nn->xmax=xmax;
    nn->ymin=ymin;
    nn->ymax=ymax;
    nn->npoints=is*js;
/*
    xbuf= malloc(sizeof(double)*nn->npoints);
    ybuf= malloc(sizeof(double)*nn->npoints);
    zbuf= malloc(sizeof(double)*nn->npoints);
*/
    alloc1darr(double,xbuf,nn->npoints);
    alloc1darr(double,ybuf,nn->npoints);
    alloc1darr(double,zbuf,nn->npoints);

  if(myid==0)
  {
    for(i=0;i<is;i++)
    {
      for(j=0;j<js;j++)
      {
        fscanf(myfile,"%lf %lf %lf",&x,&y,&z);
        //bi->arr[i][j]=temp;
        xbuf[i*js+j]=x;
        ybuf[i*js+j]=y;
        zbuf[i*js+j]=z;
      }
    }
    fclose(myfile);
  }

    nn->points = (point*) malloc(nn->npoints * sizeof(point));
    //NOW RECONSTRUCT arr from buf
    MPI_Bcast(xbuf,is*js,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(ybuf,is*js,MPI_DOUBLE,0,cpugrid);
    MPI_Bcast(zbuf,is*js,MPI_DOUBLE,0,cpugrid);

    point* p;
    for(i=0;i<is;i++)
    {
      for(j=0;j<js;j++)
      {
        p = &(nn->points)[i*js+j];

        p->x=xbuf[i*js+j];
        p->y=ybuf[i*js+j];
        p->z=zbuf[i*js+j];
      }
    }

/*    
    free(xbuf);
    free(ybuf);
    free(zbuf);
*/
    free1darr(xbuf);
    free1darr(ybuf);
    free1darr(zbuf);
    //Build delaunay triangulation
    nn->d = delaunay_build(nn->npoints,nn->points, 0, NULL, 0, NULL);
    //build hash interpolator
    //nn->nn = nnhpi_create(nn->d, 1);
    //nn->interpolator=nnhpi_create(nn->d,1); //naturla neigh, sibson-rule
    nn->interpolator=lpi_build(nn->d); // linear

 return 0;
}

