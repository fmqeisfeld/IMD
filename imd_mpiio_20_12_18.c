/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2012 Institute for Theoretical and Applied Physics
* University of Stuttgart, D-70550 Stuttgart
*
* Copyright 2013-2018 Institute for Functional Matter and Quantum Technologies
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_mpiio.c -- dimension independent IO routines 
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

#ifdef MPIIO
/******************************************************************************
*
*  new function for parallel MPI output
*
*  filter function for write_config_select
*  writes data of all atoms for checkpoints
*
******************************************************************************/

/******************************************************************************
*
* Betroffene Dateien: imd_param.c imd_mpi_util.c imd_io.c imd_forces_covalent.c
*                     imd_main_3d.c 
*
* Parameters:
* 
* parallel_output 4 # parallel output flag (0 == serial (default),
*                   # 1 == parallel, 2 == serial, but parallel picture writes)
*                   # 3 == mpiio, impliziert binary_output  1 
*                   # 4 == mpiio, wie 3 == mpiio, 
*                   # aber Anfangskonfiguration wird nicht geschrieben 
* 
******************************************************************************/
#define DATAREP "native"
#define OBSERVABLES 7

void write_atoms_config_mpiio(void)
{
  if(myid==0)
    printf("writing atoms using mpiio\n");
  int i, k, n;
  long mpi_atoms=0;
  char *atomsbuffer;
  double mpidata[OBSERVABLES];
  long startidx=0;
  int observables=OBSERVABLES;
  int byte=8;

  char *mpioutfilename, *filename;
  int atombuffersize, atombufferpos=0;
  MPI_Pack_size (OBSERVABLES, MPI_DOUBLE, MPI_COMM_WORLD, &atombuffersize);
  MPI_File fh;
  MPI_Offset disp=1024, myoffset;
  MPI_Status iostatus;

  long atoms_total;

  /* Outputfilename */

  if ((filename = malloc(256))==NULL) {
    error("Cannot malloc name buffer\n");
  };
  sprintf(filename,"%s.%05d.%s",outfilename,steps/checkpt_int,"mpiio");
  mpioutfilename=filename;

  /* Zählen der Atome pro CPU */
  for (k=0; k<NCELLS; k++) {
    cell *p;
    p = CELLPTR(k);
    mpi_atoms += p->n;
  }

  MPI_Allreduce(&mpi_atoms,&atoms_total,1,MPI_LONG,MPI_SUM,cpugrid);
  atomsbuffer = (char *) malloc(mpi_atoms * atombuffersize);
  
  for (k=0; k<NCELLS; k++) {
    cell *p;
    p = CELLPTR(k);
    for (i=0; i<p->n; i++) {

      /* Speichern der Atome in MPI Buffer */
        n=0;
        mpidata[n++] = MASSE(p,i);
        mpidata[n++] = ORT(p,i,X);
        mpidata[n++] = ORT(p,i,Y);
        mpidata[n++] = ORT(p,i,Z);
	mpidata[n++] = (IMPULS(p,i,X) / MASSE(p,i));
	mpidata[n++] = (IMPULS(p,i,Y) / MASSE(p,i));
	mpidata[n++] = (IMPULS(p,i,Z) / MASSE(p,i));

	/* Packen des MPI Buffers */
	MPI_Pack (mpidata, OBSERVABLES, MPI_DOUBLE, atomsbuffer, 
                  observables*byte*mpi_atoms, &atombufferpos, MPI_COMM_WORLD);
    }
  }

  /* Bestimmen der Position für alle CPUs */
  MPI_Exscan(&mpi_atoms, &startidx, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  
  /* Öffnen des Outputfiles */
  MPI_File_open(MPI_COMM_WORLD, mpioutfilename, 
		MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  if (myid==0){

   /* Header schreiben: Filetyp*/
   MPI_File_write_at(fh, 0, "IMD", 3, MPI_CHAR, &iostatus);
   /* Länge des Headers schreiben */
   MPI_File_write_at(fh, 3, &disp, 1, MPI_INT, &iostatus);
   /* Gesamtzahl Atome */
   MPI_File_write_at(fh, 5, &atoms_total, 1, MPI_LONG, &iostatus);
   /* Zahl der Observablen */
   MPI_File_write_at(fh, 13, &observables, 1, MPI_SHORT, &iostatus);
   // SIMBOX
   //box_x
   MPI_File_write_at(fh,15,&box_x.x,1,MPI_DOUBLE,&iostatus);
   MPI_File_write_at(fh,23,&box_x.y,1,MPI_DOUBLE,&iostatus);
   MPI_File_write_at(fh,31,&box_x.z,1,MPI_DOUBLE,&iostatus);

   //box_y
   MPI_File_write_at(fh,39,&box_y.x,1,MPI_DOUBLE,&iostatus);
   MPI_File_write_at(fh,47,&box_y.y,1,MPI_DOUBLE,&iostatus);
   MPI_File_write_at(fh,55,&box_y.z,1,MPI_DOUBLE,&iostatus);

   //box_z
   MPI_File_write_at(fh,63,&box_z.x,1,MPI_DOUBLE,&iostatus);
   MPI_File_write_at(fh,71,&box_z.y,1,MPI_DOUBLE,&iostatus);
   MPI_File_write_at(fh,79,&box_z.z,1,MPI_DOUBLE,&iostatus);


   /* Typ der Observablen */
   MPI_File_write_at(fh, 87, "Mass, Position X,Y,Z, Velocity X,Y,Z", 36, 
                     MPI_CHAR, &iostatus);
  };

  /* Setzen eines Offsets für Header (disp) und Datenlänge */
  myoffset=disp+startidx*byte*OBSERVABLES;

  /* Setzen der Position für alle CPUs */
  MPI_File_set_view(fh, myoffset, MPI_PACKED, MPI_PACKED, DATAREP, 
                    MPI_INFO_NULL);

  /* Daten schreiben */
  MPI_File_write_all(fh, atomsbuffer, atombufferpos, MPI_PACKED, &iostatus);

  /* File schließen */
  MPI_File_close(&fh); 

printf("WRITTEN\n");

}
// *******************************************************************************************
// * MPI-INPUT
// * Hier ist ein anderes Kommunikationsschema noetig als bei den bisherigen Varianten,
// * denn bisher hat nur proc0 gelesen und alle anderen haben in einer schleife (recv_atoms)
// * auf pakete gewartet 
// * (bei parallel_input natuerlich nicht, aber da muss ohnehin nicht kommuniziert werden)
// * Loesung: parallel einlesen, buffer packen,anschließend eine loop ueber alle procs wobei jeder proc abwechselnd
// * seine nachrichten verschickt, während alle anderen in der recv-schleife warten
// * ACHTUNG: der msg-buffer ist begrenzt (etwa 60e6 atome). Sollen mehr verschickt werden,
// *	      führt dies zum segmentation fault
// ********************************************************************************************
void read_atoms_mpiio(str255 infilename)
{
  short observables=OBSERVABLES;
  short readobservables=0;
  long readatoms=0;
  //char readatoms[10];
  int byte=8;
  long startidx=0;
  int len;
  char mpifilename[255];
  char obstype[255];
 
  len=strlen(infilename);
  strncpy(mpifilename,"\0",len-4);
  strncpy(mpifilename,infilename,len-5);
  strcat(mpifilename,"mpiio");

  if ((0 == myid) && (0==myrank)) 
     printf("Reading atoms from %s via MPIIO.\n", infilename); //mpifilename);
  short disp=0;
  char incar[3];
  MPI_File fh;
  MPI_Status iostatus;

  ///////////////////////
  // READ HEADER       //
  ///////////////////////
  MPI_File_open(MPI_COMM_WORLD, mpifilename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  MPI_File_read_at(fh, 0, &incar, 3, MPI_CHAR, &iostatus);
  MPI_File_read_at(fh, 3, &disp, 1, MPI_INT, &iostatus);

  MPI_File_read_at(fh, 5, &readatoms, 1, MPI_LONG, &iostatus);
  MPI_File_read_at(fh, 13,&readobservables, 1, MPI_SHORT, &iostatus);
  MPI_Barrier(MPI_COMM_WORLD);

  //box_x
  MPI_File_read_at(fh,15,&box_x.x,1,MPI_DOUBLE,&iostatus);
  MPI_File_read_at(fh,23,&box_x.y,1,MPI_DOUBLE,&iostatus);
  MPI_File_read_at(fh,31,&box_x.z,1,MPI_DOUBLE,&iostatus);

  //box_y
  MPI_File_read_at(fh,39,&box_y.x,1,MPI_DOUBLE,&iostatus);
  MPI_File_read_at(fh,47,&box_y.y,1,MPI_DOUBLE,&iostatus);
  MPI_File_read_at(fh,55,&box_y.z,1,MPI_DOUBLE,&iostatus);

  //box_z
  MPI_File_read_at(fh,63,&box_z.x,1,MPI_DOUBLE,&iostatus);
  MPI_File_read_at(fh,71,&box_z.y,1,MPI_DOUBLE,&iostatus);
  MPI_File_read_at(fh,79,&box_z.z,1,MPI_DOUBLE,&iostatus);


  ////////// //////////////
  // NOW READ ATOMS
  ////////////////////////7
  int atoms_per_cpu=(int) ((double) readatoms)/((double)num_cpus);
  int atoms_per_cpu_extra; //nuer fuer letzten proc von bedeutung
  if(myid==num_cpus-1)
        atoms_per_cpu_extra=atoms_per_cpu+(((int) readatoms)-atoms_per_cpu*num_cpus);

  header_info_t info;
  cell *input;
  long addnumber = 0;
  int  p=0, k, maxc1=0, maxc2, count_atom;
  int  i, s, n, to_cpu, have_header=0, count;
  vektor   pos, axe;
  ivektor  cellc;
  real     m, d[MAX_ITEMS_CONFIG];
  minicell *to;
  msgbuf   *input_buf=NULL, *b;


  //READ BOX
  
//  if (box_from_header == 1) 
  read_box(infilename); //header mit box-grosse einlesen
  make_box();  //enthält init-cells
  box_x.x+=shiftx_front+shiftx_rear;

  // Set up 1 atom input cell
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Cannot allocate input cell.");
  input->n_max=0;
  alloc_cell(input, 1);
  natoms=0;

  
  /* allocate num_sort and num_vsort on all CPUs */
  if ((num_sort = (long *) calloc(ntypes,sizeof(long)))==NULL)
    error("cannot allocate memory for num_sort\n");
  if ((num_vsort = (long *) calloc(vtypes,sizeof(long)))==NULL)
    error("cannot allocate memory for num_vsort\n");


  //SET UP MESSAGE BUF
  /* allocate temporary input buffer */
  inbuf_size /= (sizeof(real) * (num_cpus-1));  //<-- mui importante, sonst fehler bei alloc (bei großen proben)
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


  int myofs;
  int obs=readobservables;
  double buf[OBSERVABLES];

  int nr=myid*atoms_per_cpu;
/*
      printf("proc:%d, reading %d atoms from %d to %d\n",myid,atoms_per_cpu,
	     disp+i*byte*obs+myid*atoms_per_cpu*byte*obs,
	     disp+i*byte*obs+myid*atoms_per_cpu*byte*obs+ atoms_per_cpu*byte*obs);
*/

  double* atomsbuffer;
  if(myid!=num_cpus-1)
    atomsbuffer=(double*) malloc(atoms_per_cpu*obs*sizeof(double));
  else  //letzter proc unter umstaenden mehr --> quick & dirty, geht schoener
    atomsbuffer=(double*) malloc(atoms_per_cpu_extra*obs*sizeof(double));

  myofs=disp+ myid*atoms_per_cpu*byte*obs;

  MPI_File_set_view(fh, myofs, MPI_DOUBLE, MPI_DOUBLE, DATAREP,
                    MPI_INFO_NULL);
  if(myid!=num_cpus-1)
    MPI_File_read_all(fh,atomsbuffer,atoms_per_cpu*obs,MPI_DOUBLE,&iostatus);
  else //letzter proc muss unter umstaenden mehr lesen
  {
    MPI_File_read_all(fh,atomsbuffer,atoms_per_cpu_extra*obs,MPI_DOUBLE,&iostatus);
    atoms_per_cpu=atoms_per_cpu_extra;
  }
   
  int l;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_File_close(&fh);
  for(k=0;k<num_cpus;k++)
  {
    if(myid==k)
    {
      for(i=0;i<atoms_per_cpu;i++)
      { 
    //printf("proc:%d,atomnr:%d / %d\n",myid,nr,atoms_per_cpu);
	//myofs=disp+ i*byte*obs + myid*atoms_per_cpu*byte*obs;
	//MPI_File_read_at(fh, myofs, &buf, obs, MPI_DOUBLE, &iostatus);
	myofs=i*obs;
	input->n=1;
	NUMMER(input,0)=nr;
	MASSE(input,0)=atomsbuffer[myofs];
	ORT(input,0,X)=atomsbuffer[myofs+1]+shiftx_front;    
	ORT(input,0,Y)=atomsbuffer[myofs+2];
	ORT(input,0,Z)=atomsbuffer[myofs+3];
	IMPULS(input,0,X)=atomsbuffer[myofs+4]*MASSE(input,0);
	IMPULS(input,0,Y)=atomsbuffer[myofs+5]*MASSE(input,0);
	IMPULS(input,0,Z)=atomsbuffer[myofs+6]*MASSE(input,0);

	SORTE(input,0)=0;	//erstmal nur 1 atom type
	VSORTE(input,0)=0;
	s=0;

	KRAFT(input,0,X)=0.0;
	KRAFT(input,0,Y)=0.0;
	KRAFT(input,0,Z)=0.0;

	pos.x=ORT(input,0,X);
	pos.y=ORT(input,0,Y);
	pos.z=ORT(input,0,Z);

	cellc = cell_coord(pos.x,pos.y,pos.z);
	to_cpu = cpu_coord(cellc);
	count_atom=0;

	/* to_cpu is in my input group, but not myself */
	if(myid != to_cpu)
	{
	  b = input_buf + to_cpu;
	  if (b->data != NULL) 
	  {
	    copy_atom_cell_buf(b, to_cpu, input, 0); //erstmal nur in buffer ablegen
	    if (b->n_max - b->n < atom_size)  //erst senden ab gewisser anzahl an atomen
	    {
	      MPI_Send(b->data, b->n, REAL, to_cpu, INBUF_TAG, MPI_COMM_WORLD);
	      b->n = 0;
	    }
	    count_atom = 1; // ?
	  }
	}
	else
	{
	  cellc = local_cell_coord(cellc);
	  to = PTR_VV(cell_array,cellc,cell_dim);
	  INSERT_ATOM(to, input, 0);
	  count_atom = 1; // ?
	}
	nr++; 
	if(count_atom)
	{
	  natoms++;
	  //Deg. of freedom, s = atom-typ
	  if(s<ntypes)
	  {
	    nactive+=DIM;
	  }
	  else
	  {
	      nactive += (long) (restrictions+s)->x;
	      nactive += (long) (restrictions+s)->y;
	      nactive += (long) (restrictions+s)->z;
	  }
	  num_sort [ SORTE(input,0)]++;
	  num_vsort[VSORTE(input,0)]++;
	}
      };// for i<atoms_per_cpu
      // send final msg
      for (l=0; l<num_cpus; l++) 
      {
	if(l==myid) continue; //an sich selber senden ist nicht noetig
        b = input_buf + l;
//printf("myid:%d, to:%d, buf.n:%d,buf.max:%d\n",
//        myid,l,b->n,b->n_max);
        if (b->data)
        {
          b->data[b->n++] = 0.0;
          MPI_Send(b->data, b->n, REAL, l, INBUF_TAG+1, MPI_COMM_WORLD);
        }
      }
    } // if myid==k
    else //restlichen procs warten derweil in recv-loop
    {
      recv_atoms_mpiio(k);
    }
  } // loop k ueber procs

  free(atomsbuffer); 
  free(input_buf);
  //read_atoms_cleanup_mpiio();
  read_atoms_cleanup();
}

void read_atoms_cleanup_mpiio(void)
{

}

void recv_atoms_mpiio(int src)
{
  MPI_Status status;
  int finished = 0;
  msgbuf b = {NULL,0,0};
//printf("myid:%d,waiting for msg from %d\n",myid,src);
  alloc_msgbuf(&b, inbuf_size);
  do {
    MPI_Recv(b.data, inbuf_size, REAL, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, REAL, &b.n);
    if (status.MPI_TAG==INBUF_TAG+1) { b.n--; finished=1; } /* last buffer */
    process_buffer( &b );
//printf("proc:%d recv from %d,fin:%d\n",myid,src,finished);
  } while (0==finished);
  free_msgbuf(&b);
//printf("proc:%d,recv from %d done\n",myid,src);
}

#endif
