
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2012 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* prototypes.h -- Function prototypes for IMD
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/******************
* MY FUNCTIONS    *
*******************/
#include <complex.h>


#ifdef NRB
int init_nrb(void);
int approx(double x,double x0,double deps);
int compare( const void* a, const void* b);
int nrb_build_ifromid(void);
int nrb_binarySearch( int lo, int hi, int x,int n,int** arr);
int nrb_forces(void);
void nrb_test_forces(void); //ext. kraft um schockwelle zu erzeugen

//void pack_nrb( msgbuf *b, int k, int l, int m);
void pack_nrb( msgbuf *b, int k, int l, int m, vektor v );
void unpack_nrb_max( msgbuf *b, int k, int l, int m );
void unpack_nrb( msgbuf *b, int k, int l, int m );
void copy_nrb_max( int k, int l, int m, int r, int s, int t, vektor v );
void copy_nrb( int k, int l, int m, int r, int s, int t, vektor v );

void nrb_sendmomenta(void (*add_func)   (int, int, int, int, int, int),
              void (*pack_func)  (msgbuf*, int, int, int),
              void (*unpack_func)(msgbuf*, int, int, int));
void nrb_addmomenta( int k, int l, int m, int r, int s, int t );
void nrb_addmomentav( int k, int l, int m, int r, int s, int t, vektor v);
void nrb_copymomenta( int k, int l, int m, int r, int s, int t, vektor v );
//void nrb_packmomenta( msgbuf *b, int k, int l, int m);
void nrb_packmomenta( msgbuf *b, int k, int l, int m);
void nrb_packmomentav( msgbuf *b, int k, int l, int m, vektor v );
void nrb_unpackmomentaadd( msgbuf *b, int k, int l, int m );
void nrb_unpackmomentacopy( msgbuf *b, int k, int l, int m );

void nrb_inverse_send_cells(void (*copy_func)  (int, int, int, int, int, int, vektor),
                void (*pack_func)  (msgbuf*, int, int, int, vektor),
                void (*unpack_func)(msgbuf*, int, int, int));

void nrb_send_cells(void (*copy_func)  (int, int, int, int, int, int, vektor),
                void (*pack_func)  (msgbuf*, int, int, int, vektor),
                void (*unpack_func)(msgbuf*, int, int, int));


int  sendrecv_buf(msgbuf *send, int to_cpu,  //das wird normalerweise nur benutzt wenn AR aus ist,sonst isend/irecv. 
                  msgbuf *recv, int from_cpu, MPI_Status *status); //ich nutze besser erst einmal die "sichere" variante

void nrb_writerestart(int fzhlr);
void write_nrb_fun(FILE*out); //==write_atoms_config(..)

void nrb_readrestart(str255); //darunter alle routinen die fürs einlesen nötig sind
void nrb_readinputfile(str255); //gründlichere suche fürs einlesen aus pre-equilib. file noetig!

void nrb_recv_atoms(void);
void nrb_copy_atom_buf_cell(minicell *p, msgbuf *b, int start);
void nrb_copy_atom_cell_buf(msgbuf *to, int to_cpu, cell *p, int ind );
void nrb_move_atom(cell *to, cell *from, int index);
void nrb_copy_atom_cell_cell(cell *to, int i, cell *from, int j);
void nrb_recv_process_buffer(msgbuf *b);
#endif

#ifdef FILTER
void filter_init(void);
void filter_forces(void);
void filter_atoms(void);
//int filter_check_neighs(cell*p,int l,cell* pstart,int lstart);
int filter_check_neighs(cell*p,int l);
void add_filter( int k, int l, int m, int r, int s, int t );
void pack_filter( msgbuf *b, int k, int l, int m);
void unpack_filter( msgbuf *b, int k, int l, int m );
#endif 

//TTM-KRAM
#ifdef TTM
void do_FILLMESH(void); //Wenn option==1 --> zusätzlich EfromT interpolieren (nur 1 mal pro md-step noetig)
void do_COMMFLUX(void);
void do_ADV(double tau);
void update_u(void);
void do_DIFF(double tau);
void factor_xi(void); //ist jetzt losgelöst von do_DIFF
void do_cell_activation(void); // Cell activ/inactive

//Interpol. tables for electronic EOS
int read_bc_interp(struct bicubinterp *bi,const char *fname);
void free_bc_interp(struct bicubinterp *bi);
double bcinterp(struct bicubinterp *bi,double xp,double yp);
int nn_read_table(nn_interp *nn, const char *fname);

int free_lin_interp(struct lninterp *bi);
int read_lin_interp(struct lninterp *bi,const char *fname,int option); // fuer scattered und gridded data je nach option
double lininterp(struct lninterp *bi,double xp,double yp,int option);  // (wird bisher nicht genutzt)

//TRICUB INTERPOL STUFF
int read_tricub_interp(struct tricubinterp *tri,const char *fname);
double do_tricubinterp(struct tricubinterp *tri,double xp,double yp,double zp);
void tricub_get_coeff(double a[64], double f[8], double dfdx[8], double dfdy[8], double dfdz[8],
                        double d2fdxdy[8], double d2fdxdz[8], double d2fdydz[8], double d3fdxdydz[8]);

double tricub_eval(double a[64], double x, double y, double z);
void tricub_get_coeff_stacked(double a[64], double x[64]);
int tricub_ijk2n(int i, int j, int k);
double KappaInterpol(double rho,double Te,double Ti); //kappa from tricub-interpol. table

double NueffInterpol(double rho,double Te,double Ti);    //Nu_eff (=gammap) from from table i(FDTD needs this)
double EpsImagInterpol(double rho,double Te,double Ti);
double EpsRealInterpol(double rho,double Te,double Ti);

//EOS STUFF
double QfromT(double T,double rho);
//double CfromT(double T,double rho);
double EOS_cve_from_r_te(double r,double t);
double EOS_ee_from_r_te(double r,double t);
double EOS_te_from_r_ee(double r,double e);
double eeminfun(double x, double r, double e);
//minimiert e_from_r_te-e für EOSMODE=1 (D.h. interne energie aus tabelle)
double fminbnd(double a, double b,double (* f)(double,double,double), double tolx,double rho,double eng);
//selbes für FEG-model
double fminbnd2(double a, double b,double (* f)(double,double,double,double), double tolx,double rho,double ne,double eng);
void CFL_maxdt(void);  //CFL Condition for max. stable diffusion timestep
double chempot(double ne,double Te);
double FEG_ee_from_r_ne_te(double r,double ne, double T);
double FEG_cve_from_ne_te(double r,double ne,double T);
double FEG_eeminfun(double x, double r, double ne, double e);
double FEG_te_from_r_ne_ee(double r,double ne,double ee);

#endif


#ifdef TMM //Transfer-matrix methode
double EE(double z, double delta, double complex kl, double complex epsl, double complex Bplus, double complex Bminus);
int do_tmm(real dt);
double Runge5(double delta, double complex kl, double complex epsl, double complex Bplus, double complex Bminus);
void invmat(double complex *a,double complex *p);/* Inverse of 2x2 matrix*/
void matmul(double complex *a,double complex *b,double complex *p); /*2x2 Matrix multiplication*/
void matmul2(double complex *mat,double complex *vec,double complex *out);/* 2x2 matrix and 2x1 vector multiplication */ 
int tmm_init(void);
int tmm_get_epsilon(double lambda, double Te,double Ti, double Z, double Ne, real* re, real* img);
double complex epsl_bb(int bi);
int tmm_read_arr(char* infilename,double ***arr,int cols,int *rowsout);
double tmm_K2(double nu);
double tmm_K1(double nu);
#endif

#ifdef COLRAD
int colrad_ydot(double t, N_Vector u, N_Vector udot, void *user_data);
void do_Saha(double Te,double totalc,double ne,N_Vector y);
int colrad_GetCoeffs(N_Vector y,double It, void * user_data);
//double ExpInt(double x);
// double fak(double t, double x, double j,double s); //aux. function for genexpint
// double genexpint(double x,double ss,double j);
double EscapeFactor(double w,double tau0);
double StarkWidth(double nu,double nl,double Te,double Ti,double Tr,double Ne,double Ni);
double EinsteinCoeff(double n1,double n2,double g2,double DeltaE);

//folgende funcs sind alle für Berechnung des Transmission/Escape-Faktors!
double Lorentzian(double w,double lambda); //both params in angs!
double integrand(double w,double lambda,double tau,double phi0); // integrant von T=int[ exp(-tau*phi*phi0)*phi dx ]
double trapz(double a,double b,int n,double w,double tau); //integriert die func. "integrand"
void colrad_read_states(void);
void colrad_init(void);
void do_colrad(double dt);
void colrad_Saha_init(int i,int j,int k);
#endif
#ifdef FDTD
int init_fdtd(void);
void init_pml(void);
int fdtd_create_mpi_datatypes(void);
void fdtd_comm_ghost_cells(void);
void do_fdtd(double t);
int fdtd_softsource(double t);

//for drude dispersion
double getomegap(int i,int j);
double getgammap(int i,int j);

void SwapTTM(int i,int j,int k);
int fitDL(int i,int j,int k);
#endif

double fermi_E(double ne);
double fermi_T(double ne);
double r0(double ni);
double omega_pl(double ne);
double bMin(double X, double omega_las, double Z, double ni, double Te);
double bMax(double omega_las, double Z, double ni, double Te);
double coulomb_log(double X, double omega_las, double Z, double ni, double Te);
double ne_cr(double omega_las);
double fd_density(double natoms,double mass,double vol);
double complex epsl_bb(int bi);
double MeanCharge(double temp, double rho, double z0, double am,int i,int j,int k);
double numet(double A1, double A2, double Te, double Ti, double TF);
double numax(double A3, double vf, double Te, double ni);
double nupl(double omega_las, double Z, double ni, double ne, double Te);

//kappa und gamma nach povar
double getKappa(double Te, double Ti, double Ne,double Z);
double getGamma(double Te, double Ti, double Ne,double Z);

double nu_e_e(double Te, double EF,double Ne,double Na,double valence); 
double v_e(double Te,double EF);
double sigma_e_e(double Te, double EF,double Na,double valence);
double sigma_e_ph(double Te,double Ti,double EF,double Na,double Z);
double nueff(double Te,double Ti,double ne,double ni,double ni0,double Ce,double lod); //für drude Kombination aus Petrov & Mazhukin


char **strsplit(const char* str, const char* delim, size_t* numtokens);
void ttm_read(int number);

/*************************
* END OF MY FUNCTIONS    *
**************************/

/* the main routines - files imd.c, imd_main_3d.c */
int  main(int argc, char **argv);
int  main_loop(int simulation);

/* read parameters - file imd_param.c */
void read_command_line(int argc, char **argv);
int  read_parameters(char *paramfname, int phase);
int  getparamfile(char *paramfname, int phase);
void check_parameters_complete(void);
void broadcast_params(void);

/* read and access potential tables - file imd_potential.c */
void setup_potentials( void );
void read_pot_table ( pot_table_t*, char*, int, int );
void read_pot_table1( pot_table_t*, int, char*, FILE*, int );
void read_pot_table2( pot_table_t*, int, char*, FILE*, int );
void free_pot_table ( pot_table_t*);
#ifdef LINPOT
void make_lin_pot_table( pot_table_t, lin_pot_table_t* );
#endif

#ifdef FEFL
/* void atom_int_ec(real *pot, real *grad, int p_typ, real r2); */
void calc_fefl(void);
#endif

#ifdef PAIR
void pair_int_lj(real *pot, real *grad, int p_typ, int q_typ, real r2);
void pair_int_ljg(real *pot, real *grad, int p_typ, int q_typ, real r2);
void pair_int_morse(real *pot, real *grad, int p_typ, int q_typ, real r2);
void pair_int_buck(real *pot, real *grad, int p_typ, int q_typ, real r2);
void pair_int_mstr(real *pot, real *grad, int p_typ, int q_typ, real r2);
#ifdef EWALD
void pair_int_ewald(real *pot, real *grad, int p_typ, int q_typ, real r2);
#endif
#ifdef COULOMB
void create_coulomb_tables(void);
void pair_int_coulomb(real *pot, real *grad, real r2);
#endif
#endif
void pair_int2  (real*, real*, int*, pot_table_t*, int, int, real);
void pair_int3  (real*, real*, int*, pot_table_t*, int, int, real);
void   val_func2(real*,        int*, pot_table_t*, int, int, real);
void   val_func3(real*,        int*, pot_table_t*, int, int, real);
void deriv_func2(       real*, int*, pot_table_t*, int, int, real);
void deriv_func3(       real*, int*, pot_table_t*, int, int, real);
void init_pre_pot(void);
void create_pot_table(pot_table_t *pt);
void test_potential(pot_table_t, char*, int);
#if   defined(FOURPOINT)
void init_fourpoint(pot_table_t*, int);
#elif defined(SPLINE)
void init_spline(pot_table_t*, int, int);
#else
void init_threepoint(pot_table_t*, int);
#endif

/* read configuration - files imd_io_2/3d.c */
void read_atoms(str255 infilename);
#ifdef MPIIO
void write_atoms_config_mpiio(void);
void read_atoms_mpiio(str255 infilename);
void recv_atoms_mpiio(int src);
void read_atoms_cleanup_mpiio(void);
#endif
void read_atoms_cleanup(void);  
#ifdef MPI
void recv_atoms(void);
#endif

/* generate configuration - file imd_generate.c */
void generate_atoms(str255 infilename);
void init_cubic(void);
void generate_fcc(int maxtyp);
void generate_lav(void);
void init_hex(void);
void generate_hex(void);
void generate_SiO2(void);
#ifdef QUASI /* generate quasicrystal - file imd_qc.c */
void init_qc(void);
void generate_qc(void);
#endif

/* miscellaneous routines - files imd_misc.c, imd_time.c, imd_maxwell.c */
void usage(void);
void imderror(char *msg);
void error_str(char *msg, char *str);
void error_str_str(char *msg, char *str1, char *str2);
void warning(char *);
void warning_str(char *, char *);
void warning_str_str(char *, char *, char *);
void imd_init_timer(imd_timer *timer, int flag, char *desc, char *color);
void imd_start_timer(imd_timer *timer);
void imd_stop_timer(imd_timer *timer);
void maxwell(real TEMP);
int  endian(void);
integer SwappedInteger(integer);
float   SwappedFloat  (float  );
double  SwappedDouble (double );

/* start and stop MPI - files imd_mpi_util.c, imd_geom_mpi_*.c */
#if defined(MPI) || defined(NEB)
void init_mpi(void);
void shutdown_mpi(void);
#endif
#ifdef MPI
void alloc_msgbuf(msgbuf*, int);
void realloc_msgbuf(msgbuf*, int);
void free_msgbuf(msgbuf*);
#endif
void setup_mpi_topology(void);
void calc_cpu_dim(void);

/* manage MPI buffers and buffer cells - file imd_mpi_util.c */
#ifdef MPI
void setup_buffers(void);
void empty_mpi_buffers(void);
#ifdef SR
int  sendrecv_buf(msgbuf *send, int to_cpu, 
                  msgbuf *recv, int from_cpu, MPI_Status *status);
#else
int  irecv_buf(msgbuf *b, int from_cpu, MPI_Request *req);
int  isend_buf(msgbuf *b, int to_cpu,   MPI_Request *req);
#endif
void empty_buffer_cells(void);
void init_io(void);
#endif

/* make and maintain cells and their geometry - files imd_geom_*.c */
ivektor maximal_cell_dim(void);
vektor  back_into_box(vektor pos);
vektor  vec_prod(vektor u, vektor v);
void make_box(void);
void init_cells(void);
void make_cell_lists(void);
void check_pairs(void);
void move_atom(cell *to, cell *from, int index);
void copy_atom_cell_cell(cell *to, int i, cell *from, int j);
#ifdef VEC
void move_atom_mini(minicell *, minicell *, int);
void insert_atom(minicell *, cell *, int);
void alloc_minicell(minicell *, int);
#endif
void memalloc(void *, int, int, int, int, int, char *);
void alloc_cell(cell *thecell, int count);
#ifdef TWOD
ivektor cell_coord(real x, real y);
#else
ivektor cell_coord(real x, real y, real z);
#endif
ivektor local_cell_coord(ivektor cellc);
int     cpu_coord(ivektor cellc);
ivektor cpu_coord_v(ivektor cellc);
int     cpu_grid_coord(ivektor cellc);

/* force computation - files imd_main_*.c, imd_forces_*.c */
void calc_forces(int steps);
#ifdef EXTPOT
void init_extpot(void);
void calc_extpot(void);
void move_extpot(real factor);
#endif
#ifdef USEFCS
void init_fcs(void);
void pair_int_fcs(real *, real *, real);
void calc_forces_fcs(int);
void fcs_update_pottab(void);
void fcs_cleanup(void);
#endif
void do_forces(cell*, cell*, vektor, real*, real*, real*, real*, real*, real*, real*, real*);
#ifdef COVALENT
void do_forces2(cell*, real*, real*, real*, real*, real*, real*, real*, real*);
#endif
#ifdef EAM2
void do_forces_eam2(cell*, cell*, vektor, real*, real*, real*, real*, real*, real*, real*);
void do_embedding_energy(void);
#endif
#ifdef NBLIST
int  estimate_nblist_size(void);
void make_nblist(void);
void check_nblist(void);
void deallocate_nblist(void);
#endif
#ifdef MEAM
void init_meam(void);
#endif
#ifdef KEATING
void init_keating(void);
#endif
#ifdef STIWEB
void init_stiweb(void);
void pair_int_stiweb(real *pot, real *grad, int p_typ, int q_typ, real r2);
#endif
#ifdef TTBP
void init_ttbp(void);
#endif
#ifdef TERSOFF
void init_tersoff(void);
void pair_int_tersoff(real *pot, int p_typ, int q_typ, real r2);
#endif
#ifdef TERSOFFMOD
void init_tersoffmod(void);
#endif
#ifdef BRENNER
void init_brenner(void);
void pair_int_brenner(real *pot, int p_typ, int q_typ, real r2);
#endif

#if defined (TERNBCC) || defined (XT)
real g(real cos_theta);
real dg(real cos_theta);
#endif

#ifdef UNIAX
void gay_berne ( vektor r12, vektor e1, vektor e2, 
		 real rsqr, vektor s1, vektor w1, 
		 real *pot12, vektor *force12, 
		 vektor *torque12, vektor *torque21 );
#endif

/* communication for force computation - files imd_comm_force_2/3d.c */
#if defined(MPI) || defined(NBLIST)
#ifdef TWOD
void send_cells (void (*copy_func)  (int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int),
                 void (*unpack_func)(msgbuf*, int, int));
void sync_cells (void (*copy_func)  (int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int),
                 void (*unpack_func)(msgbuf*, int, int));
void send_forces(void (*add_func)   (int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int),
                 void (*unpack_func)(msgbuf*, int, int));
void copy_cell    ( int j, int k, int l, int m );
void pack_cell    ( msgbuf *b, int j, int k );
void unpack_cell  ( msgbuf *b, int j, int k );
void add_forces   ( int j, int k, int l, int m );
void pack_forces  ( msgbuf *b, int j, int k );
void unpack_forces( msgbuf *b, int j, int k );
#else  /* 3D */
void send_cells (void (*copy_func)  (int, int, int, int, int, int, vektor),
                 void (*pack_func)  (msgbuf*, int, int, int, vektor),
                 void (*unpack_func)(msgbuf*, int, int, int));
void sync_cells (void (*copy_func)  (int, int, int, int, int, int, vektor),
                 void (*pack_func)  (msgbuf*, int, int, int, vektor),
                 void (*unpack_func)(msgbuf*, int, int, int));
void send_forces(void (*copy_func)  (int, int, int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int, int),
                 void (*unpack_func)(msgbuf*, int, int, int));
void copy_cell    ( int k, int l, int m, int r, int s, int t, vektor v );
void pack_cell    ( msgbuf *b, int k, int l, int m, vektor v );
void unpack_cell  ( msgbuf *b, int k, int l, int m );
void add_forces   ( int k, int l, int m, int r, int s, int t );
void pack_forces  ( msgbuf *b, int k, int l, int m);
void unpack_forces( msgbuf *b, int k, int l, int m );
#ifdef EAM2
void copy_dF       ( int k, int l, int m, int r, int s, int t, vektor );
void add_rho       ( int k, int l, int m, int r, int s, int t );
void pack_dF       ( msgbuf *b, int k, int l, int m, vektor );
void pack_rho      ( msgbuf *b, int k, int l, int m );
void unpack_dF     ( msgbuf *b, int k, int l, int m );
void unpack_add_rho( msgbuf *b, int k, int l, int m );
#endif
#if defined(DIPOLE) || defined(KERMODE)
void copy_pind     ( int k, int l, int m, int r, int s, int t, vektor );
void copy_Eind     ( int k, int l, int m, int r, int s, int t, vektor );
void add_dipole    ( int k, int l, int m, int r, int s, int t );
void add_field     ( int k, int l, int m, int r, int s, int t );
void pack_pind     ( msgbuf *b, int k, int l, int m, vektor );
void pack_Eind     ( msgbuf *b, int k, int l, int m, vektor );
void pack_dipole   ( msgbuf *b, int k, int l, int m );
void pack_field    ( msgbuf *b, int k, int l, int m );
void unpack_pind   ( msgbuf *b, int k, int l, int m );
void unpack_Eind   ( msgbuf *b, int k, int l, int m );
void unpack_add_dipole( msgbuf *b, int k, int l, int m );
void unpack_add_field( msgbuf *b, int k, int l, int m );
#endif
#ifdef KERMODE
void yukawa_smooth_function(real ke_r,real rcut, real dcut,real *fc,real *dfc);
#endif
#ifdef CNA
void add_mark       ( int k, int l, int m, int r, int s, int t );
void pack_mark      ( msgbuf *b, int k, int l, int m );
void unpack_add_mark( msgbuf *b, int k, int l, int m );
#endif
#endif /* 3D  */
#endif /* MPI or NBLIST */

/* integrators - file imd_integrate.c */
void move_atoms_nve(void);
void move_atoms_mik(void);
void move_atoms_nvt(void);
void calc_dyn_pressure(void);
void move_atoms_npt_iso(void);
void move_atoms_npt_axial(void);
void move_atoms_mc(void);
void move_atoms_frac(void);
void move_atoms_sllod(void);
void move_atoms_nvx(void);
void move_atoms_msd(void);
void move_atoms_stm(void);
void move_atoms_ftg(void);
void move_atoms_finnis(void);
void move_atoms_ttm(void);
#ifdef SHOCK
void calc_pxavg(void);
#endif

/* fix distribution on cells - files imd_main_*.c, imd_mpi_util.c */
void do_boundaries(void);
void fix_cells(void);
#ifdef MPI
void copy_atoms_buf(msgbuf *to, msgbuf *from);
void copy_one_atom (msgbuf *to, int to_cpu, minicell *from, int index,int del);
void copy_atom_cell_buf(msgbuf *to, int to_cpu, cell *from, int index );
void copy_atom_buf_cell(minicell *to, msgbuf *b, int start);
void process_buffer(msgbuf *b);
void send_atoms(void);
#endif
/* write properties - file imd_io_*.c */
void write_eng_file(int steps);
void write_eng_file_header(void);
#ifdef EXTPOT
void write_fext(int steps);
void write_fext_header(void);
#endif
#ifdef RELAX
void write_ssdef(int steps);
void write_ssdef_header(void);
#endif
/* write configurations - files imd_io.c, imd_io_*.c */
void read_box(str255);
int  read_header(header_info_t *, str255);
#ifdef MPI
void broadcast_header(header_info_t *);
#endif
void flush_outbuf(FILE *out, int *len, int tag);
void write_itr_file(int fzhlr, int steps,char *suffix);
void write_config(int fzhlr, int steps);
#ifdef RELAX
void write_ssconfig(int steps);
#endif
void write_config_select(int fzhlr, char *suffix,
  void (*write_atoms_fun)(FILE *out), void (*write_header_fun)(FILE *out));
void write_atoms_config(FILE *out);
void write_header_config(FILE *out);
void write_atoms_pic(FILE *out); 
void write_header_pic(FILE *out); 
void write_atoms_pos(FILE *out); 
void write_header_pos(FILE *out); 
#ifdef DISLOC
void write_atoms_dem(FILE *out);
void write_header_dem(FILE *out);
void write_config_dem(int nr);
void write_atoms_dsp(FILE *out);
void write_header_dsp(FILE *out);
void write_config_dsp(int nr);
#endif
#ifdef CNA
void do_cna(void);
void do_cna_func(cell *p, cell* q, vektor pbc);
void domino(int start, int end, int listlength, int *max_chain, int *chain);
void write_atoms_cna(FILE *out);
void write_header_cna(FILE *out);
void write_atoms_crist(FILE *out);
void write_header_crist(FILE *out);
void init_cna(void);
void write_statistics(void);
void sort_pair_types(void);
#endif
#ifdef EFILTER
void write_atoms_ef(FILE *out);
void write_header_ef(FILE *out);
void write_config_ef(int nr);
#endif
#ifdef NNBR
void write_atoms_nb(FILE *out);
void write_header_nb(FILE *out);
void write_config_nb(int nr);
#endif
#ifdef WRITEF
void write_atoms_wf(FILE *out);
void write_header_wf(FILE *out);
void write_config_wf(int nr);
#endif
#ifdef STRESS_TENS
void write_atoms_press(FILE *out);
void write_header_press(FILE *out);
void calc_tot_presstens(void);
void write_config_press(int nr);
#endif
#ifdef REFPOS
void init_refpos(void);
#endif
#ifdef ZAPP
void init_zapp(void);
void zapp(void);
#endif
#ifdef BEND
void init_bend(void);
void update_bend(void);
#endif
#ifdef AVPOS
void write_atoms_avp(FILE *out);
void write_header_avp(FILE *out);
void write_avpos_itr_file(int fzhlr, int steps);
void write_config_avpos(int nr);
#endif
#ifdef FORCE
void write_atoms_force(FILE *out);
void write_header_force(FILE *out);
void write_config_force(int nr);
#endif
#ifdef MSQD
void write_atoms_sqd(FILE *out);
void write_header_sqd(FILE *out);
void write_config_sqd(int nr);
#endif
void reduce_displacement(vektor *d);

#ifdef SOCKET_IO
void init_socket(void);
int  connect_visualization(void);
void close_socket(void);
void check_socket(void);
void vis_init(void);
void vis_write_config_quit(void);
void vis_send_msg(char *msg);
void vis_check_atoms_flags(void);
void vis_write_atoms_buf(int *len, int tag);
void vis_write_atoms_fun(void);
void vis_init_atoms(void);
void vis_write_atoms(void);
void vis_change_params(void);
void vis_change_params_deform(integer flag);
void vis_restart_simulation(void);
#endif

#ifdef BBOOST
void  init_bboost(void);
int bb_minimize(int);
#endif


/* write distributions - file imd_distrib.c */
void make_distrib_density(void); 
void write_distrib_density(int, int);
void make_write_distrib_select(int, void (*fun)(float*, cell*, int), 
			       int, int, char*, char*);
void write_distrib_header(FILE*, int, int, char*);
void write_distrib(int);
void dist_Ekin_fun       (float*, cell*, int);
void dist_Epot_fun       (float*, cell*, int);
void dist_Ekin_long_fun  (float*, cell*, int);
void dist_Ekin_trans_fun (float*, cell*, int);
void dist_Ekin_comp_fun (float*, cell*, int);
void dist_shock_shear_fun(float*, cell*, int);
void dist_shear_aniso_fun(float*, cell*, int);
void dist_press_fun      (float*, cell*, int);
void dist_presstens_fun  (float*, cell*, int);
void dist_presstens_xx_fun  (float*, cell*, int);
void dist_presstens_yy_fun  (float*, cell*, int);
#ifndef TWOD
void dist_presstens_zz_fun  (float*, cell*, int);
void dist_presstens_yz_fun  (float*, cell*, int);
void dist_presstens_zx_fun  (float*, cell*, int);
#endif
void dist_presstens_xy_fun  (float*, cell*, int);
void dist_pressoff_fun   (float*, cell*, int);
void dist_pressxy_fun   (float*, cell*, int);
void dist_pressyz_fun   (float*, cell*, int);
void dist_presszx_fun   (float*, cell*, int);
void dist_vxavg_fun   (float*, cell*, int);

#ifdef ATDIST
void   init_atdist(void);
void update_atdist(void);
void  write_atdist(void);
void write_atoms_atdist_pos(FILE*);
void write_header_atdist_pos(FILE*);
void write_config_atdist_pos(int nr);
#endif

#ifdef DIFFPAT
void   init_diffpat(void);
void update_diffpat(int );
void  write_diffpat(void);
#endif

/* write pictures - files imd_io_*.c */
void write_pictures(int steps);
void write_pictures_bitmap(int steps);

/* shear, deform of load sample - files imd_deform.c, imd_load.c */
#ifdef HOMDEF
void shear_sample(void);
void expand_sample(void);
#ifdef TWOD
void lin_deform(vektor,vektor,real);
#else
void lin_deform(vektor,vektor,vektor,real);
#endif
void relax_pressure(void);
#endif
#ifdef DEFORM
void deform_sample(void);
#endif
#ifdef FRAC
void load_sample(void);
#endif

/* support for neighbor tables - files imd_alloc.c, imd_forces_covalent.c */
#if defined(COVALENT) || defined(NNBR_TABLE)
void do_neightab(cell *p, cell *q, vektor pbc);
neightab *alloc_neightab(neightab *neigh, int count);
void increase_neightab(neightab *neigh, int count);
#endif
#ifdef NNBR_TABLE
void do_neightab_complete();
void do_neightab2(cell *p, cell *q, vektor pbc);
#endif

/*Angle distribution analysis imd_ada.c */
#ifdef ADA
void do_ada(void);
void init_ada(void);
void pack_hopsToDefect( msgbuf *b, int k, int l, int m, vektor v);
void unpack_hopsToDefect( msgbuf *b, int k, int l, int m );
void copy_hopsToDefect( int k, int l, int m, int r, int s, int t, vektor v);
void pack_AdaAndNumber(msgbuf *b, int k, int l, int m, vektor v);
void unpack_AdaAndNumber(msgbuf *b, int k, int l, int m);
void copy_AdaAndNumber(int k, int l, int m, int r, int s, int t, vektor v);
void buildHopsToDefect();
void filterSurface();

/*Output for ada imd_io.c */
void write_header_ada(FILE *out);
void write_atoms_ada(FILE *out);
#endif

/* Nye Tensor Analysis imd_nyeTensorAnalysis_3d.c*/
#ifdef NYETENSOR
void init_NyeTensor();
void removeNyeTensorData();
void calculateNyeTensorData();
#endif

#ifdef LOADBALANCE
/* Communication routines for direct communication scheme*/
void lb_copyCellDataToSend(msgbuf*, ivektor*, int,
		void (*pack_func)(msgbuf*, int, int, int, vektor), int);
void lb_copyForcesDataToSend(msgbuf*, ivektor*, int,
		void (*pack_func)(msgbuf*, int, int, int));
void sync_cells_direct (void (*copy_func)  (int, int, int, int, int, int, vektor),
                 void (*pack_func)  (msgbuf*, int, int, int, vektor),
                 void (*unpack_func)(msgbuf*, int, int, int), int);
void lb_unpackCellDataFromBuffer(msgbuf*, int originCPU, void (*unpack_func)(msgbuf*, int, int, int));
void lb_unpackForcesDataFromBuffer(msgbuf*, void (*unpack_func)(msgbuf*, int, int, int));
void lb_initDirect(void);

void init_loadBalance(void);
cell* lb_accessCell(cell*, int, int, int, ivektor);

int balanceLoad(real (*getLoad)(void), vektor (*getCog)(void), int, int);
real lb_getLoad(void);
real lb_getVolume(void);
int lb_identifyCellType(int, int, int);
int lb_syncBufferCellAffinity(void);
void lb_relocateParticles(ivektor, ivektor, cell*);

void lb_makeNormal(int, int, int, int, lb_domainInfo*);
int lb_getTetraederVolumeIndexed(int, int, int, int, lb_domainInfo*);
real lb_getTetrahedronVolume(int, int, int, int, lb_domainInfo*);
void lb_makeNormals(lb_domainInfo*);
real lb_computeWarpage(int, int, int, int, lb_domainInfo*);
real lb_getDistanceToPlane(int, int, int, int, lb_domainInfo*);
int lb_isPointInDomain(vektor, lb_domainInfo*);
real lb_angleCos(int, int, int, lb_domainInfo*);
vektor lb_getCellCenter(int, int, int);
vektor lb_getCenterOfGravity();
void lb_updateDomain(lb_domainInfo*);
void lb_moveCorner(int, real, real, real);
void lb_moveAllCorners(real*, int);
int lb_isGeometryChangeValid(void);
int lb_isGeometryValid(lb_domainInfo*);

void write_lb_status(int);
void write_config_lb(FILE*);
void write_header_lb(FILE*);
void write_lb_file(int, int);

int lb_countAtoms(void);
void lb_computeVariance(void);

void lb_processCellDataBuffer(msgbuf*,
		void (*copy_func)(int, int, int, int, int, int, vektor),
		void (*unpack_func)(msgbuf*, int, int, int), int, int);
void lb_processForceDataBuffer(msgbuf*,
		void (*unpack_func)(msgbuf*, int, int, int), int);

void lb_moveCornersReset(real*, int iteration);

void lb_balanceOneAxisOrthogonal(int, int, int*, int*);
void balanceOrtho(void);

#endif

#ifdef BBOOST
void do_bb_neightab(cell *p, cell *q, vektor pbc);
bb_neightab *alloc_bb_neightab(bb_neightab *bb_neigh, int count);
#endif

#ifdef EWALD
/* support for computation of Coulomb forces */
void do_forces_ewald(int);
void do_forces_ewald_real(void);
void do_forces_ewald_fourier(void);
void init_ewald(void);
#endif
#if defined(EWALD) || defined(COULOMB)
real erfc1(real x);
#endif

#ifdef EPITAX
/* support for molecular beam epitaxy - file imd_epitax.c */
void create_atom(int type, real mass, real temp);
void delete_atoms(void);
real substrate_level(void);
void calc_poteng_min(void);
void check_boxheight(void);
ivektor cell_map(ivektor cellc);
#endif

/* support for dislocations - file imd_io.c */
#ifdef DISLOC
void reset_Epot_ref(void);
void update_ort_ref(void);
#endif

/* support for average over positions */
#ifdef AVPOS
void update_avpos(void);
void add_positions(void);
#ifdef STRESS_TENS
void update_avpress(void);
void add_presstensors(void);
#endif
#endif

/* support for correlation functions - file imd_correl.c */
#if defined(CORRELATE) || defined(MSQD)
void init_correl(int, int);
void correlate(int istep, int refstep, unsigned seqnum);
void write_msqd(int);
#endif
#ifdef CORRELATE
void write_add_corr(int it, int steps, unsigned seqnum);
#endif

/* support for heat transport - file imd_transport.c */
#ifdef NVX
void write_temp_dist(int steps);
#endif
#ifdef HC
void do_heat_cond(int steps);
void write_heat_current(int steps);
#endif

#ifdef CG
void reset_cg(void);
void cg_step(int steps);
void write_cgconfig(int steps);
real fonedim (real) ;
real fonedim_sd (real) ;
void cg_calcgamma(void);
void set_hg(void);
void calc_fnorm(void);
void calc_fnorm_g_h(void);
void move_atoms_cg(real);
void move_atoms_sd(real);
int linmin();
int mnbrak(real *ax, real *bx, real *cx, real *fa, real *fb, real *fc);
int brent(real ax, real bx, real cx, real fa, real fb, real fc,real *alphamin);
#endif

#ifdef ACG
void acg_step(int steps);
int findalpha();
#endif

#ifdef GLOK
void update_glok(void);
void reset_glok(void);
#endif

#ifdef NMOLDYN
void init_nmoldyn(void);
void write_nmoldyn(int);
#endif

#ifdef DSF
void update_dsf(void);
void write_dsf(void);
#endif

/* Laser heating - file imd_laser.c */
#ifdef LASER
void init_laser(void);
real laser_calc_depth(cell *, int);
void laser_rescale_dummy(void);
void laser_rescale_1(void);
void laser_rescale_2(void);
void laser_rescale_3(void);
double get_surface();
double calc_laser_atom_vol(double, int, int, int *);
#ifdef LASERYZ
double laser_intensity_profile_laguerre_00(double, double, double);
double laser_intensity_profile_laguerre_01(double, double, double);
double laser_intensity_profile_laguerre_02(double, double, double);
double laser_intensity_profile_laguerre_03(double, double, double);
double laser_intensity_profile_laguerre_04(double, double, double);

double laser_intensity_profile_laguerre_10(double, double, double);
double laser_intensity_profile_laguerre_11(double, double, double);
double laser_intensity_profile_laguerre_12(double, double, double);
double laser_intensity_profile_laguerre_13(double, double, double);

double laser_intensity_profile_laguerre_20(double, double, double);
double laser_intensity_profile_laguerre_21(double, double, double);
double laser_intensity_profile_laguerre_22(double, double, double);
double laser_intensity_profile_laguerre_23(double, double, double);

double laser_intensity_profile_laguerre_30(double, double, double);
double laser_intensity_profile_laguerre_31(double, double, double);
double laser_intensity_profile_laguerre_32(double, double, double);
double laser_intensity_profile_laguerre_33(double, double, double);

double laser_intensity_profile_hermite_00(double, double, double);
double laser_intensity_profile_hermite_01(double, double, double);
double laser_intensity_profile_hermite_02(double, double, double);
double laser_intensity_profile_hermite_03(double, double, double);

double laser_intensity_profile_hermite_10(double, double, double);
double laser_intensity_profile_hermite_11(double, double, double);
double laser_intensity_profile_hermite_12(double, double, double);
double laser_intensity_profile_hermite_13(double, double, double);

double laser_intensity_profile_hermite_20(double, double, double);
double laser_intensity_profile_hermite_21(double, double, double);
double laser_intensity_profile_hermite_22(double, double, double);
double laser_intensity_profile_hermite_23(double, double, double);

double laser_intensity_profile_hermite_30(double, double, double);
double laser_intensity_profile_hermite_31(double, double, double);
double laser_intensity_profile_hermite_32(double, double, double);
double laser_intensity_profile_hermite_33(double, double, double);

#endif

#ifdef TTM
real ttm_calc_depth(int,int,int);
void laser_rescale_ttm(void);
#endif /*TTM*/
#endif /*LASER*/

#ifdef CYCLE
void init_cycle(void);
#endif 

#ifdef TTM
/* Two temperature model TTM - file imd_ttm.c */
void init_ttm(void);
void ttm_create_mpi_datatypes(void);
void ttm_fill_ghost_layers(void);
void ttm_writeout(int);
void calc_ttm(void);
void update_fd(void);
/* TODO allow variable K */
void ttm_overwrite(void);
#endif

#ifdef SM
/* Streitz and Mintmire model SM - file imd_sm.c */
void init_sm(void);
void do_electronegativity(void);
void do_v_real(void);
void do_cg(void);
void do_charge_update(void);
void charge_update_sm(void);
void calc_sm_pot(void);
void calc_sm_chi(void);
void copy_sm_charge(int, int, int, int, int, int, vektor);
void pack_sm_charge(msgbuf *, int, int, int, vektor);
void unpack_sm_charge( msgbuf *, int, int, int);
void add_sm_pot(int, int, int, int, int, int);
void pack_sm_pot(msgbuf *, int, int, int);
void unpack_add_sm_pot(msgbuf *, int, int, int);
void add_sm_chi(int, int, int, int, int, int);
void pack_sm_chi(msgbuf *, int, int, int);
void unpack_add_sm_chi(msgbuf *, int, int, int);
#endif

#ifdef BG
int get_free_mem(void);
#endif

#ifdef CBE
void mk_pt(void);
#endif

real RealGetElm( real *p, int i);
void RealSetElm( real *p, int i, real val); 
int  IntGetElm( int *p, int i);
void IntSetElm( int *p, int i, int val);

#ifdef FBC
void init_fbc(void);
void update_fbc();
#endif
#ifdef BEND
void init_bfbc(void);
void update_bfbc();
#endif

#ifdef RIGID
void calc_superforces(void);
#endif

#ifdef RELAX
void check_relaxed(void);
#endif

#ifdef TEMPCONTROL
void increment_temperature(void);
#endif

void check_write(void);
int  check_stop(void);
int  check_walltime(void);
void close_files(void);

#ifdef NEB
void read_atoms_neb(str255);
void calc_forces_neb(void);
void write_neb_eng_file(int);
void constrain_move(void);
#endif
