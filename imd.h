/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd.h -- Header file for all modules of IMD
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/* C stuff */
#ifdef MEMALIGN
#define  _XOPEN_SOURCE 600
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>

//MYMOD
#include <stdbool.h>
#include "nn_interpol/nn.h"
#include "nn_interpol/delaunay.h"
//ENDOF MYMOD
//
//
#ifdef CBE
#define USE_WALLTIME
#endif

/* support for timers */
#ifndef MPI
#if defined(USE_RUSAGE) || defined(USE_WALLTIME)
#include <sys/time.h>
#include <sys/resource.h>
#else
#include <sys/times.h>
#include <sys/types.h>
#endif
#endif

/* Machine specific headers */
#if defined(MPI) || defined(NEB)
#include <mpi.h>
#ifdef MPE
#include <mpe.h>
#endif
#endif
#ifdef OMP
#include <omp.h>
#endif

/* FFT for diffraction patterns */
#ifdef DIFFPAT
#include <fftw3.h>
#endif

/* IMD version */
#include "version.h"

/* Configuration */
#include "config.h"

/* Data types */
#include "types.h"

/* Some constants */
#include "constants.h"

/* Some makros */
#include "makros.h"

/* Function Prototypes */
#include "prototypes.h"

/* Global Variables */
#include "globals.h"

/************
 * MY MOD   *
 ************/
#include <complex.h>

#ifdef COLRAD

#include <gsl/gsl_sf_expint.h>       //FUER EXPONENTIAL INTEGRAL
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sunmatrix/sunmatrix_dense.h>
//#include <sundials/sundials_types.h> /* definition of type realtype */
#include "sunnonlinsol/sunnonlinsol_newton.h"
#define LAPACK  /* LAPACK SOLVER MULTO BENE*/
#ifdef LAPACK
#include <sunlinsol/sunlinsol_lapackdense.h>
#else
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#endif

#endif //ifdef COLRAD
