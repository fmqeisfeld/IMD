###########################################################################
#
#  Defaults for some variables
#
###########################################################################

MV		= mv      # program to move imd binary to ${BIN_DIR}
LIBS		+= ./nn_interpol/libnn.a -lm
MPI_FLAGS	+=# -DMPI
OMP_FLAGS	+= #-DOMP
OMPI_FLAGS	+= #-DMPI -DOMP
PACX_FLAGS	+= #-DMPI -DPACX
DEBUG_FLAGS	+= #-DDEBUG # -Wall # very noisy

# directory of the openkim API
KIM_DIR 	= /data/daniel/openkim/openkim-api

#ifeq (icc,${IMDSYS})
#  CC_SERIAL     = icc
#  CC_MPI        = mpicc
#  MPICH_CC      = icc
#  MPICH_CLINKER = icc
#  OMPI_MPICC    = icc #mpicc
#  BIN_DIR       = ${HOME}/bin/
#  OPT_FLAGS     += -g #-O3 -Wno-deprecated-declarations -funroll-loops -m64 #-vec-report 2 -ffast-math #-fopt-info-vec #-O0 -m64 -Wno-unused
#  OMPI_FLAGS    += #-openmp
#  DEBUG_FLAGS   += -g
#  PROF_FLAGS    += #-g3 -pg
#  LFLAGS        += #-openmp # -static
#  export        OMPI_MPICC MPICH_CC # MPICH_CLINKER
#endif



#ifeq (gcc,${IMDSYS})
  CC_SERIAL     = gcc
  CC_MPI        = mpicc
  MPICH_CC      = gcc
  MPICH_CLINKER = gcc
  OMPI_MPICC    = gcc #mpicc
  BIN_DIR       = ${HOME}/bin/
  OPT_FLAGS     += -g #-O2 -Wno-deprecated-declarations -funroll-loops -march=corei7-avx -mtune=corei7-avx -mavx2 -ftree-vectorize -m64 -ffast-math\
		 -ftree-vectorizer-verbose=2 -fopt-info-vec #-ffast-math #-fopt-info-vec #-O0 -m64 -Wno-unused
  OMPI_FLAGS    += #-openmp
  DEBUG_FLAGS   += -g
  PROF_FLAGS    += -g3 -pg
  LFLAGS        += # -static
  export        OMPI_MPICC MPICH_CC # MPICH_CLINKER
#endif


###########################################################################
#
#  Parallelization method
#
###########################################################################

# default is serial
PARALLEL = SERIAL
# MPI
ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
PARALLEL = MPI
PP_FLAGS += -DMPI
endif
# OpenMP
ifneq (,$(strip $(findstring omp,${MAKETARGET})))
PARALLEL = OMP
PP_FLAGS += -DOMP
endif
# MPI + OpenMP
ifneq (,$(strip $(findstring ompi,${MAKETARGET})))
PARALLEL = OMPI
PP_FLAGS += -DMPI -DOMPI
endif
# PACX
ifneq (,$(strip $(findstring pacx,${MAKETARGET})))
PARALLEL = PACX
PP_FLAGS += -DMPI -DPACX
endif
# MPIIO
ifneq (,$(strip $(findstring mpiio,${MAKETARGET})))
MPI_FLAGS += -DMPIIO
endif

###########################################################################
#
#  Compiler, flags, libraries
#
###########################################################################

# compiler; if empty, we issue an error later
CC = ${CC_${PARALLEL}}

# optimization flags
OPT_FLAGS   += ${${PARALLEL}_FLAGS} ${OPT_${PARALLEL}_FLAGS}
DEBUG_FLAGS += ${${PARALLEL}_FLAGS} ${DEBUG_${PARALLEL}_FLAGS}

# libraries
LIBS += ${${PARALLEL}_LIBS}

# optimization or debug
CFLAGS := ${FLAGS}
PP_FLAGS += ${FLAGS}
ifneq (,$(findstring debug,${MAKETARGET}))
CFLAGS += ${DEBUG_FLAGS}
else
CFLAGS += ${OPT_FLAGS}
endif

# profiling support
ifneq (,$(findstring prof,${MAKETARGET}))
CFLAGS += ${PROF_FLAGS}
LIBS   += ${PROF_LIBS}
endif

# MPE logging
ifneq (,$(findstring mpel,${MAKETARGET}))
PP_FLAGS += -DMPE
MPI_LIBS += -llmpe -lmpe
endif

# MPE tracing
ifneq (,$(findstring mpet,${MAKETARGET}))
PP_FLAGS += -DMPE
MPI_LIBS += -ltmpe -lmpe
endif

ifneq (,$(findstring debug,${MAKETARGET}))
PP_FLAGS += -DDEBUG
endif

###########################################################################
#
# IMD sources
#
###########################################################################

IMDHEADERS      = config.h globals.h imd.h makros.h potaccess.h \
                  prototypes.h types.h

SOURCES         = imd_maxwell.c imd_integrate.c imd_misc.c \
	          imd_param.c imd_alloc.c imd_io.c imd_io_3d.c \
                  imd_potential.c imd_time.c imd_generate.c \
                  imd_distrib.c imd_main_3d.c
SOURCES2D       = ${SOURCES} imd_geom_2d.c imd_pictures_2d.c
SOURCES3D       = ${SOURCES} imd_geom_3d.c imd_pictures_3d.c
RISCSOURCES2D   = imd_main_risc_2d.c
RISCSOURCES3D   = imd_main_risc_3d.c

MPISOURCES      = imd_mpi_util.c
MPISOURCES2D    = ${MPISOURCES} imd_main_mpi_2d.c imd_geom_mpi_2d.c \
		  imd_comm_force_2d.c
MPISOURCES3D    = ${MPISOURCES} imd_main_mpi_3d.c imd_geom_mpi_3d.c \
		  imd_comm_force_3d.c imd_fix_cells_3d.c imd_mpiio.c

NBLSOURCES3D    = imd_geom_mpi_3d.c imd_comm_force_3d.c imd_fix_cells_3d.c imd_mpiio.c

PAIRSOURCES     = imd_forces.c
EAM2SOURCES     = imd_forces_eam2.c
MEAMSOURCES     = imd_forces_meam.c
CGSOURCES	= imd_cg.c
COVALENTSOURCES = imd_forces_covalent.c
UNIAXSOURCES    = imd_forces_uniax.c imd_gay_berne.c
EWALDSOURCES    = imd_forces_ewald.c

CNASOURCES      = imd_cna.c

EPITAXSOURCES   = imd_epitax.c

FEFLSOURCES	= imd_fefl.c

DEFORMSOURCES   = imd_deform.c

SOCKHEADERS     = sockets.h sockutil.h socket_io.h
SOCKSOURCES     = socket_io.c sockutil.c

QUASISOURCES    = imd_qc.c

CORRSOURCES     = imd_correl.c

TRANSSOURCES    = imd_transport.c

LASERSOURCES	= imd_laser.c

TTMSOURCES	= imd_ttm.c imd_interpol.c fminbnd3.c

FDTDSOURCES	= imd_fdtd.c

TMMSOURCES	= imd_tmm.c

COLRADSOURCES   = imd_colrad.c

NRBSOURCES	= imd_nrb.c

SMSOURCES	= imd_sm.c

FCSSOURCES      = imd_forces_fcs.c

#BBOOSTSOURCES	= imd_bboost.c imd_bb_core1.c imd_bb_core2.c
BBOOSTSOURCES	= imd_bboost.c

#########################################################
#
# IMD Configuration rules
#
#########################################################

HEADERS := ${IMDHEADERS}

# twod or not twod
ifneq (,$(strip $(findstring 2d,${MAKETARGET})))

  # 2d, serial or mpi
  ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
    SOURCES	:= ${SOURCES2D} ${MPISOURCES2D}
    else
    SOURCES	:= ${SOURCES2D} ${RISCSOURCES2D}
  endif
  PP_FLAGS  += -DTWOD

else

  # 3d, vec, serial, mpi
  ifneq (,$(strip $(findstring vec,${MAKETARGET})))

    SOURCES	:= ${SOURCES3D} ${NBLSOURCES3D}
    ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
      SOURCES	+= ${MPISOURCES}
    endif

  else

    ifneq (,$(strip $(findstring cbe,${MAKETARGET})))

      SOURCES	:= ${SOURCES3D} ${NBLSOURCES3D}
      ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
        SOURCES	+= ${MPISOURCES}
      endif

    else

      ifneq (,$(strip $(findstring nbl,${MAKETARGET})))

        SOURCES	:= ${SOURCES3D} ${NBLSOURCES3D}
        ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
          SOURCES += ${MPISOURCES}
        endif

      else

	ifneq (,$(strip $(findstring kim, ${MAKETARGET})))
	  SOURCES := ${SOURCES3D} ${NBLSOURCES3D}
	  ifneq (,$(strip $(findstring mpi, ${MAKETARGET})))
	    SOURCES += ${MPISOURCES}
          endif
	else
          ifneq (,$(strip $(findstring mpi,${MAKETARGET})))
            SOURCES := ${SOURCES3D} ${MPISOURCES3D}
          else
            SOURCES := ${SOURCES3D} ${RISCSOURCES3D}
          endif
        endif
      endif
    endif
  endif
endif

SRCMAIN = imd.c
ifneq (,$(strip $(findstring py,${MAKETARGET})))
  SRCMAIN  = py_imd.c
  CFLAGS  += -fPIC
endif
ifneq (,$(strip $(findstring jvis,${MAKETARGET})))
  SRCMAIN   = jvis_imd.c
  PP_FLAGS += -DNVE -DNVT -DNPT -DNPT_iso -DREFPOS
endif
SOURCES += ${SRCMAIN}

###  INTERACTION  #######################################

# pair interaction is the default
ifneq (,$(strip $(findstring vec,${MAKETARGET})))
  FORCESOURCES = imd_main_vec_3d.c
else
  ifneq (,$(strip $(findstring cbe,${MAKETARGET})))
    FORCESOURCES = imd_forces_cbe.c
    PP_FLAGS += -DCBE
  else
    ifneq (,$(strip $(findstring nbl,${MAKETARGET})))
      FORCESOURCES = imd_forces_nbl.c
      PP_FLAGS += -DNBL
      ifneq (,$(strip $(findstring kim, ${MAKETARGET})))
	ERROR = "KIM is not compatible with any NBL force routine"
      endif
    else
      FORCESOURCES = ${PAIRSOURCES}
    endif
  endif
endif

# PAIR
ifneq (,$(strip $(findstring pair,${MAKETARGET})))
  ifneq (,$(strip $(findstring kim,${MAKETARGET})))
    ERROR = "KIM is not compatible with PAIR force routines"
  endif
  PP_FLAGS  += -DPAIR
endif

# EAM2 or EAM  -  this is now the same
ifneq (,$(strip $(findstring eam,${MAKETARGET})))
  ifneq (,$(strip $(findstring kim,${MAKETARGET})))
    ERROR = "KIM is not compatible with EAM force routines"
  endif
  # EEAM
  ifneq (,$(strip $(findstring eeam,${MAKETARGET})))
    PP_FLAGS  += -DEEAM
  endif
  # MEAM
  ifneq (,$(strip $(findstring meam,${MAKETARGET})))
    PP_FLAGS  += -DMEAM
    FORCESOURCES += ${MEAMSOURCES}  ${COVALENTSOURCES}
  else
    PP_FLAGS  += -DEAM2
    ifneq (,$(strip $(findstring asympot,${MAKETARGET})))
      PP_FLAGS += -DASYMPOT
      FORCESOURCES += ${EAM2SOURCES}
    else
      ifeq (,$(strip $(findstring nbl,${MAKETARGET})))
        ifeq (,$(strip $(findstring vec,${MAKETARGET})))
          FORCESOURCES += ${EAM2SOURCES}
        endif
      endif
    endif
  endif
endif

# ADP
ifneq (,$(strip $(findstring adp,${MAKETARGET})))
  ifneq (,$(strip $(findstring kim,${MAKETARGET})))
    ERROR = "KIM is not compatible with ADP force routines"
  endif
PP_FLAGS  += -DADP
endif

# TTBP
ifneq (,$(strip $(findstring ttbp,${MAKETARGET})))
PP_FLAGS  += -DTTBP
FORCESOURCES += ${COVALENTSOURCES}
endif

# TTBPXT
ifneq (,$(strip $(findstring ttbpxt,${MAKETARGET})))
PP_FLAGS  += -DTTBP -DXT
FORCESOURCES += ${COVALENTSOURCES}
endif

# STIWEB
ifneq (,$(strip $(findstring stiweb,${MAKETARGET})))
  ifneq (,$(strip $(findstring kim,${MAKETARGET})))
    ERROR = "KIM is not compatible with STIWEB force routines"
  endif
PP_FLAGS  += -DSTIWEB
FORCESOURCES += ${COVALENTSOURCES}
endif

# TERNBCC
ifneq (,$(strip $(findstring ternbcc,${MAKETARGET})))
PP_FLAGS  += -DSTIWEB -DTERNBCC
FORCESOURCES += ${COVALENTSOURCES}
endif

# TERSOFF
ifneq (,$(strip $(findstring tersoffmod,${MAKETARGET})))
  ifneq (,$(strip $(findstring kim,${MAKETARGET})))
    ERROR = "KIM is not compatible with TERSOFFMOD force routines"
  endif
PP_FLAGS  += -DTERSOFFMOD
FORCESOURCES += ${COVALENTSOURCES}
ifneq (,$(strip $(findstring tersoffmod2,${MAKETARGET})))
PP_FLAGS  += -DTERSOFFMOD2
endif
else
ifneq (,$(strip $(findstring tersoff,${MAKETARGET})))
  ifneq (,$(strip $(findstring kim,${MAKETARGET})))
    ERROR = "KIM is not compatible with TERSOFF force routines"
  endif
PP_FLAGS  += -DTERSOFF
FORCESOURCES += ${COVALENTSOURCES}
ifneq (,$(strip $(findstring tersoff2,${MAKETARGET})))
PP_FLAGS  += -DTERSOFF2
endif
endif
endif

# BRENNER
ifneq (,$(strip $(findstring brenner,${MAKETARGET})))
PP_FLAGS  += -DBRENNER
FORCESOURCES += ${COVALENTSOURCES}
ifneq (,$(strip $(findstring brenner2,${MAKETARGET})))
PP_FLAGS  += -DBRENNER2
endif
endif

# KEATING
ifneq (,$(strip $(findstring keating,${MAKETARGET})))
PP_FLAGS  += -DKEATING
FORCESOURCES += ${COVALENTSOURCES}
endif

# LASER & LASERYZ
ifneq (,$(strip $(findstring laseryz,${MAKETARGET})))
        ifneq (,$(strip $(findstring laser,${MAKETARGET})))
		CFLAGS += -DLASERYZ -DLASER
		SOURCES += ${LASERSOURCES}
	else
	       	CFLAGS += -DLASERYZ -DLASER
		SOURCES += ${LASERSOURCES}
	endif
else
	ifneq (,$(strip $(findstring laser,${MAKETARGET})))
		PP_FLAGS += -DLASER
		SOURCES += ${LASERSOURCES}
	endif
endif

# PDECAY
ifneq (,$(findstring pdecay,${MAKETARGET}))
PP_FLAGS += -DPDECAY
endif

# TWO TEMPERATURE MODEL TTM
ifneq (,$(strip $(findstring ttm,${MAKETARGET})))
PP_FLAGS += -DTTM
SOURCES += ${TTMSOURCES}
endif


#MYMOD
# FDTD
ifneq (,$(strip $(findstring fdtd,${MAKETARGET})))
PP_FLAGS += -DFDTD
SOURCES += ${FDTDSOURCES}
endif
ifneq (,$(strip $(findstring dirichlet,${MAKETARGET})))
PP_FLAGS += -DDIRICHLET
endif

#TMM
ifneq (,$(strip $(findstring tmm,${MAKETARGET})))
PP_FLAGS += -DTMM
SOURCES += ${TMMSOURCES}
endif

#COLRAD
ifneq (,$(strip $(findstring colrad,${MAKETARGET})))
PP_FLAGS += -DCOLRAD
SOURCES += ${COLRADSOURCES}
#LIBS          += -L${KIM_DIR}/KIM_API -lkim
#CFLAGS        += -I${KIM_DIR}/KIM_API

LIBS          += -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsollapackdense -lgsl -lgslcblas
#CFLAGS        += -fopenmp

endif

# SIMD --> d.h. vektorisierung der loop in imd_forces.c
ifneq (,$(strip $(findstring simd,${MAKETARGET})))
PP_FLAGS += -DSIMD
SOURCES += imd_forces_simd.c
endif


#local order parameter
ifneq (,$(strip $(findstring lod,${MAKETARGET})))
PP_FLAGS += -DLOD
endif

#non-reflecting boundary conditions
ifneq (,$(strip $(findstring nrb,${MAKETARGET})))
PP_FLAGS += -DNRB
#PP_FLAGS += -DNBL
#PP_FLAGS += -DNBLIST
PP_FLAGS += -DREFPOS
SOURCES += ${NRBSOURCES}
endif

#endof MYMOD

# STREITZ MINTMIRE MODEL SM
ifneq (,$(strip $(findstring sm,${MAKETARGET})))
  PP_FLAGS += -DSM
  SOURCES += ${SMSOURCES}
  ifeq (,$(strip $(findstring nbl,${MAKETARGET})))
    SOURCES += ${EWALDSOURCES}
  endif
endif

# UNIAX
ifneq (,$(strip $(findstring uniax,${MAKETARGET})))
PP_FLAGS  += -DUNIAX
FORCESOURCES = ${UNIAXSOURCES}
endif

# EWALD
ifneq (,$(strip $(findstring ewald,${MAKETARGET})))
PP_FLAGS  += -DEWALD
FORCESOURCES  += ${EWALDSOURCES}
endif

# FCS
ifneq (,$(strip $(findstring fcs,${MAKETARGET})))
PP_FLAGS      += -DUSEFCS
FORCESOURCES  += ${FCSSOURCES}
LIBS          += $(shell pkg-config --libs scafacos-fcs)
CFLAGS        += $(shell pkg-config --cflags scafacos-fcs)
endif

# VARCHG
ifneq (,$(strip $(findstring varchg,${MAKETARGET})))
CFLAGS  += -DVARCHG
endif

# DIPOLE
ifneq (,$(strip $(findstring dipole,${MAKETARGET})))
CFLAGS  += -DDIPOLE
endif

ifneq (,$(strip $(findstring kermode,${MAKETARGET})))
CFLAGS  += -DKERMODE
endif

# BUCK
ifneq (,$(strip $(findstring buck,${MAKETARGET})))
CFLAGS  += -DBUCK
endif

# MORSE
ifneq (,$(strip $(findstring morse,${MAKETARGET})))
CFLAGS  += -DMORSE
endif

# EXTF
ifneq (,$(strip $(findstring extf,${MAKETARGET})))
CFLAGS  += -DEXTF
endif

# LJ - computed Lennard-Jones (for vector versions only)
ifneq (,$(strip $(findstring lj,${MAKETARGET})))
PP_FLAGS  += -DLJ
endif

ifneq (,$(strip $(findstring kim,${MAKETARGET})))
FORCESOURCES = imd_forces_kim.c
PP_FLAGS += -DKIM
LIBS          += -L${KIM_DIR}/KIM_API -lkim
CFLAGS        += -I${KIM_DIR}/KIM_API
endif

SOURCES += ${FORCESOURCES}

###  ENSEMBLES  #########################################

ifneq (,$(findstring nve,${MAKETARGET}))
PP_FLAGS += -DNVE
endif

ifneq (,$(findstring mik,${MAKETARGET}))
PP_FLAGS += -DMIK
endif


# CG
ifneq (,$(strip $(findstring cg,${MAKETARGET})))
SOURCES += ${CGSOURCES}
IMDHEADERS += cg_util.h
PP_FLAGS  += -DCG
endif
ifneq (,$(strip $(findstring acg,${MAKETARGET})))
PP_FLAGS  += -DACG
endif

ifneq (,$(findstring nvt,${MAKETARGET}))
PP_FLAGS += -DNVT
endif

ifneq (,$(findstring npt_iso,${MAKETARGET}))
PP_FLAGS += -DNPT -DNPT_iso
endif

ifneq (,$(findstring npt_axial,${MAKETARGET}))
PP_FLAGS += -DNPT -DNPT_axial
endif

ifneq (,$(findstring frac,${MAKETARGET}))
PP_FLAGS += -DFRAC
# SOURCES += ${FRACSOURCES}
endif

ifneq (,$(findstring damp,${MAKETARGET}))
PP_FLAGS += -DDAMP
endif

ifneq (,$(findstring ftg,${MAKETARGET}))
PP_FLAGS += -DFTG
endif

ifneq (,$(findstring finnis,${MAKETARGET}))
PP_FLAGS += -DFINNIS
endif

ifneq (,$(findstring stm,${MAKETARGET}))
PP_FLAGS += -DSTM
endif

###  OPTIONS  ############################################

# vector mode
ifneq (,$(findstring vec,${MAKETARGET}))
PP_FLAGS += -DVEC

ifneq (,$(findstring vec2,${MAKETARGET}))
PP_FLAGS += -DVEC2
else
ifneq (,$(findstring vec3,${MAKETARGET}))
PP_FLAGS += -DVEC3
endif
endif

endif

# nudged elastic band method
ifneq (,$(findstring neb,${MAKETARGET}))
SOURCES  += imd_neb.c
PP_FLAGS += -DNEB
CC        = ${CC_MPI}
endif

# common-neighbour analysis
ifneq (,$(findstring cna,${MAKETARGET}))
SOURCES += ${CNASOURCES}
SOURCES += ${COVALENTSOURCES}
PP_FLAGS += -DCNA
endif

# socket interface
ifneq (,$(strip $(findstring sock,${MAKETARGET})))
HEADERS += ${SOCKHEADERS}
SOURCES += ${SOCKSOURCES}
PP_FLAGS  += -DSOCKET_IO
endif

# timing
ifneq (,$(findstring timing,${MAKETARGET}))
PP_FLAGS += -DTIMING
endif

ifneq (,$(findstring and,${MAKETARGET}))
PP_FLAGS += -DAND
endif

ifneq (,$(findstring ber,${MAKETARGET}))
PP_FLAGS += -DBER
endif

ifneq (,$(findstring fbc,${MAKETARGET}))
PP_FLAGS += -DFBC
endif

ifneq (,$(findstring bend,${MAKETARGET}))
PP_FLAGS += -DBEND
PP_FLAGS += -DFBC
endif

ifneq (,$(findstring bboost,${MAKETARGET}))
PP_FLAGS += -DBBOOST
SOURCES += ${BBOOSTSOURCES}
endif

ifneq (,$(findstring debugLo,${MAKETARGET}))
PP_FLAGS += -DdebugLo
endif


ifneq (,$(findstring zapp,${MAKETARGET}))
PP_FLAGS += -DZAPP
endif

ifneq (,$(findstring flagatoms,${MAKETARGET}))
PP_FLAGS += -DFLAGEDATOMS
endif


ifneq (,$(findstring rigid,${MAKETARGET}))
PP_FLAGS += -DRIGID
endif

ifneq (,$(findstring sendrec,${MAKETARGET}))
PP_FLAGS += -DSR
endif

ifneq (,$(findstring einstein,${MAKETARGET}))
PP_FLAGS += -DEINSTEIN
endif

ifneq (,$(findstring fnorm,${MAKETARGET}))
PP_FLAGS += -DFNORM
endif

ifneq (,$(findstring getsaddle,${MAKETARGET}))
PP_FLAGS += -DGETSADDLE
endif
ifneq (,$(findstring getmin,${MAKETARGET}))
PP_FLAGS += -DGETMIN
endif


ifneq (,$(findstring relaxinfo,${MAKETARGET}))
PP_FLAGS += -DRELAXINFO
endif

ifneq (,$(findstring norhoh,${MAKETARGET}))
PP_FLAGS += -DNORHOH
endif

ifneq (,$(findstring glok,${MAKETARGET}))
PP_FLAGS += -DGLOK
endif

ifneq (,$(findstring adaptglok,${MAKETARGET}))
PP_FLAGS += -DADAPTGLOK
endif

ifneq (,$(findstring mix,${MAKETARGET}))
PP_FLAGS += -DMIX
endif

ifneq (,$(findstring fire,${MAKETARGET}))
PP_FLAGS += -DMIK
PP_FLAGS += -DMIX
PP_FLAGS += -DADAPTGLOK
PP_FLAGS += -DGLOK
PP_FLAGS += -DRELAXINFO
PP_FLAGS += -DFNORM
endif

ifneq (,$(findstring efilter,${MAKETARGET}))
PP_FLAGS += -DEFILTER
endif

ifneq (,$(findstring clone,${MAKETARGET}))
PP_FLAGS += -DCLONE
endif

ifneq (,$(findstring nnbr,${MAKETARGET}))
PP_FLAGS += -DNNBR
endif

ifneq (,$(findstring writef,${MAKETARGET}))
PP_FLAGS += -DWRITEF
endif

ifneq (,$(findstring deform,${MAKETARGET}))
PP_FLAGS += -DDEFORM
ifeq (,$(findstring homdef,${MAKETARGET}))
SOURCES += ${DEFORMSOURCES}
endif
endif

ifneq (,$(findstring cycle,${MAKETARGET}))
PP_FLAGS += -DCYCLE
PP_FLAGS += -DHOMDEF
SOURCES += ${DEFORMSOURCES}
endif
ifneq (,$(findstring homdef,${MAKETARGET}))
PP_FLAGS += -DHOMDEF
SOURCES += ${DEFORMSOURCES}
endif

ifneq (,$(findstring shock,${MAKETARGET}))
PP_FLAGS += -DSHOCK
endif

ifneq (,$(findstring stress,${MAKETARGET}))
PP_FLAGS += -DSTRESS_TENS
endif

ifneq (,$(findstring quasi,${MAKETARGET}))
PP_FLAGS += -DQUASI
SOURCES += ${QUASISOURCES}
endif

ifneq (,$(findstring disloc,${MAKETARGET}))
PP_FLAGS += -DDISLOC
endif

ifneq (,$(findstring sllod,${MAKETARGET}))
PP_FLAGS += -DSLLOD -DNVT
endif

ifneq (,$(findstring avpos,${MAKETARGET}))
PP_FLAGS += -DAVPOS
endif

ifneq (,$(findstring force,${MAKETARGET}))
PP_FLAGS += -DFORCE
endif

####################
# MY MOD
###################
ifneq (,$(findstring multijump,${MAKETARGET}))
PP_FLAGS += -DMULTIJUMP
endif

ifneq (,$(findstring filter,${MAKETARGET}))
PP_FLAGS += -DFILTER
SOURCES += imd_filter.c
endif

ifneq (,$(findstring filterxy,${MAKETARGET}))
PP_FLAGS += -DFILTERY
endif

ifneq (,$(findstring filterxyz,${MAKETARGET}))
PP_FLAGS += -DFILTERZ
endif




################
# END OF MODS
################

ifneq (,$(findstring nmoldyn,${MAKETARGET}))
PP_FLAGS += -DNMOLDYN
endif

ifneq (,$(findstring dsf,${MAKETARGET}))
PP_FLAGS += -DDSF
LFLAGS += -lfftw3
endif

ifneq (,$(findstring atdist,${MAKETARGET}))
PP_FLAGS += -DATDIST
endif

ifneq (,$(findstring diffpat,${MAKETARGET}))
PP_FLAGS += -DDIFFPAT -I ${FFTW_DIR}/include
ifneq (,$(findstring omp,${MAKETARGET}))
LIBS   += -L ${FFTW_DIR}/lib -lfftw3f_threads -lfftw3f -lgomp
else
LIBS   += -L ${FFTW_DIR}/lib -lfftw3f
endif
endif

ifneq (,$(findstring ordpar,${MAKETARGET}))
PP_FLAGS += -DORDPAR
endif

# EPITAX
ifneq (,$(strip $(findstring epitax,${MAKETARGET})))
PP_FLAGS  += -DEPITAX
SOURCES += ${EPITAXSOURCES}
endif

# FEFL
ifneq (,$(findstring fefl,${MAKETARGET}))
PP_FLAGS += -DFEFL
SOURCES += ${FEFLSOURCES}
endif

# Correlation
ifneq (,$(strip $(findstring corr,${MAKETARGET})))
PP_FLAGS  += -DCORRELATE
SOURCES += ${CORRSOURCES}
endif

# heat transport
ifneq (,$(findstring nvx,${MAKETARGET}))
  PP_FLAGS += -DNVX
  SOURCES  += ${TRANSSOURCES}
else
  ifneq (,$(strip $(findstring hc,${MAKETARGET})))
    PP_FLAGS += -DHC
    SOURCES  += ${TRANSSOURCES}
  endif
endif

# mean square displacement
ifneq (,$(strip $(findstring msqd,${MAKETARGET})))
PP_FLAGS  += -DMSQD
SOURCES += ${CORRSOURCES}
endif

# MONOLJ Case
ifneq (,$(findstring monolj,${MAKETARGET}))
PP_FLAGS += -DMONOLJ
endif

# Single precision
ifneq (,$(findstring single,${MAKETARGET}))
PP_FLAGS += -DSINGLE
endif

# monoatomic system (performance tweak)
ifneq (,$(findstring mono,${MAKETARGET}))
PP_FLAGS += -DMONO
endif

# high precision output (checkpoints)
ifneq (,$(findstring hpo,${MAKETARGET}))
PP_FLAGS += -DHPO
endif

# use 4-point-interpolation
ifneq (,$(findstring 4point,${MAKETARGET}))
PP_FLAGS += -DFOURPOINT
endif

# use spline interpolation
ifneq (,$(findstring spline,${MAKETARGET}))
PP_FLAGS += -DSPLINE
endif

# use papi
ifneq (,$(findstring papi,${MAKETARGET}))
PP_FLAGS += -DPAPI ${PAPI_INC}
LIBS   += ${PAPI_LIBS}
endif

# use external potential for indenters
ifneq (,$(findstring extpot,${MAKETARGET}))
CFLAGS += -DEXTPOT
SOURCES += imd_extpot.c
endif

# directory for python modules
ifeq (,${IMDPYDIR})
  PYDIR = ${HOME}/python/imd
else
  PYDIR = ${IMDPYDIR}
endif

# Added by Frank Pister
ifneq (,$(findstring cbe,${MAKETARGET}))
SOURCES += spu.c imd_cbe_calc_ppu.c imd_cbe_util.c
endif

ifneq (,$(findstring cbe2,${MAKETARGET}))
PP_FLAGS  += -DCBE2
endif

ifneq (,$(findstring ppu,${MAKETARGET}))
PP_FLAGS += -DON_PPU
endif

#Angle distribution & slip vector analysis
ifneq (,$(findstring ada,${MAKETARGET}))
PP_FLAGS += -DADA
SOURCES += imd_ada.c
ifeq (,$(findstring ${COVALENTSOURCES},${SOURCES}))
	SOURCES += ${COVALENTSOURCES}
endif
endif

#Viscous damping
ifneq (,$(findstring viscous,${MAKETARGET}))
PP_FLAGS += -DVISCOUS
endif
#Langevin thermostat
ifneq (,$(findstring langevin,${MAKETARGET}))
PP_FLAGS += -DLANGEVIN -DVISCOUS
endif

ifneq (,$(findstring nye,${MAKETARGET}))
PP_FLAGS += -DNYETENSOR
SOURCES += imd_nyeTensorAnalysis_3d.c
ifeq (,$(findstring ${COVALENTSOURCES},${SOURCES}))
	SOURCES += ${COVALENTSOURCES}
endif
ifeq (,$(findstring imd_ada.c,${SOURCES}))
	SOURCES += imd_ada.c
	PP_FLAGS += -DADA
endif
endif

#Load Balancer, only works for MPI
ifneq (,$(findstring loadbalance,${MAKETARGET}))
ifneq (,$(findstring mpi,${MAKETARGET}))
	PP_FLAGS += -DLOADBALANCE
	SOURCES += imd_loadBalance.c imd_loadBalance_direct.c
endif
endif

# processor affinity
ifneq (,$(findstring aff,${MAKETARGET}))
PP_FLAGS += -DAFF
endif


# Substitute .o for .c to get the names of the object files
OBJECTS := $(subst .c,.o,${SOURCES})


###########################################################################
#
#	 Rules
#
###########################################################################

# all objects depend on headers
${OBJECTS}: ${HEADERS}

# How to compile *.c files
# special rules for force computation
imd_forces.o: imd_forces.c
	${CC} ${CFLAGS} ${PP_FLAGS} ${RCD_FLAGS} -c imd_forces.c

imd_forces_nbl.o: imd_forces_nbl.c
	${CC} ${CFLAGS} ${PP_FLAGS} ${RCD_FLAGS} ${NOALIAS} -c imd_forces_nbl.c

# Uncommented by Frank Pister
imd_forces_cbe.o: imd_forces_cbe.c
	${CC}   ${CFLAGS} ${PP_FLAGS} -c imd_forces_cbe.c

# rule to make imd_forces_cbe.o

imd_forces_eam2.o: imd_forces_eam2.c
	${CC} ${CFLAGS} ${PP_FLAGS} ${RCD_FLAGS} -c imd_forces_eam2.c

imd_forces_covalent.o: imd_forces_covalent.c
	${CC} ${CFLAGS} ${PP_FLAGS} ${RCD_FLAGS} -c imd_forces_covalent.c

# generic compilation rule
.c.o:
	${CC} ${CFLAGS} ${PP_FLAGS} -c $<

# How to link
ifneq (,$(strip $(findstring py,${MAKETARGET})))
${MAKETARGET}: ${OBJECTS}
	@echo ${PYDIR} ${IMDPYDIR}
	swig ${PP_FLAGS} -python -module $@ -outdir ${PYDIR} imd.i
	${CC} ${CFLAGS} ${PP_FLAGS} -I/usr/include/python -c imd_wrap.c
	${CC} ${LFLAGS} -fPIC -shared -o _$@.so ${OBJECTS} imd_wrap.o ${LIBS}
	${MV} _$@.so ${PYDIR}; rm -f _$@.so
else
${MAKETARGET}: ${OBJECTS}
	${CC} ${LFLAGS} -o $@ ${OBJECTS} ${LIBS}
	${MV} $@ ~/bin/; rm -f $@
endif

# First recursion only set the MAKETARGET Variable
.DEFAULT:
ifneq (,${CC})
	./version.sh
	@echo \#define COMPILE_TARGET \"$@\" >> version.h
	${MAKE} MAKETARGET='$@' STAGE2
else
ifneq (,${IMDSYS})
	@echo "IMDSYS variable ${IMDSYS} is not recognized"
else
	@echo "IMDSYS variable is not set"
endif
endif

# Second recursion sets MAKETARGET variable and compiles
# An empty MAKETARGET variable would create an infinite recursion, so we check
STAGE2:
ifneq (,${ERROR})
	@echo -e "\nError: ${ERROR}\n"
else
ifneq (,${MAKETARGET})
	${MAKE} MAKETARGET='${MAKETARGET}' ${MAKETARGET}
else
	@echo 'No TARGET specified.'
endif
endif

###########################################################################
#
#	 Misc. TARGETs
#
###########################################################################

clean: clean_spu
	rm -f *.o *.u *~ \#* *.V *.T *.O *.il

# Remove SPU-related temp. files (e.g. assembler/timing files)
clean_spu:
	rm -f  spu spu1 spu2   spu.map   spu.s  imd_cbe_util_spu.s imd_cbe_calc_spu.s     spu.s.timing imd_cbe_util_spu.s.timing imd_cbe_calc_spu.s.timing *.lst



help:
	@echo "Usage: gmake imd[_<parallel>][_<option>[_<option>...]]"

socktest:
	gcc -o ${BINDIR}/socktest sockutil.c socktest.c









# SPU specific targets for CBE
# The variables SPUCXX, SPUCC must be set to the SPU C++ and C compilers
# e.g. to spu-c++ and spu-cc respectively.
# SPUCXXFLAGS and SPUCFLAGS are the corresponding flags passed to the
# compiler which may contain optimization or warning flags, for instance.
# SPUCXX and SPUCC may also be set to use the IBM  compilers.
# EMBEDSPU must contain the name of the ppu-embedspu executable which
# generates object files containing code for the SPU and which may
# be linked with the main program.



# The SPU binary containing the main program as well as some tools
# and an spu version of calc_wp
# The SPU object file which is linked to the main programm
spu.o: spu1
	${EMBEDSPU} hndle_cbe_calc spu1 spu.o



# 2 SPUlets created in a different manner

# Compile first translate C to assembler, then compile assembler sources
spu1:  spu.s  imd_cbe_util_spu.s  imd_cbe_calc_spu.s
	${SPUCC} -o spu1   ${SPUCFLAGS} ${SPULDFLAGS}     spu.s  imd_cbe_calc_spu.s  imd_cbe_util_spu.s

# Compile the C sources without generating assembler
spu2: spu.c imd_cbe_util.c imd_cbe_calc_spu.c imd_cbe.h  config.h
	${SPUCC} -o spu2 ${SPUCFLAGS} ${SPULDFLAGS} ${PP_FLAGS}   spu.c imd_cbe_util.c imd_cbe_calc_spu.c


# Generate sssembler code for some SPU sources
# This assembler code may be analyzed using the spu-timing tool
imd_cbe_calc_spu.s: imd_cbe_calc_spu.c imd_cbe.h  config.h
	${SPUCC} -o imd_cbe_calc_spu.s  ${SPUCFLAGS} ${PP_FLAGS}   -S imd_cbe_calc_spu.c

imd_cbe_util_spu.s: imd_cbe_util.c imd_cbe.h  config.h
	${SPUCC} -o imd_cbe_util_spu.s  ${SPUCFLAGS} ${PP_FLAGS}    -S imd_cbe_util.c

spu.s: spu.c imd_cbe.h  config.h
	${SPUCC} -o spu.s ${SPUCFLAGS} ${PP_FLAGS}  -S  spu.c


# Timing of SPU assembly code (for profiling purposes)
spu.s.timing: spu.s
	${SPUTIMING} spu.s

imd_cbe_calc_spu.s.timing: imd_cbe_calc_spu.s
	${SPUTIMING} imd_cbe_calc_spu.s

imd_cbe_util_spu.s.timing: imd_cbe_util_spu.s
	${SPUTIMING} imd_cbe_util_spu.s


spu_timing: spu.s.timing  imd_cbe_calc_spu.s.timing  imd_cbe_util_spu.s.timing
