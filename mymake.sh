#!/bin/bash


# make clean
# make -j8 imd_mpi_eam_ttm_tmm_nbl_colrad


# #CUSTOM MAKE MIT OPENMP ABER NUR FUER COLRAD
 # optflags=" -funroll-loops -march=corei7-avx -mtune=corei7-avx -mavx2 -ftree-vectorize -m64 -ffast-math"

 #INTEL LINK OPTINOS
 libintel=" -L/opt/intel/mkl/lib/intel64-Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl"
 cflagintel=" -m64 -I/opt/intel//include "
 make clean
 mpicc  $cflagintel -O3 -DMPI  -DNBL -DEAM2 -DTTM -DTMM -DCOLRAD -c imd_maxwell.c imd_misc.c imd_param.c imd_alloc.c imd_io.c imd_io_3d.c imd_potential.c imd_time.c imd_generate.c imd_distrib.c imd_main_3d.c imd_geom_3d.c imd_pictures_3d.c imd_geom_mpi_3d.c imd_comm_force_3d.c imd_fix_cells_3d.c imd_mpiio.c imd_mpi_util.c imd.c imd_ttm.c imd_interpol.c fminbnd3.c imd_tmm.c imd_colrad.c imd_forces_nbl.c imd_integrate.c -fopenmp
 #mpicc  -o imd_mpi_eam_ttm_tmm_nbl_colrad imd_maxwell.o imd_integrate.o imd_misc.o imd_param.o imd_alloc.o imd_io.o imd_io_3d.o imd_potential.o imd_time.o imd_generate.o imd_distrib.o imd_main_3d.o imd_geom_3d.o imd_pictures_3d.o imd_geom_mpi_3d.o imd_comm_force_3d.o imd_fix_cells_3d.o imd_mpiio.o imd_mpi_util.o imd.o imd_ttm.o imd_interpol.o fminbnd3.o imd_tmm.o imd_colrad.o imd_forces_nbl.o ./nn_interpol/libnn.a -lm    -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsollapackdense -lgsl -lgslcblas -fopenmp
 mpicc  -o imd_mpi_eam_ttm_tmm_nbl_colrad imd_maxwell.o imd_integrate.o imd_misc.o imd_param.o imd_alloc.o \
 			imd_io.o imd_io_3d.o imd_potential.o imd_time.o imd_generate.o imd_distrib.o imd_main_3d.o imd_geom_3d.o \
 			imd_pictures_3d.o imd_geom_mpi_3d.o imd_comm_force_3d.o imd_fix_cells_3d.o imd_mpiio.o imd_mpi_util.o\
 			 imd.o imd_ttm.o imd_interpol.o fminbnd3.o imd_tmm.o imd_colrad.o imd_forces_nbl.o \
 			 ./nn_interpol/libnn.a \
 			 -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsollapackdense -lgsl -lgslcblas\
 			 -L/opt/intel/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
 mv imd_mpi_eam_ttm_tmm_nbl_colrad ~/bin/





#profile
 # make clean
 # mpicc  -pg -g -DMPI  -DNBL -DEAM2 -DTTM -DTMM -DCOLRAD -c imd_maxwell.c imd_misc.c imd_param.c imd_alloc.c imd_io.c imd_io_3d.c imd_potential.c imd_time.c imd_generate.c imd_distrib.c imd_main_3d.c imd_geom_3d.c imd_pictures_3d.c imd_geom_mpi_3d.c imd_comm_force_3d.c imd_fix_cells_3d.c imd_mpiio.c imd_mpi_util.c imd.c imd_ttm.c imd_interpol.c fminbnd3.c imd_tmm.c imd_colrad.c imd_forces_nbl.c imd_integrate.c -fopenmp
 # mpicc  -pg -g -o imd_mpi_eam_ttm_tmm_nbl_colrad imd_maxwell.o imd_integrate.o imd_misc.o imd_param.o imd_alloc.o imd_io.o imd_io_3d.o imd_potential.o imd_time.o imd_generate.o imd_distrib.o imd_main_3d.o imd_geom_3d.o imd_pictures_3d.o imd_geom_mpi_3d.o imd_comm_force_3d.o imd_fix_cells_3d.o imd_mpiio.o imd_mpi_util.o imd.o imd_ttm.o imd_interpol.o fminbnd3.o imd_tmm.o imd_colrad.o imd_forces_nbl.o ./nn_interpol/libnn.a -lm    -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsollapackdense -lgsl -lgslcblas -fopenmp
 # mv imd_mpi_eam_ttm_tmm_nbl_colrad ~/bin/