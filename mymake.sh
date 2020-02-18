#!/bin/bash

# NORMAL
# make clean
#make -j8 imd_mpi_eam_ttm_tmm_nbl_loadbalance_stress
# make -j8 imd_mpi_eam_ttm_tmm_colrad_nbl_loadbalance
#make -j8 imd_mpi_eam_nve_nbl_nrb_stress_filter_loadbalance
#make -j8 imd_mpi_eam_nve_nrb_stress_loadbalance_ar_filter
#make -j8 imd_mpi_eam_nve_nrb_stress_loadbalance_nbl_filter
# make -j8 imd_mpi_eam_nve_stress_loadbalance_nbl_filter

# make -j8 imd_mpi_eam_nve_nbl_nrb_loadbalance
# make -j8 imd_mpi_eam_ttm_tmm_nbl_loadbalance
#make -j8 imd_mpi_eam_ttm_tmm_nbl
#make -j8 imd_mpi_eam_nve_nbl_nrb_stress


# make -j8 imd_mpi_eam_ttm_tmm_colrad_nbl
# mv imd_mpi_eam_ttm_tmm_colrad ~/bin/
# exit 1


#standart
 #  make clean
 #  mpicc  -g3 -DMPI -DNBL -DEAM2 -DLOADBALANCE -DTTM -DTMM  -c \
 #  				 imd_maxwell.c imd_misc.c imd_param.c imd_alloc.c imd_io.c imd_io_3d.c \
 #  				 imd_loadBalance.c imd_loadBalance_direct.c\
 #  				 imd_potential.c imd_time.c imd_generate.c imd_distrib.c imd_main_3d.c imd_geom_3d.c \
 #  				 imd_pictures_3d.c imd_geom_mpi_3d.c imd_comm_force_3d.c imd_fix_cells_3d.c imd_mpiio.c imd_mpi_util.c imd.c imd_ttm.c \
 #  				 imd_interpol.c fminbnd3.c imd_tmm.c imd_forces_nbl.c imd_integrate.c

 #  mpicc  -g3 -o imd_mpi_eam_ttm_tmm_nbl_loadbalance imd_maxwell.o imd_integrate.o imd_misc.o imd_param.o imd_alloc.o \
 #  				imd_io.o imd_io_3d.o imd_potential.o imd_time.o imd_generate.o imd_distrib.o imd_main_3d.o imd_geom_3d.o imd_pictures_3d.o \
 #  				imd_geom_mpi_3d.o imd_comm_force_3d.o imd_fix_cells_3d.o imd_mpiio.o imd_mpi_util.o imd.o imd_ttm.o imd_interpol.o \
 #  				fminbnd3.o imd_tmm.o imd_forces_nbl.o imd_loadBalance.o imd_loadBalance_direct.o ./nn_interpol/libnn.a -lm  -lgsl -lgslcblas
 #  mv imd_mpi_eam_ttm_tmm_nbl_loadbalance ~/bin/
 # exit 1

#profiling . colrad OHNE OPENMP !
 #  make clean
 #  mpicc  -pg -g -DMPI -DNBL -DEAM2 -DTTM -DTMM -DCOLRAD -c imd_maxwell.c imd_misc.c imd_param.c imd_alloc.c imd_io.c imd_io_3d.c imd_potential.c imd_time.c imd_generate.c imd_distrib.c imd_main_3d.c imd_geom_3d.c imd_pictures_3d.c imd_geom_mpi_3d.c imd_comm_force_3d.c imd_fix_cells_3d.c imd_mpiio.c imd_mpi_util.c imd.c imd_ttm.c imd_interpol.c fminbnd3.c imd_tmm.c imd_colrad.c imd_forces_nbl.c imd_integrate.c
 #  mpicc  -pg -g -o imd_mpi_eam_ttm_tmm_nbl_colrad imd_maxwell.o imd_integrate.o imd_misc.o imd_param.o imd_alloc.o imd_io.o imd_io_3d.o imd_potential.o imd_time.o imd_generate.o imd_distrib.o imd_main_3d.o imd_geom_3d.o imd_pictures_3d.o imd_geom_mpi_3d.o imd_comm_force_3d.o imd_fix_cells_3d.o imd_mpiio.o imd_mpi_util.o imd.o imd_ttm.o imd_interpol.o fminbnd3.o imd_tmm.o imd_colrad.o imd_forces_nbl.o ./nn_interpol/libnn.a -lm    -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsollapackdense -lgsl -lgslcblas
 #  mv imd_mpi_eam_ttm_tmm_nbl_colrad ~/bin/
 # exit 1



# #CUSTOM MAKE MIT OPENMP ABER NUR FUER COLRAD
 # optflags=" -funroll-loops -march=corei7-avx -mtune=corei7-avx -mavx2 -ftree-vectorize -m64 -ffast-math"

gcc simpson_omp.c -pg -g -fopenmp -lm -o a.out
exit 0
 #
 # 
 # #MIT OPENMP
 # make clean
 AVXFLAGS=" -m64 -mavx2 -march=broadwell "
 CFLAGS="-Wno-unused-variable -I/user/eisfeld/sundials-4.1.0/instdir/include  \
 		  -I/usr/local/include/gsl"  
 LFLAGS="-L/user/eisfeld/sundials-4.1.0/instdir/lib64 -lsundials_cvode -lsundials_nvecserial \
 		  -lsundials_sunlinsollapackdense  -L/usr/local/lib -lgsl -L/usr/local/lib -lgslcblas"

 # LFLAGS=" /user/eisfeld/sundials-4.1.0/instdir/lib64/libsundials_cvode.a \
 # 		  /user/eisfeld/sundials-4.1.0/instdir/lib64/libsundials_nvecserial.a \
 # 		  /user/eisfeld/sundials-4.1.0/instdir/lib64/libsundials_sunlinsollapackdense.a \
 # 		  -lgsl -lgslcblas -llapack"
 OPTFLAGS="-O2" #-ffast-math"
 mpicc   $CFLAGS -DMPI  -DNBL -DEAM2 -DTTM -DTMM -DCOLRAD -DLOADBALANCE -c -Wno-unused-variable -Wno-unused-result -fopenmp\
 		 imd_maxwell.c imd_misc.c imd_param.c imd_alloc.c imd_io.c imd_io_3d.c imd_potential.c\
 		 imd_time.c imd_generate.c imd_distrib.c imd_main_3d.c imd_geom_3d.c imd_pictures_3d.c \
 		 imd_geom_mpi_3d.c imd_comm_force_3d.c \
 		 imd_fix_cells_3d.c imd_mpiio.c imd_mpi_util.c imd.c imd_ttm.c imd_interpol.c fminbnd3.c \
 		 imd_tmm.c imd_colrad.c imd_forces_nbl.c imd_integrate.c \
 		 imd_loadBalance.c imd_loadBalance_direct.c
 mpicc   -o imd_mpi_eam_ttm_tmm_nbl_colrad_loadbalance imd_maxwell.o imd_integrate.o imd_misc.o imd_param.o imd_alloc.o imd_io.o imd_io_3d.o \
 		 imd_loadBalance.o imd_loadBalance_direct.o\
		 imd_potential.o imd_time.o imd_generate.o imd_distrib.o imd_main_3d.o imd_geom_3d.o imd_pictures_3d.o imd_geom_mpi_3d.o \
		 imd_comm_force_3d.o imd_fix_cells_3d.o imd_mpiio.o imd_mpi_util.o imd.o imd_ttm.o imd_interpol.o fminbnd3.o imd_tmm.o imd_colrad.o \
		 imd_forces_nbl.o ./nn_interpol/libnn.a -lm  $LFLAGS -fopenmp
 mv imd_mpi_eam_ttm_tmm_nbl_colrad_loadbalance ~/bin/

# export LD_LIBRARY_PATH=/usr/local/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/user/eisfeld/sundials-4.1.0/instdir/lib64

 # 
 #
 #
 #

 # MIT INTEL MKL LINK OPTINOS
 # libintel=" -L/opt/intel/mkl/lib/intel64-Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl"
 # cflagintel=" -m64 -I/opt/intel//include " 
 # make clean 
 # mpicc  $cflagintel -O3 -DMPI  -DNBL -DEAM2 -DTTM -DTMM -DCOLRAD -c imd_maxwell.c imd_misc.c imd_param.c imd_alloc.c imd_io.c imd_io_3d.c imd_potential.c imd_time.c imd_generate.c imd_distrib.c imd_main_3d.c imd_geom_3d.c imd_pictures_3d.c imd_geom_mpi_3d.c imd_comm_force_3d.c imd_fix_cells_3d.c imd_mpiio.c imd_mpi_util.c imd.c imd_ttm.c imd_interpol.c fminbnd3.c imd_tmm.c imd_colrad.c imd_forces_nbl.c imd_integrate.c -fopenmp 
 # mpicc  -o imd_mpi_eam_ttm_tmm_nbl_colrad imd_maxwell.o imd_integrate.o imd_misc.o imd_param.o imd_alloc.o \
 # 			imd_io.o imd_io_3d.o imd_potential.o imd_time.o imd_generate.o imd_distrib.o imd_main_3d.o imd_geom_3d.o \
 # 			imd_pictures_3d.o imd_geom_mpi_3d.o imd_comm_force_3d.o imd_fix_cells_3d.o imd_mpiio.o imd_mpi_util.o\
 # 			imd.o imd_ttm.o imd_interpol.o fminbnd3.o imd_tmm.o imd_colrad.o imd_forces_nbl.o \
 # 			 ./nn_interpol/libnn.a \
 # 			 -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsollapackdense -lgsl -lgslcblas\
 # 			 -L/opt/intel/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
 # mv imd_mpi_eam_ttm_tmm_nbl_colrad ~/bin/





