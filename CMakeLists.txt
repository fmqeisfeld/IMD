cmake_minimum_required(VERSION 2.10)


# set the project name
project(MYIMD C)
include_directories(~/imd_FDTD/src/ ~/imd_FDTD/src/nn_interpol /usr/lib64/mpi/gcc/openmpi/include/)
link_directories(./nn_interpol)							


set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} -Wall -DMPI -DEAM -DTTM -DTMM -DNBL -DDOUBLE -lm")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -DMPI -DEAM -DTTM -DTMM -DNBL -DDOUBLE -lm")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -DMPI -DEAM -DTTM -DTMM -DNBL -DDOUBLE -lm")

# set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ~/imd_FDTD/src/nn_interpol)
# set(EXTRA_LINK_LIBS -lm)

set( SOURCES imd.c imd_maxwell.c imd_integrate.c imd_misc.c 
	        imd_param.c imd_alloc.c imd_io.c imd_io_3d.c 
            imd_potential.c imd_time.c imd_generate.c 
            imd_distrib.c imd_main_3d.c imd_geom_3d.c )

set( MPISOURCES imd_mpi_util.c imd_main_mpi_3d.c imd_geom_mpi_3d.c 
		  	   imd_comm_force_3d.c imd_fix_cells_3d.c imd_mpiio.c )
set( NBLSOURCES imd_forces_nbl.c )		  	   
set( PAIRSOURCES imd_forces.c )
set( EAM2SOURCES imd_forces_eam2.c )
set( TTMSOURCES imd_ttm.c imd_interpol.c fminbnd3.c )
set( FDTDSOURCES imd_fdtd.c )
set( TMMSOURCES imd_tmm.c )



# add the executable
add_executable(MYIMD ${SOURCES} ${MPISOURCES} ${NBLSOURCES} ${EAM2SOURCES} ${TTMSOURCES} ${TMMSOURCES})
# target_compile_options(MYIMD -DMPI -DEAM -DTTM -DTMM -DNBL)
# target_link_libraries(MYIMD mpi)