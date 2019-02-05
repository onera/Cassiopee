# --------------------------
# Construction des Tests
# --------------------------

if (MPI_C_COMPILER)
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif() 

if (MPI_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

if (MPI_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif() 

unset(LINK_LIBRARIES)

if (PDM_ENABLE_STATIC)
 if (PDM_IN_PDMA) 
  list(APPEND LINK_LIBRARIES
       pdma_static) 
 endif()
 list(APPEND LINK_LIBRARIES
      pdm_static
      pdm_mpi_static
      ${COMMON_LIBRARIES})
 if (PDM_IN_PDMA) 
  list(APPEND LINK_LIBRARIES
       pdma_static) 
 endif()
else ()
 list(APPEND LINK_LIBRARIES
      pdm_shared
      pdm_mpi_shared ${COMMON_LIBRARIES})
 if (PDM_IN_PDMA) 
  list(APPEND LINK_LIBRARIES
       pdma_shared) 
 endif()
endif()

get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)                          


# Test pdm_t_part_cube
# --------------

include(TestCreate)

if (PDM_HAVE_PARMETIS OR PDM_HAVE_PTSCOTCH)
  test_c_create (pdm_t_part_dcube ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_global_mean ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_writer_ensight_3d ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_plane_agglo ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_cellface_orient ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_gnum ${MPIEXEC_MAX_NUMPROCS})
  test_c_create (pdm_t_gnum_with_merge ${MPIEXEC_MAX_NUMPROCS})
endif()
test_c_create (pdm_t_cellface_orient_poly 1)
