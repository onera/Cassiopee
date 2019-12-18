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
  list(APPEND LINK_LIBRARIES pdma_static pdm_static pdma_static ${COMMON_LIBRARIES})
else ()
  list(APPEND LINK_LIBRARIES pdma_shared pdm_shared ${COMMON_LIBRARIES})
endif()

if (OPENMP_FOUND)
  list(APPEND LINK_LIBRARIES ${OpenMP_C_FLAGS})
endif()

get_target_property(PDMA_INC pdma INCLUDE_DIRECTORIES)
get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)

# Test 2D (plusieurs sous-domaines par proc)
# ------------------------------------------

include(TestCreate)

test_c_create (pdm_t_plane_meshes_intersect ${MPIEXEC_MAX_NUMPROCS})
if (PDM_HAVE_PARMETIS OR PDM_HAVE_PTSCOTCH)
 test_c_create (pdm_t_cellface_orient_renum ${MPIEXEC_MAX_NUMPROCS})
 test_c_create (pdm_t_cell_center_dist ${MPIEXEC_MAX_NUMPROCS})
 test_c_create (pdm_t_graph_multicolor       2)
endif()
