
# Liste de tous les sources
# -------------------------

#
# TODO: On peut mieux faire en separant les sources des differents plugins

if (PDM_ENABLE_ANISO_AGGLO)
file(GLOB_RECURSE SOURCES 
     *.[ch]
     *.cxx
     *.hxx)
else ()
file(GLOB_RECURSE SOURCES 
     *.[ch])
endif ()

file(GLOB_RECURSE SOURCESF 
     *.[fF] 
     *.[fF]90)

add_library(pdma OBJECT ${SOURCES})
add_dependencies (pdma pdm)

add_library(pdmaf OBJECT ${SOURCESF})
add_dependencies (pdmaf pdmf)

get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)                          
target_include_directories(pdma PRIVATE ${PDM_INC})

get_target_property(PDM_INC pdmf INCLUDE_DIRECTORIES)                          
target_include_directories(pdmaf PRIVATE ${PDM_INC})

target_include_directories(pdma
                           INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE ${CMAKE_BINARY_DIR}
                           PRIVATE ${CMAKE_BINARY_DIR}/pdm
                           PRIVATE ${CMAKE_BINARY_DIR}/pdm/src
                           PUBLIC  ${CMAKE_SOURCE_DIR}
                           PUBLIC  ${CMAKE_SOURCE_DIR}/pdm
                           PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src
                           PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins
                           PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins/renum
                           PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins/coarse_mesh
                           PUBLIC  ${CMAKE_SOURCE_DIR}/src/overlay)

target_include_directories(pdmaf
                           INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE ${CMAKE_BINARY_DIR}
                           PRIVATE ${CMAKE_BINARY_DIR}/pdm
                           PRIVATE ${CMAKE_BINARY_DIR}/pdm/src
                           PUBLIC  ${CMAKE_SOURCE_DIR}
                           PUBLIC  ${CMAKE_SOURCE_DIR}/pdm
                           PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src
                           PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins
                           PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins/renum
                           PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins/coarse_mesh
                           PUBLIC  ${CMAKE_SOURCE_DIR}/src/overlay)

if (PTSCOTCH_FOUND)
    target_include_directories(pdma 
                               PRIVATE ${PTSCOTCH_INCLUDE_DIRS})
endif (PTSCOTCH_FOUND)

if (PARMETIS_FOUND)
    target_include_directories(pdma 
                               PRIVATE ${PARMETIS_INCLUDE_DIRS})
endif (PARMETIS_FOUND)


# Liste des header
# ----------------

unset ( LIST_HEADER)
set (LIST_HEADER "")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/overlay/pdm_overlay.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/overlay/pdm_edges_intersect.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/overlay/pdm_overlayf.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/overlay/pdm_poly_clipp.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/plugins/pdm_plugin.h")

if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)

  set (LIST_HEADERF "")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/overlay/mod_overlay.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/mod_overlay.mod")

endif ()

# Flags pour MPI
# --------------

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
   set_target_properties(pdm PROPERTIES
                              COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()

if ((NOT MPI_Fortran_COMPILER) AND MPI_Fortran_COMPILE_FLAGS)
    set_target_properties(pdm PROPERTIES
                              COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
endif()

# Options du link
# ---------------

unset(LINK_LIBRARIES)
list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})
#if (${PARMETIS_FOUND})
#  list(APPEND LINK_LIBRARIES ${PARMETIS_LIBRARIES})
#  list(APPEND LINK_LIBRARIES ${METIS_LIBRARIES})
#endif ()
#if (${PTSCOTCH_FOUND})
#  list(APPEND LINK_LIBRARIES ${PTSCOTCH_LIBRARIES})
#endif ()

if (PDM_ENABLE_STATIC)
  set_property(TARGET ${pdma} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdma_static STATIC  $<TARGET_OBJECTS:pdma>)
  set_target_properties(pdma_static PROPERTIES OUTPUT_NAME "pdma")
  target_link_libraries(pdma_static pdm_static ${LINK_LIBRARIES})
  set_target_properties(pdma_static PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADER}")
  add_dependencies (pdma_static pdm_static)

  set_target_properties(pdma_static PROPERTIES VERSION ${PDM_VERSION})
 
  target_include_directories(pdma_static
                             INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                             PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
                             PRIVATE ${CMAKE_BINARY_DIR}
                             PRIVATE ${CMAKE_BINARY_DIR}/pdm
                             PRIVATE ${CMAKE_BINARY_DIR}/pdm/src
                             PUBLIC  ${CMAKE_SOURCE_DIR}
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/pario
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/ppart
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/mpi_wrapper
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/mesh
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/struct
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins/renum
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins/coarse_mesh
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/overlay)

endif()

if (PDM_ENABLE_SHARED)
  set_property(TARGET ${pdma} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdma_shared SHARED  $<TARGET_OBJECTS:pdma>)
  
  set_target_properties(pdma_shared PROPERTIES OUTPUT_NAME "pdma")
  add_dependencies (pdma_shared pdm_shared)

  set_target_properties(pdma_shared PROPERTIES VERSION ${PDM_VERSION})

  target_link_libraries(pdma_shared pdm_shared ${LINK_LIBRARIES})
  set_target_properties(pdma_shared PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADER}")
  target_include_directories(pdma_shared
                             INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                             PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
                             PRIVATE ${CMAKE_BINARY_DIR}
                             PRIVATE ${CMAKE_BINARY_DIR}/pdm
                             PRIVATE ${CMAKE_BINARY_DIR}/pdm/src
                             PUBLIC  ${CMAKE_SOURCE_DIR}
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/pario
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/ppart
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/mpi_wrapper
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/mesh
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/struct
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins/renum
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins/coarse_mesh
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/overlay)

endif()

if (PDM_ENABLE_STATIC_Fortran_INTERFACE)
  set_property(TARGET ${pdmaf} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdmaf_static STATIC  $<TARGET_OBJECTS:pdmaf>)
  set_target_properties(pdmaf_static PROPERTIES OUTPUT_NAME "pdmaf")
  target_link_libraries(pdmaf_static pdmf_static ${LINK_LIBRARIES})
  set_target_properties(pdma_static PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADERF}")
  add_dependencies (pdmaf_static pdmf_static)

  set_target_properties(pdmaf_static PROPERTIES VERSION ${PDM_VERSION})
 
  target_include_directories(pdmaf_static
                             INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                             PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
                             PRIVATE ${CMAKE_BINARY_DIR}
                             PRIVATE ${CMAKE_BINARY_DIR}/pdm
                             PRIVATE ${CMAKE_BINARY_DIR}/pdm/src
                             PUBLIC  ${CMAKE_SOURCE_DIR}
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/pario
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/ppart
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/mpi_wrapper
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/mesh
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/struct
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/overlay
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins/coarse_mesh
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins/renum)

endif()

if (PDM_ENABLE_SHARED_Fortran_INTERFACE)
  set_property(TARGET ${pdmaf} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdmaf_shared SHARED  $<TARGET_OBJECTS:pdmaf>)
  
  set_target_properties(pdmaf_shared PROPERTIES OUTPUT_NAME "pdmaf")
  add_dependencies (pdmaf_shared pdmf_shared)

  set_target_properties(pdmaf_shared PROPERTIES VERSION ${PDM_VERSION})

  target_link_libraries(pdmaf_shared pdmf_shared ${LINK_LIBRARIES})
  set_target_properties(pdmaf_shared PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADER}")
  target_include_directories(pdmaf_shared
                             INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                             PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
                             PRIVATE ${CMAKE_BINARY_DIR}
                             PRIVATE ${CMAKE_BINARY_DIR}/pdm
                             PRIVATE ${CMAKE_BINARY_DIR}/pdm/src
                             PUBLIC  ${CMAKE_SOURCE_DIR}
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/pario
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/ppart
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/mpi_wrapper
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/mesh
                             PUBLIC  ${CMAKE_SOURCE_DIR}/pdm/src/struct
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/overlay
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins/coarse_mesh
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins
                             PUBLIC  ${CMAKE_SOURCE_DIR}/src/plugins/renum)

endif()

message(STATUS "Add library pdma - Succes")

# Spécification de l'installation
# -------------------------------

set (_targets "")

if (PDM_ENABLE_SHARED)
  list (APPEND _targets pdma_shared)
endif()

if (PDM_ENABLE_STATIC)
  list (APPEND _targets pdma_static)
endif()

if (PDM_ENABLE_STATIC_Fortran_INTERFACE)
  list (APPEND _targets pdmaf_static)
endif ()

if (PDM_ENABLE_SHARED_Fortran_INTERFACE)
  list (APPEND _targets pdmaf_shared)
endif ()

install(TARGETS ${_targets}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        ARCHIVE DESTINATION lib)

