
add_subdirectory(io)
add_subdirectory(mpi_wrapper)


# Liste de tous les sources
# -------------------------

file(GLOB_RECURSE SOURCES
     mesh/*.[ch]
     meshgen/*.[ch]
     pario/*.[ch]
     ppart/*.[ch]
     ppart/*.cxx
     struct/*.[ch]
     util/*.[ch]
     ext_wrapper/*.[ch]
     mpi_wrapper/pdm_mpi_node_first_rank.[hc]
     pdm.h
     pdm_priv.h)

file(GLOB_RECURSE SOURCESF
     *.[fF]
     *.[fF]90)

add_library(pdm OBJECT ${SOURCES})

if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)
  add_library(pdmf OBJECT ${SOURCESF})
endif()

target_include_directories(pdm
                           INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE  ${CMAKE_BINARY_DIR}
                           PRIVATE  ${PDM_BINARY_DIR}
                           PUBLIC  ${CMAKE_SOURCE_DIR}
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/pario
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/ppart
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/io
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/mpi_wrapper
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/ext_wrapper
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/mesh
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/meshgen
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/struct
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/util)
if (PDM_IN_PDMA)
  target_include_directories(pdm
                             PUBLIC ${PDMA_SOURCE_DIR}/src/plugins/renum
                             PUBLIC ${PDMA_SOURCE_DIR}/src/overlay)
endif()

if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)
  target_include_directories(pdmf
                           INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_BINARY_DIR}
                           PRIVATE  ${PDM_BINARY_DIR}
                           PUBLIC    ${CMAKE_SOURCE_DIR}
                           PUBLIC    ${CMAKE_CURRENT_SOURCE_DIR}
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/pario
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/ppart
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/io
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/mpi_wrapper
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/ext_wrapper
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/mesh
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/meshgen
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/struct
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/util)
    if (PDM_IN_PDMA)
      target_include_directories(pdmf
                                 PUBLIC ${PDMA_SOURCE_DIR}/src/plugins/renum
                                 PUBLIC ${PDMA_SOURCE_DIR}/src/overlay)
    endif()
endif()

if (PTSCOTCH_FOUND)
    target_include_directories(pdm
                              PRIVATE ${PTSCOTCH_INCLUDE_DIRS})
endif (PTSCOTCH_FOUND)

if (PARMETIS_FOUND)
    target_include_directories(pdm
                              PRIVATE ${PARMETIS_INCLUDE_DIRS})
endif (PARMETIS_FOUND)

# Liste des header
# ----------------

set (LIST_HEADER "")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/pdm.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/pdm_priv.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/pario/pdm_io.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_part_coarse_mesh.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_part.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/pario/pdm_io_tab.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/pario/pdm_writer.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_mesh_nodal.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/pario/pdm_writer_priv.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/util/pdm_mpi_node_first_rank.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/util/pdm_timer.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/util/pdm_fortran_to_c_string.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/meshgen/pdm_dcube_gen.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_gnum.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_part_to_block.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_block_to_part.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_hash_tab.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_handles.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/struct/pdm_global_mean.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mesh/pdm_cellface_orient.h")
if (PDM_IN_PDMA)
  list(APPEND LIST_HEADER "${CMAKE_BINARY_DIR}/pdm_config.h")
else()
  list(APPEND LIST_HEADER "${PDM_BINARY_DIR}/pdm_config.h")
endif()

if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)

  set (LIST_HEADERF "")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/pario/mod_pdm_io.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/pario/mod_pdm_writer.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/mod_pdm.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/ppart/mod_pdm_part.F90")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/ppart/pdm_partf.h")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_SOURCE_DIR}/pdmf.h")
    if (PDM_IN_PDMA)
    list(APPEND LIST_HEADERF "${CMAKE_BINARY_DIR}/pdm_configf.h")
  else()
    list(APPEND LIST_HEADERF "${PDM_BINARY_DIR}/pdm_configf.h")
  endif()
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_io.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_writer.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_gnum_location.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_part.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_mesh_nodal.mod")
  list(APPEND LIST_HEADERF "${CMAKE_CURRENT_BINARY_DIR}/pdm_dist_cloud_surf.mod")

endif ()

if (OPENMP_FOUND)
   set_target_properties(pdm PROPERTIES
                         COMPILE_FLAGS ${OpenMP_C_FLAGS})
endif()

# Options du link
# ---------------

unset(LINK_LIBRARIES)
list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})

if (PDM_ENABLE_SHARED)
  set_property(TARGET ${pdm} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdm_shared SHARED  $<TARGET_OBJECTS:pdm>)
  set_target_properties(pdm_shared PROPERTIES OUTPUT_NAME "pdm")
  set_target_properties(pdm_shared PROPERTIES VERSION ${PDM_VERSION})
  if (OPENMP_FOUND)
    set_target_properties(pdm_shared PROPERTIES
                          COMPILE_FLAGS ${OpenMP_C_FLAGS}
                          LINK_FLAGS ${OpenMP_C_FLAGS})
  endif()
  target_link_libraries(pdm_shared ${LINK_LIBRARIES} pdm_io_shared pdm_mpi_shared)
  set_target_properties(pdm_shared PROPERTIES
                            PUBLIC_HEADER "${LIST_HEADER}")

  get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)
  target_include_directories(pdm_shared PUBLIC ${PDM_INC})

endif()

if (PDM_ENABLE_STATIC)
  set_property(TARGET ${pdm} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdm_static STATIC  $<TARGET_OBJECTS:pdm>)
  set_target_properties(pdm_static PROPERTIES OUTPUT_NAME "pdm")
  set_target_properties(pdm_static PROPERTIES VERSION ${PDM_VERSION})

  if (OPENMP_FOUND)
    set_target_properties(pdm_static PROPERTIES
                          COMPILE_FLAGS ${OpenMP_C_FLAGS}
                          LINK_FLAGS ${OpenMP_C_FLAGS})
  endif()
  target_link_libraries(pdm_static ${LINK_LIBRARIES} pdm_io_static pdm_mpi_static)
  set_target_properties(pdm_static PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADER}")
  get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)
  target_include_directories(pdm_static PUBLIC ${PDM_INC})

endif()

if (PDM_ENABLE_SHARED_Fortran_INTERFACE)
  set_property(TARGET ${pdmf} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdmf_shared SHARED  $<TARGET_OBJECTS:pdmf>)
  set_target_properties(pdmf_shared PROPERTIES OUTPUT_NAME "pdmf")
  set_target_properties(pdmf_shared PROPERTIES VERSION ${PDM_VERSION})
  target_link_libraries(pdmf_shared pdm_shared)
  set_target_properties(pdmf_shared PROPERTIES
                            PUBLIC_HEADER "${LIST_HEADERF}")

  get_target_property(PDM_INC pdmf INCLUDE_DIRECTORIES)
  target_include_directories(pdmf_shared PUBLIC ${PDM_INC})

endif()

if (PDM_ENABLE_STATIC_Fortran_INTERFACE)
  set_property(TARGET ${pdmf} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdmf_static STATIC  $<TARGET_OBJECTS:pdmf>)
  set_target_properties(pdmf_static PROPERTIES OUTPUT_NAME "pdmf")
  set_target_properties(pdmf_static PROPERTIES VERSION ${PDM_VERSION})

  target_link_libraries(pdmf_static pdm_static)
  set_target_properties(pdmf_static PROPERTIES
                             PUBLIC_HEADER "${LIST_HEADERF}")

  get_target_property(PDM_INC pdmf INCLUDE_DIRECTORIES)
  target_include_directories(pdmf_static PUBLIC ${PDM_INC})

endif()

# Spécification de l'installation
# -------------------------------

set (_targets "")

if (PDM_ENABLE_SHARED)
  list (APPEND _targets pdm_shared)
endif()

if (PDM_ENABLE_STATIC)
  list (APPEND _targets pdm_static)
endif()

if (PDM_ENABLE_SHARED_Fortran_INTERFACE)
  list (APPEND _targets pdmf_shared)
endif()

if (PDM_ENABLE_STATIC_Fortran_INTERFACE)
  list (APPEND _targets pdmf_static)
endif()

install(TARGETS ${_targets}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        ARCHIVE DESTINATION lib)
