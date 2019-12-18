
# Liste de tous les sources
# -------------------------

if (MPI_C_COMPILER)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

if (MPI_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif(MPI_Fortran_COMPILER)


file(GLOB_RECURSE SOURCES
     pdm_mpi.c
     pdm_mpi_ext_dependencies.c
     ../pdm_mpi.h
     ../pdm_mpi_ext_dependencies.h
     ../../pdm.h
     ../../pdm_priv.h)


add_library(pdm_mpi OBJECT ${SOURCES})


# Flags pour MPI
# --------------

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
   set_target_properties(pdm_mpi PROPERTIES
                              COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()

target_include_directories(pdm_mpi
                           INTERFACE ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE   ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE  ${PDM_BINARY_DIR}
                           PRIVATE  ${CMAKE_BINARY_DIR}
                           PUBLIC  ${CMAKE_SOURCE_DIR}
                           PUBLIC  ${PDM_SOURCE_DIR}/src/io
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/..
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/../..)

if (PTSCOTCH_FOUND)
    target_include_directories(pdm_mpi
                              PRIVATE ${PTSCOTCH_INCLUDE_DIRS})
endif (PTSCOTCH_FOUND)

if (PARMETIS_FOUND)
    target_include_directories(pdm_mpi
                              PRIVATE ${PARMETIS_INCLUDE_DIRS})
endif (PARMETIS_FOUND)

set (LIST_HEADER "")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/../pdm_mpi.h")
list(APPEND LIST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/../pdm_mpi_ext_dependencies.h")

# Options du link
# ---------------

unset(LINK_LIBRARIES)
list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})
if (${PARMETIS_FOUND})
  list(APPEND LINK_LIBRARIES ${PARMETIS_LIBRARIES})
  list(APPEND LINK_LIBRARIES ${METIS_LIBRARIES})
endif ()
if (${PTSCOTCH_FOUND})
  list(APPEND LINK_LIBRARIES ${PTSCOTCH_LIBRARIES})
endif ()

if (PDM_ENABLE_SHARED)
  set_property(TARGET ${pdm_mpi} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdm_mpi_shared SHARED  $<TARGET_OBJECTS:pdm_mpi>)
  set_target_properties(pdm_mpi_shared PROPERTIES OUTPUT_NAME "pdm_mpi")
  set_target_properties(pdm_mpi_shared PROPERTIES VERSION ${PDM_VERSION})
  target_link_libraries(pdm_mpi_shared ${LINK_LIBRARIES} pdm_io_shared)
  set_target_properties(pdm_mpi_shared PROPERTIES
                            PUBLIC_HEADER "${LIST_HEADER}")

  get_target_property(PDM_INC pdm_mpi INCLUDE_DIRECTORIES)
  if (NOT MPI_C_COMPILER)
    if (MPI_C_LIBRARIES)
      target_link_libraries(pdm_mpi_shared ${MPI_C_LIBRARIES})
    endif()
  endif ()
  target_include_directories(pdm_mpi_shared PUBLIC ${PDM_INC})
endif()

if (PDM_ENABLE_STATIC)
  set_property(TARGET ${pdm_mpi} PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(pdm_mpi_static STATIC  $<TARGET_OBJECTS:pdm_mpi>)
  set_target_properties(pdm_mpi_static PROPERTIES OUTPUT_NAME "pdm_mpi")
  set_target_properties(pdm_mpi_static PROPERTIES VERSION ${PDM_VERSION})
  target_link_libraries(pdm_mpi_static ${LINK_LIBRARIES} pdm_io_static)
  set_target_properties(pdm_mpi_static PROPERTIES
                            PUBLIC_HEADER "${LIST_HEADER}")

  get_target_property(PDM_INC pdm_mpi INCLUDE_DIRECTORIES)
  target_include_directories(pdm_mpi_static PUBLIC ${PDM_INC})
endif()

# Spécification de l'installation
# -------------------------------

set (_targets "")

if (PDM_ENABLE_SHARED)
  list (APPEND _targets pdm_mpi_shared)
endif()

if (PDM_ENABLE_STATIC)
  list (APPEND _targets pdm_mpi_static)
endif()

install(TARGETS ${_targets}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        ARCHIVE DESTINATION lib)
