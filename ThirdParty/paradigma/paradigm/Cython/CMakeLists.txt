
if (MPI_C_COMPILER)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

if (MPI_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif(MPI_Fortran_COMPILER)

if (PDM_ENABLE_LONG_G_NUM)
   set(PYX_NPY_G_NUM NPY.int64_t)
   set(PYX_NPY_G_NUM_DTYPE NPY.int64)
   set(PDM_G_NUM_NPY_INT NPY.NPY_INT64)
   set(PYX_G_NUM long)
   set(PYX_L_NUM int)
else ()
   set(PYX_NPY_G_NUM NPY.int32_t)
   set(PYX_NPY_G_NUM_DTYPE NPY.int32)
   set(PDM_G_NUM_NPY_INT NPY.NPY_INT32)
   set(PYX_G_NUM int)
   set(PYX_L_NUM int)
endif ()

configure_file(Pypdm.pyx.in Pypdm.pyx)
configure_file(pdm_part_coarse.pxi.in pdm_part_coarse.pxi)

file ( COPY . DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" FILES_MATCHING PATTERN "*.pxi")
file ( COPY . DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" FILES_MATCHING PATTERN "*.py")

if (PDM_IN_PDMA)
  file ( COPY ${PDMA_SOURCE_DIR}/Cython DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" FILES_MATCHING PATTERN "*.pxi")
  file ( COPY ${PDMA_SOURCE_DIR}/Cython DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" FILES_MATCHING PATTERN "*.py")
endif()

set (include_dirs "")
list (APPEND include_dirs ${Mpi4Py_INCLUDE_DIR})
list (APPEND include_dirs ${NumPy_INCLUDE_DIR})
list (APPEND include_dirs ${CMAKE_CURRENT_BINARY_DIR})
list (APPEND include_dirs ${CMAKE_BINARY_DIR})
if (PDM_IN_PDMA)
  list (APPEND include_dirs ${PDMA_SOURCE_DIR}/src/mesh)
  list (APPEND include_dirs ${PDMA_SOURCE_DIR}/src/overlay)
  list (APPEND include_dirs ${PDMA_SOURCE_DIR}/src/plugins)
else ()
  set(COMMENT_NOT_PDMA \#)
endif()

include(UseCython)
if (PDM_IN_PDMA)
cython_add_module (Pypdm
                   PYX_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/Pypdm.pyx"
                   OTHER_SOURCES $<TARGET_OBJECTS:pdm> $<TARGET_OBJECTS:pdma>
                   INCLUDE_DIRECTORIES ${include_dirs})
else ()
cython_add_module (Pypdm
                   PYX_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/Pypdm.pyx"
                   OTHER_SOURCES $<TARGET_OBJECTS:pdm>
                   INCLUDE_DIRECTORIES ${include_dirs})
endif()

get_target_property(PDM_INC pdm INCLUDE_DIRECTORIES)
target_include_directories(Pypdm PUBLIC ${PDM_INC})

unset(LINK_LIBRARIES)
list(APPEND LINK_LIBRARIES ${COMMON_LIBRARIES})
if (${PARMETIS_FOUND})
  list(APPEND LINK_LIBRARIES ${PARMETIS_LIBRARIES})
  list(APPEND LINK_LIBRARIES ${METIS_LIBRARIES})
endif ()
if (${PTSCOTCH_FOUND})
  list(APPEND LINK_LIBRARIES ${PTSCOTCH_LIBRARIES})
  list(APPEND LINK_LIBRARIES ${SCOTCH_LIBRARIES})
endif ()
if (PDM_IN_PDMA)
  list(APPEND LINK_LIBRARIES pdma_shared)
endif()

target_link_libraries(Pypdm ${LINK_LIBRARIES} pdm_mpi_shared)

if ((NOT MPI_C_COMPILER) AND MPI_C_COMPILE_FLAGS)
   set_target_properties(Pypdm PROPERTIES
                              COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS})
endif()

set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/Pypdm/__init__.py")
set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}timestamp")

add_custom_command(OUTPUT ${OUTPUT}
                       COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} --quiet build
                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                       DEPENDS ${DEPS})

set (SITE_PACKAGES_PDM "${CMAKE_PDM_INSTALL_PYTHON_DIR}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/Pypdm")

add_custom_target(target ALL DEPENDS ${OUTPUT})


if ((CMAKE_C_COMPILER_ID STREQUAL "GNU") OR
    (CMAKE_C_COMPILER_ID STREQUAL "Intel") OR
    (CMAKE_C_COMPILER_ID STREQUAL "AppleClang")  OR
    (CMAKE_C_COMPILER_ID STREQUAL "PGI"))
  set_target_properties (Pypdm PROPERTIES COMPILE_FLAGS "-w")
endif()

install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} --quiet install --prefix=${CMAKE_PDM_INSTALL_PYTHON_DIR} --install-data=${SITE_PACKAGES_PDM} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")
