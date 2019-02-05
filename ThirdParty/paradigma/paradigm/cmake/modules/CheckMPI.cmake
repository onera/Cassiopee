  
set(WRAPPER_COMPILER_SCRIPT  "checkWrapperCompiler.sh")

if (MPI_C_COMPILER)
  message(STATUS "Checking compatibility between MPI_C_COMPILER and CMAKE_C_COMPILER")
  if (NOT (MPI_C_COMPILER STREQUAL CMAKE_C_COMPILER))
    execute_process(COMMAND bash ${CMAKE_MODULE_PATH}/${WRAPPER_COMPILER_SCRIPT} CC ${MPI_C_COMPILER} ${CMAKE_C_COMPILER} 
    RESULT_VARIABLE WRAPPER_COMPILER_MATCH_C)
    if ( NOT ${WRAPPER_COMPILER_MATCH_C}  )
      message(FATAL_ERROR "${MPI_C_COMPILER} is not compatible with ${CMAKE_C_COMPILER}")
    endif()
  endif()
  message(STATUS "Checking compatibility between MPI_C_COMPILER and CMAKE_C_COMPILER - done")
endif ()

if (MPI_CXX_COMPILER)
  message(STATUS "Checking compatibility between MPI_CXX_COMPILER and CMAKE_CXX_COMPILER")
  if (NOT (MPI_CXX_COMPILER STREQUAL CMAKE_CXX_COMPILER))
    execute_process(COMMAND bash ${CMAKE_MODULE_PATH}/${WRAPPER_COMPILER_SCRIPT} CXX ${MPI_CXX_COMPILER} ${CMAKE_CXX_COMPILER} 
    RESULT_VARIABLE WRAPPER_COMPILER_MATCH_CXX)
    if ( NOT ${WRAPPER_COMPILER_MATCH_CXX}  )
      message(FATAL_ERROR "${MPI_CXX_COMPILER} is not compatible with ${CMAKE_CXX_COMPILER}")
    endif()
  endif()
  message(STATUS "Checking compatibility between MPI_CXX_COMPILER and CMAKE_CXX_COMPILER - done")
endif ()

if (MPI_Fortran_COMPILER)
  message(STATUS "Checking compatibility between MPI_Fortran_COMPILER and CMAKE_Fortran_COMPILER")
  if (NOT (MPI_Fortran_COMPILER STREQUAL CMAKE_Fortran_COMPILER))
    execute_process(COMMAND bash ${CMAKE_MODULE_PATH}/${WRAPPER_COMPILER_SCRIPT} FC ${MPI_Fortran_COMPILER} ${CMAKE_Fortran_COMPILER} 
    RESULT_VARIABLE WRAPPER_COMPILER_MATCH_CXX)
    if ( NOT ${WRAPPER_COMPILER_MATCH_CXX}  )
      message(FATAL_ERROR "${MPI_Fortran_COMPILER} is not compatible with ${CMAKE_Fortran_COMPILER}")
    endif()
  endif()
  message(STATUS "Checking compatibility between MPI_Fortran_COMPILER and CMAKE_Fortran_COMPILER - done")
endif ()
