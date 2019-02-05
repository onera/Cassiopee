# - Try to find PTSCOTCH
# Once done this will define
#
#  PTSCOTCH_FOUND        - system has found PTSCOTCH
#  PTSCOTCH_INCLUDE_DIRS - include directories for PTSCOTCH
#  PTSCOTCH_LIBARIES     - libraries for PTSCOTCH
#  PTSCOTCH_VERSION      - version for PTSCOTCH
#

# mettre un if si on ne peut pas compiler et mettre un message

set(PTSCOTCH_FOUND FALSE)

# Check for header file
find_path(PTSCOTCH_INCLUDE_DIRS ptscotch.h
  HINTS ${PTSCOTCH_DIR}/include $ENV{PTSCOTCH_DIR}/include
  DOC "Directory where the Pt-PTSCOTCH header is located"
  )

# Check for scotch
find_library(PTSCOTCH_LIBRARY
  NAMES scotch
  HINTS ${PTSCOTCH_DIR}/lib $ENV{PTSCOTCH_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The PTSCOTCH library"
  )
find_library(PTSCOTCH_LIBRARY
  NAMES scotch
  DOC "The PTSCOTCH library"
  )

# Check for scotcherr
find_library(PTSCOTCHERR_LIBRARY
  NAMES scotcherr
  HINTS ${PTSCOTCH_DIR}/lib $ENV{PTSCOTCH_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The PTSCOTCH-ERROR library"
  )
find_library(PTSCOTCHERR_LIBRARY
  NAMES scotcherr
  DOC "The PTSCOTCH-ERROR library"
  )

# Check for ptscotch
find_library(PTPTSCOTCH_LIBRARY
  NAMES ptscotch
  HINTS ${PTSCOTCH_DIR}/lib $ENV{PTSCOTCH_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The Pt-PTSCOTCH library"
  )
find_library(PTPTSCOTCH_LIBRARY
  NAMES ptscotch
  DOC "The Pt-PTSCOTCH library"
  )

# Check for ptscotcherr
find_library(PTPTSCOTCHERR_LIBRARY
  NAMES ptscotcherr
  HINTS ${PTSCOTCH_DIR}/lib $ENV{PTSCOTCH_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The Pt-PTSCOTCH-ERROR library"
  )
find_library(PTPTSCOTCHERR_LIBRARY
  NAMES ptscotcherr
  DOC "The Pt-PTSCOTCH-ERROR library"
  )

set(PTSCOTCH_LIBRARIES ${PTPTSCOTCH_LIBRARY} ${PTPTSCOTCHERR_LIBRARY} ${PTSCOTCH_LIBRARY} ${PTSCOTCHERR_LIBRARY} CACHE STRING "SOCTCH parallel libraries")
set(PTSCOTCH_SEQ_LIBRARIES ${PTSCOTCH_LIBRARY} ${PTSCOTCHERR_LIBRARY} CACHE STRING "SOCTCH sequential libraries")

if (PTSCOTCH_INCLUDE_DIRS AND PTSCOTCH_LIBRARIES)

  mark_as_advanced(PTSCOTCH_INCLUDE_DIRS PTSCOTCH_LIBRARIES PTSCOTCH_SEQ_LIBRARIES)

  if (PTPTSCOTCH_LIBRARY)
    string(REGEX REPLACE "(^.*)/lib/libptscotch.*$" "\\1" PTSCOTCH_LIBRARY_PATH ${PTPTSCOTCH_LIBRARY} )
  endif (PTPTSCOTCH_LIBRARY)

  # Set flags for building test program
  set(CMAKE_REQUIRED_INCLUDES ${PTSCOTCH_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${PTSCOTCH_LIBRARIES})

  # Add MPI variables if MPI has been found
  if (MPI_C_LIBRARIES)
    set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${MPI_C_INCLUDE_PATH})
    set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${MPI_C_LIBRARIES})
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${MPI_C_COMPILE_FLAGS}")
  endif()
    
  set(CMAKE_REQUIRED_FLAGS  "${CMAKE_REQUIRED_FLAGS} ${CMAKE_C_FLAGS}")

  set(PTSCOTCH_CONFIG_TEST_VERSION_C
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/scotch_config_test_version.c")
 
  file(WRITE ${PTSCOTCH_CONFIG_TEST_VERSION_C} "
#include <stdint.h>
#include <stdio.h>
#include <mpi.h>
#include <ptscotch.h>

int main() {
  printf(\"%i.%i.%i\\n\", SCOTCH_VERSION,
	                  SCOTCH_RELEASE,
	                  SCOTCH_PATCHLEVEL);
  return 0;
}
")

try_run(
    PTSCOTCH_CONFIG_TEST_VERSION_EXITCODE
    PTSCOTCH_CONFIG_TEST_VERSION_COMPILED
    ${CMAKE_CURRENT_BINARY_DIR}
    ${PTSCOTCH_CONFIG_TEST_VERSION_C}
    CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
      "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}"
    COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE OUTPUT
)


if (   (NOT PTSCOTCH_CONFIG_TEST_VERSION_COMPILED)
    OR (NOT (PTSCOTCH_CONFIG_TEST_VERSION_EXITCODE EQUAL 0))) 
    message(WARNING "Unable to determine PTSCOTCH version")
    set(PTSCOTCH_VERSION_OK TRUE)
    set(PTSCOTCH_VERSION "??.??.??" CACHE TYPE STRING)
    set(PTSCOTCH_TEST_COMPILE TRUE)
else ()
    set(PTSCOTCH_VERSION ${OUTPUT} CACHE TYPE STRING)
    mark_as_advanced(PTSCOTCH_VERSION)

   if (PTSCOTCH_FIND_VERSION)
  # Check if version found is >= required version
      if (NOT "${PTSCOTCH_VERSION}" VERSION_LESS "${PTSCOTCH_FIND_VERSION}")
  	set(PTSCOTCH_VERSION_OK TRUE)
      endif()
   else()
  # No specific version requested
       set(PTSCOTCH_VERSION_OK TRUE)
   endif()
   mark_as_advanced(PTSCOTCH_VERSION_OK)
endif ()

unset(PTPTSCOTCH_LIBRARY CACHE)
unset(PTPTSCOTCHERR_LIBRARY CACHE)
unset(PTSCOTCH_LIBRARY CACHE)
unset(PTSCOTCHERR_LIBRARY CACHE)
unset(PTSCOTCH_LIBRARY_PATH CACHE)

endif ()
#
# Standard package handling
#
INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTSCOTCH
                                  "PTSCOTCH could not be found. Be sure to set PTSCOTCH_DIR."
                                  PTSCOTCH_LIBRARIES
                                  PTSCOTCH_SEQ_LIBRARIES
                                  PTSCOTCH_INCLUDE_DIRS
                                  PTSCOTCH_VERSION
                                  PTSCOTCH_VERSION_OK)
