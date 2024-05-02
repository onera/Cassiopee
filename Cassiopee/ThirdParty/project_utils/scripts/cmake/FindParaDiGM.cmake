#
# CMake module for ParaDIGM (E. Quémerais ONERA )
#
#

# set(PARADIGM_PATHS $ENV{PARADIGM_HOME})

# message ("-- PARADIGM: Searching..." )
# message (${PARADIGM_PATHS})
# message ($ENV${PARADIGM_HOME})
# message (${ENV${PARADIGM_HOME}})

find_path(PARADIGM_INCLUDE_DIR
          NAMES    pdm.h
          PATH_SUFFIXES include
          HINTS "$ENV{PARADIGM_HOME}"
)

find_library(PARADIGM_LIBRARY
             NAMES libpdm.so
             HINTS "$ENV{PARADIGM_HOME}/lib"
)
find_library(PARADIGMA_LIBRARY
             NAMES libpdma.so
             HINTS "$ENV{PARADIGM_HOME}/lib"
)

find_library(PARADIGM_MPI_LIBRARY
             NAMES libpdm_mpi.so
             HINTS "$ENV{PARADIGM_HOME}/lib"
)
# message (${PARADIGM_INCLUDE_DIR})
# message (${PARADIGM_LIBRARY})
# message (${PARADIGMA_LIBRARY})
# message (${PARADIGM_MPI_LIBRARY})

if (PARADIGM_INCLUDE_DIR AND
    PARADIGM_LIBRARY     AND
    PARADIGMA_LIBRARY    AND
    PARADIGM_MPI_LIBRARY)

  set(PARADIGM_LIBRARIES ${PARADIGM_LIBRARY} ${PARADIGMA_LIBRARY} ${PARADIGM_MPI_LIBRARY})

else()

  set(PARADIGM_INCLUDE_DIR "")
  set(PARADIGM_LIBRARIES   "")

endif()

# set(PARADIGM_LIBRARIES ${PARADIGM_LIBRARY} ${PARADIGMA_LIBRARY} ${PARADIGM_MPI_LIBRARY})
# CACHE STRING "ParaDigM library")

# message (${PARADIGM_LIBRARIES})
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARADIGM DEFAULT_MSG PARADIGM_LIBRARIES PARADIGM_INCLUDE_DIR )

# set(CMAKE_REQUIRED_INCLUDES ${PARADIGM_INCLUDE_DIR})
# set(CMAKE_REQUIRED_LIBRARIES ${PARADIGM_LIBRARIES})
# find_package_handle_standard_args(PARADIGM REQUIRED_VARS PARADIGM_INCLUDE_DIR)

