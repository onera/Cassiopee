# Find the Aurora compiler.
#
# This code sets the following variables:
#
#  NEC_FORTRAN_COMPILER
#
# See also UseAURORA.cmake
#=============================================================================

find_program( NEC_FORTRAN_COMPILER
              NAMES nfort
              HINTS "$ENV{NEC_HOME}/bin"
              )

find_program( NEC_C_COMPILER
              NAMES ncc
              HINTS "$ENV{NEC_HOME}/bin"
              )

find_program( NEC_CXX_COMPILER
              NAMES nc++
              HINTS "$ENV{NEC_HOME}/bin"
              )

find_library(NEC_VEOS_LIBRARY
             NAMES libveo.so
             HINTS "$ENV{NEC_HOME}/veos/lib64"
)

# find_library(NEC_VEOSTRACE_LIBRARY
#              NAMES libveftrace.so
#              HINTS "$ENV{NEC_HOME}/veos/lib64"
# )

find_library(NEC_VEOSINFO_LIBRARY
             NAMES libveosinfo.so
             HINTS "$ENV{NEC_HOME}/veos/lib64"
)

# find_library(NEC_VEIO_LIBRARY
#              NAMES libveio.so
#              HINTS "$ENV{NEC_HOME}/lib"
# )
# message("Found NEC_VEIO_LIBRARY=",${NEC_VEIO_LIBRARY})

find_path(NEC_VEOS_INCLUDE_DIR
          NAMES ve_offload.h
          HINTS "$ENV{NEC_HOME}/veos/include")

find_path(NEC_VHSHM_INCLUDE_DIR
             NAMES vhshm.h
             HINTS "$ENV{NEC_HOME}/include")

if (NEC_FORTRAN_COMPILER  AND
    NEC_C_COMPILER        AND
    NEC_VEOS_LIBRARY      AND
    # NEC_VEOSTRACE_LIBRARY AND
    NEC_VEOSINFO_LIBRARY  AND
    NEC_VEOS_INCLUDE_DIR  AND
    NEC_VHSHM_INCLUDE_DIR)

  set(NEC_AURORA_INCLUDE_DIR ${NEC_VEOS_INCLUDE_DIR} ${NEC_VHSHM_INCLUDE_DIR})
  # set(NEC_AURORA_INCLUDE_DIR ${NEC_VEOS_INCLUDE_DIR})
  set(NEC_AURORA_LIBRARIES ${NEC_VEOS_LIBRARY} ${NEC_VEOSINFO_LIBRARY}) #${NEC_VEOSTRACE_LIBRARY})

else()

  set(NEC_AURORA_INCLUDE_DIR "")
  set(NEC_AURORA_LIBRARIES "")

endif()

include( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( AURORA REQUIRED_VARS NEC_FORTRAN_COMPILER NEC_C_COMPILER NEC_CXX_COMPILER)

mark_as_advanced( NEC_FORTRAN_COMPILER NEC_C_COMPILER)
