#------------------------------------------------------------------------------
#
# Default flags for several compiler
#
#------------------------------------------------------------------------------

cmake_host_system_information(RESULT HOSTNAME QUERY HOSTNAME)

#------------------------------------------------------------------------------
# Fortran default flags
#------------------------------------------------------------------------------

if (NOT PASS_DEFAULT_FLAGS)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

  set (CMAKE_Fortran_FLAGS "-fPIC -Wall -pedantic -std=gnu -Wno-unused-dummy-argument -Wno-maybe-uninitialized")

  set (CMAKE_Fortran_FLAGS_RELEASE         "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG           "-O0 -g -fcheck=bounds -fbacktrace")
  set (CMAKE_Fortran_FLAGS_PROFILING       "-O3 -pg")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO  "-O3 -g")
  set (CMAKE_Fortran_FLAGS_MINSIZEREL      "-O2 -g")

  set (FORTRAN_LIBRARIES                   )
  set (FORTRAN_LIBRARIES_FLAG              )

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")

  set (CMAKE_Fortran_FLAGS "-cpp -fpic -warn -diag-disable 7712")

  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")

  include(CheckCCompilerFlag)

  set (CMAKE_Fortran_FLAGS_DEBUG          "-O0 -g -check all -check nopointer -traceback")
  set (CMAKE_Fortran_FLAGS_PROFILING      "${CMAKE_Fortran_FLAGS_RELEASE} -pg")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -g")
  set (CMAKE_Fortran_FLAGS_MINSIZEREL     "-O2 -g")

  find_library(FORTRAN_LIBRARIES ifcore)
  mark_as_advanced (FORTRAN_LIBRARIES)
  set (FORTRAN_LIBRARIES_FLAG          )

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "XL")

  # xlf
  # ---

  set (CMAKE_Fortran_FLAGS " ${CMAKE_Fortran_FLAGS} -q64 -qextname -qsuffix=cpp=f90")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -qhot -qstrict")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g -qcheck")
  set (CMAKE_Fortran_FLAGS_PROFILING       "${CMAKE_Fortran_FLAGS_RELEASE} -p")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO  "-O3 -qhot -g")
  set (CMAKE_Fortran_FLAGS_MINSIZEREL      "-O3")

  set(FORTRAN_LIBRARIES xl xlf90_r xlsmp xlopt ${FORTRAN_LIBRARIES})

  if (${HOSTNAME} STREQUAL "tanit")
    link_directories(/opt/ibmcmp/xlsmp/3.1/lib64 /opt/ibmcmp/vacpp/12.1/lib64 /opt/ibmcmp/xlf/14.1/lib64)
  endif()

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")

  # pgi
  # ---

  # Ajout des flags communs

  set (CMAKE_Fortran_FLAGS "-Mpreprocess -noswitcherror")

  set (CMAKE_Fortran_FLAGS_RELEASE         "-fast")
  set (CMAKE_Fortran_FLAGS_DEBUG           "-g -Mbounds")
  set (CMAKE_Fortran_FLAGS_PROFILING       "${CMAKE_Fortran_FLAGS_RELEASE} -Mprof=func,lines")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO  "${CMAKE_Fortran_FLAGS_RELEASE} -g")
  set (CMAKE_Fortran_FLAGS_MINSIZEREL      "-O2")

  set (FORTRAN_LIBRARIES                   )
  set (FORTRAN_LIBRARIES_FLAG    -pgf90libs)

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")

  set (CMAKE_Fortran_FLAGS "-eF -em -J.")

  set (CMAKE_Fortran_FLAGS_RELEASE         "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG           "-g")
  set (CMAKE_Fortran_FLAGS_PROFILING       "${CMAKE_Fortran_FLAGS_RELEASE} -h profile_generate")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO  "${CMAKE_Fortran_FLAGS_RELEASE} -g")
  set (CMAKE_Fortran_FLAGS_MINSIZEREL      "-O2")

  set (FORTRAN_LIBRARIES         )
  set (FORTRAN_LIBRARIES_FLAG    )

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PathScale")

  set (CMAKE_Fortran_FLAGS "-Wall -Wno-unused -cpp")

  set (CMAKE_Fortran_FLAGS_RELEASE         "-fast")
  set (CMAKE_Fortran_FLAGS_DEBUG           "-g  -ffortran-bounds-check")
  set (CMAKE_Fortran_FLAGS_PROFILING       "${CMAKE_Fortran_FLAGS_RELEASE} ")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO  "${CMAKE_Fortran_FLAGS_RELEASE} -g")
  set (CMAKE_Fortran_FLAGS_MINSIZEREL      "-O")

  set (FORTRAN_LIBRARIES         )
  set (FORTRAN_LIBRARIES_FLAG    )

else ()

  message (WARNING "Default flags are not defined for ${CMAKE_Fortran_COMPILER_ID}")

  set (CMAKE_Fortran_FLAGS "")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g")
  set (CMAKE_Fortran_FLAGS_PROFILING       "${CMAKE_Fortran_FLAGS_RELEASE}")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO  "${CMAKE_Fortran_FLAGS_RELEASE}")
  set (CMAKE_Fortran_FLAGS_MINSIZEREL      "-O")

  set (FORTRAN_LIBRARIES                   )
  set (FORTRAN_LIBRARIES_FLAG              )

endif ()

set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" CACHE STRING "Flags used by the compiler during all build types." FORCE)
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" CACHE STRING "Flags used by the compiler during release builds." FORCE)
set (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG}" CACHE STRING "Flags used by the compiler during debug builds."  FORCE)
set (CMAKE_Fortran_FLAGS_PROFILING       "${CMAKE_Fortran_FLAGS_PROFILING}" CACHE STRING "Flags used For profiling."  FORCE)
set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO  "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" CACHE STRING "Flags used by the compiler during release builds with debug info." FORCE)
set (CMAKE_Fortran_FLAGS_MINSIZEREL      "${CMAKE_Fortran_FLAGS_MINSIZEREL}" CACHE STRING "Flags used by the compiler during release builds for minimum size" FORCE)

set (FORTRAN_LIBRARIES "${FORTRAN_LIBRARIES}" CACHE STRING "Fortran libraries" FORCE)
set (FORTRAN_LIBRARIES_FLAG "${FORTRAN_LIBRARIES_FLAG}" CACHE STRING "Fortran libraries flag" FORCE)

mark_as_advanced (CMAKE_Fortran_FLAGS_PROFILING FORTRAN_LIBRARIES FORTRAN_LIBRARIES_FLAG)


#------------------------------------------------------------------------------
# C Default Flags
#------------------------------------------------------------------------------

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")

  link_libraries ("m")

  set (CMAKE_C_FLAGS "-std=gnu99 -fPIC -funsigned-char -pedantic -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused -Wfloat-equal  -Wno-unused-dummy-argument")

  set (CMAKE_C_FLAGS_RELEASE         "-O3")
  set (CMAKE_C_FLAGS_DEBUG           "-O0 -g")
  set (CMAKE_C_FLAGS_PROFILING       "-O3 -pg")
  set (CMAKE_C_FLAGS_RELWITHDEBINFO  "-O3 -g")
  set (CMAKE_C_FLAGS_MINSIZEREL      "-O2 -g")

elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")

  set (CMAKE_C_FLAGS "-std=gnu99 -restrict -fpic -funsigned-char -Wall -Wcheck -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused -wd869,3656")

  set (CMAKE_C_FLAGS_RELEASE "-O3")

  set (CMAKE_C_FLAGS_DEBUG           "-g -O0 -traceback -w2")
  set (CMAKE_C_FLAGS_PROFILING       "${CMAKE_C_FLAGS_RELEASE} -p")
  set (CMAKE_C_FLAGS_RELWITHDEBINFO  "${CMAKE_C_FLAGS_RELEASE} -g")
  set (CMAKE_C_FLAGS_MINSIZEREL      "-O2 -g")

elseif (CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
  set (CMAKE_C_FLAGS "-std=c99 -fPIC -funsigned-char -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused -Wno-empty-translation-unit -Wno-unused-function")

  set (CMAKE_C_FLAGS_RELEASE "-O3")
  set (CMAKE_C_FLAGS_DEBUG "-g -O0")
  set (CMAKE_C_FLAGS_PROFILING       "${CMAKE_C_FLAGS_RELEASE} -p")
  set (CMAKE_C_FLAGS_RELWITHDEBINFO  "-O3 -g")
  set (CMAKE_C_FLAGS_MINSIZEREL      "-O2")

elseif (CMAKE_C_COMPILER_ID MATCHES "XL")

  set (CMAKE_C_FLAGS "-qlanglvl=stdc99 -q64")
  set (CMAKE_C_FLAGS_RELEASE "-O3 -qhot")
  set (CMAKE_C_FLAGS_DEBUG "-g -qfullpath")
  set (CMAKE_C_FLAGS_PROFILING       "${CMAKE_C_FLAGS_RELEASE} -pg -qfullpath")
  set (CMAKE_C_FLAGS_RELWITHDEBINFO  "-O3 -g")
  set (CMAKE_C_FLAGS_MINSIZEREL      "-O2")

  if (${HOSTNAME} STREQUAL "tanit")
    link_directories(/opt/ibm/xlC/13.1.0/lib64)
  endif()

  link_libraries(m)

elseif (CMAKE_C_COMPILER_ID STREQUAL "PGI")

  set (CMAKE_C_FLAGS "-c99 -noswitcherror")

  set (CMAKE_C_FLAGS_RELEASE         "-fast")
  set (CMAKE_C_FLAGS_DEBUG           "-g -Mbounds")
  set (CMAKE_C_FLAGS_PROFILING       "${CMAKE_C_FLAGS_RELEASE} -Mprof=func,lines")
  set (CMAKE_C_FLAGS_RELWITHDEBINFO  "${CMAKE_C_FLAGS_RELEASE} -g")
  set (CMAKE_C_FLAGS_MINSIZEREL      "-O2")

elseif (CMAKE_C_COMPILER_ID STREQUAL "Cray")

  set (CMAKE_C_FLAGS "")

  set (CMAKE_C_FLAGS_RELEASE         "-O3")
  set (CMAKE_C_FLAGS_DEBUG           "-g")
  set (CMAKE_C_FLAGS_PROFILING       "${CMAKE_C_FLAGS_RELEASE} -h profile_generate")
  set (CMAKE_C_FLAGS_RELWITHDEBINFO  "${CMAKE_C_FLAGS_RELEASE} -g")
  set (CMAKE_C_FLAGS_MINSIZEREL      "-O2")

elseif (CMAKE_C_COMPILER_ID STREQUAL "PathScale")

  set (CMAKE_C_FLAGS "-std=gnu99 -fPIC -funsigned-char -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused -Wunused-value")

  set (CMAKE_C_FLAGS_RELEASE         "-0fast")
  set (CMAKE_C_FLAGS_DEBUG           "-g")
  set (CMAKE_C_FLAGS_PROFILING       "${CMAKE_C_FLAGS_RELEASE} ")
  set (CMAKE_C_FLAGS_RELWITHDEBINFO  "${CMAKE_C_FLAGS_RELEASE} -g")
  set (CMAKE_C_FLAGS_MINSIZEREL      "-O2")

else ()

  message (WARNING "Default flags are not defined for ${CMAKE_C_COMPILER_ID}")

  set (CMAKE_C_FLAGS "")
  set (CMAKE_C_FLAGS_RELEASE "-O")
  set (CMAKE_C_FLAGS_DEBUG   "-g")
  set (CMAKE_C_FLAGS_PROFILING       "${CMAKE_C_FLAGS_RELEASE}")
  set (CMAKE_C_FLAGS_RELWITHDEBINFO  "${CMAKE_C_FLAGS_RELEASE}")
  set (CMAKE_C_FLAGS_MINSIZEREL      "-O")

endif ()

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "Flags used by the compiler during all build types." FORCE)
set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "Flags used by the compiler during release builds." FORCE)
set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
set (CMAKE_C_FLAGS_PROFILING       "${CMAKE_C_FLAGS_PROFILING}" CACHE STRING "Flags used For profiling." FORCE)
set (CMAKE_C_FLAGS_RELWITHDEBINFO  "${CMAKE_C_FLAGS_RELWITHDEBINFO}" CACHE STRING "Flags used by the compiler during release builds with debug info." FORCE)
set (CMAKE_C_FLAGS_MINSIZEREL      "${CMAKE_C_FLAGS_MINSIZEREL}" CACHE STRING "Flags used by the compiler during release builds for minimum size" FORCE)

mark_as_advanced (CMAKE_C_FLAGS_PROFILING)


#------------------------------------------------------------------------------
# C++ Default Flags
#------------------------------------------------------------------------------

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

  set (CMAKE_CXX_FLAGS "-std=c++11 -fPIC -funsigned-char -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wunused -Wno-long-long -Wfloat-equal -Wno-unused-dummy-argument")

  set (CMAKE_CXX_FLAGS_RELEASE         "-O3")
  set (CMAKE_CXX_FLAGS_DEBUG           "-O0 -g")
  set (CMAKE_CXX_FLAGS_PROFILING       "-O3 -pg")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O3 -g")
  set (CMAKE_CXX_FLAGS_MINSIZEREL      "-O2 -g")

  set (CXX_LIBRARIES          stdc++)
  set (CXX_LIBRARIES_FLAG        )

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")

  set (CMAKE_CXX_FLAGS "-std=c++11 -fpic -funsigned-char -Wall -Wcheck -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused")

  set (CMAKE_CXX_FLAGS_RELEASE "-O3")

  set (CMAKE_CXX_FLAGS_DEBUG           "-g -O0 -traceback -w2")
  set (CMAKE_CXX_FLAGS_PROFILING       "${CMAKE_CXX_FLAGS_RELEASE} -p")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELEASE} -g")
  set (CMAKE_CXX_FLAGS_MINSIZEREL      "-O2 -g")

  set (CXX_LIBRARIES          -cxxlib)
  set (CXX_LIBRARIES_FLAG        )

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set (CMAKE_CXX_FLAGS "-std=c++11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused -Wempty-translation-unit -Wno-unused-function")
  set (CMAKE_CXX_FLAGS_RELEASE "-O3")
  set (CMAKE_CXX_FLAGS_DEBUG "-g -O0")
  set (CMAKE_CXX_FLAGS_PROFILING       "${CMAKE_CXX_FLAGS_RELEASE} -p")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O3 -g")
  set (CMAKE_CXX_FLAGS_MINSIZEREL      "-O2")

  set (CXX_LIBRARIES             )
  set (CXX_LIBRARIES_FLAG        )

elseif (CMAKE_CXX_COMPILER_ID MATCHES "XL")

  set (CMAKE_CXX_FLAGS " -q64 -qlanglvl=extended0x")
  set (CMAKE_CXX_FLAGS_RELEASE "-O3 -qhot")
  set (CMAKE_CXX_FLAGS_DEBUG "-g -qfullpath")
  set (CMAKE_CXX_FLAGS_PROFILING       "${CMAKE_CXX_FLAGS_RELEASE} -pg -qfullpath")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O3 -g")
  set (CMAKE_CXX_FLAGS_MINSIZEREL      "-O2")

  set(CXX_LIBRARIES stdc++ ibmc++ ${CXX_LIBRARIES})
  if (${HOSTNAME} STREQUAL "tanit")
    link_directories(/opt/ibm/xlC/13.1.0/lib64)
  endif()

#link_libraries (m)

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "PGI")

  set (CMAKE_CXX_FLAGS "-Xa -noswitcherror")

  set (CMAKE_CXX_FLAGS_RELEASE         "-fast")
  set (CMAKE_CXX_FLAGS_DEBUG           "-g -Mbounds")
  set (CMAKE_CXX_FLAGS_PROFILING       "${CMAKE_CXX_FLAGS_RELEASE} -Mprof=func,lines")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELEASE} -g")
  set (CMAKE_CXX_FLAGS_MINSIZEREL      "-O2")

  set (CXX_LIBRARIES             )
  set (CXX_LIBRARIES_FLAG        )

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Cray")

  set (CMAKE_CXX_FLAGS "")

  set (CMAKE_CXX_FLAGS_RELEASE         "-O3")
  set (CMAKE_CXX_FLAGS_DEBUG           "-g")
  set (CMAKE_CXX_FLAGS_PROFILING       "${CMAKE_CXX_FLAGS_RELEASE} -h profile_generate")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELEASE} -g")
  set (CMAKE_CXX_FLAGS_MINSIZEREL      "-O2")

  set (CXX_LIBRARIES             )
  set (CXX_LIBRARIES_FLAG        )

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "PathScale")

  set (CMAKE_CXX_FLAGS "-ansi -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused -Wunused-value")

  set (CMAKE_CXX_FLAGS_RELEASE         "-0fast")
  set (CMAKE_CXX_FLAGS_DEBUG           "-g")
  set (CMAKE_CXX_FLAGS_PROFILING       "${CMAKE_CXX_FLAGS_RELEASE} ")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELEASE} -g")
  set (CMAKE_CXX_FLAGS_MINSIZEREL      "-O2")

  set (CXX_LIBRARIES             )
  set (CXX_LIBRARIES_FLAG        )

else ()

  message (WARNING "Default flags are not defined for ${CMAKE_CXX_COMPILER_ID}")

  set (CMAKE_CXX_FLAGS "")
  set (CMAKE_CXX_FLAGS_RELEASE "-O")
  set (CMAKE_CXX_FLAGS_DEBUG   "-g")
  set (CMAKE_CXX_FLAGS_PROFILING       "${CMAKE_CXX_FLAGS_RELEASE}")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELEASE}")
  set (CMAKE_CXX_FLAGS_MINSIZEREL      "-O")

  set (CXX_LIBRARIES             )
  set (CXX_LIBRARIES_FLAG        )

endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "Flags used by the compiler during all build types." FORCE)
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "Flags used by the compiler during release builds." FORCE)
set (CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "Flags used by the compiler during debug builds." FORCE)
set (CMAKE_CXX_FLAGS_PROFILING       "${CMAKE_CXX_FLAGS_PROFILING}" CACHE STRING "Flags used For profiling." FORCE)
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "Flags used by the compiler during release builds with debug info." FORCE)
set (CMAKE_CXX_FLAGS_MINSIZEREL      "${CMAKE_CXX_FLAGS_MINSIZEREL}" CACHE STRING "Flags used by the compiler during release builds for minimum size" FORCE)

set (CXX_LIBRARIES "${CXX_LIBRARIES}" CACHE STRING "C++ libraries" FORCE)
set (CXX_LIBRARIES_FLAG "${CXX_LIBRARIES_FLAG}" CACHE STRING "C++ flags" FORCE)

set (PASS_DEFAULT_FLAGS 1 CACHE STRING "")
mark_as_advanced (CMAKE_CXX_FLAGS_PROFILING CXX_LIBRARIES CXX_LIBRARIES_FLAG PASS_DEFAULT_FLAGS)
endif()
