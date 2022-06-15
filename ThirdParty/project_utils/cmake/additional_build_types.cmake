cmake_host_system_information(RESULT HOSTNAME QUERY HOSTNAME)

# ----------------------------------------------------------------------------------------------------------------------
# Additionnal CMAKE_BUILD_TYPE: Profiling and Sanitize
# ----------------------------------------------------------------------------------------------------------------------

# Note on how it works (CMake conventions):
#   We set predefined flags inside CMAKE_<lang>_FLAGS_PROFILING and CMAKE_<lang>_FLAGS_SANITIZE.
#   Then, calling `cmake -D CMAKE_BUILD_TYPE=build_type` will trigger CMake to add
#   CMAKE_<lang>_FLAGS_<upper-case build_type> to the other CMAKE_<lang>_FLAGS flags
# SEE https://cmake.org/cmake/help/latest/manual/cmake-buildsystem.7.html#id37

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS_SANITIZE "-O0 -g -fcheck=bounds -fbacktrace -fsanitize=address -fno-omit-frame-pointer -Wall -Wextra")
  set (CMAKE_C_FLAGS_SANITIZE "-O0 -g -fsanitize=address -fno-omit-frame-pointer -Wall -Wextra")
  set (CMAKE_CXX_FLAGS_SANITIZE "-O0 -g -fsanitize=address -fno-omit-frame-pointer -Wall -Wextra ")
  # Other interesting flags
  # for python : export LD_PRELOAD=$(gcc -print-file-name=libasan.so)
  #-ftree-vectorize -ftree-loop-vectorize -fvect-cost-model=unlimited -mprefer-vector-width=512
  #-ftree-loop-optimize -ftree-vectorize -fopt-info -fopt-info-all
  #-fprofile-arcs -ftest-coverage --coverage

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set (CMAKE_CXX_FLAGS_SANITIZE "-g -O0 -traceback -w2")
  #set (CXX_LIBRARIES -cxxlib)

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set (CMAKE_CXX_FLAGS_SANITIZE "-O0 -g -fsanitize=address -fno-omit-frame-pointer")

else()
  message (WARNING "Unknown C++ compiler ${CMAKE_CXX_COMPILER_ID}")
endif()

set(CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_RELEASE} -march=native -p" CACHE STRING "Flags used for profiling." FORCE)
set(CMAKE_CXX_FLAGS_SANITIZE  "${CMAKE_CXX_FLAGS_SANITIZE}" CACHE STRING "Flags used by the compiler during sanitize builds" FORCE)

mark_as_advanced (CMAKE_CXX_FLAGS_PROFILING CMAKE_CXX_FLAGS_SANITIZE)

set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Profiling Sanitize."
    FORCE)
