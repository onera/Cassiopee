
#==============================================================================
#
#                Top level CMakeLists.txt file for PARADIGMA
#
#==============================================================================

cmake_minimum_required(VERSION 3.0)
cmake_policy(VERSION 3.0)

enable_testing()

#------------------------------------------------------------------------------
# Internal CMake module (shared with pk)
#------------------------------------------------------------------------------

set(PDMA_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/paradigm/cmake")
set(CMAKE_MODULE_PATH "${PDMA_CMAKE_DIR}/modules")

include(${CMAKE_MODULE_PATH}/version.cmake)
project(PDMA VERSION ${PDM_DEF_VERSION} LANGUAGES C CXX Fortran)

if (NOT PDM_IN_PDMA)
  set(PDM_IN_PDMA 1)
  set(PDM_IN_PDMA_BOOL 1)
endif ()

#------------------------------------------------------------------------------
# PDM version from PDM
#------------------------------------------------------------------------------

set(PDM_VERSION_MAJOR ${PDMA_VERSION_MAJOR})
set(PDM_VERSION_MINOR ${PDMA_VERSION_MINOR})
set(PDM_VERSION_PATCH ${PDMA_VERSION_PATCH})

set(PDM_VERSION "${PDM_VERSION_MAJOR}.${PDM_VERSION_MINOR}.${PDM_VERSION_PATCH}")

#------------------------------------------------------------------------------
# Enable languages
#------------------------------------------------------------------------------

enable_language(Fortran)
enable_language(C)
enable_language(CXX)

#------------------------------------------------------------------------------
# User option Definition
#------------------------------------------------------------------------------

# OpenMP
# ------

option(PDM_ENABLE_OPENMP "Enable OpenMP" OFF)

# Enable Aniso_Agglo plugin
# ----------------

option(PDM_ENABLE_ANISO_AGGLO "Enable Anisotropic agglomeration for coarse mesh function" ON)

if (PDM_ENABLE_ANISO_AGGLO)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
      message(FATAL_ERROR "Insufficient g++ version for Anisotropic agglomeration : set PDM_ENABLE_ANISO_AGGLO to OFF and re-run")
    endif()
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15)
      message(FATAL_ERROR "Insufficient icpc version for Anisotropic agglomeration : set PDM_ENABLE_ANISO_AGGLO to OFF and re-run")
    endif()
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.2)
      message(FATAL_ERROR "Insufficient clang version for Anisotropic agglomeration : set PDM_ENABLE_ANISO_AGGLO to OFF and re-run")
    endif()
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "XL")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.1.3)
      message(FATAL_ERROR "Insufficient xlc++ version for Anisotropic agglomeration : set PDM_ENABLE_ANISO_AGGLO to OFF and re-run")
    endif()
  endif ()
  set (PDM_HAVE_ANISO_AGGLO 1)
endif ()


# Enable long for absolute number
# -------------------------------

option(PDM_ENABLE_LONG_G_NUM "Enable long for global number" OFF)

# Pt-scotch
# ---------

option(PDM_ENABLE_PTSCOTCH "Compile with support for Pt-SCOTCH." ON)

# ParMetis
# --------

option(PDM_ENABLE_PARMETIS "Compile with support for ParMetis." ON)

# Shared libraries
# ----------------

option(PDM_ENABLE_SHARED "Build Shared Libraries" OFF)

# Static libraries
# ----------------

option(PDM_ENABLE_STATIC "Build Static Libraries" ON)

# Check MPI Wrapper
# -----------------

option(PDM_ENABLE_MPI_CHECK "Check MPI Wrapper" ON)

# Python Bindings
# ---------------

option(PDM_ENABLE_PYTHON_BINDINGS "Build Python module" OFF)

# Fortran interface
# -----------------

option(PDM_ENABLE_Fortran "Build Fortran Interface" ON)

if (PDM_ENABLE_STATIC AND PDM_ENABLE_Fortran)
   set(PDM_ENABLE_STATIC_Fortran_INTERFACE ON)
endif ()

if (PDM_ENABLE_SHARED AND PDM_ENABLE_Fortran)
   set(PDM_ENABLE_SHARED_Fortran_INTERFACE ON)
endif ()

#------------------------------------------------------------------------------
# Default build : Release
#------------------------------------------------------------------------------

if ((PDM_ENABLE_PYTHON_BINDINGS OR PDM_ENABLE_SHARED_Fortran_INTERFACE) AND (NOT PDM_ENABLE_SHARED))
   set (PDM_ENABLE_SHARED "ON" CACHE STRING "Build Shared Libraries" FORCE)
endif()

if ((NOT PDM_ENABLE_STATIC) AND (NOT PDM_ENABLE_SHARED))
   message (FATAL_ERROR "No enabled library. Please, enable shared or static library")
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING  "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel")
endif()

# ADD "-DDEBUG_CLASSE" in Debug
#------------------------------

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions (-DDEBUG_CLASSE)
endif()

# Add "-DHAVE_GETRUSAGE" if function getrusage and gettimeofday exist
#--------------------------------------------------------------------

include(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(getrusage PDM_HAVE_GETRUSAGE)
CHECK_FUNCTION_EXISTS(gettimeofday PDM_HAVE_GETTIMEOFDAY)

#------------------------------------------------------------------------------
# Default compiler flags
#------------------------------------------------------------------------------

include(${CMAKE_MODULE_PATH}/default_flags.cmake)

#------------------------------------------------------------------------------
# Libraries to add to link
#------------------------------------------------------------------------------

set(COMMON_LIBRARIES)

#------------------------------------------------------------------------------
# Check dependecies
#------------------------------------------------------------------------------

# Check MPI (TODO : check MPI wrapper function check )
# -------------

find_package(MPI REQUIRED)

mark_as_advanced (MPI_EXTRA_LIBRARY MPI_LIBRARY)

# ---------------------
# Check MPI Wrapper
# ---------------------

if (PDM_ENABLE_MPI_CHECK)
  include(${CMAKE_MODULE_PATH}/CheckMPI.cmake)
endif()

#------------------------------------------------------------------------------
# Check dependencies
#------------------------------------------------------------------------------

# Checck SCOTCH an ParMETIS
# -------------------------

if (PDM_ENABLE_PARMETIS)
  find_package(ParMETIS 4.0.3)
  if (PARMETIS_FOUND)
    set(PDM_HAVE_PARMETIS 1)
  endif()
endif()

if (PDM_ENABLE_PTSCOTCH)
  find_package(PTSCOTCH 6.0.0)
  if (PTSCOTCH_FOUND)
    set(PDM_HAVE_PTSCOTCH 1)
  endif()
endif()

# Check OpenMP (TODO : check MPI wrapper function check )
# -------------------------------------------------------

if (PDM_ENABLE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(PDM_HAVE_OPENMP 1)
  endif()
endif()

# Python bindings
# ---------------

if (PDM_ENABLE_PYTHON_BINDINGS)
  find_package(PythonInterp REQUIRED)
  find_package(PythonLibs REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(Mpi4Py REQUIRED)
  include(UseCython)
  if (NOT CMAKE_PDM_INSTALL_PYTHON_DIR)

     set (CMAKE_PDM_INSTALL_PYTHON_DIR  "${CMAKE_INSTALL_PREFIX}")
     set (CMAKE_PDM_INSTALL_PYTHON_DIR
         ${CMAKE_PDM_INSTALL_PYTHON_DIR} CACHE STRING "Install path directory for Python packages")

  endif ()

endif()

#------------------------------------------------------------------------------
# Store variable in configure file
#------------------------------------------------------------------------------

if (PDM_ENABLE_LONG_G_NUM)
   set(PDM_LONG_G_NUM 1)
   set(PDM_LONG_G_NUM_BOOL 1)
else ()
   set(PDM_LONG_G_NUM_BOOL 0)
endif ()

configure_file(pdm_config.h.in "${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h")

if (PDM_ENABLE_STATIC_Fortran_INTERFACE OR PDM_ENABLE_SHARED_Fortran_INTERFACE)
  if ("${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h" IS_NEWER_THAN  "${CMAKE_CURRENT_BINARY_DIR}/pdm_configf.h")
    include (FortranizeFile)
    fortranize ("${CMAKE_CURRENT_BINARY_DIR}/pdm_config.h")
  endif()
endif ()

#------------------------------------------------------------------------------
# Doc
#------------------------------------------------------------------------------

find_package(Doxygen)

#------------------------------------------------------------------------------
# Print summary of found and not found optional packages
#------------------------------------------------------------------------------

include(FeatureSummary)
FEATURE_SUMMARY(WHAT ALL)

#------------------------------------------------------------------------------
# Installation de PDM
#------------------------------------------------------------------------------

# Set PDM install sub-directories
set(PDMA_BIN_DIR "bin")
set(PDMA_LIB_DIR "lib")
set(PDMA_INCLUDE_DIR "include")
set(PDMA_PKGCONFIG_DIR "lib/pkgconfig")
set(PDMA_SHARE_DIR "share/cedre")
set(PDMA_MAN_DIR "share/man")
set(PDMA_DOC_DIR "${PDM_SHARE_DIR}/doc")
set(PDMA_ETC_DIR "etc")

#------------------------------------------------------------------------------
# Ajout des sources de tous les elements du projet
# En fonction des top_targets a installer
#------------------------------------------------------------------------------

add_subdirectory(paradigm)

add_subdirectory(src)

add_subdirectory(test)

add_subdirectory(doc)

#------------------------------------------------------------------------------
# Affichage du message de fin d'installation
#------------------------------------------------------------------------------

add_subdirectory(${PDMA_CMAKE_DIR}/post-install)
