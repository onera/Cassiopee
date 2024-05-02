# ------------------------------------------------------------------------------
# General CMake settings
# ------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
cmake_policy(SET CMP0074 NEW) # force find_package to take <PackageName>_ROOT variables into account

if (NOT DEFINED PROJECT_ROOT)
  if (DEFINED CASSIOPEE_ROOT)
    set(PROJECT_ROOT ${CASSIOPEE_ROOT} CACHE PATH)
  else()
    set(PROJECT_ROOT ${PROJECT_SOURCE_DIR} CACHE PATH)
  endif()
endif()

set(PROJECT_UTILS_DIR ${PROJECT_ROOT}/ThirdParty/project_utils)
set(PROJECT_UTILS_CMAKE_DIR ${PROJECT_UTILS_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_UTILS_CMAKE_DIR}/find_package") # for custom Find*.cmake files
include(${PROJECT_UTILS_CMAKE_DIR}/dependency_management.cmake) # project_add_subdirectory, project_find_package and target_install

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_Fortran_PREPROCESS ON)

# ------------------------------------------------------------------------------
# Project
# ------------------------------------------------------------------------------
project(
  KCore VERSION 3.3
  DESCRIPTION "Cassiopee KCore"
  LANGUAGES C CXX Fortran
)

## Compiler flags
### C++ standard
set(CMAKE_CXX_STANDARD          14 )
set(CMAKE_CXX_EXTENSIONS        OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON )
### fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
### Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE) # default to Release
endif()

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------

## Dependencies installed in the environment ##
project_find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)

# ------------------------------------------------------------------------------
# Compilation
# ------------------------------------------------------------------------------
set(src_dir ${PROJECT_SOURCE_DIR}/${PROJECT_NAME})
set(test_dir ${PROJECT_SOURCE_DIR}/test)
set(include_dir "${PROJECT_SOURCE_DIR}/${PROJECT_NAME}")

# TODO list files here in CMakeLists.txt once we will be using CMake only
execute_process(
  COMMAND python print_source_filenames.py
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE src_files
)
string(STRIP "${src_files}" src_files) # remove final empty line of previous command
add_library(${PROJECT_NAME} SHARED "${src_files}")

target_include_directories(${PROJECT_NAME}
  PUBLIC
    $<BUILD_INTERFACE:${include_dir}>
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
  PRIVATE
    $<BUILD_INTERFACE:${include_dir}/Metis>
)
target_compile_definitions(${PROJECT_NAME}
  PUBLIC
    INTEGER_E=INTEGER*4
    REAL_E=REAL*8
)
target_link_libraries(${PROJECT_NAME}
  PUBLIC
    Python::Python
    Python::NumPy
)

# ------------------------------------------------------------------------------
# Installation
# ------------------------------------------------------------------------------
target_install(${PROJECT_NAME})


# ------------------------------------------------------------------------------
# Python and wrapping
# ------------------------------------------------------------------------------
include(${PROJECT_UTILS_CMAKE_DIR}/python_wrapping.cmake)
compile_install_cpython_module(${PROJECT_NAME} "${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/kcore.cpp")
install_python_modules(${PROJECT_NAME})

# TODO DEL: will not be useful once using CMake only
install(
  FILES
    "${PROJECT_SOURCE_DIR}/Dist.py"
    "${PROJECT_SOURCE_DIR}/config.py"
  DESTINATION "${SITE_PACKAGES_OUTPUT_DIRECTORY}/${PROJECT_NAME}"
  COMPONENT   "python installation utils"
)
install(
  FILES
    "${PROJECT_SOURCE_DIR}/config.py"
    "${PROJECT_SOURCE_DIR}/installBase.py"
    "${PROJECT_SOURCE_DIR}/print_source_filenames.py"
  DESTINATION "${SITE_PACKAGES_OUTPUT_DIRECTORY}"
  COMPONENT   "python installation utils"
)

# ------------------------------------------------------------------------------
# Testing
# ------------------------------------------------------------------------------
include(CTest)
add_subdirectory(test)
