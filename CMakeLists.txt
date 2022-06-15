# ------------------------------------------------------------------------------
# General CMake settings
# ------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
cmake_policy(SET CMP0074 NEW) # force find_package to take <PackageName>_ROOT variables into account

if(NOT DEFINED PROJECT_ROOT)
  set(PROJECT_ROOT ${CMAKE_SOURCE_DIR} CACHE PATH "Root directory, where the submodules are populated")
endif()

set(CASSIOPEE_ROOT ${PROJECT_ROOT})

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
  Cassiopee VERSION 3.3
  DESCRIPTION "Cassiopee suite"
  LANGUAGES C CXX Fortran
)

project_find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)
project_add_subdirectory(KCore     )
project_add_subdirectory(Converter )
project_add_subdirectory(Compressor)
project_add_subdirectory(Connector )
project_add_subdirectory(Dist2Walls)
project_add_subdirectory(Generator )
project_add_subdirectory(Post      )
project_add_subdirectory(Transform )

add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME}
  INTERFACE
    KCore::KCore
    Converter::Converter

    Compressor::Compressor
    Connector::Connector
    Dist2Walls::Dist2Walls
    Generator::Generator
    Post::Post
    Transform::Transform
)
