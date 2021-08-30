cmake_minimum_required(VERSION 3.12)

find_path(ParMETIS_INCLUDE_DIR
  NAMES
    parmetis.h
  PATHS ENV
    PARMETIS_ROOT # note: cmake is only searching ParMETIS_ROOT (case-sensitive)
  PATH_SUFFIXES
    include
  DOC "ParMETIS include directory")
mark_as_advanced(ParMETIS_INCLUDE_DIR)

find_library(ParMETIS_LIBRARY
  NAMES
    parmetis
  PATHS ENV
    PARMETIS_ROOT
  PATH_SUFFIXES
    lib
  DOC "ParMETIS library")
mark_as_advanced(ParMETIS_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ParMETIS
  REQUIRED_VARS ParMETIS_LIBRARY ParMETIS_INCLUDE_DIR
)

if(ParMETIS_FOUND AND NOT TARGET ParMETIS::ParMETIS)
  add_library(ParMETIS::ParMETIS UNKNOWN IMPORTED)
  set_target_properties(ParMETIS::ParMETIS PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
    IMPORTED_LOCATION "${ParMETIS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${ParMETIS_INCLUDE_DIR}"
  )
endif()
