#############################################################################
# This file was found and adapted from original LLNL tools called perf-dump.
#
# Copyright (c) 2013-2014, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# This file is part of perf-dump.
# Written by Todd Gamblin, tgamblin@llnl.gov, All rights reserved.
# LLNL-CODE-647187
#
#For details, see https://scalability-llnl.github.io/perf-dump
#
#############################################################################
#
# Try to find HWLOC headers and libraries.
#
# Usage of this module as follows:
#
#     find_package(HWLOC)
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#
#  HWLOC_ROOT         Set this environment variable to the root installation of
#                    libHWLOC if the module has problems finding the
#                    proper installation path.
#
# Variables defined by this module:
#
#  HWLOC_FOUND              System has HWLOC libraries and headers
#  HWLOC_LIBRARY            The HWLOC library
#  HWLOC_INCLUDE_DIR        The location of HWLOC headers

find_library(HWLOC_LIBRARY
    NAMES libhwloc.so libhwloc.a hwloc
    HINTS ENV HWLOC_ROOT
    PATH_SUFFIXES lib lib64
)

find_path(HWLOC_INCLUDE_DIR
    NAMES hwloc.h
    HINTS ENV HWLOC_ROOT
    PATH_SUFFIXES include
)

if (HWLOC_INCLUDE_DIR AND HWLOC_LIBRARY)
  set(HWLOC_LIBRARIES ${HWLOC_LIBRARY})
else()

  set(HWLOC_INCLUDE_DIR "")
  set(HWLOC_LIBRARIES   "")
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HWLOC DEFAULT_MSG HWLOC_LIBRARIES HWLOC_INCLUDE_DIR)

mark_as_advanced(HWLOC_INCLUDE_DIR HWLOC_LIBRARIES)
