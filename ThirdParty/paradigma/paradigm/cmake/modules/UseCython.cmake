# Define a function to create Cython modules.
#
# For more information on the Cython project, see http://cython.org/.
# "Cython is a language that makes writing C extensions for the Python language
# as easy as Python itself."
#
# This file defines a CMake function to build a Cython Python module.
# To use it, first include this file.
#
#   include( UseCython )
#
# Then call cython_add_module to create a module.
#
#   cython_add_module( <module_name> <src1> <src2> ... <srcN> )
#
# To create a standalone executable, the function
#
#   cython_add_standalone_executable( <executable_name> [MAIN_MODULE src1] <src1> <src2> ... <srcN> )
#
# To avoid dependence on Python, set the PYTHON_LIBRARY cache variable to point
# to a static library.  If a MAIN_MODULE source is specified,
# the "if __name__ == '__main__':" from that module is used as the C main() method
# for the executable.  If MAIN_MODULE, the source with the same basename as
# <executable_name> is assumed to be the MAIN_MODULE.
#
# Where <module_name> is the name of the resulting Python module and
# <src1> <src2> ... are source files to be compiled into the module, e.g. *.pyx,
# *.py, *.c, *.cxx, etc.  A CMake target is created with name <module_name>.  This can
# be used for target_link_libraries(), etc.
#
# The sample paths set with the CMake include_directories() command will be used
# for include directories to search for *.pxd when running the Cython complire.
#
# Cache variables that effect the behavior include:
#
#  CYTHON_ANNOTATE
#  CYTHON_NO_DOCSTRINGS
#  CYTHON_FLAGS
#
# Source file properties that effect the build process are
#
#  CYTHON_IS_CXX
#
# If this is set of a *.pyx file with CMake set_source_files_properties()
# command, the file will be compiled as a C++ file.
#
# See also FindCython.cmake

#=============================================================================
# Copyright 2011 Kitware, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#=============================================================================

# Configuration options.
set( CYTHON_ANNOTATE OFF
  CACHE BOOL "Create an annotated .html file when compiling *.pyx." )
set( CYTHON_NO_DOCSTRINGS OFF
  CACHE BOOL "Strip docstrings from the compiled module." )
set( CYTHON_FLAGS "" CACHE STRING
  "Extra flags to the cython compiler." )
mark_as_advanced( CYTHON_ANNOTATE CYTHON_NO_DOCSTRINGS CYTHON_FLAGS )

find_package( Cython REQUIRED)
find_package( PythonLibs REQUIRED )

set( CYTHON_CXX_EXTENSION "cxx" )
set( CYTHON_C_EXTENSION "c" )

# Create a *.c or *.cxx file from a *.pyx file.
# Input the generated file basename.  The generate file will put into the variable
# placed in the "generated_file" argument. Finally all the *.py and *.pyx files.
include( CMakeParseArguments )
function( COMPILE_PYX _name generated_file)
  set(options)
  set(oneValueArgs)
  set(multiValueArgs PYX_SOURCES INCLUDE_DIRECTORIES)
  cmake_parse_arguments( COMPILE_PYX "${options}" "${oneValueArgs}"
                        "${multiValueArgs}" ${ARGN} )

  set( cxx_arg "" )
  set( extension ${CYTHON_C_EXTENSION} )
  set( pyx_lang "C" )
  set( comment "Compiling Cython C source for ${_name}..." )

  # Determining generated file name.
  set( _generated_file "${CMAKE_CURRENT_BINARY_DIR}/${_name}.${extension}" )
  set_source_files_properties( ${_generated_file} PROPERTIES GENERATED TRUE )
  set( ${generated_file} ${_generated_file} PARENT_SCOPE )

  set( include_directory_arg "" )
  if (COMPILE_PYX_INCLUDE_DIRECTORIES)
    list( REMOVE_DUPLICATES COMPILE_PYX_INCLUDE_DIRECTORIES)
    foreach( _include_dir ${COMPILE_PYX_INCLUDE_DIRECTORIES} )
      set( include_directory_arg ${include_directory_arg} "-I" "${_include_dir}" )
    endforeach()
  endif()

  set( pyx_locations "" )
  if (NOT COMPILE_PYX_PYX_SOURCES)
    message (FATAL_ERROR "No PYX_SOURCES to compile")
  endif()

  foreach( pyx_file ${COMPILE_PYX_PYX_SOURCES} )

    get_filename_component( pyx_file_basename "${pyx_file}" NAME_WE )

    # Determine if it is a C or C++ file.
    get_source_file_property( property_is_cxx ${pyx_file} CYTHON_IS_CXX )
    if( ${property_is_cxx} )
      set( cxx_arg "--cplus" )
      set( extension ${CYTHON_CXX_EXTENSION} )
      set( pyx_lang "CXX" )
      set( comment "Compiling Cython CXX source for ${_name}..." )
    endif()

    get_source_file_property( pyx_location ${pyx_file} LOCATION )
    list( APPEND pyx_locations "${pyx_location}" )

  endforeach()

  # Set additional flags.
  if( CYTHON_ANNOTATE )
    set( annotate_arg "--annotate" )
  endif()

  if( CYTHON_NO_DOCSTRINGS )
    set( no_docstrings_arg "--no-docstrings" )
  endif()

  if( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR
        "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo" )
      set( cython_debug_arg "--gdb" )
  endif()

  if( "${PYTHONLIBS_VERSION_STRING}" MATCHES "^2." )
    set( version_arg "-2" )
  elseif( "${PYTHONLIBS_VERSION_STRING}" MATCHES "^3." )
    set( version_arg "-3" )
  else()
    set( version_arg )
  endif()

  add_custom_command( OUTPUT ${_generated_file}
    COMMAND ${CYTHON_EXECUTABLE}
    ARGS ${cxx_arg} ${include_directory_arg} ${version_arg}
    ${annotate_arg} ${no_docstrings_arg} ${cython_debug_arg} ${CYTHON_FLAGS}
    --output-file  ${_generated_file} ${pyx_locations}
    DEPENDS ${pyx_locations}
    IMPLICIT_DEPENDS ${pyx_lang}
    COMMENT ${comment}
   )

endfunction()


# cython_add_module( <name> src1 src2 ... srcN )
# Build the Cython Python module.
function( CYTHON_ADD_MODULE _name )
  set(options)
  set(oneValueArgs)
  set(multiValueArgs PYX_SOURCES OTHER_SOURCES INCLUDE_DIRECTORIES)
  cmake_parse_arguments( CYTHON_ADD_MODULE "${options}" "${oneValueArgs}"
                        "${multiValueArgs}" ${ARGN} )

  if (NOT CYTHON_ADD_MODULE_PYX_SOURCES)
    message (${CYTHON_ADD_MODULE_PYX_SOURCES})
    message (FATAL_ERROR "No PYX_SOURCES to compile")
  endif()

  compile_pyx( ${_name} generated_file
                PYX_SOURCES ${CYTHON_ADD_MODULE_PYX_SOURCES}
                INCLUDE_DIRECTORIES ${CYTHON_ADD_MODULE_INCLUDE_DIRECTORIES})

  include_directories( ${PYTHON_INCLUDE_DIRS} )
  python_add_module( ${_name} ${generated_file} ${CYTHON_ADD_MODULE_OTHER_SOURCES})
  target_include_directories(${_name} PRIVATE ${CYTHON_ADD_MODULE_INCLUDE_DIRECTORIES})

  if( APPLE )
    set_target_properties( ${_name} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup" )
  else()
    target_link_libraries( ${_name} ${PYTHON_LIBRARIES} )
  endif()
endfunction()


include( CMakeParseArguments )
# cython_add_standalone_executable( _name [MAIN_MODULE src3.py] src1 src2 ... srcN )
# Creates a standalone executable the given sources.
function( cython_add_standalone_executable _name )
  set( pyx_module_sources "" )
  set( other_module_sources "" )
  set( main_module "" )
  cmake_parse_arguments( cython_arguments "" "MAIN_MODULE" "" ${ARGN} )
  include_directories( ${PYTHON_INCLUDE_DIRS} )
  foreach( _file ${cython_arguments_UNPARSED_ARGUMENTS} )
    if( ${_file} MATCHES ".*\\.py[x]?$" )
      get_filename_component( _file_we ${_file} NAME_WE )
      if( "${_file_we}" STREQUAL "${_name}" )
        set( main_module "${_file}" )
      elseif( NOT "${_file}" STREQUAL "${cython_arguments_MAIN_MODULE}" )
        set( PYTHON_MODULE_${_file_we}_static_BUILD_SHARED OFF )
        compile_pyx( "${_file_we}_static" generated_file "${_file}" )
        list( APPEND pyx_module_sources "${generated_file}" )
      endif()
    else()
      list( APPEND other_module_sources ${_file} )
    endif()
  endforeach()

  if( cython_arguments_MAIN_MODULE )
    set( main_module ${cython_arguments_MAIN_MODULE} )
  endif()
  if( NOT main_module )
    message( FATAL_ERROR "main module not found." )
  endif()
  get_filename_component( main_module_we "${main_module}" NAME_WE )
  set( CYTHON_FLAGS ${CYTHON_FLAGS} --embed )
  compile_pyx( "${main_module_we}_static" generated_file ${main_module} )
  add_executable( ${_name} ${generated_file} ${pyx_module_sources} ${other_module_sources} )
  target_link_libraries( ${_name} ${PYTHON_LIBRARIES} ${pyx_module_libs} )
endfunction()
