# Find the Cython compiler.
#
# This code sets the following variables:
#
#  CYTHON_EXECUTABLE
#  CYTHON_VERSION
#
# See also UseCython.cmake

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

# Use the Cython executable that lives next to the Python executable
# if it is a local installation.
find_package( PythonInterp )
if( PYTHONINTERP_FOUND )
  get_filename_component( _python_path ${PYTHON_EXECUTABLE} PATH )
  find_program( CYTHON_EXECUTABLE
    NAMES cython cython.bat cython3
    HINTS ${_python_path}
    )
else()
  find_program( CYTHON_EXECUTABLE
    NAMES cython cython.bat cython3
    )
endif()

if (CYTHON_EXECUTABLE)

    EXECUTE_PROCESS(COMMAND
      ${CYTHON_EXECUTABLE} --version
      OUTPUT_VARIABLE CYTHON_VERSION
      ERROR_VARIABLE CYTHON_VERSION_ERROR
      )

    if (CYTHON_VERSION_ERROR MATCHES "^Cython version")
        string (STRIP ${CYTHON_VERSION_ERROR} CYTHON_VERSION_STRIP)
        string(REGEX REPLACE "Cython version "  "" CYTHON_VERSION ${CYTHON_VERSION_STRIP})
        include( FindPackageHandleStandardArgs )
        FIND_PACKAGE_HANDLE_STANDARD_ARGS(Cython
                                     FOUND_VAR CYTHON_FOUND
                                     REQUIRED_VARS CYTHON_EXECUTABLE CYTHON_VERSION
                                     VERSION_VAR CYTHON_VERSION)

       if (CYTHON_FIND_VERSION)

           if (${CYTHON_VERSION} VERSION_LESS $Cython_FIND_VERSION})
  	          MESSAGE(FATAL_ERROR
	              "Cython version " ${Cython_VERSION}
	              " is less than required version " ${Cython_FIND_VERSION}
	              )
    
           endif()
        endif()

        mark_as_advanced( CYTHON_EXECUTABLE CYTHON_VERSION)

   else ()

       if (CYTHON_FIND_REQUIRED)
           MESSAGE(SEND_ERROR
              "Required Cython python module not found")

       endif ()
       
   endif ()
   
else ()

  if (CYTHON_FIND_REQUIRED)
     MESSAGE(SEND_ERROR
            "Required Cython python module not found")
  endif()

endif ()


