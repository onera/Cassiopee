#  Mpi4Py_FOUND       = Set to TRUE if Mpi4Py is found
#  Mpi4Py_VERSION     = Mpi4Py version number
#  Mpi4Py_INCLUDE_DIR = Path to Mpi4Py include files
#  Target Mpi4Py::Mpi4Py

if(NOT Python_EXECUTABLE)
  if(Mpi4Py_FIND_REQUIRED)
    message(SEND_ERROR
      "Python executable not found, so required Mpi4Py module not found"
      )
  endif()
endif()

# Continue processing if python executable is known
# Retrieve the Mpi4Py version
execute_process(COMMAND
  ${Python_EXECUTABLE} -c "import mpi4py; print(mpi4py.__version__)"
  OUTPUT_VARIABLE Mpi4Py_VERSION
  ERROR_VARIABLE  Mpi4Py_VERSION_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

# If Mpi4Py_VERSION_ERROR does not exist, then we know Mpi4Py exists;
# now look for the Mpi4Py include directory
if(NOT Mpi4Py_VERSION_ERROR)
  execute_process(COMMAND
    ${Python_EXECUTABLE} -c "import mpi4py; print(mpi4py.get_include())"
    OUTPUT_VARIABLE Mpi4Py_INCLUDE_DIR
    ERROR_VARIABLE  Mpi4Py_INCLUDE_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Handle the QUIETLY and REQUIRED arguments and set Mpi4Py_FOUND to
  # TRUE if all listed variables are TRUE
  include(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(Mpi4Py
    FOUND_VAR     Mpi4Py_FOUND
    REQUIRED_VARS Mpi4Py_VERSION Mpi4Py_INCLUDE_DIR
    VERSION_VAR   Mpi4Py_VERSION)

  # Version checking: If a version check is requested, compare
  # Mpi4Py_VERSION to the requested version
  if(Mpi4Py_FIND_VERSION)
    if(${Mpi4Py_VERSION} VERSION_LESS ${Mpi4Py_FIND_VERSION})
      message(FATAL_ERROR
        "Mpi4Py version " ${Mpi4Py_VERSION}
        " is less than required version " ${Mpi4Py_FIND_VERSION}
      )
    endif()
  endif()

  if (Mpi4Py_FOUND AND NOT TARGET Mpi4Py::Mpi4Py)
    add_library(Mpi4Py::Mpi4Py INTERFACE IMPORTED)
    set_target_properties(Mpi4Py::Mpi4Py
      PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${Mpi4Py_INCLUDE_DIR}")
  endif ()

# A Mpi4Py version error means that Mpi4Py was not found
else()
  if(Mpi4Py_FIND_REQUIRED)
    message(SEND_ERROR "Required Mpi4Py python module not found")
  endif()
endif()
