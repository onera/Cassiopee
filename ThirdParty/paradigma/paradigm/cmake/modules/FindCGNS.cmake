#
# Find the native CGNS includes and library
#
# CGNS_FOUND       - True if cgns was found.
# CGNS_INCLUDE_DIRS- where to find cgns.h, etc.
# CGNS_LIBRARIES   - List of fully qualified libraries to link against when using CGNS.

FIND_PATH(CGNS_INCLUDE_DIRS cgnslib.h
  HINTS ${CGNS}/include $ENV{CGNS_DIR}/include
  DOC "Directory where the CGNS header is located"
  )

if(CGNS_USE_STATIC_LIBRARIES)
  FIND_LIBRARY(CGNS_LIBRARIES
    NAMES libcgns.a
    HINTS ${CGNS}/lib $ENV{CGNS_DIR}/lib
    NO_DEFAULT_PATH
    DOC "The CGNS library"
    )
  FIND_LIBRARY(CGNS_LIBRARIES cgns
    DOC "The CGNS library"
    )
else()
  FIND_LIBRARY(CGNS_LIBRARIES cgns
    HINTS ${CGNS}/lib $ENV{CGNS_DIR}/lib
    NO_DEFAULT_PATH
    DOC "The CGNS library"
    )
  FIND_LIBRARY(CGNS_LIBRARIES cgns
    DOC "The CGNS library"
    )
endif()

SET( CGNS_FOUND FALSE )

IF(CGNS_INCLUDE_DIRS)
  IF(CGNS_LIBRARIES)
    SET( CGNS_FOUND "YES" )
  ENDIF(CGNS_LIBRARIES)
ENDIF(CGNS_INCLUDE_DIRS)

IF(CGNS_FIND_REQUIRED AND NOT CGNS_FOUND)
  message(SEND_ERROR "Unable to find the requested CGNS libraries.")
ENDIF(CGNS_FIND_REQUIRED AND NOT CGNS_FOUND)

# handle the QUIETLY and REQUIRED arguments and set CGNS_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CGNS DEFAULT_MSG CGNS_LIBRARIES CGNS_INCLUDE_DIRS)

MARK_AS_ADVANCED(
  CGNS_INCLUDE_DIRS
  CGNS_LIBRARIES
)
