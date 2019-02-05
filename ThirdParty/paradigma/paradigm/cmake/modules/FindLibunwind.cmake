# Find the libunwind library
#
# LIBUNWIND_FOUND - True if libunwind was found.
# LIBUNWIND_LIBRARIES - The libraries needed to use libunwind
# LIBUNWIND_INCLUDE_DIRS - Location of unwind.h and libunwind.h

set(LIBUNWIND_FOUND FALSE)

find_path(LIBUNWIND_INCLUDE_DIRS libunwind.h
  HINTS ${LIBUNWIND_DIR}/include $ENV{LIBUNWIND_DIR}/include
  DOC "Directory where the Libunwind header is located"
  )

find_library(LIBUNWIND_LIBRARIES unwind
  HINTS ${LIBUNWIND_DIR}/lib $ENV{LIBUNWIND_DIR}/lib
  NO_DEFAULT_PATH
  DOC "The Libunwind header library"
  )

find_library(LIBUNWIND_LIBRARIES unwind
  DOC "The Libunwind header library"
  )

if (LIBUNWIND_LIBRARIES)
  string(REGEX REPLACE "(^.*)/libunwind.*$" "\\1" LIBUNWIND_LIBRARY_PATH ${LIBUNWIND_LIBRARIES} )
endif (LIBUNWIND_LIBRARIES)

IF(LIBUNWIND_FIND_REQUIRED AND NOT LIBUNWIND_FOUND)
  message(SEND_ERROR "Unable to find the requested Libunwind libraries.")
ENDIF(LIBUNWIND_FIND_REQUIRED AND NOT LIBUNWIND_FOUND)
 
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIBUNWIND DEFAULT_MSG LIBUNWIND_LIBRARIES LIBUNWIND_INCLUDE_DIRS LIBUNWIND_LIBRARY_PATH)

MARK_AS_ADVANCED(
  LIBUNWIND_INCLUDE_DIRS
  LIBUNWIND_LIBRARIES
  LIBUNWIND_LIBRARY_PATH
)

