cmake_host_system_information(RESULT HOSTNAME QUERY HOSTNAME)

option(${PROJECT_NAME}_USE_DEFAULT_FLAGS "Enable default compiler flags" OFF)
if (${PROJECT_NAME}_USE_DEFAULT_FLAGS)
  if ($ENV{CXXFLAGS})
    message(WARNING "CXXFLAGS is not empty but ${PROJECT_NAME}_USE_DEFAULT_FLAGS is ON. Your C++ compiler flags will be overriden")
  endif()

  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wpointer-arith -Wcast-align -fmax-errors=4 ")

  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set (CMAKE_CXX_FLAGS "-Wall -Wcheck -Wpointer-arith")

  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wpointer-arith -Wno-missing-braces") # missing-braces disabled, else GCC warning
  endif()
endif()
