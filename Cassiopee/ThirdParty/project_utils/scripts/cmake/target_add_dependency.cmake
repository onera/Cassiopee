cmake_minimum_required(VERSION 3.12)

# create global variables storing a target dependencies
# these variables are useful later when installing the target
# note: both append_to_target_dependency_list and append_to_target_thirdparty_dependency_list
#       append to the same
#         ${target}_DEPENDENCIES_FIND_PACKAGE_STRING
#       but append_to_target_dependency_list appends to
#         ${target}_DEPENDENCIES_STRING
#       and append_to_target_thirdparty_dependency_list appends to
#         ${target}_THIRDPARTY_DEPENDENCIES_STRING
macro(append_to_target_dependency_list target)
  set(dependency ${ARGN})
  list(APPEND ${target}_DEPENDENCIES_FIND_PACKAGE_STRING "find_package(${dependency})\n")
  list(APPEND ${target}_DEPENDENCIES_STRING "\"${dependency}\",")
endmacro()
macro(append_to_target_thirdparty_dependency_list target) # TODO rename target->project
  set(dependency ${ARGN})
  list(APPEND ${target}_DEPENDENCIES_FIND_PACKAGE_STRING "find_package(${dependency})\n")
  list(APPEND ${target}_THIRDPARTY_DEPENDENCIES_STRING "\"${dependency}\",")
endmacro()

# target_add_dependency
# add the subdirectory ${dependency} located in ${project_root}/external/
# the string ${target}_DEPENDENCIES_FIND_PACKAGE_STRING is appended the corresponding find_package() command
#   the idea is that we will be able to use this string
#   when adding dependencies to the ${target}Config.cmake file further down the installation process
macro(target_add_dependency target dependency)
  append_to_target_dependency_list(${target} ${dependency})
  if(NOT TARGET ${dependency})
    add_subdirectory(${PROJECT_ROOT}/external/${dependency} ${CMAKE_BINARY_DIR}/external/${dependency})
  endif()
endmacro()

# target_add_dependencies
# same as target_add_dependency except for several dependencies
macro(target_add_dependencies target)
  set(dependencies ${ARGN})
  foreach(dep IN ITEMS ${dependencies})
    target_add_dependency(${target} ${dep})
  endforeach()
endmacro()

# target_add_thirdparty_dependency
# same as target_add_dependency except we call find_package instead of add_directory
macro(target_add_thirdparty_dependency target)
  append_to_target_thirdparty_dependency_list(${target} ${ARGN})
  find_package(${ARGN})
endmacro()
# same but better name # TODO remove old name
macro(project_find_package)
  append_to_target_thirdparty_dependency_list(${PROJECT_NAME} ${ARGV})
  find_package(${ARGV})
endmacro()

# find_and_link_optional_dependency
# call target_add_thirdparty_dependency and if dependency is found, add it to the target_link_libraries()
macro(find_and_link_optional_dependency target package_name)
  target_add_thirdparty_dependency(${target} ${package_name})
  if (${package_name}_FOUND)
    target_link_libraries(${target} ${ARGN})
  endif()
endmacro()

macro(target_install target)
  if(NOT DEFINED PROJECT_ROOT)
    set(PROJECT_ROOT ${CMAKE_SOURCE_DIR} CACHE PATH "Root directory, where the submodules are populated")
  endif()
  set(PROJECT_UTILS_CMAKE_DIR ${PROJECT_ROOT}/external/project_utils/scripts/cmake)

  install(TARGETS ${target} EXPORT ${target}Targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
  )
  install(EXPORT ${target}Targets
    FILE ${target}Targets.cmake
    NAMESPACE ${target}::
    DESTINATION lib/cmake/${target}
  )
  install(DIRECTORY ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}
    DESTINATION include
    FILES_MATCHING
      PATTERN "*.h"
      PATTERN "*.hpp"
      PATTERN "*.hxx"
      PATTERN "*.cxx"
  )

  set(TARGET_NAME ${target}) # WARNING Seems not used but actually used in target_config.cmake.in
  string(REPLACE ";" " " TARGET_DEPENDENCIES_FIND_PACKAGE_STRING "${${target}_DEPENDENCIES_FIND_PACKAGE_STRING}") # Same, used below
  configure_file(
    ${PROJECT_UTILS_CMAKE_DIR}/target_config.cmake.in
    ${target}Config.cmake
    @ONLY
  )
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${target}Config.cmake"
    DESTINATION lib/cmake/${target}
  )

  add_library(${target}::${target} ALIAS ${target})
endmacro()
