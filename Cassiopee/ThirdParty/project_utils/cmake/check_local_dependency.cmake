# check_local_dependency(sub_repo_name [REQUIRED])
#     Check if sub_repo_name is non-empty at ${PROJECT_ROOT}/external/
#     If it is non-empty, ${sub_repo_name}_FOUND is set to ON, else it is set to OFF
macro(check_local_dependency sub_repo_name)
  set(sub_repo_path "${PROJECT_ROOT}/external/${sub_repo_name}")
  file(GLOB sub_repo_files ${sub_repo_path}/*)
  list(LENGTH sub_repo_files sub_repo_nb_files)

  if (${ARGC} GREATER 2)
    message(FATAL_ERROR "Error: incorrect use of check_local_dependency(sub_repo_name [REQUIRED])")
  endif()
  if ((${ARGC} EQUAL 2) AND NOT ("${ARGN}" STREQUAL "REQUIRED"))
    message(FATAL_ERROR "Error: incorrect use of check_local_dependency(sub_repo_name [REQUIRED])")
  endif()

  if (sub_repo_nb_files EQUAL 0)
    if ((${ARGC} EQUAL 2) AND ("${ARGN}" STREQUAL "REQUIRED"))
      message(FATAL_ERROR
        "${PROJECT_ROOT}/external/${sub_repo_name} is empty. Maybe you forgot to initialize it with \"git submodule update --init\""
      )
    else()
      message(WARNING
       "Did not find optional submodule external/${sub_repo_name}."
      )
      set(${sub_repo_name}_FOUND OFF)
    endif()
  else()
    set(${sub_repo_name}_FOUND ON)
  endif()
endmacro()
