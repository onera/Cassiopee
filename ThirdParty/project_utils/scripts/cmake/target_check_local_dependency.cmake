
macro(target_check_local_dependency sub_repo_name mandatory)
  set(sub_repo_path "${PROJECT_ROOT}/external/${sub_repo_name}")
  file(GLOB sub_repo_files ${sub_repo_path}/*)
  list(LENGTH sub_repo_files sub_repo_nb_files)
  if(sub_repo_nb_files EQUAL 0)
    if(mandatory)
      message(FATAL_ERROR
        "${PROJECT_ROOT}/external/${sub_repo_name} is empty. Maybe you forgot to initialize it with \"git submodule update --init\""
      )
    else()
      message(WARNING
       "Did not find optional submodule ${PROJECT_ROOT}/external/${sub_repo_name}."
      )
      set(${sub_repo_name}_FOUND OFF)
    endif()
  else()
    set(${sub_repo_name}_FOUND ON)
  endif()
endmacro()