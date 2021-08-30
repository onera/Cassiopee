cmake_minimum_required(VERSION 3.12)

# https://stackoverflow.com/a/52136398/1583122
function(target_link_libraries_system visibility target)
  set(libs ${ARGN})
  foreach(lib IN LISTS ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${lib})
  endforeach(lib)
endfunction()
