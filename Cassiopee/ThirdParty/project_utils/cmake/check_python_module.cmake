# check_python_module
# also see https://cmake.org/pipermail/cmake/2011-January/041666.html
function(check_python_module module)
  if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
    set(is_required TRUE)
  endif()
  execute_process(COMMAND "${Python_EXECUTABLE}" "-c" "import ${module}"
    RESULT_VARIABLE module_status
  )
  if (NOT ${module_status} EQUAL 0)
    if (is_required)
      message(FATAL_ERROR " Could not find \"${module}\" python required module")
    else()
      message(WARNING "Could not find \"${module}\" python optional module")
    endif()
  endif()
endfunction()
