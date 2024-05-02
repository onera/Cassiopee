function(create_mpi_doctest)
  include(CTest)
  set(options)
  set(one_value_args)
  set(multi_value_args TESTED_TARGET LABEL SOURCES SERIAL_RUN N_PROC)
  cmake_parse_arguments(ARGS "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  set(n_proc ${ARGS_N_PROC})
  set(tested_target ${ARGS_TESTED_TARGET})

  set(executable_for_tests ${tested_target}_${ARGS_LABEL})
  add_executable(${executable_for_tests} ${ARGS_SOURCES})

  target_link_libraries(${executable_for_tests}
    PUBLIC
      ${tested_target}::${tested_target}
    PRIVATE
      doctest::doctest
  )

  install(TARGETS ${executable_for_tests} RUNTIME DESTINATION bin)
  if(${ARGS_SERIAL_RUN})
    add_test(
      NAME ${executable_for_tests}
      COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${executable_for_tests}
    )
  else()
    add_test(
      NAME ${executable_for_tests}
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${n_proc}
              ${MPIEXEC_PREFLAGS}
              ${CMAKE_CURRENT_BINARY_DIR}/${executable_for_tests}
              ${MPIEXEC_POSTFLAGS}
    )
  endif()

  set_tests_properties(${executable_for_tests} PROPERTIES LABELS "${ARGS_LABEL}")
  set_tests_properties(${executable_for_tests} PROPERTIES PROCESSORS ${n_proc})
  set_tests_properties(${executable_for_tests} PROPERTIES SERIAL_RUN ${ARGS_SERIAL_RUN})

  # Fails in non-slurm
  # set_tests_properties(${target_name} PROPERTIES PROCESSOR_AFFINITY true)
endfunction()


# --------------------------------------------------------------------------------
function(create_mpi_pytest name n_proc)
  set(options)
  set(one_value_args)
  set(multi_value_args SOURCES LABELS SERIAL_RUN)
  cmake_parse_arguments(ARGS "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  set(output_python_file ${CMAKE_CURRENT_BINARY_DIR}/${name}.py)
  add_custom_command(OUTPUT  "${output_python_file}"
                     DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${name}.py"
                     COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                     "${CMAKE_CURRENT_SOURCE_DIR}/${name}.py"
                     "${output_python_file}"
                     COMMENT "Copying ${name} to the binary directory")

  string(REPLACE / __ target_name t_${name})
  add_custom_target(${target_name} ALL DEPENDS "${output_python_file}" )

  # WORKING_DIRECTORY
  add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${n_proc}
            ${MPIEXEC_PREFLAGS}
            ${Python_EXECUTABLE} -m pytest -Wignore -r a -v -s ${output_python_file}
            ${MPIEXEC_POSTFLAGS})

  # > Set properties for the current test
  # pytest test/python_unit_tests.py --html=test.html --self-contained-html
  if(NOT DEFINED PROJECT_ROOT)
    set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})
  endif()

  set_tests_properties(${name} PROPERTIES LABELS "${ARGS_LABELS}")
  set_tests_properties("${name}" PROPERTIES
                       ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}:${PROJECT_BINARY_DIR}/test:${CMAKE_BINARY_DIR}/external/pytest-mpi-check:$ENV{PYTHONPATH}
                       DEPENDS t_${name})

  set_property(TEST "${name}" APPEND PROPERTY
                       ENVIRONMENT LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}:$ENV{LD_LIBRARY_PATH})
  set_property(TEST "${name}" APPEND PROPERTY
                       ENVIRONMENT PYTEST_PLUGINS=pytest_mpi_check)
  set_tests_properties(${name} PROPERTIES PROCESSORS n_proc)
  if(${ARGS_SERIAL_RUN})
    set_tests_properties(${name} PROPERTIES RUN_SERIAL true)
  endif()
endfunction()
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
function(mpi_pytest_directory_python_create name tested_dir n_proc)
  set(options)
  set(one_value_args)
  set(multi_value_args SOURCES LABELS SERIAL_RUN)
  cmake_parse_arguments(ARGS "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  # --------------------------------------------------------------------------------
  file(GLOB_RECURSE py_test_files
       CONFIGURE_DEPENDS
       "${PROJECT_SOURCE_DIR}/${tested_dir}/test_*.py")

  foreach(py_test_file ${py_test_files})
    get_filename_component( output_python_file_name ${py_test_file} NAME_WE   )
    get_filename_component( output_python_dir       ${py_test_file} DIRECTORY )
    file(RELATIVE_PATH output_python_dir_rel ${CMAKE_CURRENT_SOURCE_DIR} ${output_python_dir})

    set(output_python_file ${CMAKE_CURRENT_BINARY_DIR}/${output_python_dir_rel}/${output_python_file_name}.py)
    add_custom_command(OUTPUT  "${output_python_file}"
                       DEPENDS "${py_test_file}"
                       COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                       "${py_test_file}"
                       "${output_python_file}"
                       COMMENT "Copying ${output_python_dir_rel}/${output_python_file_name}.py to the binary directory")

    set(rel_output_python_file ${output_python_dir_rel}/${output_python_file_name}.py)
    string(REPLACE "/" "_" target_name "${rel_output_python_file}")
    add_custom_target(t_${target_name} ALL DEPENDS "${output_python_file}")
  endforeach()
  # --------------------------------------------------------------------------------

  # --------------------------------------------------------------------------------
  file(GLOB_RECURSE __conftest_files
       CONFIGURE_DEPENDS
       "${PROJECT_SOURCE_DIR}/${tested_dir}/*conftest.py")

  foreach(__conftest_file ${__conftest_files})
    get_filename_component( output_python_file_name ${__conftest_file} NAME_WE   )
    get_filename_component( output_python_dir       ${__conftest_file} DIRECTORY )
    file(RELATIVE_PATH output_python_dir_rel ${CMAKE_CURRENT_SOURCE_DIR} ${output_python_dir})

    set(output_python_file ${CMAKE_CURRENT_BINARY_DIR}/${output_python_dir_rel}/${output_python_file_name}.py)
    add_custom_command(OUTPUT  "${output_python_file}"
                       DEPENDS "${__conftest_file}"
                       COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                       "${__conftest_file}"
                       "${output_python_file}"
                       COMMENT "Copying ${output_python_dir_rel}/${output_python_file_name}.py to the binary directory")

    set(rel_output_python_file ${output_python_dir_rel}/${output_python_file_name}.py)
    string(REPLACE "/" "_" target_name "${rel_output_python_file}")
    add_custom_target(t_${target_name} ALL DEPENDS "${output_python_file}")
  endforeach()
  # --------------------------------------------------------------------------------

  # --------------------------------------------------------------------------------
  file(GLOB_RECURSE __pytestini_files
       CONFIGURE_DEPENDS
       "${PROJECT_SOURCE_DIR}/${tested_dir}/*pytest.ini")

  foreach(__pytestini_file ${__pytestini_files})
    get_filename_component( output_pytestini_file_name ${__pytestini_file} NAME_WE   )
    get_filename_component( output_pytestini_dir       ${__pytestini_file} DIRECTORY )
    file(RELATIVE_PATH output_pytestini_dir_rel ${CMAKE_CURRENT_SOURCE_DIR} ${output_pytestini_dir})

    set(output_pytestini_file ${CMAKE_CURRENT_BINARY_DIR}/${output_pytestini_dir_rel}/${output_pytestini_file_name}.ini)
    add_custom_command(OUTPUT  "${output_pytestini_file}"
                       DEPENDS "${__pytestini_file}"
                       COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                       "${__pytestini_file}"
                       "${output_pytestini_file}"
                       COMMENT "Copying ${output_pytestini_dir_rel}/${output_pytestini_file_name}.ini to the binary directory")
    set(rel_output_pytestini_file ${output_pytestini_dir_rel}/${output_pytestini_file_name}.ini)
    string(REPLACE "/" "_" target_name "${rel_output_pytestini_file}")
    add_custom_target(t_${target_name} ALL DEPENDS "${output_pytestini_file}")
  endforeach()
  # --------------------------------------------------------------------------------


  # -r : display a short test summary info, with a == all except passed (i.e. report failed, skipped, error)
  # -s : no capture (print statements output to stdout)
  # -v : verbose
  # -Wignore : ignore warnings
  add_test (${name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${n_proc}
            ${MPIEXEC_PREFLAGS}
            ${Python_EXECUTABLE} -m pytest ${PROJECT_BINARY_DIR}/${tested_dir} -Wignore -ra -v -s --with-mpi
            ${MPIEXEC_POSTFLAGS})

  # > Set properties for the current test
  set_tests_properties(${name} PROPERTIES LABELS "${ARGS_LABELS}")
  set_tests_properties(${name} PROPERTIES
                       ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}:${PROJECT_BINARY_DIR}/${tested_dir}:${CMAKE_BINARY_DIR}/external/pytest-mpi-check:$ENV{PYTHONPATH}
                       DEPENDS t_${name})
  set_property(TEST ${name} APPEND PROPERTY
                       ENVIRONMENT LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}:$ENV{LD_LIBRARY_PATH})
  set_property(TEST ${name} APPEND PROPERTY
                       ENVIRONMENT PYTEST_PLUGINS=pytest_mpi_check)
  set_tests_properties(${name} PROPERTIES PROCESSORS n_proc)
  if(${ARGS_SERIAL_RUN})
    set_tests_properties(${name} PROPERTIES RUN_SERIAL true)
  endif()
  # > Not working if not launch with srun ...
  # set_tests_properties(${name} PROPERTIES PROCESSOR_AFFINITY true)
endfunction()
# --------------------------------------------------------------------------------



## --------------------------------------------------------------------------------
#function(seq_test_python_create name)
#  # configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${name}.py ${CMAKE_CURRENT_BINARY_DIR}/${name}.py COPYONLY)
#  set(output_python_file ${CMAKE_CURRENT_BINARY_DIR}/${name}.py)
#  add_custom_command(OUTPUT  "${output_python_file}"
#                     DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${name}.py"
#                     COMMAND "${CMAKE_COMMAND}" -E copy_if_different
#                     "${CMAKE_CURRENT_SOURCE_DIR}/${name}.py"
#                     "${output_python_file}"
#                     COMMENT "Copying ${name} to the binary directory")
#
#  add_custom_target(tpyseq_${name} ALL DEPENDS "${output_python_file}")
#
#
#  add_test (${name}
#             ${Python_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/${name}.py)
#  set_tests_properties("${name}" PROPERTIES
#                       ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/mod:$ENV{PYTHONPATH}
#                       DEPENDS tpyseq_${name})
#
#endfunction()
## --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
function(mpi_test_create target_file name tested_target n_proc )
  set(options)
  set(one_value_args)
  set(multi_value_args SOURCES INCLUDES LIBRARIES LABELS SERIAL_RUN)
  cmake_parse_arguments(ARGS "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  add_executable(${name} ${target_file} ${ARGS_SOURCES})

  target_include_directories(${name} PRIVATE ${ARGS_INCLUDES})
  target_link_libraries(${name} ${ARGS_LIBRARIES})
  target_link_libraries(${name} ${tested_target}::${tested_target})

  install(TARGETS ${name} RUNTIME DESTINATION bin)
  add_test (NAME ${name}
            COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${n_proc}
                    ${MPIEXEC_PREFLAGS}
                    ${CMAKE_CURRENT_BINARY_DIR}/${name}
                    ${MPIEXEC_POSTFLAGS})
  # add_test (NAME ${name}
  #           COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${name})

  # > Set properties for the current test
  set_tests_properties(${name} PROPERTIES LABELS "${ARGS_LABELS}")
  set_tests_properties(${name} PROPERTIES PROCESSORS nproc)
  if(${ARGS_SERIAL_RUN})
    set_tests_properties(${name} PROPERTIES RUN_SERIAL true)
  endif()
  # > Fail in non slurm
  # set_tests_properties(${name} PROPERTIES PROCESSOR_AFFINITY true)

  # > Specific environement :
  # set_tests_properties(${name} PROPERTIES ENVIRONMENT I_MPI_DEBUG=5)

endfunction()
# --------------------------------------------------------------------------------
