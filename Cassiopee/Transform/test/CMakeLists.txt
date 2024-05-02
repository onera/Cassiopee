include(${PROJECT_UTILS_CMAKE_DIR}/doctest_pytest.cmake)

add_test(
  NAME t1
  COMMAND python test/t1.py
)

set(ld_library_path "${PROJECT_BINARY_DIR}")
set(pythonpath "${PROJECT_BINARY_DIR}:${PROJECT_SOURCE_DIR}") # binary for compiled (warping) modules, source for regular .py files
set_tests_properties(
  t1
  PROPERTIES
    LABELS "sanity_test"
    ENVIRONMENT "LD_LIBRARY_PATH=${ld_library_path}:$ENV{LD_LIBRARY_PATH};PYTHONPATH=${pythonpath}:$ENV{PYTHONPATH}"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
