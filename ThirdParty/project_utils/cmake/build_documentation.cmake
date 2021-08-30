# SEE https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/
macro(build_documentation)
# 0. Doxygen
  find_package(Doxygen REQUIRED)

  file(GLOB_RECURSE HEADERS ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/*.hpp)

  set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/${PROJECT_NAME})
  set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen)
  set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
  set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile)
  # replace @DOXYGEN_INPUT_DIR@ and @DOXYGEN_OUTPUT_DIR@ values in DOXYFILE_IN and output it to DOXYFILE_OUT
  configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

  set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)

  file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) # Doxygen won't create this for us
  add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                     DEPENDS ${HEADERS}
                     COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                     MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                     COMMENT "Generating Doxygen documentation")

  add_custom_target(${PROJECT_NAME}_doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

# 1. Sphinx
  find_package(Sphinx 3 REQUIRED)
  set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/doc)
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/doc/sphinx/html)
  set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

  set(SPHINX_CONF_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py.in)
  set(SPHINX_CONF_OUT ${CMAKE_CURRENT_BINARY_DIR}/doc/conf.py)
  # replace @PROJECT_NAME@ and @PROJECT_SOURCE_DIR@ values in SPHINX_CONF_IN and output it to SPHINX_CONF_OUT
  configure_file(${SPHINX_CONF_IN} ${SPHINX_CONF_OUT} @ONLY)

  file(GLOB_RECURSE doc_files ${CMAKE_CURRENT_SOURCE_DIR}/doc/*)
  add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                     COMMAND ${SPHINX_EXECUTABLE} -b html -c ${CMAKE_CURRENT_BINARY_DIR}/doc
                     -Dbreathe_projects.${PROJECT_NAME}=${DOXYGEN_OUTPUT_DIR}/xml # Tell Breathe where to find the Doxygen output
                     ${SPHINX_SOURCE} ${SPHINX_BUILD}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     DEPENDS
                       ${doc_files}
                       #${CMAKE_CURRENT_SOURCE_DIR}/index.rst
                       ${DOXYGEN_INDEX_FILE}
                     MAIN_DEPENDENCY ${SPHINX_CONF_OUT}
                     COMMENT "Generating Sphinx documentation, using Breathe to recover xml files from Doxygen")

  add_custom_target(${PROJECT_NAME}_sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})

# 2. Install
  install(DIRECTORY ${SPHINX_BUILD}
          DESTINATION doc/${PROJECT_NAME})
endmacro()
