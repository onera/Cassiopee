find_package(Doxygen 1.8.19 REQUIRED)
find_package(Perl REQUIRED)
# find_package(PythonInterp REQUIRED)
if(CMAKE_VERSION VERSION_LESS 3.14)
  find_package(Python REQUIRED COMPONENTS Interpreter Development)
else()
  find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)
endif()
find_package(Sphinx 3 REQUIRED)
include(FindPythonModule)
find_python_module(breathe REQUIRED)

function(add_breathe_doc)
  set(options)
  set(oneValueArgs
    SOURCE_DIR
    BUILD_DIR
    CACHE_DIR
    HTML_DIR
    DOXY_FILE
    CONF_FILE
    TARGET_NAME
    ENV_PATH
    COMMENT
  )
  set(multiValueArgs)

  cmake_parse_arguments(BREATHE_DOC
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
  )

  configure_file(
    ${BREATHE_DOC_CONF_FILE}
    ${BREATHE_DOC_BUILD_DIR}/conf.py
    @ONLY
  )

  configure_file(
    ${BREATHE_DOC_DOXY_FILE}
    ${BREATHE_DOC_BUILD_DIR}/Doxyfile
    @ONLY
  )

  add_custom_target(${BREATHE_DOC_TARGET_NAME}
    COMMAND
      ${CMAKE_COMMAND} -E env PYTHONPATH=$ENV{PYTHONPATH}:${BREATHE_DOC_ENV_PATH} ${SPHINX_EXECUTABLE}
         -q
         -b html
         -c ${BREATHE_DOC_BUILD_DIR}
         -d ${BREATHE_DOC_CACHE_DIR}
         ${BREATHE_DOC_SOURCE_DIR}
         ${BREATHE_DOC_HTML_DIR}
    COMMENT
      "Building ${BREATHE_DOC_TARGET_NAME} documentation with Breathe, Sphinx and Doxygen"
    VERBATIM
  )

  message(STATUS "Added ${BREATHE_DOC_TARGET_NAME} [Breathe+Sphinx+Doxygen] target to build documentation")
endfunction()
