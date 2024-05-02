set(SITE_PACKAGES_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/")


# compile ${cpython_module_file} into a module of the same name, depending on ${target}, and install it
function(compile_install_cpython_module target cpython_module_file)
  get_filename_component(mod_name ${cpython_module_file} NAME_WE)

  add_library(${mod_name} SHARED ${cpython_module_file})
  set_target_properties(${mod_name} PROPERTIES PREFIX "") # do not prefix by lib
  target_include_directories(${mod_name} PUBLIC
    $<BUILD_INTERFACE:${include_dir}>
    $<INSTALL_INTERFACE:include/${target}>
  )

  target_link_libraries(${mod_name}
    PUBLIC
      Python::Python
      ${target}
  )

  set_target_properties(${mod_name}
    PROPERTIES
      LINKER_LANGUAGE C
  )

# Install target
  install(TARGETS ${mod_name}
          LIBRARY DESTINATION ${SITE_PACKAGES_OUTPUT_DIRECTORY}/${target})
endfunction()


# Take all the .pybind.cpp in the source of ${project_name}
# and create a hierarchy of python modules mimicking the folder structure of ${project_name}'s source tree
# WARNING:
#    in Python, the module located at
#        ${project_name}/my_sub_folder/my_pybind_module.pybind.cpp
#    can by imported by
#        import c${project_name}.my_sub_folder.my_pybind_module
#    => NOTICE the "c" prefix before ${project_name}.
#       The "c" is needed because it is possible to have a pure python module at
#          ${project_name}/my_sub_folder/my_py_module.py
#       and then a desambiguation is needed
function(compile_install_pybind_modules project_name)
  file(GLOB_RECURSE pybind_files CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${project_name}/*.pybind.cpp")
  foreach(pybind_file ${pybind_files})
    # Relative path
    get_filename_component(pybind_dir ${pybind_file} DIRECTORY)
    file(RELATIVE_PATH pybind_dir_rel ${CMAKE_CURRENT_SOURCE_DIR}/${project_name} ${pybind_dir})
    set(pybind_dir_rel c${project_name}/${pybind_dir_rel})

    # Create target
    get_filename_component(mod_name ${pybind_file} NAME_WE)
    pybind11_add_module(${mod_name} ${pybind_file})
    target_link_libraries(${mod_name} PUBLIC ${project_name}::${project_name})
    set_target_properties(${mod_name} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${pybind_dir_rel}")

    # Install target
    install(TARGETS ${mod_name}
            LIBRARY DESTINATION ${SITE_PACKAGES_OUTPUT_DIRECTORY}/${pybind_dir_rel})
  endforeach()
endfunction()


function(compile_install_cython_modules project_name)
  include(UseCython)
  file(GLOB_RECURSE pyx_files CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${project_name}/*.pyx")
  foreach(pyx_file ${pyx_files})
    # Relative path
    get_filename_component(pyx_dir ${pyx_file} DIRECTORY )
    file(RELATIVE_PATH pyx_dir_rel ${CMAKE_CURRENT_SOURCE_DIR} ${pyx_dir})

    # Create target
    get_filename_component(mod_name ${pyx_file} NAME_WE) # If same name : problem
    cython_add_one_file_module(${mod_name} ${pyx_file})
    target_link_libraries(${mod_name} PUBLIC ${project_name}::${project_name})
    set_target_properties(${mod_name} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${pyx_dir_rel}")
    set_source_files_properties(${pyx_file} PROPERTIES CYTHON_IS_CXX TRUE)

    # Install target
    install(TARGETS ${mod_name}
            LIBRARY DESTINATION ${SITE_PACKAGES_OUTPUT_DIRECTORY}/${pyx_dir_rel})
  endforeach()
endfunction()


function(compile_install_swig_modules project_name)
  set(UseSWIG_TARGET_NAME_PREFERENCE STANDARD) # SWIG: use standard target name (no underscore)
  find_package(SWIG 3.0.12 REQUIRED)
  include(${PROJECT_UTILS_CMAKE_DIR}/find_package/UseSWIG-fixed.cmake) # fix to circunvent properties not propagated to dependencies

  file(GLOB_RECURSE swig_files CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${project_name}/*.i")
  foreach(swig_file ${swig_files})
    # Relative path
    get_filename_component(swig_dir ${swig_file} DIRECTORY )
    file(RELATIVE_PATH swig_dir_rel ${CMAKE_CURRENT_SOURCE_DIR} ${swig_dir})
    set(swig_dir_rel ${swig_dir_rel}_swig)

    # Create target
    get_filename_component(mod_name ${swig_file} NAME_WE) # If same name : problem
    set(swig_mod_name ${mod_name}_swig)
    set_property(SOURCE ${swig_file} PROPERTY CPLUSPLUS ON)
    set(CMAKE_SWIG_FLAGS -module ${swig_mod_name}) # prevent the need to prefix by _ when importing the module in python
    swig_add_library(
      ${swig_mod_name}
      TYPE MODULE
      LANGUAGE python
      OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/${swig_dir_rel}"
      SOURCES
        ${swig_file}
    )
    target_include_directories(${swig_mod_name}
      PUBLIC
        $<BUILD_INTERFACE:${include_dir}>
        $<INSTALL_INTERFACE:include/${project_name}>
    )
    target_link_libraries(${swig_mod_name}
      PUBLIC
        ${project_name}
    )
    set_target_properties(${swig_mod_name} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${swig_dir_rel}")

    # Install
    install(TARGETS ${swig_mod_name}
            LIBRARY DESTINATION ${SITE_PACKAGES_OUTPUT_DIRECTORY}/${swig_dir_rel})
    install(FILES       "${CMAKE_CURRENT_BINARY_DIR}/${swig_dir_rel}/${swig_mod_name}.py"
            DESTINATION "${SITE_PACKAGES_OUTPUT_DIRECTORY}/${swig_dir_rel}"
            COMPONENT   "Python file generated by SWIG")
  endforeach()
endfunction()


function(install_python_modules project_name)
  file(GLOB_RECURSE py_files CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${project_name}/*.py")

  foreach (py_file ${py_files})
    file(RELATIVE_PATH py_file_rel  ${CMAKE_CURRENT_SOURCE_DIR} ${py_file})

    get_filename_component(py_dir_rel "${py_file_rel}" DIRECTORY)
    install(FILES       "${py_file_rel}"
            DESTINATION "${SITE_PACKAGES_OUTPUT_DIRECTORY}/${py_dir_rel}"
            COMPONENT   "python")
  endforeach ()
endfunction()
